#!/usr/bin/env python3
"""
Calculate sensory-axon -> ExN-dendrite path overlap from close-contact NPZ files.

Expected server layout
----------------------
misc/
├── batch_calculate_sensory_exn_overlap.py
├── all_possible_sensory_exn_pairs_with_overlap.csv
├── overlap_results/
│   └── close_contacts_<sensory_id>_<exn_id>.npz
└── path_calculation_log/              # created automatically

Usage
-----
Test one exact file:
    python batch_calculate_sensory_exn_overlap.py \
        close_contacts_720575940862185831_720575940903188925.npz

Process every NPZ in misc/overlap_results:
    python batch_calculate_sensory_exn_overlap.py

Useful options:
    python batch_calculate_sensory_exn_overlap.py --workers 8
    python batch_calculate_sensory_exn_overlap.py --csv another_file.csv
    python batch_calculate_sensory_exn_overlap.py --no-update-csv

The script calculates only the required direction:
    root_id1 axon -> root_id2 dendrite

The overlap calculation reproduces the final axon-length calculation in
visualiza_overlap.ipynb: it measures the induced skeleton-edge length of the
union of neuron-1 axon vertices participating in close pairs, with the same
half-edge boundary correction.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_NPZ_DIR = SCRIPT_DIR / "overlap_results"
DEFAULT_CSV = SCRIPT_DIR / "all_possible_sensory_exn_pairs_with_overlap.csv"
DEFAULT_LOG_DIR = SCRIPT_DIR / "path_calculation_log"
DEFAULT_WORKERS = min(8, os.cpu_count() or 1)
DEFAULT_CHECKPOINT_EVERY = 1000

OVERLAP_NM_COL = "sensory_axon_exn_dendrite_overlap_nm"
OVERLAP_UM_COL = "sensory_axon_exn_dendrite_overlap_um"
N_CLOSE_PAIRS_COL = "n_close_vertex_pairs_1axon_2dend"
N_CLOSE_AXON_VERTICES_COL = "n_close_sensory_axon_vertices"
PAIR_CAPPED_COL = "close_pairs_possibly_capped"
STATUS_COL = "path_overlap_status"
ERROR_COL = "path_overlap_error"
SOURCE_FILE_COL = "path_overlap_npz"

FILENAME_RE = re.compile(r"^close_contacts_(\d+)_(\d+)\.npz$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate root_id1 axon -> root_id2 dendrite path overlap from "
            "one NPZ file or every NPZ in overlap_results."
        )
    )
    parser.add_argument(
        "npz_file",
        nargs="?",
        help=(
            "Optional NPZ filename or path. If omitted, all "
            "close_contacts_*.npz files are processed."
        ),
    )
    parser.add_argument(
        "--npz-dir",
        type=Path,
        default=DEFAULT_NPZ_DIR,
        help=f"NPZ directory (default: {DEFAULT_NPZ_DIR})",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help=f"CSV to update (default: {DEFAULT_CSV})",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help=f"Log directory (default: {DEFAULT_LOG_DIR})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Worker processes for batch mode (default: {DEFAULT_WORKERS})",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=DEFAULT_CHECKPOINT_EVERY,
        help=(
            "Write a batch checkpoint after this many completed files "
            f"(default: {DEFAULT_CHECKPOINT_EVERY})."
        ),
    )
    parser.add_argument(
        "--no-update-csv",
        action="store_true",
        help="Calculate and log results without modifying the CSV.",
    )
    return parser.parse_args()


def configure_logging(log_dir: Path, mode: str) -> tuple[logging.Logger, Path]:
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = log_dir / f"path_calculation_{mode}_{timestamp}.log"

    logger = logging.getLogger("path_overlap")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger, log_path


def resolve_npz_path(npz_argument: str, npz_dir: Path) -> Path:
    supplied = Path(npz_argument).expanduser()

    if supplied.is_file():
        return supplied.resolve()

    candidate = npz_dir.expanduser().resolve() / supplied.name
    if candidate.is_file():
        return candidate

    raise FileNotFoundError(
        f"Could not find '{npz_argument}' directly or inside '{npz_dir}'."
    )


def build_adjacency(n_vertices: int, edges: np.ndarray) -> list[list[int]]:
    adjacency: list[list[int]] = [[] for _ in range(n_vertices)]
    for u, v in edges:
        u_int = int(u)
        v_int = int(v)
        adjacency[u_int].append(v_int)
        adjacency[v_int].append(u_int)
    return adjacency


def induced_path_length_with_half_edges_nm(
    coords_nm: np.ndarray,
    edges: np.ndarray,
    selected_vertices: np.ndarray,
) -> float:
    """
    Match induced_geodesic_diameter() length logic from visualiza_overlap.ipynb.

    Despite that notebook function's name, its returned length is:
      1. the sum of every skeleton edge fully inside the selected vertex set;
      2. plus half of the shortest outward edge at every induced endpoint;
      3. for a singleton, half of up to its two shortest outward edges.
    """
    coords_nm = np.asarray(coords_nm, dtype=float)
    edges = np.asarray(edges, dtype=np.int64)
    selected = np.unique(np.asarray(selected_vertices, dtype=np.int64))

    if selected.size == 0:
        return 0.0

    n_vertices = len(coords_nm)
    if selected.min() < 0 or selected.max() >= n_vertices:
        raise IndexError("A selected skeleton vertex index is out of range.")

    if edges.size == 0:
        return 0.0
    if edges.ndim != 2 or edges.shape[1] != 2:
        raise ValueError(f"skel1_edges must have shape (N, 2), got {edges.shape}.")
    if edges.min() < 0 or edges.max() >= n_vertices:
        raise IndexError("A skeleton edge vertex index is out of range.")

    allowed = np.zeros(n_vertices, dtype=bool)
    allowed[selected] = True

    u = edges[:, 0]
    v = edges[:, 1]
    edge_lengths = np.linalg.norm(coords_nm[u] - coords_nm[v], axis=1)

    internal = allowed[u] & allowed[v]
    internal_length = float(edge_lengths[internal].sum())

    induced_degree = np.zeros(n_vertices, dtype=np.int64)
    np.add.at(induced_degree, u[internal], 1)
    np.add.at(induced_degree, v[internal], 1)

    adjacency = build_adjacency(n_vertices, edges)

    def edge_length(a: int, b: int) -> float:
        return float(np.linalg.norm(coords_nm[a] - coords_nm[b]))

    half_extension = 0.0

    if selected.size == 1:
        node = int(selected[0])
        outward_lengths = sorted(
            edge_length(node, neighbor)
            for neighbor in adjacency[node]
            if not allowed[neighbor]
        )
        if len(outward_lengths) >= 2:
            half_extension = 0.5 * (outward_lengths[0] + outward_lengths[1])
        elif len(outward_lengths) == 1:
            half_extension = 0.5 * outward_lengths[0]
    else:
        endpoints = selected[induced_degree[selected] == 1]
        for node_value in endpoints:
            node = int(node_value)
            outward_lengths = [
                edge_length(node, neighbor)
                for neighbor in adjacency[node]
                if not allowed[neighbor]
            ]
            if outward_lengths:
                half_extension += 0.5 * min(outward_lengths)

    return float(internal_length + half_extension)


def calculate_one_npz(npz_path: str | Path) -> dict:
    """Worker-safe calculation for one NPZ."""
    path = Path(npz_path)
    result = {
        "root_id1": None,
        "root_id2": None,
        OVERLAP_NM_COL: np.nan,
        OVERLAP_UM_COL: np.nan,
        N_CLOSE_PAIRS_COL: None,
        N_CLOSE_AXON_VERTICES_COL: None,
        PAIR_CAPPED_COL: False,
        STATUS_COL: "error",
        ERROR_COL: "",
        SOURCE_FILE_COL: path.name,
    }

    try:
        with np.load(path, allow_pickle=False) as data:
            required = {
                "root_id1",
                "root_id2",
                "skel1_vertices",
                "skel1_edges",
                "global_pairs_1axon_2dend",
            }
            missing = sorted(required.difference(data.files))
            if missing:
                raise KeyError(f"Missing NPZ keys: {missing}")

            root_id1 = int(np.asarray(data["root_id1"]).ravel()[0])
            root_id2 = int(np.asarray(data["root_id2"]).ravel()[0])
            coords1_nm = np.asarray(data["skel1_vertices"], dtype=float)
            edges1 = np.asarray(data["skel1_edges"], dtype=np.int64)
            pairs12 = np.asarray(
                data["global_pairs_1axon_2dend"], dtype=np.int64
            )

            result["root_id1"] = root_id1
            result["root_id2"] = root_id2

            filename_match = FILENAME_RE.match(path.name)
            if filename_match:
                filename_id1 = int(filename_match.group(1))
                filename_id2 = int(filename_match.group(2))
                if (filename_id1, filename_id2) != (root_id1, root_id2):
                    raise ValueError(
                        "Filename IDs do not match root_id1/root_id2 stored in NPZ."
                    )

            if pairs12.size == 0:
                pairs12 = np.empty((0, 2), dtype=np.int64)
            elif pairs12.ndim != 2 or pairs12.shape[1] != 2:
                raise ValueError(
                    "global_pairs_1axon_2dend must have shape (N, 2); "
                    f"got {pairs12.shape}."
                )

            close_axon_vertices = (
                np.unique(pairs12[:, 0])
                if len(pairs12)
                else np.empty(0, dtype=np.int64)
            )

            overlap_nm = induced_path_length_with_half_edges_nm(
                coords_nm=coords1_nm,
                edges=edges1,
                selected_vertices=close_axon_vertices,
            )

            result[OVERLAP_NM_COL] = overlap_nm
            result[OVERLAP_UM_COL] = overlap_nm / 1000.0
            result[N_CLOSE_PAIRS_COL] = int(len(pairs12))
            result[N_CLOSE_AXON_VERTICES_COL] = int(len(close_axon_vertices))

            # The original generator capped stored pairs at 5,000.
            # Exactly 5,000 therefore means the file may have been downsampled.
            result[PAIR_CAPPED_COL] = bool(len(pairs12) >= 5000)
            result[STATUS_COL] = "ok"

    except Exception as exc:  # keep batch processing alive
        result[ERROR_COL] = f"{type(exc).__name__}: {exc}"

    return result


def detect_id_columns(df: pd.DataFrame) -> tuple[str, str]:
    candidate_pairs = [
        ("pre_id", "post_id"),
        ("pre_pt_root_id", "post_pt_root_id"),
        ("sensory_id", "exn_id"),
        ("root_id1", "root_id2"),
    ]

    for id1_col, id2_col in candidate_pairs:
        if id1_col in df.columns and id2_col in df.columns:
            return id1_col, id2_col

    raise KeyError(
        "Could not identify the neuron-ID columns in the CSV. Expected one "
        "of: (pre_id, post_id), (pre_pt_root_id, post_pt_root_id), "
        "(sensory_id, exn_id), or (root_id1, root_id2)."
    )


def read_pair_csv(csv_path: Path) -> tuple[pd.DataFrame, str, str]:
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    id1_col, id2_col = detect_id_columns(df)

    df[id1_col] = pd.to_numeric(df[id1_col], errors="raise").astype("int64")
    df[id2_col] = pd.to_numeric(df[id2_col], errors="raise").astype("int64")

    duplicate_count = int(df.duplicated([id1_col, id2_col]).sum())
    if duplicate_count:
        raise ValueError(
            f"CSV contains {duplicate_count} duplicate sensory-ExN ID pairs."
        )

    return df, id1_col, id2_col


def ensure_result_columns(df: pd.DataFrame) -> pd.DataFrame:
    defaults = {
        OVERLAP_NM_COL: np.nan,
        OVERLAP_UM_COL: np.nan,
        N_CLOSE_PAIRS_COL: pd.NA,
        N_CLOSE_AXON_VERTICES_COL: pd.NA,
        PAIR_CAPPED_COL: pd.NA,
        STATUS_COL: pd.NA,
        ERROR_COL: pd.NA,
        SOURCE_FILE_COL: pd.NA,
    }
    for column, default in defaults.items():
        if column not in df.columns:
            df[column] = default
    return df


def apply_results_to_csv_dataframe(
    df: pd.DataFrame,
    id1_col: str,
    id2_col: str,
    results: Iterable[dict],
    logger: logging.Logger,
) -> tuple[pd.DataFrame, int]:
    df = ensure_result_columns(df)
    pair_to_index = {
        (int(id1), int(id2)): index
        for index, id1, id2 in zip(df.index, df[id1_col], df[id2_col])
    }

    matched = 0
    for result in results:
        if result["root_id1"] is None or result["root_id2"] is None:
            continue

        key = (int(result["root_id1"]), int(result["root_id2"]))
        row_index = pair_to_index.get(key)
        if row_index is None:
            logger.warning(
                "NPZ pair %s -> %s is not present in the CSV.", key[0], key[1]
            )
            continue

        for column in (
            OVERLAP_NM_COL,
            OVERLAP_UM_COL,
            N_CLOSE_PAIRS_COL,
            N_CLOSE_AXON_VERTICES_COL,
            PAIR_CAPPED_COL,
            STATUS_COL,
            ERROR_COL,
            SOURCE_FILE_COL,
        ):
            df.at[row_index, column] = result[column]
        matched += 1

    return df, matched


def atomic_write_csv(df: pd.DataFrame, csv_path: Path) -> None:
    csv_path = csv_path.expanduser().resolve()
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = csv_path.with_name(csv_path.name + ".tmp")
    df.to_csv(temporary_path, index=False)
    os.replace(temporary_path, csv_path)


def write_results_checkpoint(
    results: list[dict], log_dir: Path, timestamp: str
) -> Path:
    checkpoint = log_dir / f"path_overlap_checkpoint_{timestamp}.csv"
    pd.DataFrame(results).to_csv(checkpoint, index=False)
    return checkpoint


def log_one_result(result: dict, logger: logging.Logger) -> None:
    if result[STATUS_COL] == "ok":
        logger.info(
            "root_id1=%s | root_id2=%s | overlap=%.6f um | "
            "close_pairs=%s | close_axon_vertices=%s | possibly_capped=%s | file=%s",
            result["root_id1"],
            result["root_id2"],
            result[OVERLAP_UM_COL],
            result[N_CLOSE_PAIRS_COL],
            result[N_CLOSE_AXON_VERTICES_COL],
            result[PAIR_CAPPED_COL],
            result[SOURCE_FILE_COL],
        )
    else:
        logger.error(
            "Failed file=%s | root_id1=%s | root_id2=%s | %s",
            result[SOURCE_FILE_COL],
            result["root_id1"],
            result["root_id2"],
            result[ERROR_COL],
        )


def run_single(args: argparse.Namespace, logger: logging.Logger) -> int:
    npz_path = resolve_npz_path(args.npz_file, args.npz_dir)
    logger.info("Single-file mode")
    logger.info("NPZ: %s", npz_path)

    result = calculate_one_npz(npz_path)
    log_one_result(result, logger)

    if result[STATUS_COL] != "ok":
        return 1

    print("\nResult")
    print("------")
    print(f"Sensory neuron (root_id1): {result['root_id1']}")
    print(f"ExN (root_id2):            {result['root_id2']}")
    print(
        "1 axon -> 2 dendrite overlap: "
        f"{result[OVERLAP_UM_COL]:.2f} um "
        f"({result[OVERLAP_NM_COL]:.2f} nm)"
    )
    print(f"Close vertex pairs:         {result[N_CLOSE_PAIRS_COL]}")
    print(f"Unique close axon vertices: {result[N_CLOSE_AXON_VERTICES_COL]}")
    print(f"Possibly capped at 5000:    {result[PAIR_CAPPED_COL]}")

    if not args.no_update_csv:
        df, id1_col, id2_col = read_pair_csv(args.csv)
        df, matched = apply_results_to_csv_dataframe(
            df, id1_col, id2_col, [result], logger
        )
        if matched == 0:
            logger.warning("CSV was not modified because the pair was not found.")
        else:
            atomic_write_csv(df, args.csv)
            logger.info("Updated CSV atomically: %s", args.csv.resolve())

    return 0


def run_batch(args: argparse.Namespace, logger: logging.Logger) -> int:
    npz_dir = args.npz_dir.expanduser().resolve()
    if not npz_dir.is_dir():
        logger.error("NPZ directory not found: %s", npz_dir)
        return 1

    npz_paths = sorted(npz_dir.glob("close_contacts_*.npz"))
    if not npz_paths:
        logger.error("No close_contacts_*.npz files found in %s", npz_dir)
        return 1

    workers = max(1, int(args.workers))
    logger.info("Batch mode")
    logger.info("NPZ directory: %s", npz_dir)
    logger.info("Found %d NPZ files", len(npz_paths))
    logger.info("Workers: %d", workers)

    results: list[dict] = []
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    with ProcessPoolExecutor(max_workers=workers) as executor:
        iterator = executor.map(
            calculate_one_npz,
            (str(path) for path in npz_paths),
            chunksize=20,
        )

        for completed, result in enumerate(iterator, start=1):
            results.append(result)
            log_one_result(result, logger)

            if completed % 100 == 0 or completed == len(npz_paths):
                ok_count = sum(r[STATUS_COL] == "ok" for r in results)
                logger.info(
                    "Progress: %d/%d | successful=%d | errors=%d",
                    completed,
                    len(npz_paths),
                    ok_count,
                    completed - ok_count,
                )

            if (
                args.checkpoint_every > 0
                and completed % args.checkpoint_every == 0
            ):
                checkpoint = write_results_checkpoint(
                    results, args.log_dir, timestamp
                )
                logger.info("Wrote checkpoint: %s", checkpoint)

    final_results_path = args.log_dir / f"path_overlap_results_{timestamp}.csv"
    pd.DataFrame(results).to_csv(final_results_path, index=False)
    logger.info("Wrote complete result table: %s", final_results_path)

    if not args.no_update_csv:
        df, id1_col, id2_col = read_pair_csv(args.csv)
        df, matched = apply_results_to_csv_dataframe(
            df, id1_col, id2_col, results, logger
        )
        atomic_write_csv(df, args.csv)
        logger.info(
            "Updated %d CSV rows atomically in %s", matched, args.csv.resolve()
        )

    ok_count = sum(r[STATUS_COL] == "ok" for r in results)
    capped_count = sum(
        r[STATUS_COL] == "ok" and bool(r[PAIR_CAPPED_COL]) for r in results
    )
    logger.info(
        "Finished | files=%d | successful=%d | errors=%d | possibly_capped=%d",
        len(results),
        ok_count,
        len(results) - ok_count,
        capped_count,
    )

    return 0 if ok_count == len(results) else 2


def main() -> int:
    args = parse_args()
    args.npz_dir = args.npz_dir.expanduser().resolve()
    args.csv = args.csv.expanduser().resolve()
    args.log_dir = args.log_dir.expanduser().resolve()

    mode = "single" if args.npz_file else "batch"
    logger, log_path = configure_logging(args.log_dir, mode)

    logger.info("Script: %s", Path(__file__).resolve())
    logger.info("Log: %s", log_path)
    logger.info("CSV: %s", args.csv)

    try:
        if args.npz_file:
            return run_single(args, logger)
        return run_batch(args, logger)
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        return 130
    except Exception as exc:
        logger.exception("Fatal error: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
