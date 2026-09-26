#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fast batch calculation of sensory axon -> ExN dendrite path overlap.

Designed for:
    ~/sensory_exn_path/

Files:
    all_possible_sensory_exn_pairs.csv
    overlap_results/close_contacts_*.npz

Usage:
    python batch_calculate_sensory_exn_overlap_fast.py

Optional:
    python batch_calculate_sensory_exn_overlap_fast.py --resume

Python 3.7 compatible.
"""

import argparse
import logging
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


BASE_DIR = Path(__file__).resolve().parent

NPZ_DIR = BASE_DIR / "overlap_results"

INPUT_CSV = BASE_DIR / "all_possible_sensory_exn_pairs.csv"
OUTPUT_CSV = BASE_DIR / "all_possible_sensory_exn_pairs_with_overlap.csv"

CHECKPOINT_CSV = BASE_DIR / "all_possible_sensory_exn_pairs_checkpoint.csv"

LOG_DIR = BASE_DIR / "path_calculation_log"
LOG_DIR.mkdir(exist_ok=True)

CHECKPOINT_EVERY = 500


def setup_logger():
    logfile = LOG_DIR / (
        "fast_overlap_%s.log" %
        time.strftime("%Y%m%d_%H%M%S")
    )

    logger = logging.getLogger("overlap")
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(str(logfile))
    sh = logging.StreamHandler()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s"
    )

    fh.setFormatter(fmt)
    sh.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(sh)

    return logger


def calculate_overlap_nm(vertices, edges, close_vertices):
    """
    Calculate sensory axon path length within close-contact regions.
    """

    close_vertices = np.unique(close_vertices)

    if len(close_vertices) == 0:
        return 0.0

    vertices = np.asarray(vertices, dtype=float)
    edges = np.asarray(edges, dtype=np.int64)

    mask = np.zeros(len(vertices), dtype=bool)
    mask[close_vertices] = True

    u = edges[:, 0]
    v = edges[:, 1]

    lengths = np.linalg.norm(
        vertices[u] - vertices[v],
        axis=1,
    )

    internal = mask[u] & mask[v]

    total = float(lengths[internal].sum())

    degree = np.zeros(len(vertices), dtype=np.int64)

    np.add.at(degree, u[internal], 1)
    np.add.at(degree, v[internal], 1)

    shortest_external = np.full(
        len(vertices),
        np.inf,
        dtype=float,
    )

    crossing = mask[u] ^ mask[v]

    for idx in np.where(crossing)[0]:

        close_v = u[idx] if mask[u[idx]] else v[idx]

        if lengths[idx] < shortest_external[close_v]:
            shortest_external[close_v] = lengths[idx]

    endpoints = close_vertices[
        degree[close_vertices] == 1
    ]

    ext = shortest_external[endpoints]
    ext = ext[np.isfinite(ext)]

    return total + 0.5 * float(ext.sum())


def process_npz(path):

    with np.load(path, allow_pickle=False) as data:

        pre_id = int(data["root_id1"][0])
        post_id = int(data["root_id2"][0])

        pairs = data["global_pairs_1axon_2dend"]

        if len(pairs) == 0:
            overlap_nm = 0.0
            n_vertices = 0

        else:
            close_vertices = np.unique(
                pairs[:, 0]
            )

            overlap_nm = calculate_overlap_nm(
                data["skel1_vertices"],
                data["skel1_edges"],
                close_vertices,
            )

            n_vertices = len(close_vertices)

        return {
            "pre_id": pre_id,
            "post_id": post_id,
            "sensory_axon_exn_dendrite_overlap_nm": overlap_nm,
            "sensory_axon_exn_dendrite_overlap_um": overlap_nm / 1000.0,
            "n_close_vertex_pairs_1axon_2dend": len(pairs),
            "n_close_sensory_axon_vertices": n_vertices,
            "path_overlap_status": "ok",
            "path_overlap_error": "",
            "path_overlap_npz": path.name,
        }


def load_table():

    if CHECKPOINT_CSV.exists():
        print(
            "Using checkpoint:",
            CHECKPOINT_CSV
        )

        df = pd.read_csv(
            CHECKPOINT_CSV,
            dtype={
                "pre_id": "Int64",
                "post_id": "Int64",
            },
        )

    else:
        df = pd.read_csv(
            INPUT_CSV,
            dtype={
                "pre_id": "Int64",
                "post_id": "Int64",
            },
        )

    string_cols = [
        "path_overlap_status",
        "path_overlap_error",
        "path_overlap_npz",
    ]

    numeric_cols = [
        "sensory_axon_exn_dendrite_overlap_nm",
        "sensory_axon_exn_dendrite_overlap_um",
        "n_close_vertex_pairs_1axon_2dend",
        "n_close_sensory_axon_vertices",
    ]

    for c in string_cols:
        if c not in df.columns:
            df[c] = ""
        df[c] = df[c].astype(object)

    for c in numeric_cols:
        if c not in df.columns:
            df[c] = np.nan

    return df


def save_checkpoint(df):

    tmp = CHECKPOINT_CSV.with_suffix(".tmp.csv")

    df.to_csv(
        tmp,
        index=False,
    )

    shutil.move(
        tmp,
        CHECKPOINT_CSV,
    )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--resume",
        action="store_true",
        help="continue from checkpoint",
    )

    args = parser.parse_args()

    logger = setup_logger()

    df = load_table()

    done = set()

    if args.resume:

        completed = df[
            df["path_overlap_status"] == "ok"
        ]

        done = set(
            zip(
                completed.pre_id.astype(int),
                completed.post_id.astype(int),
            )
        )

        logger.info(
            "Resume mode: %d completed pairs",
            len(done),
        )

    index_map = {
        (int(row.pre_id), int(row.post_id)): i
        for i, row in df.iterrows()
    }

    files = sorted(
        NPZ_DIR.glob(
            "close_contacts_*.npz"
        )
    )

    logger.info(
        "Found %d NPZ files",
        len(files),
    )

    counter = 0

    for path in tqdm(files):

        try:

            result = process_npz(path)

            key = (
                result["pre_id"],
                result["post_id"],
            )

            if key in done:
                continue

            if key in index_map:

                idx = index_map[key]

                for k, v in result.items():
                    df.at[idx, k] = v

            counter += 1

            if counter % CHECKPOINT_EVERY == 0:

                logger.info(
                    "Checkpoint after %d files",
                    counter,
                )

                save_checkpoint(df)

        except Exception as e:

            logger.exception(
                "Failed %s: %s",
                path.name,
                e,
            )

    save_checkpoint(df)

    shutil.copy(
        CHECKPOINT_CSV,
        OUTPUT_CSV,
    )

    logger.info(
        "Finished. Output: %s",
        OUTPUT_CSV,
    )


if __name__ == "__main__":
    main()
