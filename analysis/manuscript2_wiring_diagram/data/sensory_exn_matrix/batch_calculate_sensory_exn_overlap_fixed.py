#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import time
import shutil
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

BASE_DIR = Path(__file__).resolve().parent
NPZ_DIR = BASE_DIR / "overlap_results"

INPUT_CSV = BASE_DIR / "all_possible_sensory_exn_pairs.csv"
OUTPUT_CSV = BASE_DIR / "all_possible_sensory_exn_pairs_with_overlap.csv"

LOG_DIR = BASE_DIR / "path_calculation_log"
LOG_DIR.mkdir(exist_ok=True)

PRE = "pre_id"
POST = "post_id"


def logger_setup():
    log = LOG_DIR / ("overlap_%s.log" % time.strftime("%Y%m%d_%H%M%S"))
    logger = logging.getLogger("overlap")
    logger.setLevel(logging.INFO)
    h = logging.FileHandler(log)
    s = logging.StreamHandler()
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    h.setFormatter(fmt)
    s.setFormatter(fmt)
    logger.addHandler(h)
    logger.addHandler(s)
    return logger


def overlap_length(vertices, edges, close_vertices):
    close_vertices = np.unique(close_vertices)
    if len(close_vertices) == 0:
        return 0.0

    vertices = np.asarray(vertices, float)
    edges = np.asarray(edges, int)

    mask = np.zeros(len(vertices), dtype=bool)
    mask[close_vertices] = True

    u = edges[:, 0]
    v = edges[:, 1]

    lengths = np.linalg.norm(vertices[u] - vertices[v], axis=1)

    internal = mask[u] & mask[v]
    total = lengths[internal].sum()

    degree = np.zeros(len(vertices), dtype=int)
    np.add.at(degree, u[internal], 1)
    np.add.at(degree, v[internal], 1)

    shortest = np.full(len(vertices), np.inf)

    for i in np.where(mask[u] ^ mask[v])[0]:
        cv = u[i] if mask[u[i]] else v[i]
        shortest[cv] = min(shortest[cv], lengths[i])

    endpoints = close_vertices[degree[close_vertices] == 1]
    ext = shortest[endpoints]
    ext = ext[np.isfinite(ext)]

    return float(total + 0.5 * ext.sum())


def process_file(path):
    with np.load(path, allow_pickle=False) as d:
        pre = int(d["root_id1"][0])
        post = int(d["root_id2"][0])

        pairs = d["global_pairs_1axon_2dend"]

        if len(pairs):
            close = np.unique(pairs[:, 0])
            overlap = overlap_length(
                d["skel1_vertices"],
                d["skel1_edges"],
                close
            )
        else:
            close = []
            overlap = 0.0

        return {
            PRE: pre,
            POST: post,
            "sensory_axon_exn_dendrite_overlap_nm": overlap,
            "sensory_axon_exn_dendrite_overlap_um": overlap / 1000,
            "n_close_vertex_pairs_1axon_2dend": len(pairs),
            "n_close_sensory_axon_vertices": len(close),
            "path_overlap_status": "ok",
            "path_overlap_error": "",
            "path_overlap_npz": path.name,
        }


def update_csv(result, logger):
    df = pd.read_csv(
        INPUT_CSV,
        dtype={PRE: "Int64", POST: "Int64"}
    )

    for c in ["path_overlap_status", "path_overlap_error", "path_overlap_npz"]:
        if c not in df:
            df[c] = pd.Series(dtype=object)
        df[c] = df[c].astype(object)

    for c in [
        "sensory_axon_exn_dendrite_overlap_nm",
        "sensory_axon_exn_dendrite_overlap_um",
        "n_close_vertex_pairs_1axon_2dend",
        "n_close_sensory_axon_vertices",
    ]:
        if c not in df:
            df[c] = np.nan

    idx = df.index[(df[PRE] == result[PRE]) & (df[POST] == result[POST])]

    if len(idx):
        for k, v in result.items():
            df.at[idx[0], k] = v

    tmp = OUTPUT_CSV.with_suffix(".tmp.csv")
    df.to_csv(tmp, index=False)
    shutil.move(tmp, OUTPUT_CSV)

    logger.info(
        "Updated %s: %s -> %s overlap %.3f um",
        OUTPUT_CSV.name,
        result[PRE],
        result[POST],
        result["sensory_axon_exn_dendrite_overlap_um"]
    )


def main():
    logger = logger_setup()

    if len(sys.argv) > 1:
        files = [NPZ_DIR / sys.argv[1]]
    else:
        files = sorted(NPZ_DIR.glob("close_contacts_*.npz"))

    logger.info("Processing %d NPZ files", len(files))

    for f in tqdm(files):
        try:
            update_csv(process_file(f), logger)
        except Exception as e:
            logger.exception("Failed %s: %s", f.name, e)


if __name__ == "__main__":
    main()
