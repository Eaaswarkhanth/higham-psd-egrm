#!/usr/bin/env python3
"""
Convert corrected NumPy kinship/eGRM matrices to RDS files for GENESIS.

Requires rpy2 and an R installation available in the active environment.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from rpy2 import robjects as ro
    from rpy2.robjects import numpy2ri, default_converter
    from rpy2.robjects.conversion import localconverter
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "ERROR: rpy2 is required for RDS conversion. Install rpy2 and make sure R is available.\n"
        f"Original import error: {exc}"
    )


def parse_range_or_list(values: list[str]) -> list[str]:
    out: list[str] = []
    for value in values:
        for token in str(value).replace(",", " ").split():
            if "-" in token:
                a, b = token.split("-", 1)
                out.extend([str(i) for i in range(int(a), int(b) + 1)])
            else:
                out.append(token[3:] if token.lower().startswith("chr") else token)
    return out


def parse_ancestry_map(values: list[str]) -> dict[str, str]:
    mapping = {}
    for value in values:
        for token in str(value).replace(",", " ").split():
            if ":" in token:
                src, dest = token.split(":", 1)
            else:
                src, dest = token, token
            mapping[src] = dest
    return mapping


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert .npy PD matrices to GENESIS-compatible .RDS matrices.")
    parser.add_argument(
        "--matrix-template",
        required=True,
        help="Template for input .npy matrices, e.g. 'results/tgtAnc-{anc}/pd_chr{chr}_higham.npy'.",
    )
    parser.add_argument("--sample-ids", required=True, help="Text file with one sample ID per line in matrix order.")
    parser.add_argument("--ancestries", nargs="+", default=["AFR", "EUR", "EAS", "AMR"], help="Ancestries or src:output pairs.")
    parser.add_argument("--chromosomes", nargs="+", default=["1-22"], help="Chromosomes, e.g. 1-22 or 1 2 3.")
    parser.add_argument("--out-dir", required=True, help="Output directory for RDS files.")
    parser.add_argument("--out-template", default="km_{label}_chr{chr}.RDS", help="Output filename template.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ids = pd.read_csv(args.sample_ids, header=None)[0].astype(str).tolist()
    n = len(ids)
    ancestry_map = parse_ancestry_map(args.ancestries)
    chromosomes = parse_range_or_list(args.chromosomes)
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    for chromosome in chromosomes:
        for anc_in, label_out in ancestry_map.items():
            matrix_path = Path(args.matrix_template.format(anc=anc_in, label=label_out, chr=chromosome)).expanduser()
            if not matrix_path.exists():
                raise FileNotFoundError(matrix_path)
            K = np.load(matrix_path)
            if K.shape != (n, n):
                raise ValueError(f"{matrix_path} has shape {K.shape}, expected {(n, n)} from sample IDs.")

            with localconverter(default_converter + numpy2ri.converter):
                rK = ro.conversion.py2rpy(K)
            ro.globalenv["km"] = rK
            ro.globalenv["ids"] = ro.StrVector(ids)
            ro.r("rownames(km) <- ids; colnames(km) <- ids")

            out_path = out_dir / args.out_template.format(anc=anc_in, label=label_out, chr=chromosome)
            ro.r(f"saveRDS(km, file='{out_path.as_posix()}')")
            print(f"[OK] wrote {out_path}")


if __name__ == "__main__":
    main()
