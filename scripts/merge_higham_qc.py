#!/usr/bin/env python3
"""Merge row-level Higham QC and PCA metric CSV files."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge Higham QC/PCA row CSVs into summary tables.")
    parser.add_argument("--rows-dir", required=True, help="Directory containing higham_qc_*.csv and higham_pca_*.csv files.")
    parser.add_argument("--out-dir", required=True, help="Output directory for merged summaries.")
    parser.add_argument("--qc-name", default="higham_qc_summary.csv", help="Merged QC filename.")
    parser.add_argument("--pca-name", default="higham_pca_metrics.csv", help="Merged PCA metrics filename.")
    return parser.parse_args()


def concat_csvs(files: list[Path]) -> pd.DataFrame | None:
    if not files:
        return None
    frames = []
    for path in files:
        df = pd.read_csv(path)
        df["Source_File"] = path.name
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    sort_cols = [c for c in ["Label", "Ancestry", "Chromosome", "Chr"] if c in out.columns]
    if sort_cols:
        out = out.sort_values(sort_cols)
    return out


def main() -> None:
    args = parse_args()
    rows_dir = Path(args.rows_dir).expanduser()
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    qc_files = sorted(rows_dir.glob("higham_qc_*.csv"))
    pca_files = sorted(rows_dir.glob("higham_pca_*.csv"))

    qc = concat_csvs(qc_files)
    if qc is not None:
        out = out_dir / args.qc_name
        qc.to_csv(out, index=False)
        print(f"[OK] wrote {out} rows={len(qc)}")
    else:
        print(f"[WARN] no QC files found in {rows_dir}")

    pca = concat_csvs(pca_files)
    if pca is not None:
        out = out_dir / args.pca_name
        pca.to_csv(out, index=False)
        print(f"[OK] wrote {out} rows={len(pca)}")
    else:
        print(f"[WARN] no PCA metric files found in {rows_dir}")


if __name__ == "__main__":
    main()
