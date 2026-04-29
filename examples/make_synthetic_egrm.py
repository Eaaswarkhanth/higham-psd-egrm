#!/usr/bin/env python3
"""
This script generates synthetic eGRM-like .npy toy matrices for testing the Higham PSD/PD workflow.
"""

import argparse
from pathlib import Path
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Generate toy symmetric .npy matrices for workflow testing.")
    p.add_argument("--out-root", default="examples/synthetic_input", help="Output root directory")
    p.add_argument("--ancestries", nargs="+", default=["AFR", "EUR", "EAS", "AMR"], help="Ancestry labels")
    p.add_argument("--chromosomes", nargs="+", default=["1", "2"], help="Chromosomes to generate")
    p.add_argument("--n", type=int, default=8, help="Matrix dimension / number of samples")
    p.add_argument("--seed", type=int, default=123, help="Random seed")
    return p.parse_args()


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    sample_ids = [f"sample_{i:03d}" for i in range(1, args.n + 1)]
    (out_root / "sample_ids.synthetic.txt").write_text("\n".join(sample_ids) + "\n")

    for anc in args.ancestries:
        for chrom in args.chromosomes:
            out_dir = out_root / f"tgtAnc-{anc}" / "merged_leave-one-out"
            out_dir.mkdir(parents=True, exist_ok=True)

            X = rng.normal(size=(args.n, args.n))
            A = 0.5 * (X + X.T)
            # Force a small negative eigenvalue so correction has something to do.
            w, V = np.linalg.eigh(A)
            w[0] = -abs(w[0]) - 0.05
            A = (V * w) @ V.T
            A = 0.5 * (A + A.T)

            out_file = out_dir / f"merged_diploid_no-chr{chrom}.npy"
            np.save(out_file, A)
            print(f"[OK] wrote {out_file}")

    print(f"[OK] wrote sample IDs: {out_root / 'sample_ids.synthetic.txt'}")


if __name__ == "__main__":
    main()
