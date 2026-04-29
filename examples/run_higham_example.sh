#!/usr/bin/env bash
set -euo pipefail

python scripts/higham_to_pd_and_pca.py \
  --input-template "/path/to/input/tgtAnc-{anc}/merged_leave-one-out/merged_diploid_no-chr{chr}.npy" \
  --out-dir results/higham \
  --ancestry AFR \
  --chromosome 1 \
  --k 20 \
  --target-eps 1e-10 \
  --max-iters 100 \
  --tol 1e-9 \
  --save-scores

python scripts/merge_higham_qc.py \
  --rows-dir results/higham/rows \
  --out-dir results/higham
