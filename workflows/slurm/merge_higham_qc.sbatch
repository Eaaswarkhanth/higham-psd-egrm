#!/usr/bin/env bash
#SBATCH --job-name=higham_merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:45:00
#SBATCH --output=logs/higham_merge_%j.out
#SBATCH --error=logs/higham_merge_%j.err

set -euo pipefail

OUT_DIR=${OUT_DIR:-"results/higham"}
PYTHON=${PYTHON:-python}
SCRIPT=${MERGE_SCRIPT:-scripts/merge_higham_qc.py}

mkdir -p logs "${OUT_DIR}"

${PYTHON} ${SCRIPT} \
  --rows-dir "${OUT_DIR}/rows" \
  --out-dir "${OUT_DIR}"
