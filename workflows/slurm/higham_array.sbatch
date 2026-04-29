#!/usr/bin/env bash
#SBATCH --job-name=higham_psd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=0-87
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --output=logs/higham_psd_%A_%a.out
#SBATCH --error=logs/higham_psd_%A_%a.err

set -euo pipefail

# User-configurable settings. You can override these at submission time, e.g.:
# sbatch --export=ALL,WORK_DIR=/path,OUT_DIR=results/higham workflows/slurm/higham_array.sbatch
WORK_DIR=${WORK_DIR:-"/path/to/input/root"}
OUT_DIR=${OUT_DIR:-"results/higham"}
INPUT_TEMPLATE=${INPUT_TEMPLATE:-"${WORK_DIR}/tgtAnc-{anc}/merged_leave-one-out/merged_diploid_no-chr{chr}.npy"}
ANCESTRIES_STR=${ANCESTRIES_STR:-"AFR EUR EAS AMR"}
CHROMOSOMES_STR=${CHROMOSOMES_STR:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"}
TOPK=${TOPK:-20}
TARGET_EPS=${TARGET_EPS:-"1e-10"}
MAX_ITERS=${MAX_ITERS:-100}
TOL=${TOL:-"1e-9"}
SAVE_SCORES=${SAVE_SCORES:-0}
PYTHON=${PYTHON:-python}
SCRIPT=${SCRIPT:-scripts/higham_to_pd_and_pca.py}

mkdir -p logs "${OUT_DIR}/rows"

read -r -a ANCESTRIES <<< "${ANCESTRIES_STR}"
read -r -a CHROMOSOMES <<< "${CHROMOSOMES_STR}"
N_CHR=${#CHROMOSOMES[@]}
N_ANC=${#ANCESTRIES[@]}
N_TASKS=$(( N_CHR * N_ANC ))

if (( SLURM_ARRAY_TASK_ID >= N_TASKS )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} exceeds N_TASKS=${N_TASKS}. Check --array range."
  exit 1
fi

ANC_IDX=$(( SLURM_ARRAY_TASK_ID / N_CHR ))
CHR_IDX=$(( SLURM_ARRAY_TASK_ID % N_CHR ))
ANC=${ANCESTRIES[$ANC_IDX]}
CHR=${CHROMOSOMES[$CHR_IDX]}

EXTRA_ARGS=()
if [[ "${SAVE_SCORES}" == "1" ]]; then
  EXTRA_ARGS+=(--save-scores)
fi

echo "[INFO] task=${SLURM_ARRAY_TASK_ID}/${N_TASKS} ancestry=${ANC} chromosome=${CHR}"
echo "[INFO] input template=${INPUT_TEMPLATE}"

${PYTHON} ${SCRIPT} \
  --input-template "${INPUT_TEMPLATE}" \
  --out-dir "${OUT_DIR}" \
  --ancestry "${ANC}" \
  --chromosome "${CHR}" \
  --k "${TOPK}" \
  --target-eps "${TARGET_EPS}" \
  --max-iters "${MAX_ITERS}" \
  --tol "${TOL}" \
  "${EXTRA_ARGS[@]}"
