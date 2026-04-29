#!/usr/bin/env bash
set -euo pipefail

# Submit array job and then merge QC/PCA metrics after successful completion.
# Override variables using environment variables, for example:
#   WORK_DIR=/path/to/root OUT_DIR=results/higham ANCESTRIES_STR="AFR EUR EAS AMR" bash workflows/slurm/run_higham_pipeline.sh

ARRAY_SCRIPT=${ARRAY_SCRIPT:-workflows/slurm/higham_array.sbatch}
MERGE_SCRIPT=${MERGE_SCRIPT:-workflows/slurm/merge_higham_qc.sbatch}

ARRAY_JOBID=$(sbatch "${ARRAY_SCRIPT}" | awk '{print $4}')
echo "Submitted Higham array job: ${ARRAY_JOBID}"

sbatch --dependency=afterok:${ARRAY_JOBID} "${MERGE_SCRIPT}"
echo "Submitted merge job with dependency afterok:${ARRAY_JOBID}"
