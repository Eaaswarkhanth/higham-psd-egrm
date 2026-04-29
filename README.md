# Higham PSD eGRM Utility

Reusable command-line utilities for correcting symmetric genomic relationship / kinship / eGRM matrices using Higham nearest-positive-semidefinite correction, followed by an optional diagonal shift to obtain a strictly positive-definite matrix.

This repository was designed for ancestry-specific eGRM matrices, but the main script is generic and can be used with any square numeric matrix.

## What the workflow does

For each matrix, the workflow:

1. Loads a square matrix from `.npy`, `.csv`, `.tsv`, or `.txt`.
2. Symmetrizes the matrix.
3. Applies Higham nearest-PSD correction.
4. Adds the minimal diagonal shift needed to make the matrix positive definite.
5. Saves both the PSD and final PD matrices.
6. Computes raw-vs-PD PCA comparison metrics.
7. Optionally saves raw and PD PC scores.
8. Optionally converts `.npy` matrices to GENESIS-compatible `.RDS` files.

## Repository layout

```text
higham-psd-egrm/
├── scripts/
│   ├── higham_to_pd_and_pca.py
│   ├── merge_higham_qc.py
│   └── npy_to_genesis_rds.py
├── workflows/
│   └── slurm/
│       ├── higham_array.sbatch
│       ├── merge_higham_qc.sbatch
│       └── run_higham_pipeline.sh
├── configs/
│   └── higham_example.yaml
├── examples/
│   ├── run_higham_example.sh
│   └── sample_ids.example.txt
├── requirements.txt
├── environment.yml
├── LICENSE
└── README.md
```

## Installation

Using conda:

```bash
conda env create -f environment.yml
conda activate higham-psd-egrm
```

Or using pip:

```bash
python -m pip install -r requirements.txt
```

`rpy2` and R are only required if you use `scripts/npy_to_genesis_rds.py`.

## Run one matrix

```bash
python scripts/higham_to_pd_and_pca.py \
  --input matrix_chr1.npy \
  --out-dir results/higham \
  --ancestry AFR \
  --chromosome 1 \
  --k 20 \
  --target-eps 1e-10 \
  --max-iters 100 \
  --tol 1e-9 \
  --save-scores
```

## Run using an input template

For ancestry-specific eGRM matrices arranged by ancestry and chromosome:

```bash
python scripts/higham_to_pd_and_pca.py \
  --input-template "/path/to/input/tgtAnc-{anc}/merged_leave-one-out/merged_diploid_no-chr{chr}.npy" \
  --out-dir results/higham \
  --ancestry AFR \
  --chromosome 1 \
  --k 20 \
  --target-eps 1e-10
```

The default public ancestry labels used in the SLURM workflow are:

```text
AFR EUR EAS AMR
```

These can be changed freely.

## Outputs

For label `AFR` and chromosome `1`, the main script writes:

```text
results/higham/tgtAnc-AFR/pd_chr1_higham_psd.npy
results/higham/tgtAnc-AFR/pd_chr1_higham.npy
results/higham/rows/higham_qc_AFR_chr1.csv
results/higham/rows/higham_pca_AFR_chr1.csv
```

If `--save-scores` is used, it also writes:

```text
results/higham/tgtAnc-AFR/scores_RAW_chr1_top20.csv
results/higham/tgtAnc-AFR/scores_PD_chr1_top20.csv
```

## Merge QC summaries

After running many ancestry/chromosome jobs:

```bash
python scripts/merge_higham_qc.py \
  --rows-dir results/higham/rows \
  --out-dir results/higham
```

This creates:

```text
results/higham/higham_qc_summary.csv
results/higham/higham_pca_metrics.csv
```

## SLURM workflow

The SLURM array is set up for 4 ancestries × 22 chromosomes = 88 jobs.

```bash
#SBATCH --array=0-87
```

The default indexing logic is:

```bash
ANCESTRIES=(AFR EUR EAS AMR)
CHROMOSOMES=(1 2 3 ... 22)
ANC_IDX=$(( SLURM_ARRAY_TASK_ID / 22 ))
CHR_IDX=$(( SLURM_ARRAY_TASK_ID % 22 ))
```

Submit the full workflow:

```bash
WORK_DIR=/path/to/input/root \
OUT_DIR=results/higham \
ANCESTRIES_STR="AFR EUR EAS AMR" \
bash workflows/slurm/run_higham_pipeline.sh
```

You can override the input template if your files use a different structure:

```bash
WORK_DIR=/path/to/input/root \
INPUT_TEMPLATE="/path/to/input/tgtAnc-{anc}/merged_leave-one-out/merged_diploid_no-chr{chr}.npy" \
OUT_DIR=results/higham \
bash workflows/slurm/run_higham_pipeline.sh
```

## Convert corrected matrices to GENESIS RDS

```bash
python scripts/npy_to_genesis_rds.py \
  --matrix-template "results/higham/tgtAnc-{anc}/pd_chr{chr}_higham.npy" \
  --sample-ids sample_ids.txt \
  --ancestries AFR EUR EAS AMR \
  --chromosomes 1-22 \
  --out-dir results/genesis_rds
```

If input and output ancestry labels differ, use `input:output` pairs:

```bash
python scripts/npy_to_genesis_rds.py \
  --matrix-template "results/higham/tgtAnc-{anc}/pd_chr{chr}_higham.npy" \
  --sample-ids sample_ids.txt \
  --ancestries AFR:AFR EUR:EUR EAS:EAS AMR:AMR \
  --chromosomes 1-22 \
  --out-dir results/genesis_rds
```

## Notes

- Do not upload private genotype, kinship, eGRM, or sample-level data to GitHub.
- Keep real data outside the repository and pass paths by command-line arguments.
- The `.gitignore` intentionally excludes large matrix files and result outputs.
- For HPC use, edit SLURM account/partition settings locally or pass site-specific options at submission time.

## Optional: create tiny synthetic `.npy` inputs for testing

Do not commit real eGRM/kinship/genomic matrices to GitHub. For testing the workflow, generate toy matrices:

```bash
python examples/make_synthetic_egrm.py \
  --out-root examples/synthetic_input \
  --ancestries AFR EUR EAS AMR \
  --chromosomes 1 2 \
  --n 8
```

Then run one test chromosome/ancestry:

```bash
python scripts/higham_to_pd_and_pca.py \
  --input-template "examples/synthetic_input/tgtAnc-{anc}/merged_leave-one-out/merged_diploid_no-chr{chr}.npy" \
  --out-dir results/higham_test \
  --ancestry AFR \
  --chromosome 1 \
  --k 4 \
  --target-eps 1e-10 \
  --save-scores
```
