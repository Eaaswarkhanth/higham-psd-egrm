#!/usr/bin/env python3
"""
higham_to_pd_and_pca.py
========================

Apply Higham nearest-positive-semidefinite (PSD) correction to a symmetric
matrix, add an optional diagonal shift to make it positive definite (PD), and
optionally compare raw vs corrected PCA scores.

Designed for ancestry-specific eGRM/kinship matrices, but the script is generic:
all input paths, labels, chromosome numbers, and output locations are supplied
by command-line arguments.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd

try:
    from scipy.linalg import eigh as scipy_eigh
    HAVE_SCIPY = True
except Exception:  # pragma: no cover
    scipy_eigh = None
    HAVE_SCIPY = False


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="[%(levelname)s] %(message)s",
        stream=sys.stderr,
    )


def set_single_threaded_blas() -> None:
    # Avoid oversubscription on HPC array jobs.
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")


def parse_chromosome(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    v = str(value).strip()
    if not v:
        return None
    if v.lower().startswith("chr"):
        v = v[3:]
    return v


def format_template(template: str, ancestry: Optional[str], chromosome: Optional[str]) -> Path:
    if ancestry is None and "{anc}" in template:
        raise ValueError("Input template contains {anc}; please provide --ancestry.")
    if chromosome is None and "{chr}" in template:
        raise ValueError("Input template contains {chr}; please provide --chromosome.")
    return Path(template.format(anc=ancestry, chr=chromosome)).expanduser()


def read_text_matrix(path: Path, delimiter: Optional[str], has_header: bool, has_index: bool) -> Tuple[np.ndarray, Optional[list[str]], Optional[list[str]]]:
    sep = delimiter
    if sep == "tab":
        sep = "\t"
    elif sep == "comma":
        sep = ","
    elif sep == "space":
        sep = r"\s+"

    header = 0 if has_header else None
    index_col = 0 if has_index else None
    df = pd.read_csv(path, sep=sep, header=header, index_col=index_col, engine="python")
    row_labels = df.index.astype(str).tolist() if has_index else None
    col_labels = df.columns.astype(str).tolist() if has_header else None
    return df.to_numpy(dtype=float), row_labels, col_labels


def load_matrix(path: Path, delimiter: Optional[str], has_header: bool, has_index: bool) -> Tuple[np.ndarray, Optional[list[str]], Optional[list[str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Input matrix not found: {path}")

    suffix = path.suffix.lower()
    if suffix == ".npy":
        return np.load(path), None, None
    if suffix in {".csv", ".tsv", ".txt"}:
        if delimiter is None:
            delimiter = "comma" if suffix == ".csv" else "tab"
        return read_text_matrix(path, delimiter, has_header, has_index)

    raise ValueError(f"Unsupported matrix format for {path}. Use .npy, .csv, .tsv, or .txt")


def save_matrix(path: Path, matrix: np.ndarray, row_labels: Optional[list[str]] = None, col_labels: Optional[list[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix == ".npy":
        np.save(path, matrix)
    elif suffix in {".csv", ".tsv", ".txt"}:
        sep = "," if suffix == ".csv" else "\t"
        df = pd.DataFrame(matrix, index=row_labels, columns=col_labels)
        df.to_csv(path, sep=sep, index=row_labels is not None)
    else:
        raise ValueError(f"Unsupported output matrix format for {path}. Use .npy, .csv, .tsv, or .txt")


def eig_sym(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if HAVE_SCIPY:
        return scipy_eigh(matrix, check_finite=False, overwrite_a=False)
    return np.linalg.eigh(matrix)


def validate_square_numeric(matrix: np.ndarray) -> None:
    if matrix.ndim != 2:
        raise ValueError(f"Matrix must be 2-dimensional; observed shape {matrix.shape}")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Matrix must be square; observed shape {matrix.shape}")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Matrix contains NaN or infinite values.")


def symmetrize(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + matrix.T)


def project_psd(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    eigvals, eigvecs = eig_sym(matrix)
    eigvals_clipped = np.maximum(eigvals, 0.0)
    psd = (eigvecs * eigvals_clipped) @ eigvecs.T
    return symmetrize(psd), eigvals_clipped


def higham_nearest_psd(matrix: np.ndarray, max_iters: int = 100, tol: float = 1e-9) -> Tuple[np.ndarray, bool, int, float]:
    """Higham nearest-PSD approximation via alternating projections."""
    A = symmetrize(matrix)
    fro_A = np.linalg.norm(A, ord="fro")
    if fro_A == 0:
        return np.zeros_like(A), True, 0, 0.0

    Y = A.copy()
    delta = np.zeros_like(A)
    rel_change = np.inf
    converged = False
    X = Y.copy()

    for iteration in range(1, max_iters + 1):
        R = Y - delta
        X, _ = project_psd(R)
        delta = X - R
        rel_change = np.linalg.norm(X - Y, ord="fro") / fro_A
        Y = X
        if rel_change < tol:
            converged = True
            break

    return symmetrize(X), converged, iteration, float(rel_change)


def topk_scores(matrix: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    eigvals, eigvecs = eig_sym(matrix)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    k = min(k, matrix.shape[0])
    eigvals_k = eigvals[:k]
    eigvecs_k = eigvecs[:, :k]
    scores = eigvecs_k * np.sqrt(np.maximum(eigvals_k, 0.0))
    return eigvals_k, eigvecs_k, scores


def per_pc_corr(raw_scores: np.ndarray, pd_scores: np.ndarray) -> np.ndarray:
    corrs = []
    for i in range(raw_scores.shape[1]):
        x = raw_scores[:, i]
        y = pd_scores[:, i]
        if np.dot(x, y) < 0:
            y = -y
        x = x - x.mean()
        y = y - y.mean()
        denom = np.linalg.norm(x) * np.linalg.norm(y)
        corrs.append(float(np.dot(x, y) / denom) if denom > 0 else np.nan)
    return np.asarray(corrs)


def principal_angle_cosines(raw_vecs: np.ndarray, pd_vecs: np.ndarray) -> np.ndarray:
    q_raw, _ = np.linalg.qr(raw_vecs)
    q_pd, _ = np.linalg.qr(pd_vecs)
    _, singular_values, _ = np.linalg.svd(q_raw.T @ q_pd, full_matrices=False)
    return singular_values


def procrustes_r(raw_scores: np.ndarray, pd_scores: np.ndarray) -> float:
    X = raw_scores - raw_scores.mean(axis=0, keepdims=True)
    Y = pd_scores - pd_scores.mean(axis=0, keepdims=True)
    _, singular_values, _ = np.linalg.svd(X.T @ Y, full_matrices=False)
    denom = np.linalg.norm(X, ord="fro") * np.linalg.norm(Y, ord="fro")
    return float(np.sum(singular_values) / denom) if denom > 0 else np.nan


def build_output_paths(out_dir: Path, ancestry: Optional[str], chromosome: Optional[str], prefix: str) -> Tuple[Path, Path, Path, Path]:
    label = ancestry or "matrix"
    chrom_part = f"chr{chromosome}" if chromosome is not None else "matrix"
    subdir = out_dir / f"tgtAnc-{label}" if ancestry is not None else out_dir
    rows_dir = out_dir / "rows"
    psd_path = subdir / f"{prefix}_{chrom_part}_higham_psd.npy"
    pd_path = subdir / f"{prefix}_{chrom_part}_higham.npy"
    qc_path = rows_dir / f"higham_qc_{label}_{chrom_part}.csv"
    pca_path = rows_dir / f"higham_pca_{label}_{chrom_part}.csv"
    return psd_path, pd_path, qc_path, pca_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Higham nearest-PSD correction, optional PD shift, and PCA QC for symmetric matrices."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input", help="Path to one input matrix (.npy, .csv, .tsv, or .txt).")
    source.add_argument(
        "--input-template",
        help="Input template with optional {anc} and {chr}, e.g. '/path/tgtAnc-{anc}/matrix_chr{chr}.npy'.",
    )
    parser.add_argument("--ancestry", help="Ancestry/label used with --input-template and output names, e.g. AFR.")
    parser.add_argument("--chromosome", help="Chromosome used with --input-template and output names, e.g. 1 or chr1.")
    parser.add_argument("--out-dir", required=True, help="Output directory.")
    parser.add_argument("--prefix", default="pd", help="Output matrix prefix. Default: pd")
    parser.add_argument("--delimiter", default=None, help="Delimiter for text matrices: tab, comma, space, or a literal delimiter.")
    parser.add_argument("--has-header", action="store_true", help="Input text matrix has a header row.")
    parser.add_argument("--has-index", action="store_true", help="Input text matrix has row labels in the first column.")
    parser.add_argument("--k", type=int, default=20, help="Number of top PCs for raw-vs-PD comparison. Default: 20")
    parser.add_argument("--target-eps", type=float, default=1e-10, help="Minimum eigenvalue target for final PD matrix.")
    parser.add_argument("--max-iters", type=int, default=100, help="Maximum Higham iterations.")
    parser.add_argument("--tol", type=float, default=1e-9, help="Higham convergence tolerance.")
    parser.add_argument("--save-scores", action="store_true", help="Save raw and PD top-k score CSV files.")
    parser.add_argument("--save-text", action="store_true", help="Also save PSD and PD matrices as TSV text files.")
    parser.add_argument("--allow-asymmetry", type=float, default=1e-8, help="Warn if max |A-A'| exceeds this before symmetrizing.")
    parser.add_argument("--verbose", action="store_true", help="Print debug logging.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)
    set_single_threaded_blas()

    chromosome = parse_chromosome(args.chromosome)
    input_path = Path(args.input).expanduser() if args.input else format_template(args.input_template, args.ancestry, chromosome)
    out_dir = Path(args.out_dir).expanduser()

    logging.info("Input matrix: %s", input_path)
    matrix, row_labels, col_labels = load_matrix(input_path, args.delimiter, args.has_header, args.has_index)
    validate_square_numeric(matrix)

    asym = float(np.max(np.abs(matrix - matrix.T)))
    if asym > args.allow_asymmetry:
        logging.warning("Input is not perfectly symmetric: max |A-A'| = %.3e. Symmetrizing before correction.", asym)
    raw = symmetrize(matrix.astype(float, copy=False))
    n = raw.shape[0]

    raw_eigvals, _ = eig_sym(raw)
    min_raw = float(np.min(raw_eigvals))

    logging.info("Running Higham nearest-PSD correction: n=%s, min raw eigenvalue=%.3e", n, min_raw)
    psd, converged, iterations, rel_change = higham_nearest_psd(raw, max_iters=args.max_iters, tol=args.tol)

    psd_eigvals, _ = eig_sym(psd)
    min_psd = float(np.min(psd_eigvals))

    shift = 0.0
    if min_psd < args.target_eps:
        shift = float(args.target_eps - min_psd)
    pd_matrix = symmetrize(psd + shift * np.eye(n))

    pd_eigvals, _ = eig_sym(pd_matrix)
    min_pd = float(np.min(pd_eigvals))

    psd_path, pd_path, qc_path, pca_path = build_output_paths(out_dir, args.ancestry, chromosome, args.prefix)
    for p in [psd_path.parent, pd_path.parent, qc_path.parent]:
        p.mkdir(parents=True, exist_ok=True)

    save_matrix(psd_path, psd)
    save_matrix(pd_path, pd_matrix)
    if args.save_text:
        save_matrix(psd_path.with_suffix(".tsv"), psd, row_labels=row_labels, col_labels=col_labels)
        save_matrix(pd_path.with_suffix(".tsv"), pd_matrix, row_labels=row_labels, col_labels=col_labels)

    qc = pd.DataFrame([{
        "Label": args.ancestry or "matrix",
        "Chromosome": chromosome or "NA",
        "N": n,
        "Input": str(input_path),
        "Output_PSD": str(psd_path),
        "Output_PD": str(pd_path),
        "Max_Asymmetry_Input": asym,
        "MinEig_Raw": min_raw,
        "MinEig_HighamPSD": min_psd,
        "MinEig_FinalPD": min_pd,
        "Shift_Added": shift,
        "Target_Eps": args.target_eps,
        "Iterations": iterations,
        "Converged": bool(converged),
        "RelChange": rel_change,
        "Status": "OK",
    }])
    qc.to_csv(qc_path, index=False)

    k = min(int(args.k), n)
    raw_eval, raw_vec, raw_scores = topk_scores(raw, k)
    pd_eval, pd_vec, pd_scores = topk_scores(pd_matrix, k)
    corrs = per_pc_corr(raw_scores, pd_scores)
    cosines = principal_angle_cosines(raw_vec, pd_vec)
    proc = procrustes_r(raw_scores, pd_scores)

    metrics = {
        "Label": args.ancestry or "matrix",
        "Chromosome": chromosome or "NA",
        "k": k,
        "Procrustes_R": proc,
        "MinCorr_PC1toK": float(np.nanmin(corrs)) if len(corrs) else np.nan,
        "MaxCorr_PC1toK": float(np.nanmax(corrs)) if len(corrs) else np.nan,
        "Iterations": iterations,
        "Converged": bool(converged),
        "Shift_Added": shift,
        "Target_Eps": args.target_eps,
    }
    for i, value in enumerate(corrs, start=1):
        metrics[f"PC{i}_corr"] = float(value)
    for i, value in enumerate(cosines, start=1):
        metrics[f"Subspace_cos{i}"] = float(value)
    for i, value in enumerate(raw_eval, start=1):
        metrics[f"Raw_Eig_PC{i}"] = float(value)
    for i, value in enumerate(pd_eval, start=1):
        metrics[f"PD_Eig_PC{i}"] = float(value)
    pd.DataFrame([metrics]).to_csv(pca_path, index=False)

    if args.save_scores:
        cols = [f"PC{i}" for i in range(1, k + 1)]
        score_dir = pd_path.parent
        pd.DataFrame(raw_scores, columns=cols).to_csv(score_dir / f"scores_RAW_chr{chromosome or 'matrix'}_top{k}.csv", index=False)
        pd.DataFrame(pd_scores, columns=cols).to_csv(score_dir / f"scores_PD_chr{chromosome or 'matrix'}_top{k}.csv", index=False)

    logging.info(
        "Done: label=%s chr=%s minPSD=%.3e shift=%.3e minPD=%.3e converged=%s iterations=%s Procrustes_R=%.6f",
        args.ancestry or "matrix", chromosome or "NA", min_psd, shift, min_pd, converged, iterations, proc,
    )


if __name__ == "__main__":
    main()
