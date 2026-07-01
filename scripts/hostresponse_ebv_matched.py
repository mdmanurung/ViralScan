#!/usr/bin/env python3
"""Prepare and run the EBV matched-cell host-response manuscript analysis."""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from viralscan.anellovirus import merged_name_map
from viralscan.scripts.hostresponse import DEFAULT_SEEDS, _safe_name, run_hostresponse
from viralscan.virus_grouping import group_genes_by_virus

EBV_NAME = "Epstein-Barr virus"
DEFAULT_RUN_DIR = Path("/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/out_full_depth_wl/lcl_5lines/SRR12682296")
DEFAULT_PAPER_BARCODES = Path(
    "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data/geo_GSE158275/"
    "GSM4796271_LCL_777_B958_UMI_barcodes.tsv.gz"
)
DEFAULT_OUTPUT_DIR = Path("results/hostresponse_ebv_matched")


@dataclass(frozen=True)
class PreparedInputs:
    host_h5ad: Path
    ebv_h5ad: Path
    analysis_txt: Path
    matched_barcodes: Path
    n_cells: int
    ebv_counts: np.ndarray

    def n_positive_at_threshold(self, threshold: int) -> int:
        return int((self.ebv_counts >= threshold).sum())


def normalize_barcode(barcode: str) -> str:
    """Normalize 10x barcode suffixes so paper and tool outputs can be intersected."""
    barcode = str(barcode).strip()
    if barcode.endswith("-1"):
        return barcode[:-2]
    return barcode


def _read_paper_barcodes(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    barcodes: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip().split("\t")[0]
            if value:
                barcodes.append(normalize_barcode(value))
    return barcodes


def _as_dense_vector(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel()
    arr = np.asarray(matrix)
    if arr.ndim == 1:
        return arr
    return arr.sum(axis=1)


def _normalize_adata_barcodes(adata_obj: ad.AnnData) -> ad.AnnData:
    normalized = [normalize_barcode(x) for x in adata_obj.obs_names]
    out = adata_obj.copy()
    out.obs_names = normalized
    if len(set(normalized)) != len(normalized):
        out = out[~out.obs_names.duplicated()].copy()
    return out


def _h5ad_path(run_dir: Path, filename: str) -> Path:
    direct = run_dir / filename
    # ANALYSIS_OK[cache]: resolves an INPUT file's location between two candidate
    # layouts (flat vs kb-python/counts_unfiltered); not an output cache. Missing
    # input raises FileNotFoundError below, so there is no stale-output risk.
    if direct.exists():
        return direct
    nested = run_dir / "kb-python" / "counts_unfiltered" / filename
    # ANALYSIS_OK[cache]: same input-location resolution as above; not a cache.
    if nested.exists():
        return nested
    raise FileNotFoundError(direct)


def _viral_genes(gene_ids: list[str]) -> set[str]:
    grouped, detected = group_genes_by_virus(gene_ids, merged_name_map())
    viral: set[str] = set()
    for virus in detected:
        viral.update(grouped.get(virus, []))
    return viral


def prepare_matched_inputs(
    run_dir: Path,
    paper_barcodes: Path,
    output_dir: Path,
) -> PreparedInputs:
    """Create leakage-free host and EBV-burden h5ad files for the matched paper cells."""
    run_dir = Path(run_dir)
    paper_barcodes = Path(paper_barcodes)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    host_full = _normalize_adata_barcodes(ad.read_h5ad(_h5ad_path(run_dir, "adata.h5ad")))
    virus_full = _normalize_adata_barcodes(ad.read_h5ad(_h5ad_path(run_dir, "adata_multimap.h5ad")))

    paper_order = _read_paper_barcodes(paper_barcodes)
    shared = sorted(
        {barcode for barcode in paper_order if barcode in host_full.obs_names and barcode in virus_full.obs_names}
    )
    if not shared:
        raise ValueError("No shared barcodes between paper list, adata.h5ad, and adata_multimap.h5ad")

    gene_ids = virus_full.var_names.tolist()
    grouped, _ = group_genes_by_virus(gene_ids, merged_name_map())
    ebv_genes = grouped.get(EBV_NAME, [])
    if not ebv_genes:
        raise ValueError("No EBV genes found in adata_multimap.h5ad")

    host_viral_genes = _viral_genes(host_full.var_names.tolist())
    host_genes = [gene for gene in host_full.var_names if gene not in host_viral_genes]
    if not host_genes:
        raise ValueError("Host predictor matrix has no non-viral genes after filtering")

    host = host_full[shared, host_genes].copy()
    ebv_counts = _as_dense_vector(virus_full[shared, ebv_genes].X).astype(np.float32)
    ebv = ad.AnnData(
        X=ebv_counts.reshape(-1, 1),
        obs=virus_full[shared].obs.copy(),
        var=pd.DataFrame(index=[EBV_NAME]),
    )

    matched_path = output_dir / "matched_barcodes.tsv"
    host_path = output_dir / "host_only_matched.h5ad"
    ebv_path = output_dir / "ebv_burden_matched.h5ad"
    analysis_path = output_dir / "analysis.txt"

    matched_path.write_text("\n".join(shared) + "\n", encoding="utf-8")
    host.write_h5ad(host_path)
    ebv.write_h5ad(ebv_path)
    analysis_path.write_text(f"{EBV_NAME}\n", encoding="utf-8")

    return PreparedInputs(
        host_h5ad=host_path,
        ebv_h5ad=ebv_path,
        analysis_txt=analysis_path,
        matched_barcodes=matched_path,
        n_cells=len(shared),
        ebv_counts=ebv_counts,
    )


def write_hostresponse_summary(
    output_dir: Path,
    n_cells: int,
    detection_threshold: int,
    stab_min_prob: float = 0.6,
    top_n: int = 10,
) -> Path:
    """Write a compact manuscript-oriented summary from hostresponse outputs."""
    output_dir = Path(output_dir)
    metrics_path = output_dir / "hostresponse_metrics.csv"
    stability_path = output_dir / f"{_safe_name(EBV_NAME)}_stability.csv"
    if not metrics_path.exists():
        raise FileNotFoundError(metrics_path)
    if not stability_path.exists():
        raise FileNotFoundError(stability_path)

    metrics = pd.read_csv(metrics_path)
    if metrics.empty:
        raise ValueError(f"{metrics_path} contains no hostresponse metrics")
    row = metrics.iloc[0].to_dict()
    stability = pd.read_csv(stability_path)
    stability = stability.assign(abs_weight=stability["weight_mean"].abs())
    stability["stable"] = stability["stab_prob"] >= stab_min_prob
    stable = stability[stability["stab_prob"] >= stab_min_prob].sort_values(
        ["stab_prob", "abs_weight"], ascending=[False, False]
    )
    ranked = stability.sort_values(
        ["stable", "abs_weight", "stab_prob"], ascending=[False, False, False]
    )

    n_positive = int(row.get("n_positive", 0))
    summary = dict(row)
    summary.update(
        {
            "n_cells": int(n_cells),
            "n_negative": int(n_cells) - n_positive,
            "detection_threshold_umi": int(detection_threshold),
            "stab_min_prob": float(stab_min_prob),
            "n_stable_genes": int(len(stable)),
            "top_stable_genes": ";".join(stable["gene"].head(top_n).astype(str).tolist()) or pd.NA,
            "top_ranked_genes": ";".join(ranked["gene"].head(top_n).astype(str).tolist()) or pd.NA,
        }
    )

    summary_path = output_dir / "hostresponse_summary.tsv"
    pd.DataFrame([summary]).to_csv(summary_path, sep="\t", index=False)
    return summary_path


def run_matched_hostresponse(
    run_dir: Path = DEFAULT_RUN_DIR,
    paper_barcodes: Path = DEFAULT_PAPER_BARCODES,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    detection_threshold: int = 10,
    use_hvg: bool = True,
    n_stab_iter: int = 100,
    stab_min_prob: float = 0.6,
    n_seeds: int = 6,
) -> Path:
    prepared = prepare_matched_inputs(run_dir, paper_barcodes, output_dir)
    run_hostresponse(
        virus_h5ad=str(prepared.ebv_h5ad),
        host_h5ad=str(prepared.host_h5ad),
        viral_accessions_file=str(prepared.analysis_txt),
        out_dir=str(output_dir),
        use_hvg=use_hvg,
        seeds=DEFAULT_SEEDS[:n_seeds],
        n_stab_iter=n_stab_iter,
        stab_min_prob=stab_min_prob,
        detection_threshold=detection_threshold,
        do_enrichment=False,
    )
    return write_hostresponse_summary(output_dir, prepared.n_cells, detection_threshold, stab_min_prob)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--paper-barcodes", type=Path, default=DEFAULT_PAPER_BARCODES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--detection-threshold", type=int, default=10)
    parser.add_argument("--n-stab-iter", type=int, default=100)
    parser.add_argument("--stab-min-prob", type=float, default=0.6)
    parser.add_argument("--n-seeds", type=int, default=6)
    parser.add_argument("--use-hvg", action=argparse.BooleanOptionalAction, default=True)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    summary_path = run_matched_hostresponse(
        run_dir=args.run_dir,
        paper_barcodes=args.paper_barcodes,
        output_dir=args.output_dir,
        detection_threshold=args.detection_threshold,
        use_hvg=args.use_hvg,
        n_stab_iter=args.n_stab_iter,
        stab_min_prob=args.stab_min_prob,
        n_seeds=args.n_seeds,
    )
    print(summary_path)


if __name__ == "__main__":
    main()
