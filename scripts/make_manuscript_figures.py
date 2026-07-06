#!/usr/bin/env python3
"""Generate manuscript Figures 1-2 from benchmark and host-response outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_MATCHED = Path("results/matched_barcode_comparison.tsv")
DEFAULT_PER_GENE = Path("results/matched_barcode_comparison_per_gene.tsv")
DEFAULT_HOST_SUMMARY = Path("results/hostresponse_ebv_matched/hostresponse_summary.tsv")
DEFAULT_OUTPUT_DIR = Path("docs/figures")


def _require(path: Path) -> Path:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def _metric(table: pd.DataFrame, section: str, key: str) -> float:
    hit = table[(table["section"] == section) & (table["key"] == key)]
    if hit.empty:
        raise KeyError(f"Missing metric {section}/{key}")
    return float(hit.iloc[0]["value"])


def _save(fig: plt.Figure, output_dir: Path, stem: str) -> list[Path]:
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.pdf"]
    for path in paths:
        fig.savefig(path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return paths


def _figure1(output_dir: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.axis("off")
    boxes = [
        ("FASTQ", "paired reads\nbarcode + UMI"),
        ("Combined reference", "host + virus\ntranscriptome"),
        ("Pseudoalignment", "kallisto/bustools\ncell-gene counts"),
        ("Multimap EM", "allocate ambiguous\nhost-virus reads"),
        ("Outputs", "viral burden\nhost response"),
    ]
    x = np.linspace(0.08, 0.92, len(boxes))
    for idx, ((title, body), xpos) in enumerate(zip(boxes, x, strict=True)):
        rect = plt.Rectangle((xpos - 0.085, 0.38), 0.17, 0.28, facecolor="#f2f6f8", edgecolor="#314d5b", lw=1.4)
        ax.add_patch(rect)
        ax.text(xpos, 0.58, title, ha="center", va="center", fontsize=11, fontweight="bold", color="#1d2f38")
        ax.text(xpos, 0.47, body, ha="center", va="center", fontsize=9, color="#314d5b")
        if idx < len(boxes) - 1:
            ax.annotate(
                "",
                xy=(x[idx + 1] - 0.095, 0.52),
                xytext=(xpos + 0.095, 0.52),
                arrowprops={"arrowstyle": "->", "lw": 1.5, "color": "#314d5b"},
            )
    ax.text(0.5, 0.22, "ViralScan quantifies viral UMIs per cell and preserves host predictors for downstream modeling.", ha="center", fontsize=10)
    return _save(fig, output_dir, "figure1_workflow")


def _figure2(matched: pd.DataFrame, per_gene: pd.DataFrame, host_summary: pd.DataFrame, output_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), gridspec_kw={"width_ratios": [1.05, 1.05, 1.0]})

    datasets = ["HHV-6B\nCAR-T", "EBV\nLCL", "HSV-1\nfibroblast"]
    published = [0.2, 2.2, 16.0]
    viralscan = [0.152, 2.73, 15.6]
    xpos = np.arange(len(datasets))
    width = 0.36
    axes[0].bar(xpos - width / 2, published, width, label="Published", color="#819a7a")
    axes[0].bar(xpos + width / 2, viralscan, width, label="ViralScan", color="#3c6e71")
    axes[0].set_xticks(xpos, datasets)
    axes[0].set_ylabel("Infected or lytic cells (%)")
    axes[0].set_title("Published benchmarks")
    axes[0].legend(frameon=False, fontsize=8)

    star_ge1 = _metric(matched, "tier_star", "pct_ebv_ge1umi")
    vs_ge1 = _metric(matched, "tier_vs", "pct_ebv_ge1umi")
    star_ge10 = _metric(matched, "tier_star", "pct_ebv_ge10umi")
    vs_ge10 = _metric(matched, "tier_vs", "pct_ebv_ge10umi")
    categories = ["EBV >=1 UMI", "EBV >=10 UMI"]
    axes[1].bar(np.arange(2) - width / 2, [star_ge1, star_ge10], width, label="STARsolo", color="#9a6b4f")
    axes[1].bar(np.arange(2) + width / 2, [vs_ge1, vs_ge10], width, label="ViralScan", color="#4f7cac")
    axes[1].set_xticks(np.arange(2), categories)
    axes[1].set_ylabel("Matched cells (%)")
    axes[1].set_title("EBV matched-cell anchor")
    axes[1].legend(frameon=False, fontsize=8)

    top = per_gene.assign(total=lambda df: df["star_total_umi"] + df["vs_total_umi"]).nlargest(5, "total")
    labels = top["gene_id"].str.replace("EPSTEIN_HHV4_", "", regex=False).str.replace("EPSTEIN_", "", regex=False)
    y = np.arange(len(top))
    axes[2].barh(y - 0.18, top["star_total_umi"], 0.36, label="STARsolo", color="#9a6b4f")
    axes[2].barh(y + 0.18, top["vs_total_umi"], 0.36, label="ViralScan", color="#4f7cac")
    axes[2].set_yticks(y, labels)
    axes[2].set_xscale("symlog", linthresh=1)
    axes[2].set_xlabel("Total EBV UMI")
    axes[2].set_title("Gene attribution")
    axes[2].legend(frameon=False, fontsize=8)

    if not host_summary.empty:
        row = host_summary.iloc[0]
        fig.text(
            0.5,
            -0.02,
            f"Host-response matched anchor: n={int(row['n_cells'])}, EBV >=10 UMI positives={int(row['n_positive'])}, "
            f"raw AUROC={float(row.get('auc_mean', np.nan)):.2f}, balanced accuracy={float(row.get('balanced_acc_mean', np.nan)):.2f}.",
            ha="center",
            fontsize=9,
        )
        # The raw >=10 corrected-UMI label tracks sequencing depth; report the
        # depth-alone baseline and depth-controlled estimates alongside the raw AUROC
        # so the figure does not present a confounded number without caveat.
        fig.text(
            0.5,
            -0.08,
            "Raw label is depth-confounded (depth alone AUROC 0.97, same design); "
            "depth-controlled AUROC 0.64 (CPM label) to 0.72 (depth-matched).",
            ha="center",
            fontsize=8,
            style="italic",
            color="#7a2f2f",
        )
    fig.tight_layout()
    return _save(fig, output_dir, "figure2_benchmark")


def make_figures(
    matched_comparison: Path = DEFAULT_MATCHED,
    per_gene_comparison: Path = DEFAULT_PER_GENE,
    hostresponse_summary: Path = DEFAULT_HOST_SUMMARY,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[Path]:
    matched_path = _require(Path(matched_comparison))
    per_gene_path = _require(Path(per_gene_comparison))
    host_path = _require(Path(hostresponse_summary))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    matched = pd.read_csv(matched_path, sep="\t")
    per_gene = pd.read_csv(per_gene_path, sep="\t")
    host = pd.read_csv(host_path, sep="\t")
    outputs: list[Path] = []
    outputs.extend(_figure1(output_dir))
    outputs.extend(_figure2(matched, per_gene, host, output_dir))
    return outputs


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matched-comparison", type=Path, default=DEFAULT_MATCHED)
    parser.add_argument("--per-gene-comparison", type=Path, default=DEFAULT_PER_GENE)
    parser.add_argument("--hostresponse-summary", type=Path, default=DEFAULT_HOST_SUMMARY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    for path in make_figures(
        matched_comparison=args.matched_comparison,
        per_gene_comparison=args.per_gene_comparison,
        hostresponse_summary=args.hostresponse_summary,
        output_dir=args.output_dir,
    ):
        print(path)


if __name__ == "__main__":
    main()
