#!/usr/bin/env python3
"""Summarize the broad covid_viralscan viral survey results.

Produces a Markdown summary covering:
  1. Top-N viral panel hits ranked by combined UMI across all samples.
  2. Specificity control: SARS-CoV-2 (NC_045512.2) vs SARS-CoV-1 (sarsp* / NC_004718.3).
  3. Cell-barcode overlap: how many viral+ barcodes are CellRanger-called cells.

Usage::

    PYTHONPATH=src python covid_viralscan/scripts/summarize_survey.py \\
        --results-dir covid_viralscan/results \\
        --samples LUM-SJ-x213-g LUM-SJ-x216-g \\
        --cellranger-outs /path/to/202502341a_count_v2/outs \\
        --output covid_viralscan/results/SURVEY_SUMMARY.md
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

# ── Accession / pattern helpers ──────────────────────────────────────────────

# SARS-CoV-2: whole-genome GTF produces a single gene_id "NC_045512.2_gene1".
_SARSCOV2_PREFIX = "NC_045512"

# SARS-CoV-1: Serratus RefSeq annotation uses locus_tag "sarsp{N}" as gene_id.
_SARSCOV1_PREFIX = "sarsp"


def _is_sarscov2(virus_name: str) -> bool:
    return virus_name.startswith(_SARSCOV2_PREFIX)


def _is_sarscov1(virus_name: str) -> bool:
    return virus_name.startswith(_SARSCOV1_PREFIX) or "NC_004718" in virus_name


# ── Name prettifier ──────────────────────────────────────────────────────────

try:
    from viralscan.constants import VIRUS_NAME_MAP  # type: ignore[import]
    from viralscan.anellovirus import merged_name_map  # type: ignore[import]

    _NAME_MAP: dict[str, str] = merged_name_map()
except ImportError:
    _NAME_MAP = {}

# Override for the two SARS targets that won't be in the map
_DISPLAY_NAMES: dict[str, str] = {
    _SARSCOV2_PREFIX: "SARS-CoV-2 (NC_045512.2)",
    _SARSCOV1_PREFIX: "SARS-CoV-1 (NC_004718.3, negative control)",
}


def _display(virus_name: str) -> str:
    if _is_sarscov2(virus_name):
        return "SARS-CoV-2 (NC_045512.2)"
    if _is_sarscov1(virus_name):
        return "SARS-CoV-1 / sarsp (negative control)"
    # Try the merged name map (covers anellovirus genera + bundled 195-virus panel)
    for key in sorted(_NAME_MAP, key=len, reverse=True):
        if virus_name == key or (
            virus_name.startswith(key)
            and len(virus_name) > len(key)
            and (virus_name[len(key)] == "_" or virus_name[len(key)].isdigit())
        ):
            return _NAME_MAP[key]
    return virus_name


# ── I/O helpers ──────────────────────────────────────────────────────────────


def _read_viral_summary(results_dir: Path, sample: str) -> pd.DataFrame:
    """Read viral_summary.tsv for one sample; add a 'sample' column."""
    path = results_dir / sample / "results" / "viral_summary.tsv"
    if not path.exists():
        sys.exit(
            f"ERROR: missing {path}\n"
            "Run Stage 3 (viralscan quant) before this script."
        )
    df = pd.read_csv(path, sep="\t")
    df["sample"] = sample
    return df


def _read_per_cell(results_dir: Path, sample: str) -> pd.DataFrame:
    """Read per_cell_viral.tsv for one sample; add a 'sample' column."""
    path = results_dir / sample / "results" / "per_cell_viral.tsv"
    if not path.exists():
        sys.exit(f"ERROR: missing {path}")
    df = pd.read_csv(path, sep="\t")
    df["sample"] = sample
    return df


def _load_cellranger_barcodes(outs_dir: Path) -> set[str]:
    """Load called-cell barcodes from CellRanger filtered_feature_bc_matrix."""
    bc_file = outs_dir / "filtered_feature_bc_matrix" / "barcodes.tsv.gz"
    if not bc_file.exists():
        # Fallback: uncompressed
        bc_file_plain = outs_dir / "filtered_feature_bc_matrix" / "barcodes.tsv"
        if bc_file_plain.exists():
            return {
                line.strip().split("-")[0]  # strip GEM-group suffix
                for line in bc_file_plain.read_text().splitlines()
                if line.strip()
            }
        print(
            f"WARNING: CellRanger barcodes not found at {bc_file} — "
            "cell-overlap section will be skipped.",
            file=sys.stderr,
        )
        return set()
    barcodes: set[str] = set()
    with gzip.open(bc_file, "rt") as fh:
        for line in fh:
            bc = line.strip()
            if bc:
                barcodes.add(bc.split("-")[0])  # strip GEM-group suffix (-1)
    return barcodes


# ── Core analysis ─────────────────────────────────────────────────────────────


def _ranked_panel(
    summaries: list[pd.DataFrame],
    top_n: int = 50,
) -> pd.DataFrame:
    """Rank all viruses by combined total UMI across all samples."""
    combined = pd.concat(summaries, ignore_index=True)
    # Aggregate across samples
    agg = (
        combined.groupby("virus_name")
        .agg(
            combined_umi=("total_umi", "sum"),
            combined_infected_cells=("infected_cells", "sum"),
            n_samples=("sample", "nunique"),
        )
        .reset_index()
        .sort_values("combined_umi", ascending=False)
        .head(top_n)
    )
    agg["display_name"] = agg["virus_name"].apply(_display)
    return agg


def _sars_specificity(summaries: list[pd.DataFrame]) -> dict[str, dict]:
    """Extract SARS-CoV-2 and SARS-CoV-1 rows from per-sample summaries."""
    combined = pd.concat(summaries, ignore_index=True)
    result: dict[str, dict] = {}
    for label, fn in [("SARS-CoV-2 (NC_045512.2)", _is_sarscov2),
                      ("SARS-CoV-1 / sarsp (NC_004718.3)", _is_sarscov1)]:
        rows = combined[combined["virus_name"].apply(fn)]
        result[label] = {
            "total_umi": int(rows["total_umi"].sum()),
            "infected_cells": int(rows["infected_cells"].sum()),
            "n_rows": len(rows),
        }
    return result


def _cell_overlap(
    per_cell_dfs: list[pd.DataFrame],
    cellranger_barcodes: set[str],
    samples: list[str],
) -> dict[str, dict]:
    """Compute per-sample intersection of viral+ barcodes with CellRanger cells."""
    result: dict[str, dict] = {}
    for sample, df in zip(samples, per_cell_dfs):
        # Strip GEM-group suffix from ViralScan barcodes for comparison
        vs_barcodes = {bc.split("-")[0] for bc in df["barcode"].unique()}
        overlap = vs_barcodes & cellranger_barcodes
        result[sample] = {
            "vs_viral_barcodes": len(vs_barcodes),
            "cr_cells": len(cellranger_barcodes),
            "overlap": len(overlap),
            "pct_viral_in_cr": (
                round(100.0 * len(overlap) / len(vs_barcodes), 1)
                if vs_barcodes
                else 0.0
            ),
        }
    return result


# ── Markdown report ───────────────────────────────────────────────────────────


def _write_report(
    output_path: Path,
    samples: list[str],
    ranked: pd.DataFrame,
    sars: dict[str, dict],
    overlap: dict[str, dict],
    summaries: list[pd.DataFrame],
) -> None:
    lines: list[str] = []
    a = lines.append

    a("# covid_viralscan — Broad Viral Survey Summary")
    a("")
    a(f"**Samples:** {', '.join(samples)}")
    a(f"**Reference:** Combined human + SARS-CoV-2 + Serratus/anellovirus panel")
    a(f"**Method:** `--multimap-method host-conservative`, 10x 5′ v3")
    a("")

    # ── Section 1: ranked panel ──────────────────────────────────────────────
    a("## 1 — Broad viral panel (top-50 by combined UMI)")
    a("")
    a("| Rank | Virus | Gene ID | Combined UMI | Infected cells | Samples |")
    a("|------|-------|---------|-------------|----------------|---------|")
    for i, row in enumerate(ranked.itertuples(), 1):
        # Only show display_name if different from virus_name
        display = row.display_name
        gene_id = row.virus_name if display != row.virus_name else ""
        a(
            f"| {i} | {display} | {gene_id} | "
            f"{row.combined_umi:,} | {row.combined_infected_cells:,} | "
            f"{row.n_samples} |"
        )
    a("")

    # ── Section 2: per-sample breakdown ─────────────────────────────────────
    a("## 2 — Per-sample viral summaries")
    a("")
    for df in summaries:
        sample = df["sample"].iloc[0]
        total_cells = df["total_cells"].iloc[0] if len(df) else "N/A"
        total_umi = df["total_umi"].sum()
        n_viruses = len(df)
        a(f"### {sample}")
        a(f"- Total barcodes tested: {total_cells:,}")
        a(f"- Viruses detected (≥ threshold): {n_viruses}")
        a(f"- Combined viral UMI: {total_umi:,}")
        a("")

    # ── Section 3: SARS specificity ─────────────────────────────────────────
    a("## 3 — SARS-CoV-2 vs SARS-CoV-1 specificity control")
    a("")
    a("| Target | Total UMI | Infected cells | Interpretation |")
    a("|--------|-----------|----------------|----------------|")
    cov2_umi = sars.get("SARS-CoV-2 (NC_045512.2)", {}).get("total_umi", 0)
    cov1_umi = sars.get("SARS-CoV-1 / sarsp (NC_004718.3)", {}).get("total_umi", 0)
    for label, stats in sars.items():
        interp = ""
        if "CoV-2" in label:
            interp = "Primary target" if stats["total_umi"] > 0 else "Not detected"
        else:
            interp = (
                "⚠ Unexpected signal — check cross-mapping"
                if stats["total_umi"] > 0
                else "✓ Negative control clean"
            )
        a(
            f"| {label} | {stats['total_umi']:,} | "
            f"{stats['infected_cells']:,} | {interp} |"
        )
    a("")
    if cov1_umi > 0:
        a(
            "> **Warning:** SARS-CoV-1 signal detected. Inspect with "
            "`viralscan evidence` to rule out cross-mapping from the "
            "~79% genome-identity overlap with SARS-CoV-2."
        )
        a("")

    # ── Section 4: cell-barcode overlap ──────────────────────────────────────
    a("## 4 — Viral+ barcodes vs CellRanger called cells")
    a("")
    if all(v["cr_cells"] == 0 for v in overlap.values()):
        a("*CellRanger barcodes not available — overlap not computed.*")
    else:
        a("| Sample | ViralScan viral+ barcodes | CellRanger cells | Overlap | % viral+ in cells |")
        a("|--------|--------------------------|------------------|---------|-------------------|")
        for sample, stats in overlap.items():
            a(
                f"| {sample} | {stats['vs_viral_barcodes']:,} | "
                f"{stats['cr_cells']:,} | {stats['overlap']:,} | "
                f"{stats['pct_viral_in_cr']}% |"
            )
        a("")
        a(
            "> High % (>70%) = viral signal concentrated in called cells → "
            "confident cell-level detection.  "
            "Low % = signal may be ambient RNA or empty droplets."
        )
    a("")

    # ── Section 5: next steps ────────────────────────────────────────────────
    a("## 5 — Next steps")
    a("")
    a(
        "- If SARS-CoV-2 detected: run `viralscan evidence` to validate reads, "
        "then cross-ref with clinical metadata."
    )
    a(
        "- If unexpected viruses rank highly: check with `viralscan evidence --blast` "
        "to confirm sequence identity."
    )
    a(
        "- Full anellovirus genus breakdown: the panel covers ~2,022 clareaulab "
        "representatives; any Alphatorquevirus/Betatorquevirus/Gammatorquevirus "
        "hits represent anellovirus reactivation."
    )

    output_path.write_text("\n".join(lines) + "\n")
    print(f"Written: {output_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--results-dir",
        required=True,
        type=Path,
        metavar="DIR",
        help="Directory containing per-sample result subdirs (covid_viralscan/results/).",
    )
    p.add_argument(
        "--samples",
        nargs="+",
        required=True,
        metavar="SAMPLE",
        help="Sample names (must match subdirectory names under results-dir).",
    )
    p.add_argument(
        "--cellranger-outs",
        type=Path,
        default=None,
        metavar="DIR",
        help=(
            "Path to CellRanger `outs/` directory containing "
            "filtered_feature_bc_matrix/barcodes.tsv.gz. "
            "If omitted, cell-overlap section is skipped."
        ),
    )
    p.add_argument(
        "--output",
        required=True,
        type=Path,
        metavar="FILE",
        help="Output path for SURVEY_SUMMARY.md.",
    )
    p.add_argument(
        "--top-n",
        type=int,
        default=50,
        metavar="N",
        help="Number of top viruses to show in the ranked panel table (default: 50).",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    print(f"Loading results for {len(args.samples)} sample(s): {args.samples}")
    summaries = [_read_viral_summary(args.results_dir, s) for s in args.samples]
    per_cell_dfs = [_read_per_cell(args.results_dir, s) for s in args.samples]

    # CellRanger barcodes (optional)
    cr_barcodes: set[str] = set()
    if args.cellranger_outs is not None:
        cr_barcodes = _load_cellranger_barcodes(args.cellranger_outs)
        print(f"Loaded {len(cr_barcodes):,} CellRanger called-cell barcodes.")

    ranked = _ranked_panel(summaries, top_n=args.top_n)
    sars = _sars_specificity(summaries)
    overlap = _cell_overlap(per_cell_dfs, cr_barcodes, args.samples)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    _write_report(
        args.output,
        args.samples,
        ranked,
        sars,
        overlap,
        summaries,
    )


if __name__ == "__main__":
    main()
