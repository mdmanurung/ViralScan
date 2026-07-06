from pathlib import Path

import pandas as pd
import pytest

from scripts.make_manuscript_figures import make_figures

# Tests a repo-root analysis script (scripts/make_manuscript_figures.py), not the
# installed viralscan package; excluded from the default (hermetic) test run.
pytestmark = pytest.mark.research


def _write_inputs(tmp_path: Path) -> dict[str, Path]:
    matched = tmp_path / "matched.tsv"
    pd.DataFrame(
        [
            {
                "section": "barcode_accounting",
                "key": "shared_barcodes",
                "value": 1906,
                "note": "anchor",
            },
            {"section": "tier_star", "key": "pct_ebv_ge1umi", "value": 77.28, "note": ""},
            {"section": "tier_star", "key": "pct_ebv_ge10umi", "value": 10.23, "note": ""},
            {"section": "tier_vs", "key": "pct_ebv_ge1umi", "value": 93.91, "note": ""},
            {"section": "tier_vs", "key": "pct_ebv_ge10umi", "value": 13.48, "note": ""},
        ]
    ).to_csv(matched, sep="\t", index=False)

    per_gene = tmp_path / "per_gene.tsv"
    pd.DataFrame(
        [
            {"gene_id": "EPSTEIN_HHV4_BRLF1", "star_total_umi": 77, "vs_total_umi": 57},
            {"gene_id": "EPSTEIN_HHV4_EBNA-2", "star_total_umi": 0, "vs_total_umi": 716},
            {"gene_id": "EPSTEIN_HHV4_LMP-1", "star_total_umi": 47220, "vs_total_umi": 99},
        ]
    ).to_csv(per_gene, sep="\t", index=False)

    host_summary = tmp_path / "hostresponse_summary.tsv"
    pd.DataFrame(
        [
            {
                "virus": "Epstein-Barr virus",
                "n_cells": 1906,
                "n_positive": 257,
                "n_negative": 1649,
                "auc_mean": 0.71,
                "balanced_acc_mean": 0.66,
                "n_stable_genes": 0,
                "top_ranked_genes": "GENE1;GENE2",
            }
        ]
    ).to_csv(host_summary, sep="\t", index=False)
    return {"matched": matched, "per_gene": per_gene, "host_summary": host_summary}


def test_make_figures_writes_expected_pngs_and_pdfs(tmp_path: Path) -> None:
    inputs = _write_inputs(tmp_path)
    out_dir = tmp_path / "figures"

    outputs = make_figures(
        matched_comparison=inputs["matched"],
        per_gene_comparison=inputs["per_gene"],
        hostresponse_summary=inputs["host_summary"],
        output_dir=out_dir,
    )

    expected = {
        out_dir / "figure1_workflow.png",
        out_dir / "figure1_workflow.pdf",
        out_dir / "figure2_benchmark.png",
        out_dir / "figure2_benchmark.pdf",
    }
    assert expected.issubset(set(outputs))
    for path in expected:
        assert path.stat().st_size > 0


def test_make_figures_names_missing_required_input(tmp_path: Path) -> None:
    inputs = _write_inputs(tmp_path)
    missing = tmp_path / "missing.tsv"

    with pytest.raises(FileNotFoundError, match=str(missing)):
        make_figures(
            matched_comparison=missing,
            per_gene_comparison=inputs["per_gene"],
            hostresponse_summary=inputs["host_summary"],
            output_dir=tmp_path / "figures",
        )
