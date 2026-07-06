from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scripts.hostresponse_ebv_matched import (
    EBV_NAME,
    prepare_matched_inputs,
    write_hostresponse_summary,
)

# Tests a repo-root analysis script (scripts/hostresponse_ebv_matched.py), not the
# installed viralscan package; excluded from the default (hermetic) test run.
pytestmark = pytest.mark.research


def test_prepare_matched_inputs_anchors_on_paper_and_whitelist_and_excludes_viral_genes(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    output_dir = tmp_path / "out"

    obs_names = ["AAAC-1", "BBBC-1", "CCCD-1", "DDDE-1"]
    var_names = ["ACTB", "MS4A1", "EPSTEIN_HHV4_EBNA-2", "EPSTEIN_HHV4_BRLF1", "HUM_HERP1_gene1"]
    counts = np.array(
        [
            [10, 1, 2, 3, 100],
            [20, 2, 5, 7, 200],
            [30, 3, 0, 0, 300],
            [40, 4, 11, 13, 400],
        ],
        dtype=np.float32,
    )
    ad.AnnData(
        counts, obs=pd.DataFrame(index=obs_names), var=pd.DataFrame(index=var_names)
    ).write_h5ad(run_dir / "adata.h5ad")
    ad.AnnData(
        counts, obs=pd.DataFrame(index=obs_names), var=pd.DataFrame(index=var_names)
    ).write_h5ad(run_dir / "adata_multimap.h5ad")

    paper_barcodes = tmp_path / "barcodes.tsv"
    paper_barcodes.write_text("BBBC\nAAAC-1\nZZZZ-1\n", encoding="utf-8")

    result = prepare_matched_inputs(run_dir, paper_barcodes, output_dir)

    assert result.n_cells == 2
    assert result.n_positive_at_threshold(10) == 1
    assert (output_dir / "matched_barcodes.tsv").read_text(encoding="utf-8").splitlines() == [
        "AAAC",
        "BBBC",
    ]
    assert (output_dir / "analysis.txt").read_text(encoding="utf-8").strip() == EBV_NAME

    host = ad.read_h5ad(output_dir / "host_only_matched.h5ad")
    assert host.obs_names.tolist() == ["AAAC", "BBBC"]
    assert host.var_names.tolist() == ["ACTB", "MS4A1"]

    burden = ad.read_h5ad(output_dir / "ebv_burden_matched.h5ad")
    assert burden.obs_names.tolist() == ["AAAC", "BBBC"]
    assert burden.var_names.tolist() == [EBV_NAME]
    np.testing.assert_array_equal(np.asarray(burden.X).ravel(), np.array([5.0, 12.0]))


def test_write_hostresponse_summary_reports_stable_and_top_ranked_genes(tmp_path: Path) -> None:
    metrics = pd.DataFrame(
        [
            {
                "virus": EBV_NAME,
                "n_positive": 2,
                "balanced_acc_mean": 0.75,
                "balanced_acc_sd": 0.05,
                "auc_mean": 0.8,
                "auc_sd": 0.03,
            }
        ]
    )
    metrics.to_csv(tmp_path / "hostresponse_metrics.csv", index=False)
    pd.DataFrame(
        [
            {"gene": "GENE1", "stab_prob": 0.7, "weight_mean": 2.0},
            {"gene": "GENE2", "stab_prob": 0.3, "weight_mean": -5.0},
        ]
    ).to_csv(tmp_path / "Epstein-Barr_virus_stability.csv", index=False)

    summary_path = write_hostresponse_summary(
        tmp_path, n_cells=5, detection_threshold=10, stab_min_prob=0.6
    )
    summary = pd.read_csv(summary_path, sep="\t")

    assert summary.loc[0, "n_cells"] == 5
    assert summary.loc[0, "n_positive"] == 2
    assert summary.loc[0, "n_negative"] == 3
    assert summary.loc[0, "n_stable_genes"] == 1
    assert summary.loc[0, "top_stable_genes"] == "GENE1"
    assert summary.loc[0, "top_ranked_genes"] == "GENE1;GENE2"


def test_write_hostresponse_summary_handles_no_stable_genes(tmp_path: Path) -> None:
    pd.DataFrame(
        [{"virus": EBV_NAME, "n_positive": 1, "auc_mean": 0.55, "balanced_acc_mean": 0.5}]
    ).to_csv(tmp_path / "hostresponse_metrics.csv", index=False)
    pd.DataFrame(
        [
            {"gene": "GENE1", "stab_prob": 0.2, "weight_mean": -3.0},
            {"gene": "GENE2", "stab_prob": 0.1, "weight_mean": 4.0},
        ]
    ).to_csv(tmp_path / "Epstein-Barr_virus_stability.csv", index=False)

    summary = pd.read_csv(
        write_hostresponse_summary(tmp_path, n_cells=4, detection_threshold=10), sep="\t"
    )

    assert summary.loc[0, "n_stable_genes"] == 0
    assert pd.isna(summary.loc[0, "top_stable_genes"])
    assert summary.loc[0, "top_ranked_genes"] == "GENE2;GENE1"
