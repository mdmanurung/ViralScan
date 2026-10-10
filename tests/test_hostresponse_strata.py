"""Independent strata/provenance contracts using hand-labelled synthetic cells."""

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.scripts import hostresponse as hr


def inputs(tmp_path):
    cells = [f"cell_{i:03}" for i in range(60)]
    obs = pd.DataFrame(index=cells)
    obs["donor"] = np.repeat([f"d{i}" for i in range(6)], 10)
    obs["type"] = np.repeat(["T/CD4", "B cells"], 30)
    positive = np.tile([True] * 4 + [False] * 6, 6)
    rng = np.random.default_rng(29)
    raw = rng.poisson(20, (60, 4)).astype(np.float32)
    raw[positive, 0] += 20
    host = ad.AnnData(sp.csr_matrix(raw), obs=obs, var=pd.DataFrame(index=list("abcd")))
    virus = ad.AnnData(
        sp.csr_matrix(positive.astype(np.float32).reshape(-1, 1)),
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=["V"]),
    )
    host.write_h5ad(tmp_path / "host.h5ad")
    virus.write_h5ad(tmp_path / "virus.h5ad")
    (tmp_path / "accessions.txt").write_text("V\n")
    return dict(
        host_h5ad=str(tmp_path / "host.h5ad"),
        virus_h5ad=str(tmp_path / "virus.h5ad"),
        viral_accessions_file=str(tmp_path / "accessions.txt"),
        out_dir=str(tmp_path / "out"),
        seeds=[0],
        use_hvg=False,
        n_stab_iter=1,
        stab_min_prob=2,
        detection_threshold=1,
        control_mito=False,
        cv_mode="group",
        groups_column="donor",
        cv_folds=3,
        cell_type_column="type",
    )


def test_selected_strata_use_separate_cohorts_and_manifest_hashes(tmp_path):
    kw = inputs(tmp_path)
    prov = hr.run_hostresponse(**kw, cell_types=["T/CD4", "B cells"])
    out = Path(kw["out_dir"])
    status = pd.read_csv(out / hr.STATUS_FILENAME, sep="\t")
    assert set(status["stratum"]) == {"T/CD4", "B cells"}
    assert set(status["status"]) == {"ok"}
    assert list(status["n_positive"]) == [12, 12]
    assert list(status["n_negative"]) == [18, 18]
    assert list(status["n_groups"]) == [3, 3]
    assert not (out / "V_gene_weights.csv").exists()
    for cell_type, donors in [("T_CD4", {"d0", "d1", "d2"}), ("B_cells", {"d3", "d4", "d5"})]:
        folds = pd.read_csv(out / f"cell_type_{cell_type}/V_cv_folds.tsv", sep="\t")
        assert set(folds["group"]) == donors
        assert len(folds) == 30
    manifest = json.loads((out / hr.MANIFEST_FILENAME).read_text())
    assert manifest == prov
    assert manifest["cell_type_filter"] == ["T/CD4", "B cells"]
    assert manifest["descriptive_inference_unit"] == "cell"
    for name, digest in manifest["artifacts"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest
    assert any(name.endswith("_cv_group_auc.tsv") for name in manifest["artifacts"])
    assert prov["strata"]["T/CD4"]["viruses"]["V"]["fold_features"]


def test_viruses_status_keeps_one_entry_per_stratum(tmp_path):
    prov = hr.run_hostresponse(**inputs(tmp_path), cell_types=["T/CD4", "B cells"])
    assert prov["viruses_status"] == {"T/CD4::V": "ok", "B cells::V": "ok"}


def test_stratum_specific_support_and_stale_cleanup_preserves_user_files(tmp_path):
    kw = inputs(tmp_path)
    hr.run_hostresponse(**kw, cell_types=["T/CD4", "B cells"])
    out = Path(kw["out_dir"])
    user_file = out / "cell_type_T_CD4" / "my_notes.txt"
    user_file.write_text("keep")
    prov = hr.run_hostresponse(**kw, cell_types=["T/CD4"], min_groups=4)
    status = pd.read_csv(out / hr.STATUS_FILENAME, sep="\t")
    assert status.loc[0, "status"] == "skipped_min_groups"
    assert prov["status"] == "no_eligible_virus"
    assert not (out / "hostresponse_metrics.csv").exists()
    assert not list(out.rglob("*_gene_weights.csv"))
    assert not (out / "cell_type_B_cells").exists()
    assert user_file.read_text() == "keep"


@pytest.mark.parametrize("types", [["unknown"], ["T/CD4", "T/CD4"], [""]])
def test_invalid_stratum_requests_preserve_existing_outputs(tmp_path, types):
    kw = inputs(tmp_path)
    out = Path(kw["out_dir"])
    out.mkdir()
    sentinel = out / "V_gene_weights.csv"
    sentinel.write_text("previous")
    with pytest.raises(ValueError, match="cell type|collision"):
        hr.run_hostresponse(**kw, cell_types=types)
    assert sentinel.read_text() == "previous"


def test_sanitization_collisions_fail_before_cleanup(tmp_path):
    kw = inputs(tmp_path)
    host = ad.read_h5ad(kw["host_h5ad"])
    host.obs["type"] = np.repeat(["T/CD4", "T_CD4"], 30)
    host.write_h5ad(kw["host_h5ad"])
    out = Path(kw["out_dir"])
    out.mkdir()
    sentinel = out / "V_gene_weights.csv"
    sentinel.write_text("previous")
    with pytest.raises(ValueError, match="collision"):
        hr.run_hostresponse(**kw, cell_types=["T/CD4", "T_CD4"])
    assert sentinel.read_text() == "previous"


def test_missing_raw_group_depth_preserves_outputs(tmp_path):
    kw = inputs(tmp_path)
    host = ad.read_h5ad(kw["host_h5ad"])
    host.uns["log1p"] = {}
    host.write_h5ad(kw["host_h5ad"])
    out = Path(kw["out_dir"])
    out.mkdir()
    sentinel = out / "V_gene_weights.csv"
    sentinel.write_text("previous")
    with pytest.raises(ValueError, match="raw depth"):
        hr.run_hostresponse(**kw)
    assert sentinel.read_text() == "previous"


def test_no_matching_virus_writes_status_and_manifest_clears_stale_summary(tmp_path):
    kw = inputs(tmp_path)
    out = Path(kw["out_dir"])
    out.mkdir()
    (out / "hostresponse_metrics.csv").write_text("stale")
    Path(kw["viral_accessions_file"]).write_text("unmatched\n")
    prov = hr.run_hostresponse(**kw)
    assert prov["status"] == "no_viral_accessions_matched"
    assert not (out / "hostresponse_metrics.csv").exists()
    assert json.loads((out / hr.MANIFEST_FILENAME).read_text()) == prov
    assert pd.read_csv(out / hr.STATUS_FILENAME, sep="\t").loc[0, "status"] == prov["status"]


def test_invalid_observed_retains_every_requested_null_replicate():
    def invalid(y):
        raise hr.FoldDesignError("invalid_folds", "both classes unavailable")

    result = hr.run_structured_null(
        invalid,
        [0, 1],
        ["a", "b"],
        ["d1", "d2"],
        unit="group",
        n_permutations=3,
    )
    assert len(result.replicates) == 3
    assert list(result.replicates["permutation_id"]) == [1, 2, 3]
    assert result.summary["n_failed"] == 3
    assert result.summary["empirical_p_auc"] is None


@pytest.mark.parametrize(
    "values,observed,requested", [([0.5], 0.7, 2), ([np.nan], 0.7, 1), ([0.5], np.inf, 1)]
)
def test_empirical_p_rejects_incomplete_or_nonfinite_null(values, observed, requested):
    with pytest.raises(ValueError, match="complete finite"):
        hr.empirical_p_value(values, observed, requested)


def test_standalone_parser_preserves_defaults_and_group_options():
    parser = hr._build_parser()
    flags = ["--virus-h5ad", "v", "--host-h5ad", "h", "--viral-accessions", "a", "--output", "o"]
    defaults = parser.parse_args(flags)
    assert defaults.cv == "cell" and defaults.permutations == 0
    args = parser.parse_args(
        flags
        + [
            "--cv",
            "group",
            "--groups",
            "donor",
            "--cell-type-column",
            "type",
            "--cell-types",
            "T, B",
        ]
    )
    assert args.groups == "donor" and args.cell_types == ["T", " B"]


def test_missing_group_splitter_reports_explicit_reason(monkeypatch):
    from sklearn import model_selection

    monkeypatch.delattr(model_selection, "StratifiedGroupKFold")
    with pytest.raises(hr.FoldDesignError, match="splitter_unavailable"):
        hr.make_group_folds([0, 1, 0, 1], ["a", "a", "b", "b"], list("1234"), n_folds=2)
    # Default cell splitting remains available with the original splitter API.
    depth = np.arange(20) + 100
    assert (
        hr._balanced_split(np.arange(10), np.arange(10, 20), depth, depth[:, None], 0, 1.0)
        is not None
    )
