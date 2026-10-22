"""Optional SCTK QC annotations stay in place and never filter cells."""

import numpy as np
import pytest


def test_scautoqc_adds_flags_and_qc_umap_without_changing_matrix():
    ad = pytest.importorskip("anndata")
    pytest.importorskip("sctk")
    pytest.importorskip("leidenalg")  # QC-cluster step; not in the base test env
    from scipy import sparse

    from viralscan.scautoqc import add_scautoqc

    rng = np.random.default_rng(19)
    X = sparse.csr_matrix(rng.poisson(3, size=(80, 665)).astype(np.float32))
    adata = ad.AnnData(X=X.copy())
    symbols = (
        ["MT-ND1", "RPLP0", "RPS3"]
        + [f"GENE{i}" for i in range(659)]
        + ["HBA1", "HBEGF", "EBV_gp1"]
    )
    adata.var_names = symbols
    adata.var["gene_name"] = symbols
    adata.var["gene_role"] = ["host"] * 664 + ["viral"]
    before = adata.X.copy()
    n_obs = adata.n_obs

    add_scautoqc(adata, random_state=7)

    assert adata.n_obs == n_obs
    np.testing.assert_array_equal(adata.X.indptr, before.indptr)
    np.testing.assert_array_equal(adata.X.indices, before.indices)
    np.testing.assert_array_equal(adata.X.data, before.data)
    assert adata.var.loc["HBA1", "scautoqc_hb"]
    assert not adata.var.loc["HBEGF", "scautoqc_hb"]
    assert not adata.var.loc["EBV_gp1", "scautoqc_hb"]
    assert "scautoqc_cell_passed_qc" in adata.obs
    assert "scautoqc_cluster_passed_qc" in adata.obs
    assert "scautoqc_qc_cluster" in adata.obs
    assert adata.obsm["scautoqc_X_umap_qc"].shape == (n_obs, 2)
    assert adata.uns["scautoqc"]["cell_removal_performed"] is False


def test_scautoqc_cli_is_opt_in_and_recorded_in_run_config():
    from viralscan.menu import _build_run_config, build_parser

    base = ["run", "--r1", "R1.fq", "--r2", "R2.fq", "-o", "out", "-x", "10xv3"]
    default_args = build_parser().parse_args(base)
    assert default_args.scautoqc is False
    enabled_args = build_parser().parse_args([*base, "--scautoqc"])
    cfg = _build_run_config(enabled_args, "out/", "index.idx", "t2g.txt", None, "R1.fq", "R2.fq")
    assert cfg.scautoqc is True
    assert cfg.to_dict()["scautoqc"] is True


def _toy(ad, n_obs=40):
    from scipy import sparse

    rng = np.random.default_rng(3)
    # scanpy's top-n expression metrics need at least 500 features
    genes = ["MT-ND1", "RPLP0", "HBA1"] + [f"G{i}" for i in range(597)] + ["V1"]
    adata = ad.AnnData(X=sparse.csr_matrix(rng.poisson(3, size=(n_obs, 601)).astype(np.float32)))
    adata.var_names = genes
    adata.var["gene_name"] = list(adata.var_names)
    adata.var["gene_role"] = ["host"] * 600 + ["viral"]
    return adata


def test_flag_off_adds_no_config_key_and_no_rule(tmp_path):
    from pathlib import Path

    from viralscan.menu import _build_run_config, build_parser
    from viralscan.runconfig import RunConfig

    base = ["run", "--r1", "R1.fq", "--r2", "R2.fq", "-o", "out", "-x", "10xv3"]
    cfg = _build_run_config(build_parser().parse_args(base), "out/", "i", "t", None, "a", "b")
    assert "scautoqc" not in cfg.to_dict()
    cfg.to_yaml(tmp_path / "config.yaml")
    assert RunConfig.from_yaml(tmp_path / "config.yaml").scautoqc is False
    snakefile = Path(__file__).parent.parent / "src" / "viralscan" / "Snakefile"
    # the rule and the detection input are both guarded by the typed config value
    assert snakefile.read_text().count('config.get("scautoqc", False)') == 2


def test_missing_sctk_gives_clear_error(monkeypatch):
    ad = pytest.importorskip("anndata")
    import sys

    from viralscan.scautoqc import add_scautoqc

    monkeypatch.setitem(sys.modules, "sctk", None)  # makes `import sctk` raise ImportError
    with pytest.raises(RuntimeError, match=r"viralscan\[scautoqc\]"):
        add_scautoqc(_toy(ad))


def test_missing_sctk_stops_the_run_at_startup(monkeypatch):
    import importlib.util

    from viralscan import menu

    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)
    with pytest.raises(SystemExit):
        menu._check_scautoqc_deps()


def test_all_cells_fail_does_not_remove_or_crash(monkeypatch):
    """SCTK 0.2.0 raises TypeError in its all-fail fallback; the guard must absorb exactly that."""
    ad = pytest.importorskip("anndata")
    sctk = pytest.importorskip("sctk")

    from viralscan.scautoqc import add_scautoqc

    def all_fail(host_adata, cell_qc_key, **_):
        host_adata.obs[cell_qc_key] = False
        raise TypeError("list indices must be integers or slices, not str")

    monkeypatch.setattr(sctk, "cellwise_qc", all_fail)
    adata = _toy(ad, 5)  # <6 cells: clustering (leidenalg) is skipped
    add_scautoqc(adata)
    assert adata.n_obs == 5
    assert not adata.obs["scautoqc_cell_passed_qc"].any()


def test_unrelated_type_error_is_not_swallowed(monkeypatch):
    ad = pytest.importorskip("anndata")
    sctk = pytest.importorskip("sctk")

    from viralscan.scautoqc import add_scautoqc

    def boom(host_adata, cell_qc_key, **_):
        host_adata.obs[cell_qc_key] = False  # all-fail state, but the message is unrelated
        raise TypeError("something else")

    monkeypatch.setattr(sctk, "cellwise_qc", boom)
    with pytest.raises(TypeError, match="something else"):
        add_scautoqc(_toy(ad))


def test_flags_only_n_obs_and_cell_set_recorded():
    ad = pytest.importorskip("anndata")
    pytest.importorskip("sctk")

    from viralscan.scautoqc import add_scautoqc

    adata = _toy(ad, 5)  # <6 cells: clustering (leidenalg) is skipped
    names = adata.obs_names.copy()
    add_scautoqc(adata, random_state=1)
    assert adata.n_obs == 5 and adata.obs_names.equals(names)
    assert adata.uns["scautoqc"]["n_cells_scored"] == 5
    assert "before cell calling" in adata.uns["scautoqc"]["cell_set"]


def test_host_only_metrics_and_curated_globins_without_leiden():
    ad = pytest.importorskip("anndata")
    pytest.importorskip("sctk")

    from viralscan.scautoqc import add_scautoqc

    adata = _toy(ad, 5)
    names = list(adata.var_names)
    names[3] = "HBEGF"
    names[4] = "mt-Nd5"  # mouse-style lowercase symbol
    adata.var_names = names
    adata.var["gene_name"] = names
    host_sum = np.asarray(adata.X[:, :600].sum(axis=1)).ravel()
    all_sum = np.asarray(adata.X.sum(axis=1)).ravel()
    add_scautoqc(adata)
    assert not adata.var.loc["HBEGF", "scautoqc_hb"]
    assert adata.var.loc["HBA1", "scautoqc_hb"]
    assert adata.var.loc["mt-Nd5", "scautoqc_mito"]
    assert not adata.var.loc["V1", "scautoqc_mito"]
    np.testing.assert_allclose(adata.obs["scautoqc_n_counts"], host_sum)
    assert not np.allclose(host_sum, all_sum)


def test_biomart_failure_is_actionable(monkeypatch):
    sc = pytest.importorskip("scanpy")

    from viralscan.scautoqc import query_gene_annotations

    def down(*_a, **_k):
        raise ConnectionError("no route")

    monkeypatch.setattr(sc.queries, "biomart_annotations", down)
    with pytest.raises(RuntimeError, match="BioMart"):
        query_gene_annotations(["ENSG00000000001"])
