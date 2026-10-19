"""Optional SCTK QC annotations stay in place and never filter cells."""

import numpy as np
import pytest


def test_scautoqc_adds_flags_and_qc_umap_without_changing_matrix():
    ad = pytest.importorskip("anndata")
    pytest.importorskip("sctk")
    from scipy import sparse

    from viralscan.scautoqc import add_scautoqc

    rng = np.random.default_rng(19)
    X = sparse.csr_matrix(rng.poisson(3, size=(80, 65)).astype(np.float32))
    adata = ad.AnnData(X=X.copy())
    symbols = (
        ["MT-ND1", "RPLP0", "RPS3"] + [f"GENE{i}" for i in range(59)] + ["HBA1", "HBEGF", "EBV_gp1"]
    )
    adata.var_names = symbols
    adata.var["gene_name"] = symbols
    adata.var["gene_role"] = ["host"] * 64 + ["viral"]
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

    base = ["-s1", "R1.fq", "-s2", "R2.fq", "-x", "10xv3"]
    default_args = build_parser().parse_args(base)
    assert default_args.scautoqc is False
    enabled_args = build_parser().parse_args([*base, "--scautoqc"])
    cfg = _build_run_config(enabled_args, "out/", "index.idx", "t2g.txt", None, "R1.fq", "R2.fq")
    assert cfg.scautoqc is True
    assert cfg.to_dict()["scautoqc"] is True
