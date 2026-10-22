"""Regression for the requested median-depth shifted logarithm, including caller reach."""

import hashlib

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.runconfig import RunConfig
from viralscan.scripts import hostresponse as hr
from viralscan.scripts import umap as um


def shifted_log_reference(raw):
    """Independent shifted-log formula for nonempty library totals."""
    totals = raw.sum(axis=1)
    median = np.median(totals[totals > 0])
    scaled = np.divide(
        raw * median, totals[:, None], out=np.zeros_like(raw), where=totals[:, None] > 0
    )
    return np.log1p(scaled)


def dense(matrix):
    return matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)


@pytest.mark.parametrize("sparse", [False, True])
@pytest.mark.parametrize("dtype", [np.int32, np.float32, np.float64])
def test_host_shifted_log_formula_and_raw_preservation(sparse, dtype):
    raw = np.array([[20, 0, 20], [1, 9, 10], [100, 0, 0]], dtype=dtype)
    data = ad.AnnData(sp.csr_matrix(raw) if sparse else raw.copy())
    data.layers["counts"] = data.X.copy()
    data.raw = data.copy()
    hr._detect_and_normalize(data)
    np.testing.assert_allclose(dense(data.X), shifted_log_reference(raw.astype(float)), rtol=1e-6)
    np.testing.assert_array_equal(dense(data.layers["counts"]), raw)
    np.testing.assert_array_equal(dense(data.raw.X), raw)
    np.testing.assert_array_equal(data.obs["_raw_depth"], raw.sum(axis=1))
    assert sp.issparse(data.X) is sparse
    before = dense(data.X).copy()
    hr._detect_and_normalize(data)
    np.testing.assert_array_equal(dense(data.X), before)
    np.testing.assert_array_equal(data.obs["_raw_depth"], raw.sum(axis=1))


@pytest.mark.parametrize("sparse", [False, True])
def test_empty_rows_stay_empty_and_finite(sparse):
    raw = np.array([[0, 0], [20, 20], [1, 19], [100, 0]], dtype=float)
    data = ad.AnnData(sp.csr_matrix(raw) if sparse else raw.copy())
    data.layers["counts"] = data.X.copy()
    hr._detect_and_normalize(data)
    assert np.isfinite(dense(data.X)).all()
    np.testing.assert_array_equal(dense(data.X)[0], [0, 0])
    np.testing.assert_array_equal(dense(data.layers["counts"]), raw)
    np.testing.assert_array_equal(data.obs["_raw_depth"], raw.sum(axis=1))


@pytest.mark.parametrize("sparse", [False, True])
@pytest.mark.parametrize("viral", [False, True])
def test_both_umap_paths_and_host_graph_receive_shifted_log(tmp_path, monkeypatch, sparse, viral):
    raw = np.array([[20, 0, 1], [1, 19, 0], [100, 0, 2], [5, 35, 0]], dtype=float)
    data = ad.AnnData(sp.csr_matrix(raw) if sparse else raw.copy())
    data.var_names = ["H1", "H2", "V"]
    data.layers["counts"] = data.X.copy()
    seen = []

    def inspect_hvg(matrix, **kwargs):
        indices = [data.var_names.get_loc(gene) for gene in matrix.var_names]
        expected = shifted_log_reference(raw[:, indices])
        np.testing.assert_allclose(dense(matrix.X), expected, rtol=1e-6)
        if "viral_counts" in matrix.obs:
            np.testing.assert_array_equal(matrix.obs["viral_counts"], raw[:, 2])
        seen.append(list(matrix.var_names))
        # Exercise the real documented no-HVG exit after checking the transformed input.
        matrix.var["highly_variable"] = False

    monkeypatch.setattr(um.sc.pp, "highly_variable_genes", inspect_hvg)
    monkeypatch.setattr(um, "config", RunConfig(output=str(tmp_path), min_counts=1, min_genes=1))
    monkeypatch.setattr(um, "identity", None)
    um.umap(data, {"V": 1} if viral else {})
    assert seen == ([["H1", "H2"], ["H1", "H2", "V"]] if viral else [["H1", "H2", "V"]])
    np.testing.assert_array_equal(dense(data.X), raw)
    np.testing.assert_array_equal(dense(data.layers["counts"]), raw)


def test_fitted_host_response_uses_shifted_features_and_records_policy(tmp_path, monkeypatch):
    from .test_hostresponse import _make_host_adata, _make_virus_adata

    host = _make_host_adata()
    virus = _make_virus_adata()
    host_file, virus_file, accessions = (
        tmp_path / name for name in ["host.h5ad", "virus.h5ad", "ids.txt"]
    )
    host.write_h5ad(host_file)
    virus.write_h5ad(virus_file)
    accessions.write_text("VIRUS_A\nVIRUS_B\n")
    inputs = {
        path: hashlib.sha256(path.read_bytes()).hexdigest() for path in [host_file, virus_file]
    }
    real = hr._detect_and_normalize
    observed = []

    def check_features(matrix):
        raw = dense(matrix.X).copy()
        real(matrix)
        np.testing.assert_allclose(dense(matrix.X), shifted_log_reference(raw), rtol=1e-6)
        observed.append(matrix.n_obs)

    monkeypatch.setattr(hr, "_detect_and_normalize", check_features)
    result = hr.run_hostresponse(
        virus_h5ad=str(virus_file),
        host_h5ad=str(host_file),
        viral_accessions_file=str(accessions),
        out_dir=str(tmp_path / "out"),
        use_hvg=False,
        seeds=[0],
        n_stab_iter=2,
        detection_threshold=1,
        cv_mode="cell",
        do_enrichment=False,
    )
    assert observed == [host.n_obs]
    policy = result["raw_input_normalization"]
    assert policy["method"] == "shifted_log"
    assert policy["target_sum"] == "scanpy_default_median_cell_total"
    assert policy["preprocessed_inputs"] == "preserved"
    metrics = pd.read_csv(tmp_path / "out/hostresponse_metrics.csv")
    assert metrics["auc_mean"].between(0, 1).all()
    assert inputs == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
