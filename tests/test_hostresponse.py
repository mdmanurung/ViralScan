"""Tests for viralscan.scripts.hostresponse."""

import os
import textwrap
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from viralscan.scripts.hostresponse import (
    DEFAULT_SEEDS,
    MIN_VIRUS_CELLS,
    _balanced_split,
    _detect_and_normalize,
    _load_viral_accessions,
    _run_l2_regression,
    _run_stability_selection,
    _safe_name,
    _select_features,
    run_hostresponse,
)


# ── fixtures ──────────────────────────────────────────────────────────────────

def _make_host_adata(n_obs: int = 80, n_vars: int = 50, seed: int = 0, raw: bool = True):
    """Synthetic host gene-expression matrix."""
    rng = np.random.default_rng(seed)
    if raw:
        X = rng.negative_binomial(5, 0.3, size=(n_obs, n_vars)).astype(np.float32)
    else:
        # pre-normalized (float, relatively small values)
        X = rng.gamma(1.0, 0.5, size=(n_obs, n_vars)).astype(np.float32)
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=[f"gene_{j}" for j in range(n_vars)])
    return ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)


def _make_virus_adata(n_obs: int = 80, viruses=("VIRUS_A", "VIRUS_B"), seed: int = 1):
    """Synthetic virus count matrix; VIRUS_A is detectable, VIRUS_B is absent."""
    rng = np.random.default_rng(seed)
    counts = np.zeros((n_obs, len(viruses)), dtype=np.float32)
    # First virus: 20 positive cells
    pos_cells = rng.choice(n_obs, size=20, replace=False)
    counts[pos_cells, 0] = rng.integers(1, 10, size=20).astype(np.float32)
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=list(viruses))
    return ad.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)


@pytest.fixture()
def analysis_txt(tmp_path):
    f = tmp_path / "analysis.txt"
    f.write_text("VIRUS_A\nVIRUS_B\n")
    return str(f)


@pytest.fixture()
def host_h5ad_path(tmp_path):
    adata = _make_host_adata(raw=True)
    p = tmp_path / "host.h5ad"
    adata.write_h5ad(p)
    return str(p)


@pytest.fixture()
def virus_h5ad_path(tmp_path):
    adata = _make_virus_adata()
    p = tmp_path / "virus.h5ad"
    adata.write_h5ad(p)
    return str(p)


# ── unit tests ────────────────────────────────────────────────────────────────

class TestSafeName:
    def test_slashes(self):
        assert _safe_name("NC/001/1") == "NC_001_1"

    def test_spaces(self):
        assert _safe_name("Human herpesvirus 1") == "Human_herpesvirus_1"

    def test_colons(self):
        assert _safe_name("virus:strain") == "virus_strain"

    def test_no_change(self):
        assert _safe_name("VIRUS_A") == "VIRUS_A"


class TestLoadViralAccessions:
    def test_reads_lines(self, tmp_path):
        f = tmp_path / "analysis.txt"
        f.write_text("ACC_1\nACC_2\nACC_3\n")
        result = _load_viral_accessions(str(f))
        assert result == {"ACC_1", "ACC_2", "ACC_3"}

    def test_strips_whitespace(self, tmp_path):
        f = tmp_path / "analysis.txt"
        f.write_text("  ACC_1  \nACC_2\n")
        assert "ACC_1" in _load_viral_accessions(str(f))

    def test_skips_blank_lines(self, tmp_path):
        f = tmp_path / "analysis.txt"
        f.write_text("ACC_1\n\nACC_2\n")
        assert len(_load_viral_accessions(str(f))) == 2


class TestDetectAndNormalize:
    def test_raw_counts_get_normalized(self):
        adata = _make_host_adata(raw=True)
        _detect_and_normalize(adata)
        assert "log1p" in adata.uns
        assert "_raw_depth" in adata.obs.columns

    def test_already_normalized_not_double_normalized(self):
        adata = _make_host_adata(raw=False)
        adata.uns["log1p"] = {}
        X_before = adata.X.toarray().copy()
        _detect_and_normalize(adata)
        np.testing.assert_array_equal(adata.X.toarray(), X_before)

    def test_raw_depth_stored_before_normalization(self):
        adata = _make_host_adata(raw=True)
        expected_depth = np.asarray(adata.X.sum(axis=1)).flatten()
        _detect_and_normalize(adata)
        np.testing.assert_allclose(adata.obs["_raw_depth"].values, expected_depth)

    def test_raw_depth_is_not_normalized_values(self):
        adata = _make_host_adata(raw=True)
        _detect_and_normalize(adata)
        # After normalization, row sums ≈ 1e4; raw depth is the original integer sum
        depth = adata.obs["_raw_depth"].values
        assert depth.max() > 100  # integer counts, not normalized


class TestSelectFeatures:
    def test_hvg_returns_subset(self):
        adata = _make_host_adata(n_obs=100, n_vars=200)
        _detect_and_normalize(adata)
        X, names = _select_features(adata, use_hvg=True)
        assert X.shape[0] == 100
        assert X.shape[1] == len(names)
        assert X.shape[1] <= 200

    def test_all_genes_returns_full_matrix(self):
        adata = _make_host_adata(n_obs=80, n_vars=50)
        _detect_and_normalize(adata)
        X, names = _select_features(adata, use_hvg=False)
        assert X.shape == (80, 50)
        assert len(names) == 50

    def test_returns_dense_float32(self):
        adata = _make_host_adata()
        _detect_and_normalize(adata)
        X, _ = _select_features(adata, use_hvg=False)
        assert isinstance(X, np.ndarray)
        assert X.dtype == np.float32


class TestBalancedSplit:
    def setup_method(self):
        rng = np.random.default_rng(0)
        n = 80
        self.X = rng.standard_normal((n, 20)).astype(np.float32)
        self.depth = rng.integers(500, 5000, size=n).astype(float)
        self.pos_idx = np.arange(20)
        self.neg_idx = np.arange(20, n)

    def test_returns_four_parts(self):
        result = _balanced_split(self.pos_idx, self.neg_idx, self.depth, self.X, seed=0)
        assert result is not None
        X_train, y_train, X_test_pos, X_test_neg = result
        assert X_train.shape[1] == 20

    def test_y_train_is_balanced(self):
        _, y_train, _, _ = _balanced_split(self.pos_idx, self.neg_idx, self.depth, self.X, seed=0)
        assert y_train.sum() == (y_train == 0).sum()

    def test_returns_none_when_too_few_cells(self):
        result = _balanced_split(np.arange(2), np.arange(2, 4), self.depth, self.X, seed=0)
        assert result is None


class TestL2Regression:
    def setup_method(self):
        rng = np.random.default_rng(42)
        n = 100
        self.n_features = 30
        self.X = rng.standard_normal((n, self.n_features)).astype(np.float32)
        # First 5 features strongly predictive of virus presence
        virus_score = self.X[:, :5].sum(axis=1)
        self.virus_presence = virus_score > np.percentile(virus_score, 70)
        self.depth = rng.integers(500, 5000, size=n).astype(float)
        self.feature_names = [f"gene_{i}" for i in range(self.n_features)]

    def test_returns_dataframe_and_metrics(self):
        weights_df, metrics = _run_l2_regression(
            self.X, self.virus_presence, self.depth, DEFAULT_SEEDS[:2], self.feature_names
        )
        assert weights_df is not None
        assert metrics is not None
        assert set(weights_df.columns) >= {"gene", "weight_mean", "weight_sd"}
        assert len(weights_df) == self.n_features

    def test_metrics_contain_expected_keys(self):
        _, metrics = _run_l2_regression(
            self.X, self.virus_presence, self.depth, DEFAULT_SEEDS[:2], self.feature_names
        )
        assert "sensitivity" in metrics
        assert "specificity" in metrics
        assert "balanced_acc" in metrics

    def test_returns_none_below_min_cells(self):
        virus_presence = np.zeros(len(self.virus_presence), dtype=bool)
        virus_presence[:MIN_VIRUS_CELLS - 1] = True
        weights_df, metrics = _run_l2_regression(
            self.X, virus_presence, self.depth, DEFAULT_SEEDS[:2], self.feature_names
        )
        assert weights_df is None
        assert metrics is None

    def test_auc_in_range(self):
        _, metrics = _run_l2_regression(
            self.X, self.virus_presence, self.depth, DEFAULT_SEEDS[:2], self.feature_names
        )
        if "auc" in metrics:
            assert 0.0 <= metrics["auc"]["mean"] <= 1.0


class TestStabilitySelection:
    def setup_method(self):
        rng = np.random.default_rng(7)
        n = 100
        self.n_features = 20
        self.X = rng.standard_normal((n, self.n_features)).astype(np.float32)
        score = self.X[:, 0] + self.X[:, 1]
        self.virus_presence = score > np.percentile(score, 75)

    def test_returns_probabilities_in_unit_interval(self):
        probs = _run_stability_selection(self.X, self.virus_presence, n_iter=20, seed=0)
        assert probs.shape == (self.n_features,)
        assert probs.min() >= 0.0
        assert probs.max() <= 1.0

    def test_strong_predictors_have_higher_prob(self):
        # Features 0 and 1 are predictive; feature 19 should have lower selection probability
        probs = _run_stability_selection(self.X, self.virus_presence, n_iter=30, seed=42)
        assert probs[0] >= probs[-1], "predictive feature should have higher stability probability"


class TestRunHostresponse:
    """End-to-end integration tests over synthetic data."""

    def _setup_files(self, tmp_path):
        n_cells = 80
        n_genes = 50

        host_adata = _make_host_adata(n_obs=n_cells, n_vars=n_genes, raw=True)
        virus_adata = _make_virus_adata(n_obs=n_cells)

        host_h5ad = tmp_path / "host.h5ad"
        virus_h5ad = tmp_path / "virus.h5ad"
        host_adata.write_h5ad(host_h5ad)
        virus_adata.write_h5ad(virus_h5ad)

        analysis_txt = tmp_path / "analysis.txt"
        analysis_txt.write_text("VIRUS_A\nVIRUS_B\n")

        out_dir = tmp_path / "hostresponse"
        return str(virus_h5ad), str(host_h5ad), str(analysis_txt), str(out_dir)

    def test_creates_output_files(self, tmp_path):
        virus_h5ad, host_h5ad, analysis_txt, out_dir = self._setup_files(tmp_path)
        run_hostresponse(
            virus_h5ad=virus_h5ad,
            host_h5ad=host_h5ad,
            viral_accessions_file=analysis_txt,
            out_dir=out_dir,
            use_hvg=False,  # avoid HVG selection on tiny synthetic data
            seeds=DEFAULT_SEEDS[:2],
            n_stab_iter=10,
            stab_min_prob=0.5,
            detection_threshold=1,
        )
        assert Path(out_dir).exists()
        # VIRUS_A has 20 positive cells → should produce output
        assert (Path(out_dir) / "VIRUS_A_gene_weights.csv").exists()
        assert (Path(out_dir) / "VIRUS_A_stability.csv").exists()
        # Summary metrics file
        assert (Path(out_dir) / "hostresponse_metrics.csv").exists()

    def test_gene_weights_csv_has_expected_columns(self, tmp_path):
        virus_h5ad, host_h5ad, analysis_txt, out_dir = self._setup_files(tmp_path)
        run_hostresponse(
            virus_h5ad=virus_h5ad,
            host_h5ad=host_h5ad,
            viral_accessions_file=analysis_txt,
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:2],
            n_stab_iter=10,
            detection_threshold=1,
        )
        df = pd.read_csv(Path(out_dir) / "VIRUS_A_gene_weights.csv")
        assert {"gene", "weight_mean", "weight_sd", "virus"}.issubset(df.columns)

    def test_stability_csv_has_stable_column(self, tmp_path):
        virus_h5ad, host_h5ad, analysis_txt, out_dir = self._setup_files(tmp_path)
        run_hostresponse(
            virus_h5ad=virus_h5ad,
            host_h5ad=host_h5ad,
            viral_accessions_file=analysis_txt,
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:2],
            n_stab_iter=10,
            stab_min_prob=0.5,
            detection_threshold=1,
        )
        df = pd.read_csv(Path(out_dir) / "VIRUS_A_stability.csv")
        assert "stable" in df.columns
        assert "stab_prob" in df.columns

    def test_skips_virus_below_min_cells(self, tmp_path):
        virus_h5ad, host_h5ad, analysis_txt, out_dir = self._setup_files(tmp_path)
        run_hostresponse(
            virus_h5ad=virus_h5ad,
            host_h5ad=host_h5ad,
            viral_accessions_file=analysis_txt,
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:2],
            n_stab_iter=5,
            detection_threshold=1,
        )
        # VIRUS_B has all-zero counts → no output files
        assert not (Path(out_dir) / "VIRUS_B_gene_weights.csv").exists()

    def test_no_shared_barcodes_returns_early(self, tmp_path):
        n_cells = 80
        host_adata = _make_host_adata(n_obs=n_cells, raw=True)
        # Give host completely different barcodes
        host_adata.obs.index = [f"OTHER_{i}" for i in range(n_cells)]
        virus_adata = _make_virus_adata(n_obs=n_cells)

        host_h5ad = tmp_path / "host.h5ad"
        virus_h5ad_p = tmp_path / "virus.h5ad"
        host_adata.write_h5ad(host_h5ad)
        virus_adata.write_h5ad(virus_h5ad_p)
        analysis_txt = tmp_path / "analysis.txt"
        analysis_txt.write_text("VIRUS_A\n")
        out_dir = str(tmp_path / "out")

        run_hostresponse(
            virus_h5ad=str(virus_h5ad_p),
            host_h5ad=str(host_h5ad),
            viral_accessions_file=str(analysis_txt),
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:2],
            n_stab_iter=5,
        )
        # Should not produce any output files (returned early)
        assert not (Path(out_dir) / "VIRUS_A_gene_weights.csv").exists()

    def test_unknown_viral_accessions_returns_early(self, tmp_path):
        virus_adata = _make_virus_adata()
        virus_h5ad_p = tmp_path / "virus.h5ad"
        virus_adata.write_h5ad(virus_h5ad_p)

        host_adata = _make_host_adata()
        host_h5ad_p = tmp_path / "host.h5ad"
        host_adata.write_h5ad(host_h5ad_p)

        analysis_txt = tmp_path / "analysis.txt"
        analysis_txt.write_text("NOT_A_REAL_VIRUS\n")
        out_dir = str(tmp_path / "out")

        run_hostresponse(
            virus_h5ad=str(virus_h5ad_p),
            host_h5ad=str(host_h5ad_p),
            viral_accessions_file=str(analysis_txt),
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:1],
            n_stab_iter=5,
        )
        assert not (Path(out_dir) / "NOT_A_REAL_VIRUS_gene_weights.csv").exists()
