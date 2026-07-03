"""Tests for viralscan.scripts.hostresponse."""

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
    _depth_alone_auc,
    _depth_match_indices,
    _detect_and_normalize,
    _e_value,
    _load_viral_accessions,
    _mt_gene_mask,
    _per_gene_evalues,
    _run_l2_regression,
    _run_stability_selection,
    _safe_name,
    _select_features,
    _virus_presence_label,
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
        virus_presence[: MIN_VIRUS_CELLS - 1] = True
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


class TestHostresponsePlantedSignal:
    """Sanity check: planted predictive genes should rank highly in stability selection."""

    def test_planted_genes_rank_high_in_stability(self, tmp_path):
        rng = np.random.default_rng(99)
        n_obs = 200
        n_vars = 100
        n_planted = 5
        n_pos = 40

        X = rng.negative_binomial(5, 0.3, size=(n_obs, n_vars)).astype(np.float32)
        pos_cells = rng.choice(n_obs, size=n_pos, replace=False)
        X[np.ix_(pos_cells, np.arange(n_planted))] *= 10.0

        obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
        var = pd.DataFrame(
            index=[f"planted_{j}" if j < n_planted else f"noise_{j}" for j in range(n_vars)]
        )
        host_adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)

        virus_counts = np.zeros((n_obs, 1), dtype=np.float32)
        virus_counts[pos_cells, 0] = rng.integers(1, 20, size=n_pos).astype(np.float32)
        virus_adata = ad.AnnData(
            X=sp.csr_matrix(virus_counts),
            obs=obs.copy(),
            var=pd.DataFrame(index=["VIRUS_A"]),
        )

        host_h5ad = tmp_path / "host.h5ad"
        virus_h5ad = tmp_path / "virus.h5ad"
        host_adata.write_h5ad(host_h5ad)
        virus_adata.write_h5ad(virus_h5ad)

        analysis_txt = tmp_path / "analysis.txt"
        analysis_txt.write_text("VIRUS_A\n")
        out_dir = tmp_path / "out"

        run_hostresponse(
            virus_h5ad=str(virus_h5ad),
            host_h5ad=str(host_h5ad),
            viral_accessions_file=str(analysis_txt),
            out_dir=str(out_dir),
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:2],
            n_stab_iter=30,
            stab_min_prob=0.5,
            detection_threshold=1,
        )

        stab_csv = out_dir / "VIRUS_A_stability.csv"
        assert stab_csv.exists()
        stab_df = pd.read_csv(stab_csv).sort_values("stab_prob", ascending=False)
        top10_genes = set(stab_df.head(10)["gene"])
        planted_genes = {f"planted_{j}" for j in range(n_planted)}
        n_recovered = len(top10_genes & planted_genes)
        assert n_recovered >= 3, (
            f"Expected ≥3 planted genes in top-10 by stability probability; "
            f"got {n_recovered}. Top-10: {top10_genes}"
        )


# ── SH1.1 depth-confound diagnostics ──────────────────────────────────────────


class TestEValue:
    """Ding & VanderWeele E-value formula."""

    def test_or_one_is_evalue_one(self):
        assert _e_value(1.0) == pytest.approx(1.0)

    def test_symmetric_in_inverse(self):
        # E-value depends on the risk ratio magnitude, not its direction.
        assert _e_value(3.0) == pytest.approx(_e_value(1.0 / 3.0))

    def test_larger_or_gives_larger_evalue(self):
        assert _e_value(5.0) > _e_value(2.0) > _e_value(1.2)

    def test_invalid_or_is_nan(self):
        assert np.isnan(_e_value(0.0))
        assert np.isnan(_e_value(-1.0))
        assert np.isnan(_e_value(float("nan")))


class TestDepthDiagnostics:
    """Depth-alone baseline and per-gene E-values (SH1.1)."""

    def test_depth_alone_auc_detects_depth_signal(self):
        # Construct a case where depth IS the label: positives are the
        # high-depth cells. Depth-alone AUC should be well above chance.
        rng = np.random.default_rng(0)
        n = 200
        depth = rng.uniform(500, 5000, size=n)
        y = (depth >= np.median(depth)).astype(int)
        pos_idx = np.where(y == 1)[0]
        neg_idx = np.where(y == 0)[0]
        res = _depth_alone_auc(pos_idx, neg_idx, depth, DEFAULT_SEEDS)
        assert res is not None
        assert res["mean"] > 0.8

    def test_synthetic_depth_only_gene_is_not_robust(self):
        # A gene that is only a NOISY proxy of depth (its apparent association
        # with the depth-driven label is entirely mediated by depth) should,
        # after adjusting for depth, collapse to an odds ratio ~1 and a small
        # E-value — the guard against silently re-confounding.
        rng = np.random.default_rng(1)
        n = 600
        depth = rng.uniform(500, 5000, size=n)
        y = (depth >= np.median(depth)).astype(int)
        zlogd = (np.log1p(depth) - np.log1p(depth).mean()) / np.log1p(depth).std()
        gene = zlogd + rng.normal(0, 1.0, size=n)  # correlated with depth, not collinear
        df = _per_gene_evalues(gene.reshape(-1, 1), ["depth_proxy"], y, depth)
        assert abs(df.loc[0, "adj_OR"] - 1.0) < 0.5
        assert df.loc[0, "E_value"] < 1.8

    def test_independent_gene_signal_survives_adjustment(self):
        # A gene carrying label information NOT explained by depth should keep a
        # non-trivial odds ratio and an E-value above 1 after depth adjustment.
        rng = np.random.default_rng(2)
        n = 400
        depth = rng.uniform(500, 5000, size=n)
        y = rng.integers(0, 2, size=n)  # label independent of depth
        gene = y * 2.0 + rng.normal(0, 1.0, size=n)  # strong, depth-free signal
        df = _per_gene_evalues(gene.reshape(-1, 1), ["real_gene"], y, depth)
        assert df.loc[0, "adj_OR"] > 1.5
        assert df.loc[0, "E_value"] > 2.0


class TestDepthDiagnosticsIntegration:
    """The end-to-end run emits the depth diagnostics."""

    def _setup_files(self, tmp_path):
        host_adata = _make_host_adata(n_obs=80, n_vars=50, raw=True)
        virus_adata = _make_virus_adata(n_obs=80)
        host_h5ad = tmp_path / "host.h5ad"
        virus_h5ad = tmp_path / "virus.h5ad"
        host_adata.write_h5ad(host_h5ad)
        virus_adata.write_h5ad(virus_h5ad)
        analysis_txt = tmp_path / "analysis.txt"
        analysis_txt.write_text("VIRUS_A\nVIRUS_B\n")
        return str(virus_h5ad), str(host_h5ad), str(analysis_txt), str(tmp_path / "hostresponse")

    def test_metrics_has_depth_alone_columns(self, tmp_path):
        virus_h5ad, host_h5ad, analysis_txt, out_dir = self._setup_files(tmp_path)
        run_hostresponse(
            virus_h5ad=virus_h5ad,
            host_h5ad=host_h5ad,
            viral_accessions_file=analysis_txt,
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:3],
            n_stab_iter=10,
            stab_min_prob=0.3,
            detection_threshold=1,
        )
        df = pd.read_csv(Path(out_dir) / "hostresponse_metrics.csv")
        assert "depth_alone_auc_mean" in df.columns
        assert "n_genes_evalue_ge2" in df.columns
        assert "n_stable_genes" in df.columns


# ── SH1.2 depth-independent label + depth-matched design ──────────────────────


class TestVirusPresenceLabel:
    """Label definitions: raw vs depth-normalized (cpm/fraction)."""

    def test_raw_is_threshold(self):
        counts = np.array([0, 5, 10, 20], dtype=float)
        depth = np.array([1000, 1000, 1000, 1000], dtype=float)
        vp = _virus_presence_label(counts, depth, detection_threshold=10, label="raw")
        assert vp.tolist() == [False, False, True, True]

    def test_cpm_preserves_prevalence(self):
        rng = np.random.default_rng(0)
        counts = rng.integers(0, 40, size=200).astype(float)
        depth = rng.uniform(500, 5000, size=200)
        n_raw = int((counts >= 10).sum())
        vp = _virus_presence_label(counts, depth, detection_threshold=10, label="cpm")
        # Same number of positives as the raw label (prevalence-matched).
        assert int(vp.sum()) == n_raw

    def test_cpm_and_fraction_select_same_cells(self):
        rng = np.random.default_rng(1)
        counts = rng.integers(0, 40, size=150).astype(float)
        depth = rng.uniform(500, 5000, size=150)
        vp_cpm = _virus_presence_label(counts, depth, 10, "cpm")
        vp_frac = _virus_presence_label(counts, depth, 10, "fraction")
        np.testing.assert_array_equal(vp_cpm, vp_frac)

    def test_cpm_differs_from_raw_when_depth_varies(self):
        # A high-count-but-deep cell (low CPM) and a low-count-but-shallow cell
        # (high CPM) should swap positivity between raw and cpm labels.
        counts = np.array([12, 8, 12, 8], dtype=float)
        depth = np.array([100_000, 1000, 100_000, 1000], dtype=float)
        raw = _virus_presence_label(counts, depth, 10, "raw")
        cpm = _virus_presence_label(counts, depth, 10, "cpm")
        assert not np.array_equal(raw, cpm)

    def test_invalid_label_raises(self):
        with pytest.raises(ValueError, match="label must be one of"):
            _virus_presence_label(np.array([1.0]), np.array([1.0]), 1, "bogus")


class TestDepthMatchIndices:
    """Coarsened-exact depth matching equalizes class depth."""

    def test_matched_classes_have_equal_counts(self):
        rng = np.random.default_rng(0)
        depth = rng.uniform(500, 5000, size=400)
        # Label is depth-driven: positives are deep cells.
        vp = depth >= np.median(depth)
        idx = _depth_match_indices(vp, depth, n_bins=10)
        assert vp[idx].sum() == (~vp[idx]).sum()

    def test_matching_removes_depth_gap(self):
        # Label is depth-BIASED but not a step function of depth, so bins contain
        # both classes and matching can equalize them.
        rng = np.random.default_rng(2)
        depth = rng.uniform(500, 5000, size=1500)
        p = (depth - depth.min()) / (depth.max() - depth.min())  # deeper → more likely +
        vp = rng.random(len(depth)) < p
        gap_before = abs(np.median(depth[vp]) - np.median(depth[~vp]))
        idx = _depth_match_indices(vp, depth, n_bins=20)
        matched_pos = depth[idx][vp[idx]]
        matched_neg = depth[idx][~vp[idx]]
        gap_after = abs(np.median(matched_pos) - np.median(matched_neg))
        assert gap_after < gap_before


class TestLabelDepthMatchIntegration:
    """End-to-end: label and depth_match options run and are recorded."""

    def _setup_files(self, tmp_path):
        host_adata = _make_host_adata(n_obs=120, n_vars=50, raw=True)
        virus_adata = _make_virus_adata(n_obs=120)
        host_h5ad = tmp_path / "host.h5ad"
        virus_h5ad = tmp_path / "virus.h5ad"
        host_adata.write_h5ad(host_h5ad)
        virus_adata.write_h5ad(virus_h5ad)
        analysis_txt = tmp_path / "analysis.txt"
        analysis_txt.write_text("VIRUS_A\nVIRUS_B\n")
        return str(virus_h5ad), str(host_h5ad), str(analysis_txt), str(tmp_path / "hr")

    def test_cpm_label_recorded_in_metrics(self, tmp_path):
        virus_h5ad, host_h5ad, analysis_txt, out_dir = self._setup_files(tmp_path)
        run_hostresponse(
            virus_h5ad=virus_h5ad,
            host_h5ad=host_h5ad,
            viral_accessions_file=analysis_txt,
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:3],
            n_stab_iter=10,
            stab_min_prob=0.3,
            detection_threshold=1,
            label="cpm",
        )
        df = pd.read_csv(Path(out_dir) / "hostresponse_metrics.csv")
        assert (df["label"] == "cpm").all()

    def test_depth_match_recorded_in_metrics(self, tmp_path):
        virus_h5ad, host_h5ad, analysis_txt, out_dir = self._setup_files(tmp_path)
        run_hostresponse(
            virus_h5ad=virus_h5ad,
            host_h5ad=host_h5ad,
            viral_accessions_file=analysis_txt,
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:3],
            n_stab_iter=10,
            stab_min_prob=0.3,
            detection_threshold=1,
            depth_match=True,
        )
        df = pd.read_csv(Path(out_dir) / "hostresponse_metrics.csv")
        assert bool(df["depth_matched"].iloc[0]) is True

    def test_invalid_label_raises_at_entry(self, tmp_path):
        virus_h5ad, host_h5ad, analysis_txt, out_dir = self._setup_files(tmp_path)
        with pytest.raises(ValueError, match="label must be one of"):
            run_hostresponse(
                virus_h5ad=virus_h5ad,
                host_h5ad=host_h5ad,
                viral_accessions_file=analysis_txt,
                out_dir=out_dir,
                label="nonsense",
            )


# ── SH1.3 %mito control ───────────────────────────────────────────────────────


class TestMitoGeneMask:
    def test_matches_mt_symbol_prefix(self):
        mask = _mt_gene_mask(["ACTB", "MT-CO1", "mt-nd2", "GAPDH"])
        assert mask.tolist() == [False, True, True, False]

    def test_matches_ensembl_mt_id_with_and_without_version(self):
        mask = _mt_gene_mask(["ENSG00000198804", "ENSG00000198804.2", "ENSG00000075624"])
        assert mask.tolist() == [True, True, False]


class TestMitoControlEValues:
    def test_mito_covariate_removes_mito_driven_signal(self):
        # Label is driven by %mito; a gene that is only a noisy proxy of %mito
        # should collapse to OR~1 once %mito is a covariate.
        rng = np.random.default_rng(0)
        n = 600
        depth = rng.uniform(500, 5000, size=n)
        pct_mito = rng.uniform(1, 30, size=n)
        y = (pct_mito >= np.median(pct_mito)).astype(int)
        zpm = (pct_mito - pct_mito.mean()) / pct_mito.std()
        gene = zpm + rng.normal(0, 1.0, size=n)  # proxy of %mito
        without = _per_gene_evalues(gene.reshape(-1, 1), ["g"], y, depth)
        with_mito = _per_gene_evalues(gene.reshape(-1, 1), ["g"], y, depth, pct_mito=pct_mito)
        # Adjusting for %mito pulls the odds ratio closer to 1.
        assert abs(with_mito.loc[0, "adj_OR"] - 1.0) < abs(without.loc[0, "adj_OR"] - 1.0)

    def test_leave_one_out_avoids_self_suppression(self):
        # An MT gene tested against a %mito covariate that INCLUDES it is trivially
        # suppressed (circularity). Leave-one-out should give a less-suppressed,
        # finite odds ratio than naive inclusion.
        rng = np.random.default_rng(1)
        n = 500
        depth = rng.uniform(500, 5000, size=n)
        mt_gene_raw = rng.integers(0, 50, size=n).astype(float)
        other_mt = rng.integers(0, 200, size=n).astype(float)
        raw_total = depth
        mt_total = mt_gene_raw + other_mt
        pct_mito = mt_total / np.maximum(raw_total, 1) * 100
        y = (mt_gene_raw >= np.median(mt_gene_raw)).astype(int)
        gene = np.log1p(mt_gene_raw)
        naive = _per_gene_evalues(gene.reshape(-1, 1), ["MT-XX"], y, depth, pct_mito=pct_mito)
        loo = _per_gene_evalues(
            gene.reshape(-1, 1),
            ["MT-XX"],
            y,
            depth,
            pct_mito=pct_mito,
            mt_total_counts=mt_total,
            raw_total=raw_total,
            mt_self_counts={"MT-XX": mt_gene_raw},
        )
        assert np.isfinite(loo.loc[0, "adj_OR"])
        assert loo.loc[0, "adj_OR"] >= naive.loc[0, "adj_OR"]


class TestMitoControlIntegration:
    def _setup_files(self, tmp_path, mt: bool):
        n_obs, n_vars = 120, 40
        host = _make_host_adata(n_obs=n_obs, n_vars=n_vars, raw=True)
        if mt:
            # Rename a few genes to mitochondrial symbols.
            names = list(host.var_names)
            for i, sym in enumerate(["MT-CO1", "MT-ND2", "MT-CYB"]):
                names[i] = sym
            host.var_names = names
        virus = _make_virus_adata(n_obs=n_obs)
        host_p, virus_p = tmp_path / "host.h5ad", tmp_path / "virus.h5ad"
        host.write_h5ad(host_p)
        virus.write_h5ad(virus_p)
        atxt = tmp_path / "analysis.txt"
        atxt.write_text("VIRUS_A\nVIRUS_B\n")
        return str(virus_p), str(host_p), str(atxt), str(tmp_path / "hr")

    def test_mito_controlled_true_when_mt_genes_present(self, tmp_path):
        virus_h5ad, host_h5ad, atxt, out_dir = self._setup_files(tmp_path, mt=True)
        run_hostresponse(
            virus_h5ad=virus_h5ad,
            host_h5ad=host_h5ad,
            viral_accessions_file=atxt,
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:3],
            n_stab_iter=10,
            stab_min_prob=0.3,
            detection_threshold=1,
            control_mito=True,
        )
        df = pd.read_csv(Path(out_dir) / "hostresponse_metrics.csv")
        assert "mito_controlled" in df.columns
        assert bool(df["mito_controlled"].iloc[0]) is True

    def test_mito_controlled_false_when_no_mt_genes(self, tmp_path):
        virus_h5ad, host_h5ad, atxt, out_dir = self._setup_files(tmp_path, mt=False)
        run_hostresponse(
            virus_h5ad=virus_h5ad,
            host_h5ad=host_h5ad,
            viral_accessions_file=atxt,
            out_dir=out_dir,
            use_hvg=False,
            seeds=DEFAULT_SEEDS[:3],
            n_stab_iter=10,
            stab_min_prob=0.3,
            detection_threshold=1,
            control_mito=True,
        )
        df = pd.read_csv(Path(out_dir) / "hostresponse_metrics.csv")
        assert bool(df["mito_controlled"].iloc[0]) is False
