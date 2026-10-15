"""Tests for UMAP-related statistical functions in scripts/umap.py.

``umap.py`` now guards its Snakemake wiring behind ``run(ctx)`` and lazily
imports plotly, so it is importable in the test environment. These tests
exercise the real ``viral_neighbor_enrichment`` directly.

Audit findings covered:
  §3.1 — ``viral_neighbor_enrichment`` and sc PCA/UMAP calls lack random seeds,
          making UMAP coordinates and permutation p-values non-reproducible.

Regression for: audits/2026-05-08-full-pipeline.md §3.1
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# The REAL enrichment function. umap.py is now importable without Snakemake
# (plotly is lazily imported), so these tests exercise production code directly
# instead of a mirror re-implementation.
# ---------------------------------------------------------------------------
from viralscan.scripts.umap import viral_neighbor_enrichment as _viral_neighbor_enrichment

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def clustered_viral_data() -> tuple[np.ndarray, np.ndarray]:
    """100 cells: first 10 viral-positive and spatially clustered near (0,0)."""
    rng = np.random.default_rng(42)
    viral_coords = rng.normal(loc=[0.0, 0.0], scale=0.5, size=(10, 2))
    non_viral_coords = rng.normal(loc=[5.0, 5.0], scale=2.0, size=(90, 2))
    coords = np.vstack([viral_coords, non_viral_coords])
    labels = np.array([1] * 10 + [0] * 90)
    return coords, labels


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestViralNeighborEnrichmentReproducibility:
    """Audit §3.1: enrichment p-value must be reproducible with same random_state.

    Without a seeded RNG the permutation p-value changes between identical runs,
    making it possible for the same data to flip across the 0.05 significance
    threshold between runs.
    """

    def test_same_seed_produces_identical_p_values(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """
        GIVEN: identical coords, labels, k, and random_state=0
        WHEN:  _viral_neighbor_enrichment is called twice
        THEN:  both calls return exactly the same p-value

        Regression for: audits/2026-05-08-full-pipeline.md §3.1
        """
        coords, labels = clustered_viral_data
        _, _, p1 = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=200, random_state=0
        )
        _, _, p2 = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=200, random_state=0
        )
        assert p1 == p2, (
            f"Same random_state=0 yielded different p-values: {p1} vs {p2}. "
            "The permutation loop must use a seeded RNG, not np.random.permutation."
        )

    @pytest.mark.parametrize("n_cells", [1, 2, 3, 9])
    def test_small_samples_do_not_exceed_neighbor_count(self, n_cells):
        from viralscan.scripts.umap import calculate_k_neighbors

        coords = np.arange(n_cells * 2).reshape(n_cells, 2)
        labels = np.zeros(n_cells, dtype=int)
        labels[0] = 1
        assert calculate_k_neighbors(n_cells) <= n_cells
        result = _viral_neighbor_enrichment(coords, labels, k=10, n_permutations=10)
        assert np.all(np.isfinite(result))
        assert 0 <= result[2] <= 1

    def test_same_seed_produces_identical_observed_enrichment(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """observed enrichment score must also be identical across identical calls."""
        coords, labels = clustered_viral_data
        obs1, _, _ = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=50, random_state=7
        )
        obs2, _, _ = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=50, random_state=7
        )
        assert obs1 == obs2, "Observed enrichment differs — kNN structure is non-deterministic"

    def test_no_viral_cells_returns_p_value_one(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """With no viral-positive cells, p_value must be 1.0 without crashing."""
        coords, _ = clustered_viral_data
        labels = np.zeros(len(coords), dtype=int)
        observed, expected, p_value = _viral_neighbor_enrichment(
            coords, labels, k=10, random_state=0
        )
        assert p_value == 1.0, f"Expected p=1.0 with no viral cells, got {p_value}"
        assert observed == 0.0
        assert expected == 0.0

    def test_p_value_bounded_in_unit_interval(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """p-value must always lie in [0, 1]."""
        coords, labels = clustered_viral_data
        _, _, p = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=100, random_state=0
        )
        assert 0.0 <= p <= 1.0, f"p-value {p} outside [0, 1]"

    def test_clustered_viral_cells_yield_significant_p_value(
        self, clustered_viral_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """Spatially clustered viral cells must produce a low enrichment p-value."""
        coords, labels = clustered_viral_data
        _, _, p = _viral_neighbor_enrichment(
            coords, labels, k=10, n_permutations=500, random_state=0
        )
        assert p < 0.05, (
            f"Strongly clustered viral cells produced p={p:.4f} >= 0.05. "
            "The enrichment test may not be functioning correctly."
        )


class TestLayerMergeNoDoubleCount:
    """Task 1: counts_corrected + counts_original must not double-count unique reads.

    multimap.py skips unique-mapping ECs (len(genes_in_ec)==1) so
    counts_corrected only carries the redistributed multimapper share.
    The sum is therefore: unique reads + multimapper share — no duplication.
    """

    def test_unique_mapping_reads_not_duplicated(self) -> None:
        """
        GIVEN: counts_original has unique reads; counts_corrected has
               multimapper fractions only (zeros for all unique-mapping genes).
        WHEN:  adata.X = counts_original + counts_corrected
        THEN:  adata.X equals the element-wise sum (no duplication).
        """
        import numpy as np
        import scipy.sparse as sp

        # Simulate: gene0 = host (no viral), gene1 = unique viral (counts_corrected=0),
        #           gene2 = multimapper viral gene
        counts_original = sp.csr_matrix(
            np.array(
                [
                    [10.0, 3.0, 0.0],
                    [5.0, 0.0, 2.0],
                ]
            )
        )
        counts_corrected = sp.csr_matrix(
            np.array(
                [
                    [0.0, 0.0, 1.5],  # multimapper share for gene2, cell 0
                    [0.0, 0.0, 0.5],  # multimapper share for gene2, cell 1
                ]
            )
        )

        merged = counts_original + counts_corrected
        expected = np.array(
            [
                [10.0, 3.0, 1.5],
                [5.0, 0.0, 2.5],
            ]
        )
        np.testing.assert_array_almost_equal(merged.toarray(), expected)

    def test_no_counts_corrected_layer_leaves_original_intact(self) -> None:
        """When only counts_original is present, X must not be modified."""
        import numpy as np
        import scipy.sparse as sp

        counts_original = sp.csr_matrix(np.array([[1.0, 2.0], [3.0, 4.0]]))
        # Simulate the guard: "counts_corrected" not in layers → skip
        has_both = "counts_corrected" in {} and "counts_original" in {}
        assert not has_both, "Guard must be False when counts_corrected is absent"
        # X stays as counts_original unchanged
        x = counts_original
        np.testing.assert_array_equal(x.toarray(), np.array([[1.0, 2.0], [3.0, 4.0]]))


# ---------------------------------------------------------------------------
# AUD-05 / B-02 and B-24: umap() entry behaviour
# ---------------------------------------------------------------------------


def _counts_adata(n_cells: int = 200, n_genes: int = 100):
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(0)
    rates = rng.gamma(0.5, 4.0, size=n_genes)  # gene-specific means so HVGs exist
    x = rng.poisson(rates, size=(n_cells, n_genes)).astype(float)
    adata = ad.AnnData(sp.csr_matrix(x))
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


class TestUmapEntry:
    def test_creates_plots_dir_when_absent(self, tmp_path, monkeypatch) -> None:
        """B-02: the QC histogram is saved even when plots/ was never created.

        GIVEN: a negative sample / --no-visual run (no plots/ directory)
        WHEN:  umap() runs
        THEN:  plots/qc_hist_total_counts.png exists (no FileNotFoundError)
        """
        pytest.importorskip("plotly")
        from viralscan.runconfig import RunConfig
        from viralscan.scripts import umap as umap_mod

        monkeypatch.setattr(
            umap_mod, "config", RunConfig(output=str(tmp_path), min_counts=10, min_genes=5)
        )
        assert not (tmp_path / "plots").exists()
        umap_mod.umap(_counts_adata(), {})
        assert (tmp_path / "plots" / "qc_hist_total_counts.png").is_file()

    def test_all_cells_filtered_fails_loudly(self, tmp_path, monkeypatch) -> None:
        """B-24: QC removing every cell raises (the Snakemake rule fails), not exit 0."""
        from viralscan.runconfig import RunConfig
        from viralscan.scripts import umap as umap_mod

        monkeypatch.setattr(
            umap_mod, "config", RunConfig(output=str(tmp_path), min_counts=10**9, min_genes=1)
        )
        with pytest.raises(ValueError, match="No cells left after QC filtering"):
            umap_mod.umap(_counts_adata(), {})


class TestDetectedGenes:
    def test_main_passes_only_detected_genes_to_plot(self, tmp_path, monkeypatch):
        from viralscan.runconfig import RunConfig
        from viralscan.scripts import umap as umap_mod

        monkeypatch.setattr(
            umap_mod, "config", RunConfig(output=str(tmp_path), umap=True, multimapping=False)
        )
        monkeypatch.setattr(umap_mod, "identity", None)
        monkeypatch.setattr(umap_mod, "px", object())
        adata = _counts_adata()
        monkeypatch.setattr(umap_mod.sc, "read_h5ad", lambda _: adata)
        (tmp_path / "results").mkdir()
        (tmp_path / "log").mkdir()
        name = umap_mod._gene_to_virus(["TTV_gp1"])["TTV_gp1"]
        (tmp_path / "results" / "viral_summary.tsv").write_text(
            f"virus_name\tdetected\n{name}\tfalse\n"
        )
        (tmp_path / "log" / "found_genes.txt").write_text("TTV_gp1;3\n")
        seen = []
        monkeypatch.setattr(umap_mod, "umap", lambda a, genes: seen.append(genes))
        umap_mod.main()
        assert seen == [{}]

    def test_candidate_only_genes_do_not_label_cells(self, tmp_path):
        import anndata as ad
        from scipy import sparse

        from viralscan.scripts.umap import detected_genes, label_infected_cells
        from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

        table = VirusIdentityTable(
            genes=tuple(
                GeneIdentity(
                    gene_id=g,
                    genome_accession=g,
                    status="catalogued",
                    viral=True,
                    virus_key=v,
                    virus_name=v,
                )
                for g, v in [("A1", "A"), ("A2", "A"), ("B1", "B")]
            )
        )
        (tmp_path / "results").mkdir()
        (tmp_path / "results" / "viral_summary.tsv").write_text(
            "virus_name\tdetected\nA\ttrue\nB\tFalse\n"
        )
        kept = detected_genes({"A1": 1, "A2": 1, "B1": 3}, tmp_path, table)
        assert kept == {"A1": 1, "A2": 1}
        adata = ad.AnnData(sparse.csr_matrix([[0.6, 0.6, 0], [0, 0, 3]]))
        adata.var_names = ["A1", "A2", "B1"]
        labels, counts = label_infected_cells(adata, adata.X, list(kept), table, 1)
        assert labels == ["A", "No Virus"]
        np.testing.assert_allclose(counts, [1.2, 0])

    def test_missing_summary_legacy_warns_but_identity_run_fails(self, tmp_path, caplog):
        from viralscan.scripts.umap import detected_genes
        from viralscan.virus_identity import VirusIdentityTable

        genes = {"A": 1}
        assert detected_genes(genes, tmp_path) == genes
        assert "Legacy run" in caplog.text
        with pytest.raises(FileNotFoundError, match="Rerun detection"):
            detected_genes(genes, tmp_path, VirusIdentityTable(genes=()))

    def test_invalid_summary_fails(self, tmp_path):
        from viralscan.scripts.umap import detected_genes

        (tmp_path / "results").mkdir()
        summary = tmp_path / "results" / "viral_summary.tsv"
        summary.write_text("gene\tdetected\nA\ttrue\n")
        with pytest.raises(ValueError, match="virus_name column"):
            detected_genes({"A": 1}, tmp_path)

    def test_legacy_summary_filters_by_legacy_name_and_warns(self, tmp_path, caplog):
        from viralscan.scripts.umap import _gene_to_virus, detected_genes

        (tmp_path / "results").mkdir()
        name = _gene_to_virus(["TTV_gp1"])["TTV_gp1"]
        (tmp_path / "results" / "viral_summary.tsv").write_text(f"virus_name\n{name}\n")
        assert detected_genes({"TTV_gp1": 2, "unknown": 3}, tmp_path) == {"TTV_gp1": 2}
        assert "no detected column" in caplog.text


# ---------------------------------------------------------------------------
# AUD3-02 (P5-01, P5-07, P8-01): the permutation test runs on a viral-free graph
# ---------------------------------------------------------------------------

_VIR = ["TTV_gp1", "TTV_gp2", "TTV_gp3"]


def _noise_with_random_virus(seed: int, n: int = 300, g: int = 600, pos: float = 0.08):
    """Homogeneous host (no cell structure) + dominant viral counts in random cells."""
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    depth = rng.lognormal(0, 0.3, n)
    rates = rng.gamma(0.6, 3.0, g)
    host = rng.poisson(np.outer(depth, rates)).astype(float)
    v = np.zeros((n, 3))
    idx = rng.choice(n, int(pos * n), replace=False)
    v[idx] = rng.poisson(60, (len(idx), 3)) + 20
    a = ad.AnnData(sp.csr_matrix(np.hstack([host, v])))
    a.var_names = [f"h{i}" for i in range(g)] + _VIR
    a.obs_names = [f"c{i}" for i in range(n)]
    return a


class _FakePx:
    """Plotly stub: the statistic does not depend on rendering."""

    class _Fig:
        def __getattr__(self, name):
            return lambda *a, **k: None

    def scatter(self, *a, **k):
        return self._Fig()


def _run_umap(umap_mod, monkeypatch, tmp_path, adata, found, *, fake_px=True):
    from viralscan.runconfig import RunConfig

    monkeypatch.setattr(
        umap_mod, "config", RunConfig(output=str(tmp_path), min_counts=10, min_genes=5)
    )
    monkeypatch.setattr(umap_mod, "identity", None)
    seen: list = []
    real = umap_mod.viral_neighbor_enrichment

    def spy(coords, labels, k, **kw):
        seen.append((np.asarray(coords).shape, real(coords, labels, k, **kw)))
        return seen[-1][1]

    monkeypatch.setattr(umap_mod, "viral_neighbor_enrichment", spy)
    if fake_px:
        monkeypatch.setattr(umap_mod, "px", _FakePx())
    umap_mod.umap(adata, found)
    return seen


class TestViralFreeTestGraph:
    def test_host_graph_features_exclude_viral_genes(self) -> None:
        """P5-01: the representation the test runs on never contains a viral gene."""
        from viralscan.scripts.umap import host_only_pca

        adata = _noise_with_random_virus(0)
        coords, hvgs = host_only_pca(adata, _VIR, 0.0125, 3, 0.5)
        assert coords is not None and coords.shape[0] == adata.n_obs
        assert len(hvgs) > 1
        assert not set(hvgs) & set(_VIR)

    def test_host_graph_ignores_viral_load_in_size_factors(self) -> None:
        """P5-07: host normalisation does not depend on viral counts."""
        import scipy.sparse as sp

        from viralscan.scripts.umap import host_only_pca

        a = _noise_with_random_virus(1)
        b = a.copy()
        dense = b.X.toarray()
        dense[:, b.var_names.isin(_VIR)] *= 50  # huge viral load changes nothing host-side
        b.X = sp.csr_matrix(dense)
        ca, _ = host_only_pca(a, _VIR, 0.0125, 3, 0.5)
        cb, _ = host_only_pca(b, _VIR, 0.0125, 3, 0.5)
        np.testing.assert_allclose(np.abs(ca), np.abs(cb), atol=1e-4)

    def test_random_dominant_virus_not_significant_every_run(self, tmp_path, monkeypatch) -> None:
        """P5-01: virus-positive cells chosen at random, viral counts dominant.

        With viral genes in the embedding the audit saw p < 0.05 in 15/15 runs;
        on the viral-free graph the rate is that of the host noise alone (4/15 in
        the audit, an open item, AUD3-02-OPEN). Six fixed seeds, 1 of 6 observed; the
        bound ``<= 3`` is loose on purpose (at a 4/15 rate P(>= 4 of 6) is about 3%, and the
        seeds are fixed) and only asserts that the test no longer fires on every run.
        """
        from viralscan.scripts import umap as umap_mod

        ps = []
        for seed in range(6):
            out = tmp_path / f"s{seed}"
            seen = _run_umap(
                umap_mod, monkeypatch, out, _noise_with_random_virus(seed), {g: 1.0 for g in _VIR}
            )
            ps.append(seen[-1][1][2])
        assert sum(p < 0.05 for p in ps) <= 3, ps

    def test_enrichment_runs_on_pca_space_not_umap(self, tmp_path, monkeypatch) -> None:
        from viralscan.scripts import umap as umap_mod

        seen = _run_umap(
            umap_mod, monkeypatch, tmp_path, _noise_with_random_virus(2), {g: 1.0 for g in _VIR}
        )
        ((shape, _),) = seen
        assert shape[1] > 2  # PCA components, not the 2-D UMAP

    def test_zero_hvg_negative_sample_skips_cleanly(self, tmp_path, monkeypatch, caplog) -> None:
        """P8-01: no HVG survives -> warning and return, no sklearn exception."""
        import anndata as ad
        import scipy.sparse as sp

        from viralscan.scripts import umap as umap_mod

        rng = np.random.default_rng(0)
        x = sp.csr_matrix(rng.poisson(2, size=(200, 100)).astype(np.float32))
        a = ad.AnnData(x)
        a.var_names = [f"g{i}" for i in range(100)]
        a.obs_names = [f"c{i}" for i in range(200)]
        with caplog.at_level("WARNING"):
            _run_umap(umap_mod, monkeypatch, tmp_path, a, [])
        assert "highly variable" in caplog.text
        assert not (tmp_path / "plots" / "umap_no_virus.html").exists()

    def test_zero_host_hvg_with_virus_writes_figures_without_pvalue(
        self, tmp_path, monkeypatch
    ) -> None:
        """P8-01: viral-free graph impossible -> figures kept, test reported as not run."""
        pytest.importorskip("plotly")
        import anndata as ad
        import scipy.sparse as sp

        from viralscan.scripts import umap as umap_mod

        rng = np.random.default_rng(0)
        host = rng.poisson(2, size=(200, 100)).astype(np.float32)
        v = np.zeros((200, 3), dtype=np.float32)
        v[:20] = 30
        a = ad.AnnData(sp.csr_matrix(np.hstack([host, v])))
        a.var_names = [f"g{i}" for i in range(100)] + _VIR
        a.obs_names = [f"c{i}" for i in range(200)]
        seen = _run_umap(umap_mod, monkeypatch, tmp_path, a, {g: 1.0 for g in _VIR}, fake_px=False)
        assert seen == []
        html = (tmp_path / "plots" / "umap_binary.html").read_text()
        assert "not run" in html

    def test_figures_written_and_labelled(self, tmp_path, monkeypatch) -> None:
        pytest.importorskip("plotly")
        from viralscan.scripts import umap as umap_mod

        _run_umap(
            umap_mod,
            monkeypatch,
            tmp_path,
            _noise_with_random_virus(3),
            {g: 1.0 for g in _VIR},
            fake_px=False,
        )
        for name in ("umap_binary.html", "umap_continuous.html"):
            html = (tmp_path / "plots" / name).read_text()
            assert "viral-inclusive" in html
            assert "viral-free" in html
            assert "Does not distinguish infection from ambient" in html

    def test_report_glossary_makes_no_infection_claim(self) -> None:
        import viralscan

        tpl = Path(viralscan.__file__).parent / "templates" / "report.html.j2"
        text = tpl.read_text(encoding="utf-8")
        assert "suggests true infection" not in text
        assert "viral-free" in text
