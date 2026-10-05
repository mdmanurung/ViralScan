"""Tests for positive-control capture measurement and negative certification.

A positive control is the only thing in a run that can measure the k-mer capture
term. Depth is measurable without one, and it was ample in every real run in
this repo (see tests/test_sensitivity.py::TestClassifyLod), which is exactly why
depth alone cannot be allowed to certify a negative: the covid samples called
SARS-CoV-2 = 0 at 21.6M quantified molecules.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse

from viralscan.runconfig import RunConfig
from viralscan.scripts import detection as D

SPIKEIN = "SPIKEIN_planted"
DEPTH_PER_CELL = 50


def _adata(spikein_molecules: float = 0.0, n_cells: int = 400, seed: int = 1) -> sc.AnnData:
    rng = np.random.default_rng(seed)
    names = (
        [f"H{i}" for i in range(18)]
        + [SPIKEIN, "EPSTEIN_HHV4_EBNA-2", "NORWALK_gp1"]
        + [f"G{i}" for i in range(19)]
    )
    x = rng.poisson(DEPTH_PER_CELL, size=(n_cells, 40)).astype(np.float32)
    x[:, 18] = 0.0  # the spike-in column starts empty
    if spikein_molecules:
        x[0, 18] = spikein_molecules
    return sc.AnnData(
        sparse.csr_matrix(x),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=names),
    )


def _config(**kwargs) -> RunConfig:
    base = {"positive_control_gene": SPIKEIN, "positive_control_expected_molecules": 100.0}
    base.update(kwargs)
    return RunConfig(**base)


class TestMeasurePositiveControl:
    def test_no_control_configured(self) -> None:
        capture, detail = D.measure_positive_control(_adata(), RunConfig())
        assert capture is None
        assert detail["status"] == "not-configured"

    def test_gene_only_without_expected_count_is_not_configured(self) -> None:
        cfg = RunConfig(positive_control_gene=SPIKEIN)
        capture, detail = D.measure_positive_control(_adata(), cfg)
        assert capture is None
        assert detail["status"] == "not-configured"

    def test_half_recovery_measures_capture_one_half(self) -> None:
        capture, detail = D.measure_positive_control(_adata(50.0), _config())
        assert capture == pytest.approx(0.5)
        assert detail["status"] == "measured"
        assert detail["observed_molecules"] == pytest.approx(50.0)
        assert detail["expected_molecules"] == pytest.approx(100.0)

    def test_measured_capture_implies_a_divergence(self) -> None:
        _, detail = D.measure_positive_control(_adata(50.0), _config())
        # capture 0.5 sits between the 10% (0.90) and 15% (0.32) curve points.
        assert 0.10 < detail["implied_divergence"] < 0.15

    def test_zero_recovery_fails_and_withholds_capture(self) -> None:
        """The load-bearing case: the control is invisible, so nothing is certifiable."""
        capture, detail = D.measure_positive_control(_adata(0.0), _config())
        assert capture is None
        assert detail["status"] == "failed"
        assert "uninterpretable" in detail["detail"]

    def test_over_recovery_is_flagged_not_silently_clamped(self) -> None:
        capture, detail = D.measure_positive_control(_adata(118.0), _config())
        assert detail["status"] == "over-recovered"
        assert "not spike-in-specific" in detail["detail"]
        # Reported as no measurable loss, which is the optimistic direction and
        # is stated as such rather than hidden.
        assert capture == 1.0

    def test_gene_absent_from_reference_is_a_config_error(self) -> None:
        cfg = RunConfig(
            positive_control_gene="NOT_IN_INDEX_gp1",
            positive_control_expected_molecules=50.0,
        )
        capture, detail = D.measure_positive_control(_adata(), cfg)
        assert capture is None
        assert detail["status"] == "gene-not-in-reference"
        assert "not a column of the count matrix" in detail["detail"]

    def test_full_recovery_is_valid(self) -> None:
        capture, detail = D.measure_positive_control(_adata(100.0), _config())
        assert capture == pytest.approx(1.0)
        assert detail["status"] == "measured"

    def test_implied_divergence_is_none_when_no_loss_measurable(self) -> None:
        assert D._implied_divergence(1.0) is None
        assert D._implied_divergence(0.999) is None

    def test_implied_divergence_is_none_below_the_identifiable_floor(self) -> None:
        """Capture under ~1e-8 cannot be attributed to divergence in float64."""
        from viralscan.sensitivity import fragment_capture

        floor = fragment_capture(D.MAX_IDENTIFIABLE_DIVERGENCE, read_length=90)
        assert D._implied_divergence(floor) is None
        assert D._implied_divergence(0.0) is None
        assert D._implied_divergence(-1.0) is None

    def test_implied_divergence_decreases_as_capture_increases(self) -> None:
        """More capture means less divergence, so the inversion is monotone down."""
        prev = 1.1
        for capture in (1e-6, 0.001, 0.008, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95):
            d = D._implied_divergence(capture)
            assert d is not None, f"capture={capture} should be identifiable"
            assert d < prev, f"capture={capture} -> {d}, previous {prev}"
            prev = d

    def test_implied_divergence_recovers_known_curve_points(self) -> None:
        from viralscan.sensitivity import fragment_capture

        for d in (0.10, 0.15, 0.20, 0.25, 0.30, 0.40):
            capture = fragment_capture(d, read_length=90)
            assert D._implied_divergence(capture) == pytest.approx(d, abs=1e-3)


class TestCertificationFlip:
    """The same negative, with and without a working control."""

    STATS = {"Betatorquevirus": {"viral_molecules_total_est": 0.0}}

    def test_without_control_a_negative_is_never_certifiable(self) -> None:
        adata = _adata()
        capture, _ = D.measure_positive_control(adata, RunConfig())
        row = D.build_sensitivity_table(adata, self.STATS, RunConfig(), capture=capture).iloc[0]
        assert row["observed_molecules"] == 0
        assert row["capture_measured"] is False or row["capture_measured"] == False  # noqa: E712
        assert not row["informative_negative"]

    def test_failed_control_leaves_the_negative_uncertifiable(self) -> None:
        adata = _adata(0.0)  # planted 100, recovered 0
        cfg = _config()
        capture, _ = D.measure_positive_control(adata, cfg)
        row = D.build_sensitivity_table(adata, self.STATS, cfg, capture=capture).iloc[0]
        assert row["capture_measured"] == False  # noqa: E712
        assert not row["informative_negative"]
        assert row["depth_sufficient"] in (True, False)  # depth is reported either way

    def test_working_control_at_good_depth_certifies_the_negative(self) -> None:
        adata = _adata(50.0, n_cells=4000)
        cfg = _config()
        capture, _ = D.measure_positive_control(adata, cfg)
        assert capture == pytest.approx(0.5)
        row = D.build_sensitivity_table(adata, self.STATS, cfg, capture=capture).iloc[0]
        assert row["capture_measured"]
        assert row["depth_sufficient"]
        assert row["informative_negative"]

    def test_measured_capture_never_tightens_below_the_depth_floor(self) -> None:
        """Over-recovery must not produce a sub-floor LOD."""
        adata = _adata(118.0)
        cfg = _config()
        capture, _ = D.measure_positive_control(adata, cfg)
        row = D.build_sensitivity_table(adata, self.STATS, cfg, capture=capture).iloc[0]
        floor = D.build_sensitivity_table(adata, self.STATS, RunConfig()).iloc[0]
        assert row["lod95_per_10k"] == pytest.approx(floor["lod95_per_10k"])


class TestZeroRows:
    """MECH-B: every indexed virus gets a sensitivity row, detected or not."""

    def test_undetected_indexed_virus_gets_a_zero_row(self) -> None:
        stats = {"Betatorquevirus": {"viral_molecules_total_est": 12.0}}
        df = D.build_sensitivity_table(
            _adata(), stats, RunConfig(), index_viruses=["EBV", "Betatorquevirus"]
        )
        assert list(df["virus_name"]) == ["Betatorquevirus", "EBV"]
        assert list(df["observed_molecules"]) == [12, 0]


class TestRunConfigValidation:
    """A control with no known abundance cannot certify anything, so it is refused."""

    def _snake(self, **overrides) -> dict:
        base = {
            "index": "i",
            "transcripts": "t",
            "sample1": "s1",
            "sample2": "s2",
            "output": "o/",
            "technology": "10xv3",
            "visual": "True",
            "multimapping": "True",
            "gtf": "None",
            "fasta": "None",
            "f1": "None",
            "reference": "False",
            "umap": "False",
            "whitelist": "None",
            "emptydrops_seed": 100,
        }
        base.update(overrides)
        return base

    def test_gene_without_molecules_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="must be supplied together"):
            RunConfig.from_snakemake_config(self._snake(positive_control_gene=SPIKEIN))

    def test_molecules_without_gene_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="must be supplied together"):
            RunConfig.from_snakemake_config(self._snake(positive_control_expected_molecules=100))

    def test_non_positive_molecules_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="must be > 0"):
            RunConfig.from_snakemake_config(
                self._snake(
                    positive_control_gene=SPIKEIN,
                    positive_control_expected_molecules=0,
                )
            )

    def test_require_control_without_a_control_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="require_positive_control"):
            RunConfig.from_snakemake_config(self._snake(require_positive_control="True"))

    def test_valid_control_pair_is_accepted(self) -> None:
        cfg = RunConfig.from_snakemake_config(
            self._snake(
                positive_control_gene=SPIKEIN,
                positive_control_expected_molecules="100",
                require_positive_control="True",
            )
        )
        assert cfg.positive_control_gene == SPIKEIN
        assert cfg.positive_control_expected_molecules == 100.0
        assert cfg.require_positive_control is True

    def test_defaults_are_off_and_permissive(self) -> None:
        cfg = RunConfig.from_snakemake_config(self._snake())
        assert cfg.positive_control_gene is None
        assert cfg.positive_control_expected_molecules is None
        assert cfg.require_positive_control is False

    @pytest.mark.parametrize("blank", ["", "none", "None", "null", "  "])
    def test_blank_molecules_is_treated_as_unset(self, blank: str) -> None:
        cfg = RunConfig.from_snakemake_config(
            self._snake(positive_control_expected_molecules=blank)
        )
        assert cfg.positive_control_expected_molecules is None


class TestControlReport:
    def test_report_records_that_an_unmeasured_control_certifies_nothing(self, tmp_path) -> None:
        path = D.write_control_report({"status": "not-configured"}, None, str(tmp_path))
        payload = json.loads((tmp_path / "results" / "positive_control.json").read_text())
        assert payload["status"] == "not-configured"
        assert payload["certifies_negatives"] is False
        assert payload["capture_used_for_sensitivity"] is None
        assert path.endswith("positive_control.json")

    def test_report_records_a_measured_capture(self, tmp_path) -> None:
        detail = {
            "status": "measured",
            "gene": SPIKEIN,
            "expected_molecules": 100.0,
            "observed_molecules": 50.0,
            "capture": 0.5,
            "implied_divergence": 0.134,
        }
        D.write_control_report(detail, 0.5, str(tmp_path))
        payload = json.loads((tmp_path / "results" / "positive_control.json").read_text())
        assert payload["certifies_negatives"] is True
        assert payload["capture_used_for_sensitivity"] == 0.5
        assert payload["gene"] == SPIKEIN
