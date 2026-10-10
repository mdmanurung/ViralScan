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

from tests._exact_control import exact_control_kwargs
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
    """A legacy control: no scope, so it maps to panel_mechanics and certifies no virus."""
    base = {"positive_control_gene": SPIKEIN, "positive_control_expected_molecules": 100.0}
    base.update(kwargs)
    return RunConfig(**base)


def _scoped_config(virus: str = "Betatorquevirus", **kwargs) -> RunConfig:
    """An exact_sequence control attached to ``virus``: certifies that row only."""
    kwargs.setdefault("positive_control_scope", "exact_sequence")
    kwargs.setdefault("positive_control_virus_key", virus)
    return _config(**kwargs)


def _verified_config(tmp_path, monkeypatch, virus: str = "Betatorquevirus") -> RunConfig:
    """exact_sequence control with a checksum-bound receipt: the only kind that certifies."""
    return _scoped_config(virus, **exact_control_kwargs(tmp_path, monkeypatch, SPIKEIN, virus))


def _table(adata, stats, cfg, index_viruses=()):
    """build_sensitivity_table fed by the same measurement run() uses."""
    _, detail = D.measure_positive_control(adata, cfg)
    return D.build_sensitivity_table(
        adata, stats, cfg, control_detail=detail, index_viruses=index_viruses
    )


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
        # Exact substitution-only model (SENS-CORR-02): capture 0.5 sits between the
        # 5 % (0.71) and 10 % (0.25) points. Under the retired overlapping-window
        # formula this was 0.134; that formula overstated capture.
        assert 0.05 < detail["implied_divergence"] < 0.10
        assert "substitution-only heuristic" in detail["implied_divergence_note"]

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
        assert D.substitution_model_implied_divergence(1.0) is None
        assert D.substitution_model_implied_divergence(0.999) is None

    def test_implied_divergence_is_none_below_the_identifiable_floor(self) -> None:
        """Capture under ~1e-8 cannot be attributed to divergence in float64."""
        from viralscan.sensitivity import fragment_capture_exact

        floor = fragment_capture_exact(D.MAX_IDENTIFIABLE_DIVERGENCE, read_length=90)
        assert D.substitution_model_implied_divergence(floor) is None
        assert D.substitution_model_implied_divergence(0.0) is None
        assert D.substitution_model_implied_divergence(-1.0) is None

    def test_implied_divergence_decreases_as_capture_increases(self) -> None:
        """More capture means less divergence, so the inversion is monotone down."""
        prev = 1.1
        for capture in (1e-6, 0.001, 0.008, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95):
            d = D.substitution_model_implied_divergence(capture)
            assert d is not None, f"capture={capture} should be identifiable"
            assert d < prev, f"capture={capture} -> {d}, previous {prev}"
            prev = d

    def test_implied_divergence_recovers_known_curve_points(self) -> None:
        from viralscan.sensitivity import fragment_capture_exact

        for d in (0.10, 0.15, 0.20, 0.25, 0.30, 0.40):
            capture = fragment_capture_exact(d, read_length=90)
            assert D.substitution_model_implied_divergence(capture) == pytest.approx(d, abs=1e-3)

    def test_inversion_is_labelled_a_heuristic_and_never_a_measurement(self) -> None:
        _, detail = D.measure_positive_control(_adata(50.0), _config())
        assert "not a measured" in detail["implied_divergence_note"]
        assert detail["status"] == "measured"


class TestCertificationFlip:
    """The same negative, with and without a working control."""

    # A negative is a virus that did NOT clear the gate, i.e. an indexed virus absent
    # from virus_stats. (These tests used to put a zero-molecule entry in virus_stats,
    # which the pipeline treats as a *called* virus; SENS-CORR-03 made membership of
    # virus_stats the detection decision, so the old fixture described a detected row.)
    STATS: dict = {}
    INDEX = ["Betatorquevirus"]

    def test_without_control_a_negative_is_never_certifiable(self) -> None:
        adata = _adata()
        row = _table(adata, self.STATS, RunConfig(), self.INDEX).iloc[0]
        assert row["observed_molecules"] == 0
        assert row["capture_measured"] is False or row["capture_measured"] == False  # noqa: E712
        assert not row["informative_negative"]

    def test_failed_control_leaves_the_negative_uncertifiable(self) -> None:
        adata = _adata(0.0)  # planted 100, recovered 0
        cfg = _scoped_config()
        row = _table(adata, self.STATS, cfg, self.INDEX).iloc[0]
        assert row["capture_measured"] == False  # noqa: E712
        assert not row["informative_negative"]
        assert row["depth_sufficient"] in (True, False)  # depth is reported either way

    def test_in_scope_control_at_good_depth_certifies_the_negative(
        self, tmp_path, monkeypatch
    ) -> None:
        adata = _adata(50.0, n_cells=4000)
        cfg = _verified_config(tmp_path, monkeypatch)
        capture, _ = D.measure_positive_control(adata, cfg)
        assert capture == pytest.approx(0.5)
        row = _table(adata, self.STATS, cfg, self.INDEX).iloc[0]
        assert row["capture_measured"]
        assert row["depth_sufficient"]
        assert row["informative_negative"]
        assert row["capture"] == pytest.approx(0.5)
        assert "certifiable for that scope only" in row["notes"]
        assert "exact_sequence" in row["notes"] and "any other virus" in row["notes"]

    def test_legacy_control_without_scope_certifies_nothing(self) -> None:
        """CHANGED ON PURPOSE (SENS-CORR-01). This used to certify every virus.

        A control with no declared scope maps to panel_mechanics: it proves the
        pipeline recovers a planted molecule, not that it sees any other virus.
        """
        adata = _adata(50.0, n_cells=4000)
        cfg = _config()  # measured capture 0.5, no scope
        capture, detail = D.measure_positive_control(adata, cfg)
        assert capture == pytest.approx(0.5) and detail["status"] == "measured"
        row = _table(adata, self.STATS, cfg, self.INDEX).iloc[0]
        assert row["capture_measured"] == False  # noqa: E712
        assert pd.isna(
            row["capture"]
        )  # empty, so the depth-only floor never reads as a measurement
        assert not row["informative_negative"]
        assert "NOT measured" in row["notes"]

    def test_out_of_scope_rows_get_no_capture_and_never_borrow_the_in_scope_value(
        self, tmp_path, monkeypatch
    ) -> None:
        adata = _adata(50.0, n_cells=4000)
        cfg = _verified_config(tmp_path, monkeypatch)
        df = _table(
            adata,
            {},
            cfg,
            index_viruses=["Betatorquevirus", "Epstein-Barr virus", "Human herpesvirus 1"],
        ).set_index("virus_name")
        assert bool(df.loc["Betatorquevirus", "informative_negative"]) is True
        for other in ("Epstein-Barr virus", "Human herpesvirus 1"):  # incl. undetected-row
            assert bool(df.loc[other, "capture_measured"]) is False
            assert bool(df.loc[other, "informative_negative"]) is False
            assert pd.isna(df.loc[other, "capture"])  # empty, never a borrowed value

    def test_over_recovered_control_certifies_nothing_even_in_scope(
        self, tmp_path, monkeypatch
    ) -> None:
        """measure_positive_control clamps over-recovery to 1.0; that must not certify."""
        adata = _adata(118.0, n_cells=4000)
        cfg = _verified_config(tmp_path, monkeypatch)
        capture, detail = D.measure_positive_control(adata, cfg)
        assert capture == 1.0 and detail["status"] == "over-recovered"
        row = _table(adata, self.STATS, cfg, self.INDEX).iloc[0]
        assert bool(row["capture_measured"]) is False
        assert bool(row["informative_negative"]) is False

    def test_measured_capture_never_tightens_below_the_depth_floor(self) -> None:
        """Over-recovery must not produce a sub-floor LOD."""
        adata = _adata(118.0)
        cfg = _scoped_config()
        row = _table(adata, self.STATS, cfg, self.INDEX).iloc[0]
        floor = D.build_sensitivity_table(
            adata, self.STATS, RunConfig(), index_viruses=self.INDEX
        ).iloc[0]
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


class TestScopeConfig:
    """--positive-control-scope / --positive-control-virus-key validation (SENS-CORR-01)."""

    _snake = TestRunConfigValidation._snake

    def _cfg(self, **overrides):
        base = {
            "positive_control_gene": SPIKEIN,
            "positive_control_expected_molecules": "100",
        }
        base.update(overrides)
        return RunConfig.from_snakemake_config(self._snake(**base))

    def test_unset_scope_is_none_and_not_defaulted(self) -> None:
        """Left as None so the run fingerprint and legacy configs are unchanged."""
        cfg = self._cfg()
        assert cfg.positive_control_scope is None
        assert cfg.positive_control_virus_key is None

    def test_exact_sequence_with_key_is_accepted(self) -> None:
        cfg = self._cfg(
            positive_control_scope="exact_sequence", positive_control_virus_key="Torque teno virus"
        )
        assert cfg.positive_control_scope == "exact_sequence"
        assert cfg.positive_control_virus_key == "Torque teno virus"

    def test_panel_mechanics_without_key_is_accepted(self) -> None:
        assert self._cfg(positive_control_scope="panel_mechanics").positive_control_scope == (
            "panel_mechanics"
        )

    @pytest.mark.parametrize("blank", ["", "none", "None"])
    def test_blank_values_are_unset(self, blank: str) -> None:
        cfg = self._cfg(positive_control_scope=blank, positive_control_virus_key=blank)
        assert cfg.positive_control_scope is None and cfg.positive_control_virus_key is None

    def test_unknown_scope_rejected(self) -> None:
        with pytest.raises(ValueError, match="positive_control_scope"):
            self._cfg(positive_control_scope="global")

    def test_virus_key_scope_rejected_without_transfer_calibration(self) -> None:
        with pytest.raises(ValueError, match="calibration"):
            self._cfg(positive_control_scope="virus_key", positive_control_virus_key="HHV4")

    def test_exact_sequence_needs_a_target(self) -> None:
        with pytest.raises(ValueError, match="positive_control_virus_key"):
            self._cfg(positive_control_scope="exact_sequence")

    def test_key_without_scope_rejected(self) -> None:
        with pytest.raises(ValueError, match="positive_control_scope"):
            self._cfg(positive_control_virus_key="HHV4")

    def test_panel_mechanics_with_key_rejected(self) -> None:
        with pytest.raises(ValueError, match="certifies no virus"):
            self._cfg(positive_control_scope="panel_mechanics", positive_control_virus_key="HHV4")

    def test_scope_without_a_control_rejected(self) -> None:
        with pytest.raises(ValueError, match="require a positive control"):
            RunConfig.from_snakemake_config(
                self._snake(
                    positive_control_scope="exact_sequence", positive_control_virus_key="HHV4"
                )
            )

    def test_yaml_and_snakemake_wire_roundtrip_keep_scope(self, tmp_path) -> None:
        cfg = self._cfg(positive_control_scope="exact_sequence", positive_control_virus_key="X y")
        path = tmp_path / "config.yaml"
        cfg.to_yaml(path)
        back = RunConfig.from_yaml(path)
        assert (back.positive_control_scope, back.positive_control_virus_key) == (
            "exact_sequence",
            "X y",
        )
        wire = cfg.to_snakemake_config_args()
        assert "positive_control_scope=exact_sequence" in wire
        assert "positive_control_virus_key=X y" in wire
        assert "positive_control_scope=" in RunConfig(**{}).to_snakemake_config_args()


class TestLegacyScopeWarning:
    def test_legacy_control_warns_once_and_maps_to_panel_mechanics(
        self, caplog, monkeypatch
    ) -> None:
        monkeypatch.setattr(D, "_LEGACY_SCOPE_WARNED", False)
        with caplog.at_level("WARNING"):
            first = D.control_claim(_config())
            second = D.control_claim(_config())
        assert first == second == (D.CaptureScope.PANEL_MECHANICS, None)
        warned = [r for r in caplog.records if "panel_mechanics" in r.getMessage()]
        assert len(warned) == 1

    def test_no_control_needs_no_warning(self, caplog, monkeypatch) -> None:
        monkeypatch.setattr(D, "_LEGACY_SCOPE_WARNED", False)
        with caplog.at_level("WARNING"):
            D.control_claim(RunConfig())
        assert not caplog.records

    def test_explicit_scope_does_not_warn(self, caplog, monkeypatch) -> None:
        monkeypatch.setattr(D, "_LEGACY_SCOPE_WARNED", False)
        with caplog.at_level("WARNING"):
            claim = D.control_claim(_scoped_config("HHV4"))
        assert claim == (D.CaptureScope.EXACT_SEQUENCE, "HHV4")
        assert not caplog.records


class TestControlReport:
    """positive_control.json: ``certifies_negatives`` is per in-scope row, never panel-wide."""

    def _write(self, tmp_path, adata, cfg, viruses):
        _, detail = D.measure_positive_control(adata, cfg)
        control = D.positive_control_from(cfg, detail)
        certified = D.certified_viruses(control, viruses)
        D.write_control_report(detail, str(tmp_path), control=control, certified_targets=certified)
        return json.loads((tmp_path / "results" / "positive_control.json").read_text())

    def test_report_records_that_an_unmeasured_control_certifies_nothing(self, tmp_path) -> None:
        path = D.write_control_report({"status": "not-configured"}, str(tmp_path))
        payload = json.loads((tmp_path / "results" / "positive_control.json").read_text())
        assert payload["status"] == "not-configured"
        assert payload["certifies_negatives"] is False
        assert payload["capture_used_for_sensitivity"] is None
        assert payload["certified_targets"] == []
        assert path.endswith("positive_control.json")

    def test_in_scope_control_certifies_only_its_target(self, tmp_path, monkeypatch) -> None:
        payload = self._write(
            tmp_path,
            _adata(50.0),
            _verified_config(tmp_path, monkeypatch),
            ["Betatorquevirus", "Epstein-Barr virus"],
        )
        assert payload["certifies_negatives"] is True
        assert payload["certified_targets"] == ["Betatorquevirus"]
        assert payload["scope"] == "exact_sequence"
        assert payload["target"] == "Betatorquevirus"
        assert payload["capture_used_for_sensitivity"] == 0.5
        assert payload["gene"] == SPIKEIN

    def test_legacy_control_measures_but_certifies_nothing(self, tmp_path) -> None:
        """CHANGED ON PURPOSE (SENS-CORR-01): a measured legacy control was global."""
        payload = self._write(tmp_path, _adata(50.0), _config(), ["Betatorquevirus"])
        assert payload["status"] == "measured"
        assert payload["capture"] == 0.5
        assert payload["certifies_negatives"] is False
        assert payload["certified_targets"] == []
        assert payload["capture_used_for_sensitivity"] is None
        assert payload["scope"] == "panel_mechanics"
        assert "legacy" in payload["scope_note"]

    def test_scoped_control_whose_target_is_not_in_the_table_certifies_nothing(
        self, tmp_path
    ) -> None:
        payload = self._write(tmp_path, _adata(50.0), _scoped_config("HHV4"), ["Betatorquevirus"])
        assert payload["certifies_negatives"] is False and payload["certified_targets"] == []

    def test_over_recovered_control_certifies_nothing_in_scope(self, tmp_path, monkeypatch) -> None:
        cfg = _verified_config(tmp_path, monkeypatch)
        payload = self._write(tmp_path, _adata(118.0), cfg, ["Betatorquevirus"])
        assert payload["status"] == "over-recovered"
        assert payload["certifies_negatives"] is False

    def test_failed_control_is_still_written_with_its_diagnostics(
        self, tmp_path, monkeypatch
    ) -> None:
        cfg = _verified_config(tmp_path, monkeypatch)
        payload = self._write(tmp_path, _adata(0.0), cfg, ["Betatorquevirus"])
        assert payload["status"] == "failed" and "uninterpretable" in payload["detail"]
        assert payload["certifies_negatives"] is False

    def test_implied_divergence_is_labelled_substitution_only(self, tmp_path) -> None:
        payload = self._write(tmp_path, _adata(50.0), _scoped_config(), ["Betatorquevirus"])
        assert "substitution-only heuristic" in payload["implied_divergence_note"]
        assert "not a measured" in payload["implied_divergence_note"]
