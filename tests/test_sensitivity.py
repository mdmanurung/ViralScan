"""Tests for the detection-sensitivity module.

The numbers asserted here are the ones measured from this repo's own data, not
invented for the test:

* Depth values are the ``counts_unfiltered`` molecule sums of the three covid
  PBMC runs and the bundled EBV LCL run. These are molecules, not reads: the
  EBV run had 103,145,071 pseudoaligned reads against 55,266,624 quantified
  molecules, so the two are not interchangeable and reads would understate the
  LOD95 by 1.9x.
* The Poisson floor was checked by molecule-level downsampling of the EBV run
  (55,266,624 molecules, 1,636,934 EBV): P(detect) stayed 1.0000 down to 1,270
  downsampled reads and first hit 0 at 127 reads, i.e. the behaviour is the
  sampling floor, not capture loss.
* Capture-curve values are the exact-match k-mer model at k=31 for a 90 bp 10x
  cDNA read.
"""

from __future__ import annotations

import math

import pytest

from viralscan.sensitivity import (
    ADEQUATE_LOD_PER_10K,
    DEFAULT_INFORMATIVE_LOD_PER_10K,
    LOD95_MOLECULES,
    classify_lod,
    expected_viral_molecules,
    fragment_capture,
    lod95,
    lod95_per_10k,
    minimum_molecules_for_detection,
    negative_result_statement,
    poisson_detection_probability,
    poisson_zero_probability,
    sensitivity_record,
)

# Measured molecule depths (sum of kb counts_unfiltered/cells_x_genes.mtx).
DEPTH_EBV_LCL = 55_266_624
DEPTH_COVID_NO_FILTER = 21_613_840
DEPTH_COVID_DLIST = 8_404_326
DEPTH_COVID_STAR = 5_308_302
DEPTH_10X_5K_1K = 5_000_000
DEPTH_10X_5K_200 = 1_000_000
DEPTH_SHALLOW = 100_000


class TestPoisson:
    def test_zero_probability_basic(self) -> None:
        assert poisson_zero_probability(0) == 1.0
        assert poisson_zero_probability(3) == pytest.approx(0.049787, rel=1e-5)
        # By construction a Poisson at the LOD95 mean is missed 5% of the time.
        assert poisson_zero_probability(LOD95_MOLECULES) == pytest.approx(0.05, rel=1e-9)

    def test_detection_probability_threshold_one(self) -> None:
        # P(>=1) = 1 - e^-lambda
        assert poisson_detection_probability(3, 1) == pytest.approx(0.950213, rel=1e-5)
        assert poisson_detection_probability(0, 1) == 0.0
        assert poisson_detection_probability(0, 0) == 1.0
        assert poisson_detection_probability(5, 0) == 1.0

    def test_detection_probability_known_values(self) -> None:
        # Cross-checked against scipy.stats.poisson (see test_matches_scipy).
        assert poisson_detection_probability(2, 2) == pytest.approx(0.5939941503, rel=1e-9)
        assert poisson_detection_probability(20, 10) == pytest.approx(0.9950045877, rel=1e-9)
        assert poisson_detection_probability(3, 1) == pytest.approx(0.9502129316, rel=1e-9)

    def test_matches_scipy(self) -> None:
        """The hand-rolled Poisson tail must agree with scipy, not just look sane.

        A wrong tail would silently mis-state P(detect) in every sensitivity
        record, which is the number a reader uses to judge a negative.
        """
        scipy_stats = pytest.importorskip("scipy.stats")
        for lam in (0.0, 0.3, 1.0, 2.0, 3.0, 7.5, 20.0, 100.0, 531.0):
            for thr in (1, 2, 3, 10, 50):
                expected = 0.0 if thr <= 0 else 1.0 - scipy_stats.poisson.cdf(thr - 1, lam)
                assert poisson_detection_probability(lam, thr) == pytest.approx(
                    expected, rel=1e-9, abs=1e-12
                ), f"lambda={lam} threshold={thr}"

    def test_detection_probability_is_monotone(self) -> None:
        prev = -1.0
        for lam in range(0, 40):
            p = poisson_detection_probability(lam, 1)
            assert p >= prev
            prev = p

    def test_detection_probability_stays_in_unit_interval(self) -> None:
        for lam in (0.0, 0.5, 1.0, 7.3, 50.0, 1e4, 1e6):
            for thr in (1, 2, 10, 100):
                p = poisson_detection_probability(lam, thr)
                assert 0.0 <= p <= 1.0

    def test_large_threshold_does_not_blow_up(self) -> None:
        # A loop summing a lower tail must saturate rather than return garbage.
        assert 0.0 <= poisson_detection_probability(3.0, 5000) <= 1.0
        assert poisson_detection_probability(3.0, 5000) == pytest.approx(0.0, abs=1e-6)

    def test_minimum_molecules_for_detection(self) -> None:
        assert minimum_molecules_for_detection(0.95, 1) == pytest.approx(LOD95_MOLECULES, rel=1e-9)
        assert minimum_molecules_for_detection(0.95, 10) == pytest.approx(
            10 * LOD95_MOLECULES, rel=1e-9
        )
        assert minimum_molecules_for_detection(0.95, 0) == 0.0
        with pytest.raises(ValueError):
            minimum_molecules_for_detection(1.5, 1)


class TestExpectedAndLod:
    def test_expected_molecules(self) -> None:
        assert expected_viral_molecules(1e6, 1e-4) == pytest.approx(100.0)
        assert expected_viral_molecules(1e6, 1e-4, capture=0.5) == pytest.approx(50.0)

    def test_expected_molecules_degenerate_inputs(self) -> None:
        assert expected_viral_molecules(0, 0.1) == 0.0
        assert expected_viral_molecules(1e6, 0.0) == 0.0
        assert expected_viral_molecules(1e6, 0.1, capture=0.0) == 0.0

    @pytest.mark.parametrize("depth", [0, -1])
    def test_lod_is_infinite_without_depth(self, depth: float) -> None:
        assert math.isinf(lod95(depth))
        assert math.isinf(lod95_per_10k(depth))

    def test_lod95_matches_hand_calculation(self) -> None:
        # 2.9957 / 5,000,000 * 1e4
        assert lod95_per_10k(DEPTH_10X_5K_1K) == pytest.approx(0.0059915, rel=1e-4)
        assert lod95(DEPTH_10X_5K_1K) == pytest.approx(5.9915e-7, rel=1e-4)

    def test_lod_scales_inversely_with_depth(self) -> None:
        assert lod95_per_10k(DEPTH_10X_5K_1K) == pytest.approx(
            2 * lod95_per_10k(2 * DEPTH_10X_5K_1K), rel=1e-9
        )

    def test_lod_scales_with_threshold(self) -> None:
        assert lod95_per_10k(DEPTH_10X_5K_1K, threshold=10) == pytest.approx(
            10 * lod95_per_10k(DEPTH_10X_5K_1K, threshold=1), rel=1e-9
        )

    def test_lod_scales_inversely_with_capture(self) -> None:
        base = lod95_per_10k(DEPTH_COVID_STAR, capture=1.0)
        assert lod95_per_10k(DEPTH_COVID_STAR, capture=0.5) == pytest.approx(2 * base, rel=1e-9)

    def test_capture_divergence_cost_matches_measurement(self) -> None:
        """A 20%-divergent query needs ~17x more molecules at 90 bp, k=31."""
        assert lod95_per_10k(DEPTH_COVID_STAR, capture=fragment_capture(0.20)) == (
            pytest.approx(lod95_per_10k(DEPTH_COVID_STAR) / 0.0577178, rel=1e-3)
        )

    def test_lod_is_monotone_in_depth(self) -> None:
        depths = [
            DEPTH_SHALLOW,
            DEPTH_10X_5K_200,
            DEPTH_10X_5K_1K,
            DEPTH_COVID_STAR,
            DEPTH_COVID_DLIST,
            DEPTH_COVID_NO_FILTER,
            DEPTH_EBV_LCL,
        ]
        vals = [lod95_per_10k(d) for d in depths]
        assert vals == sorted(vals, reverse=True), "deeper runs must have lower LOD"


class TestFragmentCapture:
    """Exact-match k-mer model, k=31.

    P(fragment captured) = 1 - (1 - (1-d)^k)^(L-k+1)
    """

    def test_no_divergence_is_certain(self) -> None:
        assert fragment_capture(0.0, read_length=90) == pytest.approx(1.0)
        assert fragment_capture(0.0, read_length=150) == pytest.approx(1.0)

    @pytest.mark.parametrize(
        "divergence,expected_90bp",
        [
            (0.05, 1.0),
            (0.10, 0.9031),
            (0.15, 0.3232),
            (0.20, 0.0577),
            (0.25, 0.0080),
            (0.30, 0.0009),
        ],
    )
    def test_curve_at_90bp(self, divergence: float, expected_90bp: float) -> None:
        assert fragment_capture(divergence, read_length=90) == pytest.approx(
            expected_90bp, abs=5e-4
        )

    @pytest.mark.parametrize(
        "divergence,expected_150bp",
        [(0.10, 0.9906), (0.15, 0.5420), (0.20, 0.1121), (0.25, 0.0159), (0.30, 0.0019)],
    )
    def test_curve_at_150bp(self, divergence: float, expected_150bp: float) -> None:
        assert fragment_capture(divergence, read_length=150) == pytest.approx(
            expected_150bp, abs=5e-4
        )

    def test_longer_reads_capture_more(self) -> None:
        for d in (0.10, 0.15, 0.20, 0.25, 0.30):
            assert fragment_capture(d, 150) > fragment_capture(d, 90)

    def test_capture_is_monotone_decreasing_in_divergence(self) -> None:
        prev = 1.1
        for i in range(0, 60):
            c = fragment_capture(i / 100, 90)
            assert c <= prev
            prev = c

    def test_read_shorter_than_k_is_all_or_nothing(self) -> None:
        # A 20 bp fragment can never contain a 31-mer.
        assert fragment_capture(0.0, read_length=20) == pytest.approx(1.0)
        assert fragment_capture(0.05, read_length=20) == pytest.approx(0.95**20, rel=1e-9)

    def test_invalid_divergence_rejected(self) -> None:
        for bad in (-0.01, 1.0, 1.5):
            with pytest.raises(ValueError):
                fragment_capture(bad)


class TestClassifyLod:
    @pytest.mark.parametrize(
        "depth,expected",
        [
            (DEPTH_EBV_LCL, "informative"),
            (DEPTH_COVID_NO_FILTER, "informative"),
            (DEPTH_COVID_DLIST, "informative"),
            (DEPTH_COVID_STAR, "informative"),
            (DEPTH_10X_5K_1K, "informative"),
            (DEPTH_10X_5K_200, "adequate"),
            (DEPTH_SHALLOW, "shallow"),
        ],
    )
    def test_measured_depths_land_in_expected_bands(self, depth: float, expected: str) -> None:
        assert classify_lod(lod95_per_10k(depth)) == expected

    def test_band_edges(self) -> None:
        assert classify_lod(0.0) == "informative"
        assert classify_lod(DEFAULT_INFORMATIVE_LOD_PER_10K) == "informative"
        assert classify_lod(ADEQUATE_LOD_PER_10K) == "adequate"
        assert classify_lod(1.0) == "shallow"
        assert classify_lod(math.inf) == "insufficient-depth"

    def test_every_measured_run_is_depth_sufficient(self) -> None:
        """The load-bearing claim: depth did not limit these runs.

        The covid samples called SARS-CoV-2 = 0 at 21.6M quantified molecules.
        If depth were the binding limit the LOD would be poor; it is not, which
        is why the binding limit has to be reference capture instead.
        """
        for depth in (DEPTH_COVID_STAR, DEPTH_COVID_DLIST, DEPTH_COVID_NO_FILTER):
            assert lod95_per_10k(depth) < DEFAULT_INFORMATIVE_LOD_PER_10K


class TestSensitivityRecord:
    def test_positive_call_has_full_row(self) -> None:
        r = sensitivity_record("Epstein-Barr virus", 1015.0, DEPTH_10X_5K_1K)
        assert r.virus_name == "Epstein-Barr virus"
        assert r.observed_molecules == 1015.0
        assert r.depth_sufficient is True
        assert r.capture_measured is False
        row = r.as_row()
        assert set(row) == {
            "virus_name",
            "observed_molecules",
            "detection_threshold",
            "capture",
            "capture_measured",
            "lod95_per_10k",
            "lod95_molecules",
            "lod_interpretation",
            "depth_sufficient",
            "informative_negative",
            "expected_molecules_at_1_per_10k",
            "p_detect_at_1_per_10k",
            "p_zero_at_lod95",
            "notes",
        }

    def test_negative_is_never_certifiable_without_a_measured_capture(self) -> None:
        """The central guarantee: depth alone cannot certify a negative."""
        r = sensitivity_record("Betatorquevirus", 0.0, DEPTH_COVID_STAR)
        assert r.observed_molecules == 0
        assert r.depth_sufficient is True, "depth was ample in the real covid run"
        assert r.capture_measured is False
        assert r.informative is False
        assert any("negative, not an absence" in n for n in r.notes)

    def test_negative_becomes_certifiable_with_depth_and_measured_capture(self) -> None:
        r = sensitivity_record(
            "Betatorquevirus",
            0.0,
            DEPTH_COVID_STAR,
            capture=fragment_capture(0.05),
            capture_measured=True,
        )
        assert r.capture_measured is True
        assert r.depth_sufficient is True
        assert r.informative is True

    def test_shallow_depth_blocks_certification_even_with_capture(self) -> None:
        r = sensitivity_record(
            "Betatorquevirus",
            0.0,
            DEPTH_SHALLOW,
            capture=0.9,
            capture_measured=True,
        )
        assert r.depth_sufficient is False
        assert r.informative is False
        assert r.lod_interpretation == "shallow"

    def test_measured_capture_worsens_the_lod(self) -> None:
        optimistic = sensitivity_record("X", 0.0, DEPTH_COVID_STAR)
        measured = sensitivity_record(
            "X",
            0.0,
            DEPTH_COVID_STAR,
            capture=fragment_capture(0.20),
            capture_measured=True,
        )
        assert measured.lod95_per_10k > optimistic.lod95_per_10k
        assert measured.lod_interpretation in {"adequate", "shallow", "insufficient-depth"}

    def test_p_detect_at_yardstick_is_high_for_deep_runs(self) -> None:
        r = sensitivity_record("X", 0.0, DEPTH_COVID_STAR)
        # 1 UMI/10k host UMI at 5.3M molecules -> ~531 expected molecules.
        assert r.expected_molecules_if_present_at_1_per_10k == pytest.approx(530.8, rel=1e-3)
        assert r.p_detect_if_present_at_1_per_10k > 0.999

    def test_p_zero_at_lod95_is_five_percent_by_construction(self) -> None:
        r = sensitivity_record("X", 0.0, DEPTH_COVID_STAR)
        assert r.p_zero_at_lod95 == pytest.approx(0.05, rel=1e-6)

    def test_rejects_bad_threshold(self) -> None:
        for bad in (0, -1):
            with pytest.raises(ValueError):
                sensitivity_record("X", 0.0, DEPTH_10X_5K_1K, detection_threshold=bad)

    def test_rejects_non_positive_capture(self) -> None:
        for bad in (0.0, -0.5):
            with pytest.raises(ValueError, match="capture must be > 0"):
                sensitivity_record("X", 0.0, DEPTH_10X_5K_1K, capture=bad)

    def test_capture_above_one_is_accepted_and_clamped_to_the_depth_floor(self) -> None:
        """A control that over-recovers must not crash the run or over-tighten LOD.

        capture > 1 means more was recovered than planted, so detection is at
        least as easy as the sampling floor; the LOD is clamped there and the
        record says the value bounds loss from above only.
        """
        at_floor = sensitivity_record("X", 0.0, DEPTH_COVID_STAR, capture=1.0)
        over = sensitivity_record("X", 0.0, DEPTH_COVID_STAR, capture=1.18, capture_measured=True)
        assert over.lod95_per_10k == pytest.approx(at_floor.lod95_per_10k)
        assert over.capture == 1.18
        assert any("no capture loss was measurable" in n for n in over.notes)
        assert any("does not demonstrate zero divergence" in n for n in over.notes)

    def test_exact_unit_capture_measured_is_allowed_and_annotated(self) -> None:
        """capture=1.0 is a legitimate measurement when the control fully recovered."""
        r = sensitivity_record("X", 0.0, DEPTH_COVID_STAR, capture=1.0, capture_measured=True)
        assert r.capture_measured is True
        assert any("no capture loss was measurable" in n for n in r.notes)


class TestNegativeResultStatement:
    def test_zero_depth_is_reported_as_quantification_failure(self) -> None:
        s = negative_result_statement(0)
        assert "NOT ESTIMABLE" in s
        assert "quantification failure" in s

    def test_statement_always_names_the_binding_limit(self) -> None:
        """Without a measured capture term, the statement must say so."""
        s = negative_result_statement(DEPTH_COVID_STAR)
        assert "ASSUMED 1.0, not measured" in s
        assert "CANNOT be read as absence" in s
        assert "0.32" in s and "0.06" in s, "must quantify the divergence cliff"

    def test_statement_reports_depth_and_molecules(self) -> None:
        s = negative_result_statement(DEPTH_COVID_STAR)
        assert "0.005643" in s
        assert "3 true viral molecules" in s
        assert f"{DEPTH_COVID_STAR:,}" in s

    def test_statement_says_depth_is_sufficient_for_measured_runs(self) -> None:
        for depth in (DEPTH_COVID_STAR, DEPTH_COVID_DLIST, DEPTH_COVID_NO_FILTER):
            assert "Depth is sufficient" in negative_result_statement(depth)

    def test_statement_agrees_with_the_depth_sufficient_column(self) -> None:
        """Prose and TSV must never disagree about the same run.

        ``depth_sufficient`` is the machine-readable field a reader may filter
        on, so the sentence has to use the same ceiling, not a looser band test.
        """
        for depth in (
            DEPTH_COVID_STAR,
            DEPTH_COVID_DLIST,
            DEPTH_COVID_NO_FILTER,
            DEPTH_10X_5K_1K,
            DEPTH_10X_5K_200,
            DEPTH_SHALLOW,
            DEPTH_EBV_LCL,
        ):
            record = sensitivity_record("X", 0.0, depth)
            statement = negative_result_statement(depth)
            says_sufficient = "Depth is sufficient" in statement
            assert says_sufficient == record.depth_sufficient, (
                f"depth={depth}: column says depth_sufficient={record.depth_sufficient} "
                f"but the statement says otherwise: {statement[:120]}"
            )

    def test_statement_flags_insufficient_depth(self) -> None:
        s = negative_result_statement(DEPTH_SHALLOW)
        assert "Depth is NOT sufficient" in s
        assert "Deepen the library" in s

    def test_measured_capture_makes_the_statement_certifiable(self) -> None:
        s = negative_result_statement(
            DEPTH_COVID_STAR, capture=fragment_capture(0.05), capture_measured=True
        )
        assert "measured at" in s
        assert "certifiable" in s
        assert "CANNOT be read as absence" not in s
