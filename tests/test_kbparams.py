"""Tests for viralscan.kbparams (PLAN KBPARAMS-01): signals combine into advice, never silently."""

import pytest

from viralscan import kbparams
from viralscan.chemistry import Detection
from viralscan.tchem import Opinion, TchemCall


def _detect(name):
    det = Detection(name, "bundled on-list", "test", 26, None, 10, 1000)
    return lambda paths: [det]


def _opinion(status, strand=None, layout="10x 3' v2"):
    call = TchemCall("match", "high", layout, "error_under_3pct", None, None, strand, "test")
    return lambda tech, r1, r2, **kw: Opinion(status, f"tchem {status}", call)


PILOT_FWD = {"forward": 0.72, "reverse": 0.03, "unstranded": 0.74}
PILOT_REV = {"forward": 0.10, "reverse": 0.52, "unstranded": 0.59}


def rec(name="10xv2", op=None, pilot=None):
    return kbparams.recommend_kb_parameters(
        "r1.fq.gz",
        "r2.fq.gz",
        pilot_rates=pilot,
        detect=_detect(name),
        second_opinion=op or _opinion("agree", "forward"),
    )


def test_three_prime_with_agreeing_signals_is_accepted_and_gives_kb_args():
    r = rec(pilot=PILOT_FWD)
    assert r.status == kbparams.ACCEPTED
    assert r.kb_args() == ["-x", "10xv2", "-w", "", "--strand", "forward"]


def test_unavailable_tchem_leaves_the_technology_for_review():
    r = rec(op=lambda *a, **k: Opinion("unavailable", "tchem is not installed"), pilot=PILOT_FWD)
    assert r.technology.status == kbparams.REVIEW
    assert r.technology.value == "10xv2"


def test_ambiguous_detector_needs_a_flag_and_gives_no_args_for_x_or_w():
    r = rec(name=None, pilot=PILOT_FWD)
    assert r.technology.status == r.whitelist.status == kbparams.NEEDS_FLAG
    assert "-x" not in r.kb_args() and "-w" not in r.kb_args()


def test_tchem_disagreement_is_a_conflict_and_its_strand_is_ignored():
    r = rec(op=_opinion("disagree", "reverse"), pilot=None)
    assert r.technology.status == kbparams.CONFLICT
    assert r.technology.value is None
    assert r.strand.status == kbparams.NEEDS_FLAG  # tchem's strand is for another chemistry


def test_reverse_agreed_by_both_signals_stays_review_with_the_hhv6b_warning():
    r = rec(name="10xv2", op=_opinion("agree", "reverse"), pilot=PILOT_REV)
    assert r.strand.value == "reverse"
    assert r.strand.status == kbparams.REVIEW
    assert kbparams.NONFORWARD_NOTE in r.strand.notes


def test_strand_conflict_falls_back_to_unstranded():
    r = rec(op=_opinion("agree", "forward"), pilot=PILOT_REV)
    assert r.strand.value == "unstranded"
    assert r.strand.status == kbparams.CONFLICT


def test_single_strand_signal_is_review():
    r = rec(op=_opinion("agree", "forward"), pilot=None)
    assert (r.strand.value, r.strand.status) == ("forward", kbparams.REVIEW)


def test_no_strand_signal_needs_a_flag():
    r = rec(op=lambda *a, **k: Opinion("unavailable", "off"), pilot=None)
    assert (r.strand.value, r.strand.status) == (None, kbparams.NEEDS_FLAG)


def test_dropseq_whitelist_is_none_so_kb_builds_no_hidden_cell_caller():
    r = rec(name="dropseq", op=_opinion("agree", "forward", "Drop-seq"), pilot=PILOT_FWD)
    assert r.whitelist.value == "None"
    assert "None" in r.kb_args()


def test_overall_status_is_the_worst_parameter():
    r = rec(op=_opinion("agree", "forward"), pilot=PILOT_REV)
    assert r.status == kbparams.CONFLICT
    assert r.as_block()["status"] == kbparams.CONFLICT


@pytest.mark.parametrize("tau", [0.5, 0.8])
def test_tau_is_passed_to_the_pilot_inference(tau):
    r = kbparams.recommend_kb_parameters(
        "r1",
        "r2",
        pilot_rates={"forward": 0.6, "reverse": 0.1, "unstranded": 0.8},
        tau=tau,
        detect=_detect("10xv2"),
        second_opinion=_opinion("agree", "forward"),
    )
    assert r.evidence["strand_pilot"]["tau"] == tau
