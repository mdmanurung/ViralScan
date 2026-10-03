"""Chemistry detection and the -x check (PLAN MECH-D, DEF-02, WP1E Q6).

The rate sets are the ones measured on the four local libraries (see the
module docstring); the reads are synthetic with the same structure.
"""

from __future__ import annotations

import gzip
import random

import pytest

from viralscan import chemistry as ch
from viralscan.chemistry import ChemistryError, classify, read_structure, resolve

random.seed(0)


def _bc(n: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(n))


def _reads(n=200, umi=10, tail="", r1_len=None):
    seqs = [_bc(16) + _bc(umi) + tail for _ in range(n)]
    return [s[:r1_len] for s in seqs] if r1_len else seqs


EBV = {"10xv2": 0.973, "10xv3": 0.093, "10xv4": 0.007}
HHV6B = {"10xv2": 0.882, "10xv3": 0.097, "10xv4": 0.009}
COVID = {"10xv2": 0.0, "10xv3": 0.005, "10xv4": 0.044}
HSV1 = {"10xv2": 0.0, "10xv3": 0.001, "10xv4": 0.002}


class TestReadStructure:
    def test_trimmed_26_is_a_10bp_umi(self):
        assert read_structure(_reads(umi=10)) == (26, None, 10)

    def test_tso_at_26_is_5p_with_a_10bp_umi(self):
        seqs = _reads(umi=10, tail=ch.TSO + "G" * 111)
        assert read_structure(seqs) == (150, "5p", 10)

    def test_polyt_at_28_is_3p_with_a_12bp_umi(self):
        seqs = _reads(umi=12, tail="T" * 30)
        assert read_structure(seqs) == (58, "3p", 12)


class TestMeasuredLibraries:
    def test_ebv_10xv2_3p(self):
        d = classify(_reads(umi=10), EBV)
        assert (d.chemistry, d.basis) == ("10xv2", "bundled on-list")

    def test_hhv6b_5p_is_10xv2(self):
        d = classify(_reads(umi=10, tail=ch.TSO + "G" * 111), HHV6B)
        assert (d.chemistry, d.end) == ("10xv2", "5p")

    def test_gemx_5p_without_w_fails_closed(self):
        d = classify(_reads(umi=12), COVID)
        assert d.chemistry is None and "-w" in d.reason

    def test_gemx_5p_with_its_list_is_10xv3_geometry(self):
        d = classify(_reads(umi=12), {**COVID, "user": 0.672}, user_list="cr.txt")
        assert (d.chemistry, d.basis) == ("10xv3", "user on-list")

    def test_dropseq_is_inferred_from_a_20bp_r1(self):
        seqs = [_bc(20) for _ in range(200)]
        assert classify(seqs, HSV1).chemistry == "dropseq"


class TestResolve:
    def test_mech_d_10xv2_run_as_10xv3_is_refused(self):
        d = classify(_reads(umi=10), EBV)
        with pytest.raises(ChemistryError, match="look like 10xv2"):
            resolve("10xv3", [d])

    def test_f005_gemx_labelled_10xv3_is_refused(self):
        with pytest.raises(ChemistryError, match="no on-list matches"):
            resolve("10xv3", [classify(_reads(umi=12), COVID)])

    def test_wrong_w_list_is_refused(self):
        d = classify(_reads(umi=10), {**EBV, "user": 0.002}, user_list="other.txt")
        with pytest.raises(ChemistryError, match="does not belong"):
            resolve("10xv2", [d])

    def test_auto_takes_the_detection(self):
        assert resolve(None, [classify(_reads(umi=10), EBV)]) == "10xv2"

    def test_user_list_compares_geometry_not_name(self):
        d = classify(_reads(umi=12), {**COVID, "user": 0.672}, user_list="cr.txt")
        assert resolve("10xv4", [d]) == "10xv4"

    def test_samples_that_disagree_are_refused(self):
        a = classify(_reads(umi=10), EBV)
        b = classify([_bc(20) for _ in range(200)], HSV1)
        with pytest.raises(ChemistryError, match="disagree"):
            resolve(None, [a, b])

    def test_force_keeps_an_explicit_x(self):
        d = classify(_reads(umi=10), EBV)
        assert resolve("10xv3", [d], force=True) == "10xv3"

    def test_force_without_x_is_refused(self):
        with pytest.raises(ChemistryError, match="explicit -x"):
            resolve(None, [], force=True)


def test_match_rates_streams_lists(tmp_path):
    on = [_bc(16) for _ in range(5)]
    listed = tmp_path / "list.txt.gz"
    with gzip.open(listed, "wt") as fh:
        fh.write("\n".join(on) + "\n")
    seqs = [on[0] + "A" * 10, on[1] + "C" * 10, _bc(26), _bc(26)]
    assert ch.match_rates([seqs], {"x": str(listed)}) == [{"x": 0.5}]


def test_every_kb_chemistry_with_an_onlist_resolves_in_ngs_tools():
    pytest.importorskip("ngs_tools")
    for name in ch.DETECTABLE:
        assert ch.onlist_path(name).endswith(ch.CHEMISTRIES[name].onlist)
    assert ch.onlist_path("dropseq") is None


class TestKbWhitelistArg:
    def test_user_list_wins(self):
        assert ch.kb_whitelist_arg("dropseq", "/wl.txt") == "/wl.txt"

    def test_10x_without_w_uses_kbs_onlist(self):
        assert ch.kb_whitelist_arg("10xv3", None) == ""

    def test_sw21_no_onlist_bypasses_correction(self):
        assert ch.kb_whitelist_arg("DROPSEQ", None) == "None"

    def test_custom_geometry_is_left_to_kb(self):
        assert ch.kb_whitelist_arg("0,0,16:0,16,28:1,0,0", None) == ""
