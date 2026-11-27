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


V3 = {"10xv2": 0.093, "10xv3": 0.96, "10xv4": 0.05}


class TestR1Length:
    """CHEM-R1LEN-01: an R1 longer than barcode + UMI is fine; only an adapter or a short R1 refuses."""

    @pytest.mark.parametrize(("r1_len", "trailing"), [(26, 0), (28, 2), (150, 124)])
    def test_10xv2_with_any_r1_at_least_26(self, r1_len, trailing):
        # BC-1 (GSE154900): 10xv2 with a 28-cycle R1; the bases after the UMI are random.
        d = classify(_reads(umi=10, tail=_bc(124), r1_len=r1_len), EBV)
        assert (d.chemistry, d.basis, d.trailing_bases) == ("10xv2", "bundled on-list", trailing)
        assert resolve(None, [d]) == "10xv2" and resolve("10xv2", [d]) == "10xv2"
        assert d.as_block()["trailing_bases"] == trailing

    def test_trailing_bases_are_in_the_reason(self):
        d = classify(_reads(umi=10, tail=_bc(2)), EBV)
        assert "ignores the 2 bp" in d.reason

    @pytest.mark.parametrize(("r1_len", "trailing"), [(28, 0), (90, 62), (150, 122)])
    def test_10xv3_with_any_r1_at_least_28(self, r1_len, trailing):
        d = classify(_reads(umi=12, tail=_bc(122), r1_len=r1_len), V3)
        assert (d.chemistry, d.trailing_bases) == ("10xv3", trailing)

    def test_10xv3_polyt_after_the_umi_is_accepted(self):
        d = classify(_reads(umi=12, tail="T" * 62), V3)
        assert (d.chemistry, d.end, d.trailing_bases) == ("10xv3", "3p", 62)

    def test_10xv2_with_a_polyt_at_26_is_accepted(self):
        d = classify(_reads(umi=10, tail="T" * 62), EBV)
        assert (d.chemistry, d.end, d.trailing_bases) == ("10xv2", "3p", 62)

    def test_hhv6b_tso_at_26_still_passes(self):
        # SRR20710641: 150 bp R1, the TSO at base 26 fixes the UMI at 10.
        d = classify(_reads(umi=10, tail=ch.TSO + "G" * 111), HHV6B)
        assert (d.chemistry, d.end, d.umi_len, d.trailing_bases) == ("10xv2", "5p", 10, 124)

    def test_v2_list_with_a_tso_at_28_is_refused(self):
        # 12 bp UMI evidence from the adapter contradicts the v2 geometry.
        d = classify(_reads(umi=12, tail=ch.TSO + "G" * 111), EBV)
        assert d.chemistry is None and "implies a 12 bp UMI, not 10" in d.reason

    def test_v3_list_with_a_tso_at_26_is_refused(self):
        d = classify(_reads(umi=10, tail=ch.TSO + "G" * 111), V3)
        assert d.chemistry is None and "implies a 10 bp UMI, not 12" in d.reason

    @pytest.mark.parametrize(
        ("rates", "umi", "r1_len"), [(EBV, 10, 24), (EBV, 10, 25), (V3, 12, 26), (V3, 12, 27)]
    )
    def test_r1_shorter_than_barcode_plus_umi_is_refused(self, rates, umi, r1_len):
        d = classify(_reads(umi=umi, r1_len=r1_len), rates)
        assert d.chemistry is None and "shorter than" in d.reason
        with pytest.raises(ChemistryError, match="Chemistry check failed"):
            resolve(None, [d])


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


def test_cell_ranger_dash_one_suffix_on_the_w_list_still_matches_b27(tmp_path):
    on = [_bc(16) for _ in range(5)]
    plain, suffixed = tmp_path / "plain.txt", tmp_path / "cr.txt"
    plain.write_text("\n".join(on) + "\n")
    suffixed.write_text("\n".join(b + "-1" for b in on) + "\n")
    seqs = [on[0] + "A" * 10, on[1] + "C" * 10, _bc(26), _bc(26)]
    rates = ch.match_rates([seqs], {"plain": str(plain), "cr": str(suffixed)})
    assert rates == [{"plain": 0.5, "cr": 0.5}]


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


def test_dropseq_with_its_own_w_list_is_dropseq():
    """A 20 bp R1 gives no UMI length, so the -w branch must not refuse it."""
    seqs = [_bc(20) for _ in range(200)]
    d = classify(seqs, {**HSV1, "user": 0.9}, user_list="dropseq_wl.txt")
    assert (d.chemistry, d.basis) == ("dropseq", "user on-list")
    assert resolve("dropseq", [d]) == "dropseq"
