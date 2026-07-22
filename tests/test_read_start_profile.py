"""A1: read-start distribution + PCR-dup handling (pure logic; no samtools needed)."""

from __future__ import annotations

import pytest

from viralscan.evidence import (
    _cb_umi,
    _cigar_ref_span,
    _parse_sam_read_starts,
    add_cell_tags_to_sam,
    deduplicate_umi_sam,
)


def _sam(qname, flag, rname, pos, cigar="50M"):
    return f"{qname}\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t*\t*"


def test_cigar_ref_span_counts_only_reference_consuming_ops():
    assert _cigar_ref_span("50M") == 50
    assert _cigar_ref_span("10M5D10M") == 25  # D consumes reference
    assert _cigar_ref_span("10M5I10M") == 20  # I does NOT
    assert _cigar_ref_span("5S40M5S") == 40  # soft-clips do not
    assert _cigar_ref_span("*") == 0


def test_cb_umi_parsing():
    assert _cb_umi("ACGT_TTTT_7") == ("ACGT", "TTTT")
    assert _cb_umi("noname") is None


def test_forward_starts_and_umi_dedup():
    sam = "\n".join(
        [
            "@HD\tVN:1.6",  # header skipped
            _sam("CB1_U1_1", 0, "SARS", 100),  # start 99
            _sam("CB1_U1_2", 0, "SARS", 100),  # PCR dup of CB1/U1 -> collapsed under umi
            _sam("CB2_U2_1", 0, "SARS", 100),  # start 99 (different molecule)
            _sam("CB3_U3_1", 0, "SARS", 200),  # start 199
        ]
    )
    umi = {(r["reference"], r["position"]): r for r in _parse_sam_read_starts(sam, dedup="umi")}
    assert umi[("SARS", 99)]["n_read_starts"] == 2  # CB1/U1 + CB2/U2
    assert umi[("SARS", 199)]["n_read_starts"] == 1
    assert umi[("SARS", 99)]["n_reads"] == 3  # total deduped reads on SARS

    none = {(r["reference"], r["position"]): r for r in _parse_sam_read_starts(sam, dedup="none")}
    assert none[("SARS", 99)]["n_read_starts"] == 3  # dup NOT collapsed


def test_reverse_strand_uses_3prime_end_when_strand_aware():
    sam = _sam("CB1_U1_1", 16, "SARS", 100, "50M")  # reverse; 5' end = 99 + 50 - 1 = 148
    aware = _parse_sam_read_starts(sam, dedup="none", strand_aware=True)
    assert aware[0]["position"] == 148
    plain = _parse_sam_read_starts(sam, dedup="none", strand_aware=False)
    assert plain[0]["position"] == 99  # leftmost POS


def test_secondary_supplementary_unmapped_are_skipped():
    sam = "\n".join(
        [
            _sam("CB1_U1_1", 0x100, "SARS", 100),  # secondary
            _sam("CB2_U2_1", 0x800, "SARS", 100),  # supplementary
            _sam("CB3_U3_1", 0x4, "*", 0),  # unmapped
            _sam("CB4_U4_1", 0, "SARS", 100),  # the only primary mapped read
        ]
    )
    rows = _parse_sam_read_starts(sam, dedup="none")
    assert len(rows) == 1 and rows[0]["n_read_starts"] == 1


def test_bin_size_bins_positions():
    sam = "\n".join([_sam("CB1_U1_1", 0, "SARS", 105), _sam("CB2_U2_1", 0, "SARS", 112)])
    rows = _parse_sam_read_starts(sam, dedup="none", bin_size=10)
    # starts 104 and 111 -> bins 100 and 110
    assert {r["position"] for r in rows} == {100, 110}


def test_invalid_bin_size():
    with pytest.raises(ValueError):
        _parse_sam_read_starts("", bin_size=0)


def test_add_cell_tags_appends_cb_ub_from_read_name():
    sam = "\n".join(
        [
            "@HD\tVN:1.6",  # header preserved
            _sam("ACGT_TTTT_1", 0, "SARS", 100),  # -> CB:Z:ACGT UB:Z:TTTT
            _sam("nobarcode", 0, "SARS", 100),  # un-parseable name -> unchanged
        ]
    )
    out = add_cell_tags_to_sam(sam).splitlines()
    assert out[0] == "@HD\tVN:1.6"
    assert out[1].endswith("\tCB:Z:ACGT\tUB:Z:TTTT")
    assert "CB:Z:" not in out[2]  # no barcode parsed, no tag added


def test_add_cell_tags_skips_malformed_and_empty_barcode_records():
    malformed = "READ\t0\tSARS\t100"  # only 4 fields (< 11 mandatory) -> untouched
    empty_cb = _sam("_TTTT_1", 0, "SARS", 100)  # empty CB -> no tag (invalid SAM value)
    out = add_cell_tags_to_sam(malformed + "\n" + empty_cb).splitlines()
    assert out[0] == malformed  # passed through, not corrupted with a tag
    assert "CB:Z:" not in out[1]


def test_cb_umi_rejects_empty_components():
    assert _cb_umi("_TTTT_1") is None
    assert _cb_umi("ACGT__1") is None


def test_umi_bam_dedup_contract_includes_position_strand_and_cigar():
    sam = "\n".join(
        [
            "@HD\tVN:1.6",
            _sam("CB_U1_1", 0, "VIRUS|v", 100, "50M"),
            _sam("CB_U1_2", 0, "VIRUS|v", 100, "50M"),  # exact PCR duplicate
            _sam("CB_U1_3", 0, "VIRUS|v", 101, "50M"),  # different start
            _sam("CB_U1_4", 16, "VIRUS|v", 100, "50M"),  # different strand
            _sam("CB_U1_5", 0, "VIRUS|v", 100, "25M1D25M"),  # different CIGAR
        ]
    )
    lines = deduplicate_umi_sam(sam).splitlines()
    assert lines[0].startswith("@HD")
    assert len(lines) == 5  # header + four distinct alignment contexts
