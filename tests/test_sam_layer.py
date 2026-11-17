"""The shared SAM layer (PLAN SAM-01)."""

import shutil
import subprocess

import pytest

from viralscan.sam import cigar_ops, int_tag, tag_value, view_fields, view_lines

FIELDS = [
    "r1",
    "0",
    "chr1",
    "10",
    "60",
    "5M2N3M",
    "*",
    "0",
    "0",
    "ACGTACGT",
    "IIIIIIII",
    "NH:i:2",
    "AS:i:-7",
    "CB:Z:AAAC",
]


def test_cigar_ops():
    assert cigar_ops("5M2N3M") == [(5, "M"), (2, "N"), (3, "M")]
    assert cigar_ops("*") == []


def test_int_tag_reads_negative_and_missing():
    assert int_tag(FIELDS, "AS") == -7
    assert int_tag(FIELDS, "NM") is None
    assert int_tag(FIELDS[:11] + ["AS:i:x"], "AS") is None
    with pytest.raises(ValueError):
        int_tag(FIELDS[:11] + ["AS:i:x"], "AS", strict=True)


def test_tag_value():
    assert tag_value(FIELDS, "CB:Z:") == "AAAC"
    assert tag_value(FIELDS, "UB:Z:") is None


def test_view_failure_raises_called_process_error(tmp_path):
    missing = tmp_path / "none.bam"
    with pytest.raises(subprocess.CalledProcessError):
        list(view_lines(missing))


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed")
def test_view_fields_streams_records(tmp_path):
    sam = tmp_path / "a.sam"
    sam.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n" + "\t".join(FIELDS) + "\n")
    rows = list(view_fields(sam))
    assert rows[0][0] == "r1"
    assert int_tag(rows[0], "AS") == -7
