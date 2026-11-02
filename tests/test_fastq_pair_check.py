"""REVIEW-FIX-H: a shorter R2 must fail the run, not silently drop reads."""

import gzip
import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "fastq_pair_check",
    Path(__file__).parents[1] / "src" / "viralscan" / "scripts" / "fastq_pair_check.py",
)
fpc = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(fpc)


def _fq(path, n, gz=False):
    text = "".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(n))
    with (gzip.open if gz else open)(path, "wt") as fh:
        fh.write(text)
    return path


def test_matching_pairs_counted_plain_and_gz(tmp_path):
    assert fpc.check_pairs(_fq(tmp_path / "a.fq", 3), _fq(tmp_path / "b.fq.gz", 3, gz=True)) == 3


def test_short_r2_fails_with_counts(tmp_path):
    with pytest.raises(ValueError, match="R1 has 2 reads, R2 has 1"):
        fpc.check_pairs(_fq(tmp_path / "a.fq", 2), _fq(tmp_path / "b.fq", 1))


def test_truncated_record_fails(tmp_path):
    bad = tmp_path / "b.fq"
    bad.write_text("@r\nACGT\n+\n")
    with pytest.raises(ValueError, match="multiple of 4"):
        fpc.check_pairs(_fq(tmp_path / "a.fq", 1), bad)
