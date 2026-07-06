"""Tests for the barcode/whitelist chemistry-mismatch preflight (SH1.4)."""

import gzip

import pytest

from viralscan.whitelist_preflight import (
    check_whitelist,
    load_whitelist,
    whitelist_match_rate,
)


def _write_fastq(path, barcodes, cb_len=16, umi_len=12, tail="ACGT"):
    """Write a minimal 4-line-per-record FASTQ; each read is barcode+UMI+tail."""
    lines = []
    for i, bc in enumerate(barcodes):
        seq = bc + ("A" * umi_len) + tail
        lines += [f"@r{i}", seq, "+", "I" * len(seq)]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _bc(n, cb_len=16):
    """Deterministic 16-mer barcode from an integer."""
    alphabet = "ACGT"
    s = ""
    x = n
    for _ in range(cb_len):
        s += alphabet[x % 4]
        x //= 4
    return s


class TestLoadWhitelist:
    def test_reads_plain(self, tmp_path):
        p = tmp_path / "wl.txt"
        p.write_text("AAAA\nCCCC\n\nGGGG\n")
        assert load_whitelist(str(p)) == {"AAAA", "CCCC", "GGGG"}

    def test_reads_gzip(self, tmp_path):
        p = tmp_path / "wl.txt.gz"
        with gzip.open(p, "wt") as fh:
            fh.write("AAAA\nCCCC\n")
        assert load_whitelist(str(p)) == {"AAAA", "CCCC"}


class TestWhitelistMatchRate:
    def test_all_match(self, tmp_path):
        bcs = [_bc(i) for i in range(50)]
        r1 = _write_fastq(tmp_path / "r1.fastq", bcs)
        rate, n = whitelist_match_rate(r1, set(bcs), cb_len=16)
        assert rate == 1.0
        assert n == 50

    def test_none_match(self, tmp_path):
        bcs = [_bc(i) for i in range(50)]
        r1 = _write_fastq(tmp_path / "r1.fastq", bcs)
        rate, n = whitelist_match_rate(r1, {_bc(1000)}, cb_len=16)
        assert rate == 0.0

    def test_partial_match(self, tmp_path):
        bcs = [_bc(i) for i in range(100)]
        r1 = _write_fastq(tmp_path / "r1.fastq", bcs)
        wl = {_bc(i) for i in range(40)}  # first 40 are whitelisted
        rate, n = whitelist_match_rate(r1, wl, cb_len=16)
        assert rate == pytest.approx(0.4)

    def test_wrong_cb_len_lowers_match(self, tmp_path):
        # Reading a 16bp barcode as 14bp (wrong technology) shifts the string and
        # should drop the match rate.
        bcs = [_bc(i) for i in range(50)]
        r1 = _write_fastq(tmp_path / "r1.fastq", bcs)
        rate16, _ = whitelist_match_rate(r1, set(bcs), cb_len=16)
        rate14, _ = whitelist_match_rate(r1, set(bcs), cb_len=14)
        assert rate16 == 1.0
        assert rate14 < 1.0

    def test_n_sample_caps_reads(self, tmp_path):
        bcs = [_bc(i) for i in range(500)]
        r1 = _write_fastq(tmp_path / "r1.fastq", bcs)
        _, n = whitelist_match_rate(r1, set(bcs), cb_len=16, n_sample=100)
        assert n == 100

    def test_empty_whitelist_raises(self, tmp_path):
        r1 = _write_fastq(tmp_path / "r1.fastq", [_bc(0)])
        with pytest.raises(ValueError, match="empty"):
            whitelist_match_rate(r1, set(), cb_len=16)


class TestCheckWhitelist:
    def test_ok_when_matching(self, tmp_path):
        bcs = [_bc(i) for i in range(60)]
        r1 = _write_fastq(tmp_path / "r1.fastq", bcs)
        res = check_whitelist(r1, set(bcs), "10xv3")
        assert res.ok
        assert res.cb_len == 16
        assert "OK" in res.message

    def test_flags_mismatch(self, tmp_path):
        bcs = [_bc(i) for i in range(60)]
        r1 = _write_fastq(tmp_path / "r1.fastq", bcs)
        res = check_whitelist(r1, {_bc(9999)}, "10xv3")
        assert not res.ok
        assert "F-005" in res.message

    def test_geometry_from_technology(self, tmp_path):
        # 10xv2 has a 16bp CB too, dropseq has 12bp — cb_len should follow.
        bcs = [_bc(i) for i in range(20)]
        r1 = _write_fastq(tmp_path / "r1.fastq", bcs)
        assert check_whitelist(r1, set(bcs), "dropseq").cb_len == 12
