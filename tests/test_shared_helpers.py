"""One definition each for the shared helpers (PLAN HELPERS-01)."""

import gzip
from pathlib import Path

from viralscan import evidence, utils
from viralscan.star_run import plain_whitelist


def test_evidence_reexports_the_shared_helpers() -> None:
    assert evidence._open_maybe_gzip is utils.open_maybe_gzip
    assert evidence.to_int is utils.to_int
    assert evidence.to_float is utils.to_float


def test_open_maybe_gzip_reads_plain_and_gz_with_kwargs(tmp_path: Path) -> None:
    (tmp_path / "a.txt").write_text("x\r\ny\n")
    with gzip.open(tmp_path / "a.txt.gz", "wt", newline="") as fh:
        fh.write("x\r\ny\n")
    for name in ("a.txt", "a.txt.gz"):
        with utils.open_maybe_gzip(tmp_path / name, newline="", encoding="utf-8") as fh:
            assert fh.read() == "x\r\ny\n"


def test_plain_whitelist_name(tmp_path: Path) -> None:
    gz = tmp_path / "wl.txt.gz"
    gz.write_bytes(gzip.compress(b"AAAA\n"))
    out = plain_whitelist(str(gz), tmp_path, name="onlist.txt")
    assert out == str(tmp_path / "onlist.txt")
    assert (tmp_path / "onlist.txt").read_bytes() == b"AAAA\n"
