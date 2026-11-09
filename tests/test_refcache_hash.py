"""Ensembl cache reuse must verify the sha256 recorded at fetch time (REF-CACHE-HASH-01)."""

import gzip
import json
import urllib.request

import pytest

from viralscan.scripts import build_reference as br


class _Resp:
    def __init__(self, data):
        self._data, self._pos = data, 0

    def read(self, n=-1):
        chunk = self._data[self._pos :] if n < 0 else self._data[self._pos : self._pos + n]
        self._pos += len(chunk)
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _setup(monkeypatch, tmp_path, payload):
    monkeypatch.setattr(br, "_ensembl_release", lambda: 116)
    monkeypatch.setattr(br, "_list_ensembl_files", _fake_listing)
    calls = []

    def fake_urlopen(url, timeout=0):
        calls.append(url)
        return _Resp(payload)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    return calls


def _fake_listing(species, base, retries=3):
    if "fasta" in base:
        return ["Homo_sapiens.GRCh38.cdna.all.fa.gz"]
    return ["Homo_sapiens.GRCh38.116.gtf.gz"]


def _fetch(tmp_path):
    return br.fetch_host_cdna("human", tmp_path / "out", cache_dir=tmp_path / "cache")


def test_first_fetch_records_hash_and_reuse_does_not_redownload(monkeypatch, tmp_path):
    calls = _setup(monkeypatch, tmp_path, gzip.compress(b">t\nACGT\n"))
    _fetch(tmp_path)
    assert len(calls) == 2
    _fetch(tmp_path)
    assert len(calls) == 2


def test_corrupted_cache_fails_loudly_and_is_not_replaced(monkeypatch, tmp_path):
    calls = _setup(monkeypatch, tmp_path, gzip.compress(b">t\nACGT\n"))
    _fetch(tmp_path)
    cached = next((tmp_path / "cache").rglob("*.cdna.all.fa.gz"))
    cached.write_bytes(cached.read_bytes()[:10])
    with pytest.raises(RuntimeError, match="does not match its recorded sha256"):
        _fetch(tmp_path)
    assert len(calls) == 2
    assert len(cached.read_bytes()) == 10


def test_cache_without_receipt_fails_loudly(monkeypatch, tmp_path):
    calls = _setup(monkeypatch, tmp_path, gzip.compress(b">t\nACGT\n"))
    _fetch(tmp_path)
    next((tmp_path / "cache").rglob("*.gtf.gz.retrieval.json")).unlink()
    with pytest.raises(RuntimeError, match="no recorded sha256"):
        _fetch(tmp_path)
    assert len(calls) == 2


def test_download_is_atomic_when_the_stream_breaks(monkeypatch, tmp_path):
    class Broken(_Resp):
        def read(self, n=-1):
            raise OSError("connection reset")

    monkeypatch.setattr(urllib.request, "urlopen", lambda url, timeout=0: Broken(b""))
    monkeypatch.setattr(br.time, "sleep", lambda s: None)
    dest = tmp_path / "x.gz"
    with pytest.raises(RuntimeError, match="Failed to download"):
        br._download("http://example.invalid/x.gz", dest, retries=1)
    assert not dest.exists()
    assert not dest.with_suffix(".gz.retrieval.json").exists()


def test_receipt_hash_matches_file(monkeypatch, tmp_path):
    _setup(monkeypatch, tmp_path, b"payload")
    dest = br._download("http://example.invalid/x.gz", tmp_path / "x.gz")
    receipt = json.loads((tmp_path / "x.gz.retrieval.json").read_text())
    assert receipt["sha256"] == br.sha256_file(dest)
