"""Tests for external ViralScan data download/cache handling."""

from __future__ import annotations

import hashlib
import json
import shutil
import tarfile
import zipfile
from pathlib import Path

import pytest
import requests

from viralscan import data_fetch


def _zip_gtfs(path: Path, files: dict[str, str]) -> None:
    with zipfile.ZipFile(path, "w") as zf:
        for name, content in files.items():
            zf.writestr(name, content)


def _zip_members(path: Path, files: dict[str, str]) -> None:
    """Write a zip archive with arbitrary member names and text content."""
    with zipfile.ZipFile(path, "w") as zf:
        for name, content in files.items():
            zf.writestr(name, content)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_manifest(data_dir: Path, gtf_count: int | None = None) -> None:
    gtfs = sorted(data_dir.glob("*.gtf"))
    (data_dir / data_fetch.MANIFEST_NAME).write_text(
        json.dumps(
            {
                "doi": data_fetch.VIRAL_DATA_DOI,
                "gtf_count": len(gtfs) if gtf_count is None else gtf_count,
                "files": {gtf.name: _sha256(gtf) for gtf in gtfs},
            }
        )
    )


class TestEnsureViralData:
    def test_missing_cache_raises_fetch_instruction(self, tmp_path: Path) -> None:
        with pytest.raises(data_fetch.ViralScanDataError, match="viralscan data fetch"):
            data_fetch.ensure_viral_data(tmp_path)

    def test_existing_gtf_cache_is_returned(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "virus.gtf").write_text('NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "V";\n')
        _write_manifest(data_dir)

        assert data_fetch.ensure_viral_data(tmp_path) == data_dir

    def test_gtf_without_manifest_is_not_valid_cache(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "partial.gtf").write_text('NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "V";\n')

        with pytest.raises(data_fetch.ViralScanDataError, match="viralscan data fetch"):
            data_fetch.ensure_viral_data(tmp_path)

    def test_manifest_count_mismatch_is_not_valid_cache(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "partial.gtf").write_text('NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "V";\n')
        _write_manifest(data_dir, gtf_count=2)

        with pytest.raises(data_fetch.ViralScanDataError, match="viralscan data fetch"):
            data_fetch.ensure_viral_data(tmp_path)

    def test_legacy_manifest_without_file_hashes_is_accepted(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "legacy.gtf").write_text('NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "V";\n')
        (data_dir / data_fetch.MANIFEST_NAME).write_text(
            json.dumps({"doi": data_fetch.VIRAL_DATA_DOI, "gtf_count": 1})
        )

        assert data_fetch.ensure_viral_data(tmp_path) == data_dir


class TestZenodoFileSelection:
    def test_download_link_preferred_over_self_metadata_link(self) -> None:
        url, checksum, key = data_fetch._select_zenodo_file(
            {
                "files": [
                    {
                        "key": "viralscan-panel.zip",
                        "checksum": "md5:abc123",
                        "links": {
                            "self": "https://zenodo.org/api/records/20112332/files/file-id",
                            "download": "https://zenodo.org/api/records/20112332/files/file-id/content",
                        },
                    }
                ]
            }
        )

        assert url.endswith("/content")
        assert "api/records" in url
        assert checksum == "md5:abc123"
        assert key == "viralscan-panel.zip"

    def test_missing_download_link_falls_back_to_content_url(self) -> None:
        url, _, key = data_fetch._select_zenodo_file(
            {"files": [{"key": "viral panel.zip", "links": {"self": "metadata-url"}}]}
        )

        assert url == f"{data_fetch.VIRAL_DATA_ARCHIVE_URL}/viral%20panel.zip/content"
        assert key == "viral panel.zip"


class TestFetchViralData:
    def test_fetch_downloads_verifies_and_unpacks_gtfs(self, tmp_path: Path, monkeypatch) -> None:
        archive = tmp_path / "panel.zip"
        _zip_gtfs(
            archive,
            {
                "nested/a.gtf": 'NC_A\t.\tgene\t1\t10\t.\t+\t.\tgene_id "A";\n',
                "nested/b.gtf": 'NC_B\t.\tgene\t1\t10\t.\t+\t.\tgene_id "B";\n',
                "README.txt": "not a gtf",
            },
        )
        digest = _sha256(archive)

        def fake_download(url: str, destination: Path) -> None:
            assert url == "https://example.org/panel.zip"
            shutil.copyfile(archive, destination)

        monkeypatch.setattr(data_fetch, "_download", fake_download)

        data_dir = data_fetch.fetch_viral_data(
            cache_dir=tmp_path / "cache",
            archive_url="https://example.org/panel.zip",
            expected_sha256=digest,
        )

        assert sorted(p.name for p in data_dir.glob("*.gtf")) == ["a.gtf", "b.gtf"]
        assert (data_dir / data_fetch.MANIFEST_NAME).exists()

    def test_existing_cache_skips_download_without_force(self, tmp_path: Path, monkeypatch) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "existing.gtf").write_text('NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "V";\n')
        _write_manifest(data_dir)

        def fail_download(url: str, destination: Path) -> None:
            raise AssertionError("download should not run when cache already exists")

        monkeypatch.setattr(data_fetch, "_download", fail_download)

        assert data_fetch.fetch_viral_data(cache_dir=tmp_path) == data_dir

    def test_incomplete_cache_triggers_download(self, tmp_path: Path, monkeypatch) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "partial.gtf").write_text('NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "OLD";\n')

        archive = tmp_path / "panel.zip"
        _zip_gtfs(archive, {"complete.gtf": 'NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "NEW";\n'})

        def fake_download(url: str, destination: Path) -> None:
            shutil.copyfile(archive, destination)

        monkeypatch.setattr(data_fetch, "_download", fake_download)

        result = data_fetch.fetch_viral_data(
            cache_dir=tmp_path,
            archive_url="https://example.org/panel.zip",
            expected_sha256=_sha256(archive),
        )

        assert result == data_dir
        assert sorted(p.name for p in data_dir.glob("*.gtf")) == ["complete.gtf"]

    def test_manifest_hash_mismatch_is_not_valid_cache(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        gtf = data_dir / "virus.gtf"
        gtf.write_text('NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "V";\n')
        _write_manifest(data_dir)

        gtf.write_text("truncated\n")

        with pytest.raises(data_fetch.ViralScanDataError, match="viralscan data fetch"):
            data_fetch.ensure_viral_data(tmp_path)

    def test_viralscan_cache_env_is_used_by_default(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setenv("VIRALSCAN_CACHE", str(tmp_path / "env-cache"))
        assert data_fetch.cache_root() == tmp_path / "env-cache"

    def test_checksum_mismatch_raises(self, tmp_path: Path, monkeypatch) -> None:
        archive = tmp_path / "panel.zip"
        _zip_gtfs(archive, {"a.gtf": 'NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "A";\n'})

        def fake_download(url: str, destination: Path) -> None:
            shutil.copyfile(archive, destination)

        monkeypatch.setattr(data_fetch, "_download", fake_download)

        with pytest.raises(data_fetch.ViralScanDataError, match="Checksum mismatch"):
            data_fetch.fetch_viral_data(
                cache_dir=tmp_path / "cache",
                archive_url="https://example.org/panel.zip",
                expected_sha256="0" * 64,
            )

    def test_archive_without_gtfs_raises(self, tmp_path: Path, monkeypatch) -> None:
        archive = tmp_path / "panel.zip"
        _zip_gtfs(archive, {"README.txt": "empty panel"})

        def fake_download(url: str, destination: Path) -> None:
            shutil.copyfile(archive, destination)

        monkeypatch.setattr(data_fetch, "_download", fake_download)

        with pytest.raises(data_fetch.ViralScanDataError, match="No .gtf files"):
            data_fetch.fetch_viral_data(
                cache_dir=tmp_path / "cache",
                archive_url="https://example.org/panel.zip",
                expected_sha256=_sha256(archive),
            )

    def test_download_request_error_raises_data_error(self, tmp_path: Path, monkeypatch) -> None:
        def fail_get(*args, **kwargs):
            raise requests.Timeout("timed out")

        monkeypatch.setattr(data_fetch.requests, "get", fail_get)

        with pytest.raises(data_fetch.ViralScanDataError, match="Failed to download"):
            data_fetch.fetch_viral_data(
                cache_dir=tmp_path / "cache",
                archive_url="https://example.org/panel.zip",
            )


# ---------------------------------------------------------------------------
# C.1 — _extract_members (GTF + FASTA + TSV extraction)
# ---------------------------------------------------------------------------


class TestExtractMembers:
    _GTF = 'NC_A\t.\tgene\t1\t10\t.\t+\t.\tgene_id "A";\n'
    _FASTA = ">NC_A\nATCGATCG\n"
    _TSV = "accession\tgenus\nNC_A\tAlphatorquevirus\n"

    def test_zip_extracts_gtf_fasta_tsv(self, tmp_path: Path) -> None:
        archive = tmp_path / "panel.zip"
        _zip_members(
            archive,
            {
                "nested/a.gtf": self._GTF,
                "nested/anellovirus.fa": self._FASTA,
                "nested/anellovirus_accessions.tsv": self._TSV,
                "README.txt": "ignored",
            },
        )
        dest = tmp_path / "out"
        counts = data_fetch._extract_members(archive, dest)

        assert counts == {"gtf": 1, "fasta": 1, "tsv": 1}
        assert (dest / "a.gtf").read_text() == self._GTF
        assert (dest / "anellovirus.fa").read_text() == self._FASTA
        assert (dest / "anellovirus_accessions.tsv").read_text() == self._TSV

    def test_zip_extracts_fasta_suffix(self, tmp_path: Path) -> None:
        archive = tmp_path / "panel.zip"
        _zip_members(archive, {"a.gtf": self._GTF, "genomes.fasta": self._FASTA})
        counts = data_fetch._extract_members(archive, tmp_path / "out")
        assert counts["fasta"] == 1

    def test_tar_gz_extracts_members(self, tmp_path: Path) -> None:
        archive = tmp_path / "panel.tar.gz"
        with tarfile.open(archive, "w:gz") as tf:
            for name, content in [
                ("a.gtf", self._GTF),
                ("anellovirus.fa", self._FASTA),
                ("anellovirus_accessions.tsv", self._TSV),
            ]:
                import io
                data = content.encode()
                info = tarfile.TarInfo(name=name)
                info.size = len(data)
                tf.addfile(info, io.BytesIO(data))

        dest = tmp_path / "out"
        counts = data_fetch._extract_members(archive, dest)
        assert counts == {"gtf": 1, "fasta": 1, "tsv": 1}

    def test_non_matching_files_are_ignored(self, tmp_path: Path) -> None:
        archive = tmp_path / "panel.zip"
        _zip_members(archive, {"a.gtf": self._GTF, "notes.txt": "skip", "data.csv": "skip"})
        counts = data_fetch._extract_members(archive, tmp_path / "out")
        assert counts == {"gtf": 1, "fasta": 0, "tsv": 0}

    def test_tsv_with_different_name_is_ignored(self, tmp_path: Path) -> None:
        archive = tmp_path / "panel.zip"
        _zip_members(archive, {"a.gtf": self._GTF, "other_accessions.tsv": self._TSV})
        counts = data_fetch._extract_members(archive, tmp_path / "out")
        assert counts["tsv"] == 0


# ---------------------------------------------------------------------------
# C.2 — cache_valid with extended manifest (FASTA / TSV checksums)
# ---------------------------------------------------------------------------


class TestCacheValidExtended:
    _GTF = 'NC_A\t.\tgene\t1\t10\t.\t+\t.\tgene_id "A";\n'
    _FASTA = ">NC_A\nATCGATCG\n"
    _TSV = "accession\tgenus\nNC_A\tAlphatorquevirus\n"

    def _setup_cache(self, data_dir: Path) -> None:
        data_dir.mkdir(parents=True)
        (data_dir / "virus.gtf").write_text(self._GTF)
        (data_dir / "anellovirus.fa").write_text(self._FASTA)
        (data_dir / data_fetch._AUX_TSV_NAME).write_text(self._TSV)

    def _write_manifest(self, data_dir: Path) -> None:
        gtfs = sorted(data_dir.glob("*.gtf"))
        fa = data_dir / "anellovirus.fa"
        tsv = data_dir / data_fetch._AUX_TSV_NAME
        manifest = {
            "doi": data_fetch.VIRAL_DATA_DOI,
            "gtf_count": len(gtfs),
            "files": {g.name: _sha256(g) for g in gtfs},
            "fasta": fa.name,
            "fasta_checksum": _sha256(fa),
            "tsv": data_fetch._AUX_TSV_NAME,
            "tsv_checksum": _sha256(tsv),
        }
        (data_dir / data_fetch.MANIFEST_NAME).write_text(json.dumps(manifest))

    def test_valid_extended_manifest_accepted(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        self._setup_cache(data_dir)
        self._write_manifest(data_dir)
        assert data_fetch.cache_valid(tmp_path)

    def test_fasta_checksum_mismatch_invalidates_cache(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        self._setup_cache(data_dir)
        self._write_manifest(data_dir)
        (data_dir / "anellovirus.fa").write_text(">NC_A\nGGGGGGGG\n")
        assert not data_fetch.cache_valid(tmp_path)

    def test_tsv_checksum_mismatch_invalidates_cache(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        self._setup_cache(data_dir)
        self._write_manifest(data_dir)
        (data_dir / data_fetch._AUX_TSV_NAME).write_text("modified\n")
        assert not data_fetch.cache_valid(tmp_path)

    def test_fasta_file_missing_invalidates_cache(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        self._setup_cache(data_dir)
        self._write_manifest(data_dir)
        (data_dir / "anellovirus.fa").unlink()
        assert not data_fetch.cache_valid(tmp_path)

    def test_manifest_without_fasta_fields_still_valid(self, tmp_path: Path) -> None:
        """Old-style manifests without fasta/tsv keys are accepted (back-compat)."""
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "virus.gtf").write_text(self._GTF)
        gtfs = sorted(data_dir.glob("*.gtf"))
        manifest = {
            "doi": data_fetch.VIRAL_DATA_DOI,
            "gtf_count": len(gtfs),
            "files": {g.name: _sha256(g) for g in gtfs},
        }
        (data_dir / data_fetch.MANIFEST_NAME).write_text(json.dumps(manifest))
        assert data_fetch.cache_valid(tmp_path)


# ---------------------------------------------------------------------------
# C.3 — bundled_anellovirus_fasta accessor
# ---------------------------------------------------------------------------


class TestBundledAnellovirusFasta:
    _GTF = 'NC_A\t.\tgene\t1\t10\t.\t+\t.\tgene_id "A";\n'
    _FASTA = ">NC_A\nATCGATCG\n"

    def _write_manifest(self, data_dir: Path, fasta_name: str | None) -> None:
        (data_dir / "virus.gtf").write_text(self._GTF)
        gtfs = sorted(data_dir.glob("*.gtf"))
        manifest: dict = {
            "doi": data_fetch.VIRAL_DATA_DOI,
            "gtf_count": len(gtfs),
            "files": {g.name: _sha256(g) for g in gtfs},
        }
        if fasta_name:
            manifest["fasta"] = fasta_name
            manifest["fasta_checksum"] = _sha256(data_dir / fasta_name)
        (data_dir / data_fetch.MANIFEST_NAME).write_text(json.dumps(manifest))

    def test_returns_path_when_fasta_in_manifest(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "anellovirus.fa").write_text(self._FASTA)
        self._write_manifest(data_dir, "anellovirus.fa")

        result = data_fetch.bundled_anellovirus_fasta(tmp_path)
        assert result == data_dir / "anellovirus.fa"
        assert result.exists()

    def test_raises_when_no_manifest(self, tmp_path: Path) -> None:
        with pytest.raises(data_fetch.ViralScanDataError, match="viralscan data fetch"):
            data_fetch.bundled_anellovirus_fasta(tmp_path)

    def test_raises_when_manifest_has_no_fasta_key(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        self._write_manifest(data_dir, None)

        with pytest.raises(data_fetch.ViralScanDataError, match="does not include an anellovirus FASTA"):
            data_fetch.bundled_anellovirus_fasta(tmp_path)

    def test_raises_when_fasta_file_missing_from_disk(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "anellovirus.fa").write_text(self._FASTA)
        self._write_manifest(data_dir, "anellovirus.fa")
        (data_dir / "anellovirus.fa").unlink()

        with pytest.raises(data_fetch.ViralScanDataError, match="not found at"):
            data_fetch.bundled_anellovirus_fasta(tmp_path)

    def test_fetch_viral_data_populates_fasta_in_manifest(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        """fetch_viral_data writes fasta + tsv fields to the manifest when archive has them."""
        archive = tmp_path / "panel.zip"
        _zip_members(
            archive,
            {
                "virus.gtf": 'NC_A\t.\tgene\t1\t10\t.\t+\t.\tgene_id "A";\n',
                "anellovirus.fa": ">NC_A\nATCG\n",
                "anellovirus_accessions.tsv": "accession\tgenus\nNC_A\tAlpha\n",
            },
        )
        digest = _sha256(archive)

        def fake_download(url: str, destination: Path) -> None:
            shutil.copyfile(archive, destination)

        monkeypatch.setattr(data_fetch, "_download", fake_download)

        data_dir = data_fetch.fetch_viral_data(
            cache_dir=tmp_path / "cache",
            archive_url="https://example.org/panel.zip",
            expected_sha256=digest,
        )

        manifest = json.loads((data_dir / data_fetch.MANIFEST_NAME).read_text())
        assert manifest.get("fasta") == "anellovirus.fa"
        assert "fasta_checksum" in manifest
        assert manifest.get("tsv") == data_fetch._AUX_TSV_NAME
        assert "tsv_checksum" in manifest
        assert (data_dir / "anellovirus.fa").exists()
        assert (data_dir / data_fetch._AUX_TSV_NAME).exists()

        # bundled_anellovirus_fasta should now work
        result = data_fetch.bundled_anellovirus_fasta(tmp_path / "cache")
        assert result == data_dir / "anellovirus.fa"
