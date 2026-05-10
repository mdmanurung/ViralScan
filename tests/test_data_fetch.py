"""Tests for external ViralScan data download/cache handling."""

from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from pathlib import Path

import pytest
import requests

from viralscan import data_fetch


def _zip_gtfs(path: Path, files: dict[str, str]) -> None:
    with zipfile.ZipFile(path, "w") as zf:
        for name, content in files.items():
            zf.writestr(name, content)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_manifest(data_dir: Path, gtf_count: int) -> None:
    (data_dir / data_fetch.MANIFEST_NAME).write_text(
        json.dumps({"doi": data_fetch.VIRAL_DATA_DOI, "gtf_count": gtf_count})
    )


class TestEnsureViralData:
    def test_missing_cache_raises_fetch_instruction(self, tmp_path: Path) -> None:
        with pytest.raises(data_fetch.ViralScanDataError, match="viralscan data fetch"):
            data_fetch.ensure_viral_data(tmp_path)

    def test_existing_gtf_cache_is_returned(self, tmp_path: Path) -> None:
        data_dir = data_fetch.viral_data_dir(tmp_path)
        data_dir.mkdir(parents=True)
        (data_dir / "virus.gtf").write_text('NC\t.\tgene\t1\t10\t.\t+\t.\tgene_id "V";\n')
        _write_manifest(data_dir, gtf_count=1)

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
        _write_manifest(data_dir, gtf_count=1)

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
