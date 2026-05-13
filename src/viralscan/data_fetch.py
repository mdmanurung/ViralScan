"""Download and locate ViralScan's external viral annotation panel."""

from __future__ import annotations

import hashlib
import json
import shutil
import tarfile
import tempfile
import zipfile
from pathlib import Path
from typing import Any
from urllib.parse import quote

import requests

VIRAL_DATA_DOI = "10.5281/zenodo.20112332"
VIRAL_DATA_RECORD_ID = "20112332"
VIRAL_DATA_RECORD_URL = f"https://zenodo.org/api/records/{VIRAL_DATA_RECORD_ID}"
VIRAL_DATA_ARCHIVE_URL = f"https://zenodo.org/records/{VIRAL_DATA_RECORD_ID}/files"
MANIFEST_NAME = "manifest.json"

ARCHIVE_SUFFIXES = (".zip", ".tar", ".tar.gz", ".tgz")


class ViralScanDataError(RuntimeError):
    """Raised when ViralScan's external annotation data cannot be used."""


def cache_root(cache_dir: str | Path | None = None) -> Path:
    """Return the root cache directory for ViralScan data."""
    if cache_dir is not None:
        return Path(cache_dir).expanduser()
    return Path("~/.cache/viralscan").expanduser()


def viral_data_dir(cache_dir: str | Path | None = None) -> Path:
    """Return the directory containing fetched viral GTF annotations."""
    return cache_root(cache_dir) / "data"


def ensure_viral_data(cache_dir: str | Path | None = None) -> Path:
    """Return the viral data directory or raise with the user-facing fetch command."""
    data_dir = viral_data_dir(cache_dir)
    if cache_valid(cache_dir):
        return data_dir
    raise ViralScanDataError(
        "Viral reference annotations were not found at "
        f"{data_dir}. Run `viralscan data fetch` before using the bundled "
        "viral reference panel, or pass custom annotations with -gtf."
    )


def _checksum(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _verify_checksum(path: Path, expected: str) -> None:
    if ":" in expected:
        algorithm, digest = expected.split(":", 1)
    else:
        algorithm, digest = "sha256", expected
    actual = _checksum(path, algorithm.lower())
    if actual.lower() != digest.lower():
        raise ViralScanDataError(
            f"Checksum mismatch for {path.name}: expected {algorithm}:{digest}, "
            f"got {algorithm}:{actual}."
        )


def _download(url: str, destination: Path) -> None:
    headers = {"User-Agent": "ViralScan data fetch"}
    try:
        with requests.get(url, headers=headers, stream=True, timeout=60) as response:
            response.raise_for_status()
            with destination.open("wb") as fh:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fh.write(chunk)
    except requests.RequestException as exc:
        raise ViralScanDataError(
            f"Failed to download viral data archive from {url}: {exc}"
        ) from exc


def _select_zenodo_file(record: dict[str, Any]) -> tuple[str, str | None, str]:
    files = record.get("files")
    if not isinstance(files, list) or not files:
        raise ViralScanDataError(f"Zenodo record {VIRAL_DATA_DOI} does not list any files.")

    candidates: list[dict[str, Any]] = []
    for item in files:
        key = str(item.get("key", ""))
        if key.endswith(ARCHIVE_SUFFIXES):
            candidates.append(item)
    if not candidates:
        candidates = [files[0]]

    chosen = candidates[0]
    key = str(chosen.get("key", "viralscan-data.zip"))
    raw_links = chosen.get("links")
    links: dict[str, Any] = raw_links if isinstance(raw_links, dict) else {}
    url = str(links.get("download") or links.get("content") or "")
    if not url:
        url = f"{VIRAL_DATA_ARCHIVE_URL}/{quote(key)}/content"
    checksum = chosen.get("checksum")
    return url, str(checksum) if checksum else None, key


def _zenodo_archive() -> tuple[str, str | None, str]:
    headers = {"User-Agent": "ViralScan data fetch"}
    try:
        response = requests.get(VIRAL_DATA_RECORD_URL, headers=headers, timeout=30)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise ViralScanDataError(
            f"Failed to read Zenodo metadata for {VIRAL_DATA_DOI}: {exc}"
        ) from exc
    return _select_zenodo_file(response.json())


def cache_valid(cache_dir: str | Path | None = None) -> bool:
    """Return True when the cache manifest matches the GTF files present."""
    data_dir = viral_data_dir(cache_dir)
    manifest_path = data_dir / MANIFEST_NAME
    if not manifest_path.exists():
        return False
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    gtfs = list(data_dir.glob("*.gtf"))
    return (
        manifest.get("doi") == VIRAL_DATA_DOI
        and isinstance(manifest.get("gtf_count"), int)
        and manifest["gtf_count"] > 0
        and manifest["gtf_count"] == len(gtfs)
    )


def _extract_gtfs(archive: Path, destination: Path) -> int:
    destination.mkdir(parents=True, exist_ok=True)
    count = 0
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as zf:
            for info in zf.infolist():
                if info.is_dir() or not info.filename.endswith(".gtf"):
                    continue
                target = destination / Path(info.filename).name
                with zf.open(info) as src, target.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
                count += 1
        return count

    if tarfile.is_tarfile(archive):
        with tarfile.open(archive) as tf:
            for member in tf.getmembers():
                if not member.isfile() or not member.name.endswith(".gtf"):
                    continue
                tar_src = tf.extractfile(member)
                if tar_src is None:
                    continue
                target = destination / Path(member.name).name
                with tar_src, target.open("wb") as dst:
                    shutil.copyfileobj(tar_src, dst)
                count += 1
        return count

    raise ViralScanDataError(f"Unsupported viral data archive format: {archive.name}")


def fetch_viral_data(
    cache_dir: str | Path | None = None,
    archive_url: str | None = None,
    expected_sha256: str | None = None,
    force: bool = False,
) -> Path:
    """Download, verify, and unpack the Zenodo viral GTF archive."""
    data_dir = viral_data_dir(cache_dir)
    if not force and cache_valid(cache_dir):
        return data_dir

    root = cache_root(cache_dir)
    root.mkdir(parents=True, exist_ok=True)

    zenodo_checksum: str | None = None
    archive_name = "viralscan-data.zip"
    if archive_url is None:
        archive_url, zenodo_checksum, archive_name = _zenodo_archive()
    archive_path = root / archive_name
    _download(archive_url, archive_path)

    if zenodo_checksum:
        _verify_checksum(archive_path, zenodo_checksum)
    if expected_sha256:
        _verify_checksum(archive_path, f"sha256:{expected_sha256}")

    with tempfile.TemporaryDirectory(prefix="viralscan-data-", dir=str(root)) as tmp:
        extracted_dir = Path(tmp) / "data"
        count = _extract_gtfs(archive_path, extracted_dir)
        if count == 0:
            raise ViralScanDataError(f"No .gtf files found in {archive_path.name}.")

        data_dir.mkdir(parents=True, exist_ok=True)
        for old in data_dir.glob("*.gtf"):
            old.unlink()
        for gtf in extracted_dir.glob("*.gtf"):
            shutil.move(str(gtf), data_dir / gtf.name)

    manifest = {
        "doi": VIRAL_DATA_DOI,
        "record_id": VIRAL_DATA_RECORD_ID,
        "archive_url": archive_url,
        "archive_checksum": zenodo_checksum,
        "sha256": _checksum(archive_path, "sha256"),
        "gtf_count": len(list(data_dir.glob("*.gtf"))),
    }
    (data_dir / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2) + "\n")
    return data_dir
