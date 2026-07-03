#!/usr/bin/env python3
"""Fetch canonical FASTQs for the reference-strategy benchmark."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import urllib.request
from pathlib import Path


ENA_REPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={srr}&result=read_run"
    "&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
)


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ena_fastqs(srr: str) -> list[dict[str, str]]:
    with urllib.request.urlopen(ENA_REPORT.format(srr=srr), timeout=60) as response:
        text = response.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"expected one ENA row for {srr}, found {len(rows)}")
    row = rows[0]
    urls = row["fastq_ftp"].split(";")
    md5s = row["fastq_md5"].split(";")
    sizes = row["fastq_bytes"].split(";")
    if not (len(urls) == len(md5s) == len(sizes) == 2):
        raise ValueError(f"expected paired FASTQs for {srr}, got {row}")
    return [
        {"url": "https://" + url, "md5": md5, "bytes": size}
        for url, md5, size in zip(urls, md5s, sizes)
    ]


def fetch(srr: str, out_dir: Path, *, force: bool = False) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    records = ena_fastqs(srr)
    metadata_path = out_dir / "source_urls.tsv"
    with metadata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["mate", "url", "md5", "bytes"], delimiter="\t")
        writer.writeheader()
        for idx, record in enumerate(records, 1):
            writer.writerow({"mate": f"R{idx}", **record})

    for idx, record in enumerate(records, 1):
        out = out_dir / f"{srr}_{idx}.fastq.gz"
        if out.exists() and not force:
            observed = _md5(out)
            if observed == record["md5"]:
                print(f"{out} already present and MD5 matches")
                continue
            raise ValueError(f"{out} exists but MD5 {observed} != expected {record['md5']}")
        tmp = out.with_suffix(out.suffix + ".tmp")
        print(f"Downloading {record['url']} -> {out}", flush=True)
        urllib.request.urlretrieve(record["url"], tmp)
        observed = _md5(tmp)
        if observed != record["md5"]:
            tmp.unlink(missing_ok=True)
            raise ValueError(f"{out} MD5 {observed} != expected {record['md5']}")
        tmp.replace(out)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--srr", default="SRR12682296")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("benchmark_inputs/reference_strategy/SRR12682296"),
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    try:
        fetch(args.srr, args.out_dir, force=args.force)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
