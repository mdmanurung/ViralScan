#!/usr/bin/env python3
"""Build the deterministic viral-annotation-panel archive for Zenodo (PLAN REF-11, preparation only).

Writes four files to ``--out-dir``; nothing is uploaded and no identifier is reserved:

* ``viralscan-viral-panel-<version>.zip``  flat archive in the layout ``viralscan data fetch``
  unpacks (``*.gtf``, at most one ``*.fa``, ``anellovirus_accessions.tsv``) plus ``PANEL_MANIFEST.json``
* ``PANEL_MANIFEST.json``                   per-file sha256/size (same bytes as inside the zip)
* ``SHA256SUMS``                            archive + manifest checksums for the Zenodo upload page
* ``zenodo_metadata.draft.json``            DRAFT deposit metadata; every ``TODO:`` value is unknown here

Deterministic: members sorted by name, stored (no compression, so no zlib-version dependence), fixed
timestamp and mode, no clock or host values in any output.  The same panel directory gives the same archive
bytes.  Before building, the GTFs are checked against ``analysis/panel_expansion/gtf_corpus_manifest.tsv``
so a stale local corpus cannot be archived by accident (``--no-corpus-check`` to skip).

    PYTHONPATH=src python scripts/make_panel_archive.py --panel-dir <dir with the GTFs> \\
        --out-dir <empty output dir> --version 3.0.0rc1 --source-commit <main sha>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
AUX_TSV = "anellovirus_accessions.tsv"
MANIFEST_NAME = "PANEL_MANIFEST.json"  # not "manifest.json": that name is the fetch cache's own file
FIXED_DATE = (2020, 1, 1, 0, 0, 0)  # zip epoch floor is 1980; any constant works
TODO = "TODO: "


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collect(panel_dir: Path, fasta: Path | None, aux_tsv: Path) -> dict[str, Path]:
    """Flat name -> source path, sorted.  Rejects duplicate basenames and more than one FASTA."""
    gtfs = sorted(panel_dir.glob("*.gtf"))
    if not gtfs:
        raise SystemExit(f"no *.gtf in {panel_dir}")
    members = {p.name: p for p in gtfs}
    extra = [aux_tsv] + ([fasta] if fasta else [])
    for p in extra:
        if not p.is_file():
            raise SystemExit(f"missing input file: {p}")
        if p.name in members:
            raise SystemExit(f"duplicate archive member name: {p.name}")
        members[p.name] = p
    # fetch_viral_data records only fastas[0], so a second FASTA would be silently ignored
    if fasta and not fasta.name.endswith((".fa", ".fasta")):
        raise SystemExit(f"--fasta must end in .fa or .fasta: {fasta.name}")
    return dict(sorted(members.items()))


def check_corpus(members: dict[str, Path], manifest: Path) -> list[str]:
    """GTF names and sha256 must equal the committed corpus manifest (file<TAB>sha256<TAB>size...)."""
    want = dict(line.split("\t")[:2] for line in manifest.read_text().splitlines()[1:] if line)
    got = {n: sha256(p) for n, p in members.items() if n.endswith(".gtf")}
    return (
        [f"missing from panel dir: {n}" for n in sorted(want.keys() - got.keys())]
        + [f"not in corpus manifest: {n}" for n in sorted(got.keys() - want.keys())]
        + [f"sha256 differs: {n}" for n in sorted(want.keys() & got.keys()) if want[n] != got[n]]
    )


def manifest_json(members: dict[str, Path], version: str, source_commit: str) -> bytes:
    files = {n: {"sha256": sha256(p), "size_bytes": p.stat().st_size} for n, p in members.items()}
    doc = {
        "schema_version": 1,
        "panel_version": version,
        "source_commit": source_commit,
        "gtf_count": sum(n.endswith(".gtf") for n in files),
        "files": files,
        "generator": "scripts/make_panel_archive.py",
    }
    return (json.dumps(doc, indent=2, sort_keys=True) + "\n").encode()


def write_zip(dest: Path, members: dict[str, Path], manifest: bytes) -> None:
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_STORED) as zf:
        items = [(n, p.read_bytes()) for n, p in members.items()] + [(MANIFEST_NAME, manifest)]
        for name, data in sorted(items):
            info = zipfile.ZipInfo(name, FIXED_DATE)
            info.compress_type = zipfile.ZIP_STORED
            info.external_attr = 0o644 << 16
            info.create_system = 3  # unix, fixed so the host OS does not leak into the bytes
            zf.writestr(info, data)


def zenodo_draft(version: str, n_gtf: int, has_fasta: bool, source_commit: str) -> dict:
    """Deposit metadata.  Only facts derivable from the build are filled; the rest are TODO markers."""
    return {
        "_status": "DRAFT. Not reviewed. Replace every 'TODO:' value, then delete this key and _todo.",
        "_todo": [
            "creators: names, ORCIDs, affiliations and order (PLAN REL-15)",
            "license: the licence holder and an SPDX/Zenodo licence id (PLAN REL-15); GTFs derive from NCBI records",
            "related_identifiers: the pinned ViralScan commit/tag that consumes this panel",
            "version: confirm it matches the software release this panel is frozen for",
            "doi: Zenodo mints it at deposit creation; then update data_fetch.VIRAL_DATA_DOI/RECORD_ID",
        ],
        "metadata": {
            "upload_type": "dataset",
            "title": "ViralScan viral annotation panel (GTF annotations for human-associated viruses)",
            "description": (
                f"{n_gtf} GTF annotation files"
                + (", one FASTA" if has_fasta else "")
                + " and the anellovirus accession table used by the ViralScan bundled viral reference "
                "panel (`viralscan data fetch`). Per-file sha256 checksums are in PANEL_MANIFEST.json. "
                f"Built from ViralScan source commit {source_commit}. "
                + TODO
                + "describe provenance of the annotations (RefSeq/GenBank records, modelled CDS) in the "
                "authors' words."
            ),
            "version": version,
            "access_right": "open",
            "license": TODO + "licence id (REL-15)",
            "creators": [{"name": TODO + "Family, Given", "orcid": TODO + "ORCID", "affiliation": TODO + "affiliation"}],
            "keywords": ["viral reference", "GTF", "single-cell RNA-seq", "kallisto", "ViralScan"],
            "related_identifiers": [
                {
                    "relation": "isSupplementTo",
                    "identifier": "https://github.com/mdmanurung/ViralScan",
                    "resource_type": "software",
                    "scheme": "url",
                }
            ],
            "doi": TODO + "leave unset; Zenodo assigns it",
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-dir", type=Path, required=True, help="directory holding the *.gtf corpus")
    ap.add_argument("--out-dir", type=Path, required=True, help="created if missing; must not hold earlier outputs")
    ap.add_argument("--version", required=True, help="panel version, e.g. 3.0.0rc1")
    ap.add_argument("--source-commit", default=TODO + "pinned main sha", help="ViralScan commit the panel is built from")
    ap.add_argument("--fasta", type=Path, help="optional bundled anellovirus FASTA (.fa/.fasta)")
    ap.add_argument("--aux-tsv", type=Path, default=REPO / "src" / "viralscan" / "data" / AUX_TSV)
    ap.add_argument("--corpus-manifest", type=Path, default=REPO / "analysis" / "panel_expansion" / "gtf_corpus_manifest.tsv")
    ap.add_argument("--no-corpus-check", action="store_true")
    a = ap.parse_args(argv)

    members = collect(a.panel_dir, a.fasta, a.aux_tsv)
    if not a.no_corpus_check:
        bad = check_corpus(members, a.corpus_manifest)
        if bad:
            print("\n".join(bad), file=sys.stderr)
            print(f"{len(bad)} corpus mismatch(es); fix the panel dir or pass --no-corpus-check", file=sys.stderr)
            return 1

    a.out_dir.mkdir(parents=True, exist_ok=True)
    archive = a.out_dir / f"viralscan-viral-panel-{a.version}.zip"
    outputs = [archive, a.out_dir / MANIFEST_NAME, a.out_dir / "SHA256SUMS", a.out_dir / "zenodo_metadata.draft.json"]
    if any(p.exists() for p in outputs):
        print(f"refusing to overwrite earlier outputs in {a.out_dir}", file=sys.stderr)
        return 1

    manifest = manifest_json(members, a.version, a.source_commit)
    write_zip(archive, members, manifest)
    (a.out_dir / MANIFEST_NAME).write_bytes(manifest)
    (a.out_dir / "SHA256SUMS").write_text(
        f"{sha256(archive)}  {archive.name}\n{sha256(a.out_dir / MANIFEST_NAME)}  {MANIFEST_NAME}\n"
    )
    n_gtf = sum(n.endswith(".gtf") for n in members)
    draft = zenodo_draft(a.version, n_gtf, a.fasta is not None, a.source_commit)
    (a.out_dir / "zenodo_metadata.draft.json").write_text(json.dumps(draft, indent=2) + "\n")
    print(f"{archive}  {n_gtf} GTFs  sha256 {sha256(archive)}")
    print("DRAFT metadata has TODO placeholders; nothing was uploaded.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
