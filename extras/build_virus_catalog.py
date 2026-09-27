#!/usr/bin/env python3
"""Generate ``src/viralscan/data/virus_catalog.tsv`` — the frozen virus catalogue.

PLAN `CAT-02`/`CAT-03`/`CAT-09`. One row per reference accession, carrying the
taxonomy, segment and host fields that the rest of the pipeline needs to turn a
gene ID back into a virus.

Why a catalogue at all. ``virus_grouping.virus_name_for_gene`` resolves a gene
ID by prefix-matching a hand-maintained ``VIRUS_NAME_MAP``. That works for the
legacy panel's ``HUM_HERP6B_U67`` style IDs and fails for genome-scoped IDs, so
13.4 % of gene IDs on the covid index resolve to no virus at all. It also has no
concept of a segment: influenza A's eight segments are eight unrelated
accessions, so one strain surfaces as eight separate "viruses". Both problems go
away once every accession maps to a species name.

Cache-first, exactly like ``extras/build_anellovirus_genes.py``: each record is
read from the ``~/.cache/viralscan/ncbi`` flatfile when a valid ``.sha256``
sidecar exists, so a re-run costs no NCBI requests and the catalogue stays
re-derivable offline.

Usage
-----
    # Catalogue whatever the current reference contains
    python extras/build_virus_catalog.py --from-fasta <viral_genome.fa> \\
        --email you@example.org

    # Add explicit accessions (Tier 1 breadth)
    python extras/build_virus_catalog.py --accessions NC_045512.2,NC_001802.1 \\
        --email you@example.org

    # Report drift without rewriting the committed TSV
    python extras/build_virus_catalog.py --from-fasta ... --check
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from viralscan.run_safety import sha256_file  # noqa: E402
from viralscan.scripts.ncbi_fetch import (  # noqa: E402
    DEFAULT_CACHE_DIR,
    NCBIFetchError,
    _locus_fields,
    _source_qualifiers,
    fetch_genbank,
    genbank_cache_path,
)

DEFAULT_OUT = REPO_ROOT / "src" / "viralscan" / "data" / "virus_catalog.tsv"

COLUMNS = [
    "accession",
    "accession_version",
    "species",
    "genus",
    "family",
    "lineage",
    "molecule",
    "topology",
    "genome_length",
    "segment",
    "host",
    "isolate",
    "refseq",
    "source",
    "gb_sha256",
    "tier",
    "persistence_class",
    "risk_class",
    "inclusion_rationale",
]

#: A GenBank ORGANISM block is "<name>\n<lineage; ...>." — the lineage ranks are
#: positional, so the family is the token ending in "viridae" and the genus the
#: token ending in "virus" that sits last.
_FAMILY_RE = re.compile(r"\b([A-Z][A-Za-z]*viridae)\b")
_SUBFAMILY_RE = re.compile(r"\b([A-Z][A-Za-z]*virinae)\b")


def parse_organism(genbank_text: str) -> tuple[str, str]:
    """Return ``(organism, lineage)`` from the ORGANISM block."""
    organism = ""
    lineage_parts: list[str] = []
    in_block = False
    for line in genbank_text.splitlines():
        if line.startswith("  ORGANISM"):
            organism = line[len("  ORGANISM") :].strip()
            in_block = True
            continue
        if in_block:
            if line.startswith(" " * 12):
                lineage_parts.append(line.strip())
            else:
                break
    return organism, " ".join(lineage_parts).rstrip(".")


def derive_taxonomy(organism: str, lineage: str) -> tuple[str, str]:
    """Return ``(family, genus)`` from a GenBank lineage string.

    The lineage is ordered high rank → low rank, so the genus is the last
    element before the species. Ranks are not labelled in a flatfile, which is
    why this leans on the ICTV suffixes rather than on position alone.
    """
    ranks = [part.strip() for part in lineage.split(";") if part.strip()]
    family_match = _FAMILY_RE.search(lineage)
    family = family_match.group(1) if family_match else ""
    genus = ""
    for rank in reversed(ranks):
        if rank.endswith("viridae") or rank.endswith("virinae") or rank.endswith("virales"):
            break
        if _SUBFAMILY_RE.fullmatch(rank):
            continue
        # The last "…virus" token above the species is the genus.
        if rank.endswith("virus") or rank.endswith("viruses"):
            genus = rank
            break
    return family, genus


def species_name(organism: str, definition: str) -> str:
    """The display name a virus groups under.

    Prefers the ORGANISM name, which is the taxonomic species. Segment records
    share it, which is exactly what makes segment grouping fall out for free.
    """
    name = organism.strip()
    if not name:
        name = re.sub(
            r",?\s*(complete|partial)\s+(genome|sequence|cds).*$", "", definition, flags=re.I
        )
    # "Influenza A virus (A/Puerto Rico/8/1934(H1N1))" → "Influenza A virus"
    name = re.sub(r"\s*\(.*\)\s*$", "", name).strip()
    return name


def parse_definition(genbank_text: str) -> str:
    out: list[str] = []
    for line in genbank_text.splitlines():
        if line.startswith("DEFINITION"):
            out.append(line[len("DEFINITION") :].strip())
        elif out and line.startswith(" " * 12):
            out.append(line.strip())
        elif out:
            break
    return " ".join(out).rstrip(".")


def catalog_row(accession: str, genbank_text: str, gb_path: Path) -> dict[str, str]:
    organism, lineage = parse_organism(genbank_text)
    definition = parse_definition(genbank_text)
    family, genus = derive_taxonomy(organism, lineage)
    locus = _locus_fields(genbank_text)
    source = _source_qualifiers(genbank_text)
    version = str(locus.get("version") or accession)
    return {
        "accession": version.split(".")[0],
        "accession_version": version,
        "species": species_name(organism, definition),
        "genus": genus,
        "family": family,
        "lineage": lineage,
        "molecule": str(locus.get("molecule") or ""),
        "topology": str(locus.get("topology") or ""),
        "genome_length": str(locus.get("genome_length") or 0),
        "segment": source.get("segment", ""),
        "host": source.get("host", ""),
        "isolate": source.get("isolate", "") or source.get("strain", ""),
        "refseq": "true" if version.startswith(("NC_", "NG_", "NM_", "NR_")) else "false",
        "source": "ncbi-genbank",
        "gb_sha256": sha256_file(gb_path),
        # Curation overlay columns: filled by CAT-10, left empty by the fetcher
        # so a re-run never silently drops a human decision.
        "tier": "",
        "persistence_class": "",
        "risk_class": "",
        "inclusion_rationale": "",
    }


def accessions_from_fasta(path: Path) -> list[str]:
    """Accessions from a FASTA, skipping non-NCBI pseudocontig identifiers."""
    out: list[str] = []
    seen: set[str] = set()
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            token = line[1:].split()[0]
            # The legacy panel carries pseudocontigs like HUM_HERP6B_U67 that are
            # gene fragments, not accessions; they have no GenBank record.
            if not re.fullmatch(r"[A-Z]{1,2}_?\d{5,8}(\.\d+)?", token):
                continue
            if token not in seen:
                seen.add(token)
                out.append(token)
    return out


def load_existing(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with open(path, newline="") as handle:
        return {row["accession_version"]: row for row in csv.DictReader(handle, delimiter="\t")}


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--from-fasta", type=Path, action="append", default=[])
    parser.add_argument("--accessions", default="", help="Comma-separated accessions to add.")
    parser.add_argument("--accession-file", type=Path, help="One accession per line.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    parser.add_argument(
        "--cache-only",
        action="store_true",
        help="Never hit the network; skip accessions with no cached flatfile.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Report what would change without writing the TSV.",
    )
    args = parser.parse_args(argv)

    wanted: list[str] = []
    for fasta in args.from_fasta:
        wanted.extend(accessions_from_fasta(fasta))
    if args.accessions:
        wanted.extend(a.strip() for a in args.accessions.split(",") if a.strip())
    if args.accession_file:
        wanted.extend(
            line.strip() for line in args.accession_file.read_text().splitlines() if line.strip()
        )
    # Preserve first-seen order, drop duplicates.
    seen: set[str] = set()
    accessions = [a for a in wanted if not (a in seen or seen.add(a))]
    if not accessions:
        parser.error("no accessions: pass --from-fasta, --accessions or --accession-file")
    if not args.cache_only and not args.email:
        parser.error("NCBI needs a contact email: pass --email, set NCBI_EMAIL, or --cache-only")

    existing = load_existing(args.out)
    rows: list[dict[str, str]] = []
    skipped: list[str] = []
    failed: list[str] = []
    fetched = 0
    for accession in accessions:
        cached = genbank_cache_path(accession, args.cache_dir)
        try:
            if args.cache_only:
                if not cached.exists():
                    skipped.append(accession)
                    continue
                gb_path, text = cached, cached.read_text()
            else:
                had = cached.exists()
                gb_path, text = fetch_genbank(accession, args.email, args.api_key, args.cache_dir)
                if not had:
                    fetched += 1
        except (NCBIFetchError, OSError) as exc:
            failed.append(f"{accession}: {exc}")
            continue
        row = catalog_row(accession, text, gb_path)
        # Carry forward any curated overlay values rather than blanking them.
        prior = existing.get(row["accession_version"]) or existing.get(accession)
        if prior:
            for column in ("tier", "persistence_class", "risk_class", "inclusion_rationale"):
                if prior.get(column):
                    row[column] = prior[column]
        rows.append(row)

    rows.sort(key=lambda r: (r["family"], r["species"], r["accession_version"]))
    print(
        f"catalogued {len(rows)} accessions "
        f"({fetched} newly fetched, {len(skipped)} skipped, {len(failed)} failed)",
        file=sys.stderr,
    )
    for message in failed[:10]:
        print(f"  FAILED {message}", file=sys.stderr)
    if skipped:
        print(f"  skipped (no cached flatfile): {len(skipped)}", file=sys.stderr)

    if args.check:
        added = [r["accession_version"] for r in rows if r["accession_version"] not in existing]
        removed = sorted(set(existing) - {r["accession_version"] for r in rows})
        print(f"check: {len(added)} added, {len(removed)} removed", file=sys.stderr)
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    _report(rows, args.out)
    return 1 if failed else 0


def _report(rows: list[dict[str, str]], out_path: Path) -> None:
    families: dict[str, int] = {}
    species: set[str] = set()
    segmented: dict[str, int] = {}
    for row in rows:
        families[row["family"] or "(unclassified)"] = (
            families.get(row["family"] or "(unclassified)", 0) + 1
        )
        species.add(row["species"])
        if row["segment"]:
            segmented[row["species"]] = segmented.get(row["species"], 0) + 1
    print(f"\nwrote {out_path} ({len(rows)} accessions, {len(species)} species)")
    print(f"  families: {len(families)}")
    for family, count in sorted(families.items(), key=lambda kv: -kv[1])[:12]:
        print(f"    {count:5d}  {family}")
    if segmented:
        print(f"  segmented species: {len(segmented)}")
        for name, count in sorted(segmented.items(), key=lambda kv: -kv[1])[:8]:
            print(f"    {count:3d} segments  {name}")
    missing_species = [r["accession_version"] for r in rows if not r["species"]]
    if missing_species:
        print(f"  WARNING: {len(missing_species)} rows have no species name")


if __name__ == "__main__":
    sys.exit(main())
