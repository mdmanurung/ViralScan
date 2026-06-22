#!/usr/bin/env python3
"""Build anellovirus_accessions.tsv — one-time curation script.

Reconciles two accession sources into a single, deduplicated TSV:
  1. clareaulab/anellovirus_reference `simple_anello_metadata_V2.csv`:
     ~2,200 CD-HIT representative human Anelloviridae genomes drawn from
     NCBI Virus (~5,500 complete genomes, host == human).
     Taxonomic classification via ``infer_genus`` column.
     Citation: Lareau et al., clareaulab/anellovirus_reference (GitHub, 2025).
  2. ViralScan's existing 20 bundled RefSeq anellovirus GTFs in
     ``src/viralscan/data/``.

Usage (run from repo root, point CSV at the downloaded metadata file):

    python extras/build_anello_table.py \\
        --csv /path/to/simple_anello_metadata_V2.csv \\
        --gtf-dir src/viralscan/data \\
        --out src/viralscan/data/anellovirus_accessions.tsv

The output TSV has columns:
    accession   virus_name   genus   family   source
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# Families that mark a row as Anelloviridae-related
_ANELLO_FAMILY = "Anelloviridae"

# Genus value used when infer_genus is empty or "Unclassified"
_FALLBACK_GENUS = "Anelloviridae"

# Pattern matching filenames of the bundled anellovirus GTFs
_ANELLO_GTF_PATTERN = re.compile(
    r"(?i)(torque|alphatorque|betatorque|gammatorque|anello)",
)


def _strip_version(accession: str) -> str:
    """Return accession without trailing .N version suffix (e.g. NC_002076.2 → NC_002076)."""
    return re.sub(r"\.\d+$", "", accession)


def _genus_from_row(row: dict[str, str]) -> str:
    """Resolve genus from a clareaulab CSV row.

    Prefer ``infer_genus`` over ``Genus`` because many NCBI-deposited records
    have ``Genus = "Unclassified"`` even when the inferred phylogenetic genus
    is known.  Fall back to family level when both are absent/unclassified.
    """
    g = row.get("infer_genus", "").strip()
    if g and g.lower() not in ("unclassified", ""):
        return g
    g = row.get("Genus", "").strip()
    if g and g.lower() not in ("unclassified", ""):
        return g
    return _FALLBACK_GENUS


def _virus_name_from_row(row: dict[str, str]) -> str:
    """Human-readable virus name from a clareaulab CSV row."""
    name = row.get("Virus Name", "").strip()
    if name:
        return name
    species = row.get("Species", "").strip()
    if species:
        return species
    return row.get("Accession", "unknown")


def load_clareaulab(csv_path: Path) -> list[dict[str, str]]:
    """Read the clareaulab metadata CSV; return representative rows only."""
    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    reps: list[dict[str, str]] = []
    for row in rows:
        acc = row.get("Accession", "").strip()
        rep = row.get("cdhit_representative", "").strip()
        if not acc:
            continue
        # A row is a CD-HIT representative when its own accession == representative field
        if acc == rep:
            reps.append(
                {
                    "accession": acc,
                    "virus_name": _virus_name_from_row(row),
                    "genus": _genus_from_row(row),
                    "family": row.get("Family", _ANELLO_FAMILY).strip() or _ANELLO_FAMILY,
                    "source": "clareaulab",
                }
            )
    return reps


def load_viralscan_gtfs(gtf_dir: Path) -> list[dict[str, str]]:
    """Extract accessions from the bundled anellovirus GTFs in *gtf_dir*."""
    results: list[dict[str, str]] = []
    for gtf_path in sorted(gtf_dir.glob("*.gtf")):
        if not _ANELLO_GTF_PATTERN.search(gtf_path.name):
            continue
        # Collect distinct seqid values (column 1) — these are the accessions
        accessions: list[str] = []
        with open(gtf_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 9:
                    continue
                acc = cols[0].strip()
                if acc and acc not in accessions:
                    accessions.append(acc)
        for acc in accessions:
            # Derive a best-effort virus name from the filename:
            # first strip any trailing "_NC_XXXXXX" accession suffix, then
            # convert remaining underscores to spaces.
            stem = re.sub(r"_NC_\d+(\.\d+)?$", "", gtf_path.stem)
            stem = stem.replace("_", " ")
            results.append(
                {
                    "accession": acc,
                    "virus_name": stem.strip(),
                    "genus": "Alphatorquevirus",  # all bundled entries are Alphatorquevirus
                    "family": _ANELLO_FAMILY,
                    "source": "viralscan-refseq",
                }
            )
    return results


def reconcile(
    clareaulab_rows: list[dict[str, str]],
    viralscan_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Union + dedup accessions; prefer NC_* RefSeq when versions clash.

    Dedup key is the bare accession (version-stripped).  When the same bare
    accession appears in both sources, the ViralScan (NC_*) RefSeq record wins
    because it has been curated for quality and gene annotation.
    """
    seen: dict[str, dict[str, str]] = {}  # bare_acc → row

    # Add clareaulab first (lower priority)
    for row in clareaulab_rows:
        bare = _strip_version(row["accession"])
        if bare not in seen:
            seen[bare] = row

    # Add ViralScan RefSeq — overwrites clareaulab if same bare accession
    for row in viralscan_rows:
        bare = _strip_version(row["accession"])
        seen[bare] = row  # RefSeq wins

    # Sort: NC_* first, then by accession
    def _sort_key(row: dict[str, str]) -> tuple[int, str]:
        return (0 if row["accession"].startswith("NC_") else 1, row["accession"])

    return sorted(seen.values(), key=_sort_key)


def validate(rows: list[dict[str, str]]) -> None:
    """Assert no duplicate accessions and every row has a non-empty genus."""
    accs = [r["accession"] for r in rows]
    dupes = [a for a in set(accs) if accs.count(a) > 1]
    if dupes:
        raise ValueError(f"Duplicate accessions: {dupes}")
    for row in rows:
        if not row.get("genus"):
            raise ValueError(f"Missing genus for accession {row['accession']}")


def write_tsv(rows: list[dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["accession", "virus_name", "genus", "family", "source"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Build anellovirus_accessions.tsv from clareaulab CSV + ViralScan GTFs.",
    )
    parser.add_argument("--csv", required=True, help="Path to simple_anello_metadata_V2.csv")
    parser.add_argument(
        "--gtf-dir",
        default="src/viralscan/data",
        help="Directory containing ViralScan bundled GTFs (default: src/viralscan/data)",
    )
    parser.add_argument(
        "--out",
        default="src/viralscan/data/anellovirus_accessions.tsv",
        help="Output TSV path",
    )
    args = parser.parse_args(argv)

    csv_path = Path(args.csv)
    gtf_dir = Path(args.gtf_dir)
    out_path = Path(args.out)

    print(f"Loading clareaulab representatives from {csv_path} …")
    clareaulab_rows = load_clareaulab(csv_path)
    print(f"  {len(clareaulab_rows)} CD-HIT representative genomes")

    print(f"Loading ViralScan GTF accessions from {gtf_dir} …")
    viralscan_rows = load_viralscan_gtfs(gtf_dir)
    print(f"  {len(viralscan_rows)} accessions from bundled GTFs")

    print("Reconciling (union + dedup) …")
    unified = reconcile(clareaulab_rows, viralscan_rows)
    validate(unified)

    from collections import Counter

    genus_counts = Counter(r["genus"] for r in unified)
    source_counts = Counter(r["source"] for r in unified)
    print(f"  Total unique accessions: {len(unified)}")
    print(f"  Sources: {dict(source_counts)}")
    print("  Genus breakdown:")
    for genus, count in sorted(genus_counts.items(), key=lambda x: -x[1]):
        print(f"    {count:5d}  {genus}")

    write_tsv(unified, out_path)
    print(f"\nWrote {len(unified)} rows → {out_path}")


if __name__ == "__main__":
    main(sys.argv[1:])
