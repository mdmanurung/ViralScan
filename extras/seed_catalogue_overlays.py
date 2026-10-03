#!/usr/bin/env python3
"""Seed the curation-overlay columns of ``virus_catalog.tsv`` (PLAN ``MECH-A``).

The Virus Identity table reads four per-accession decisions from the catalogue.
Today each lives in code, keyed by a *name string*:

``common_name``
    The display name existing outputs, tests and vignettes use
    ("Epstein-Barr virus", "Human herpesvirus 1"). It is taken from the legacy
    prefix map (``constants.VIRUS_NAME_MAP`` + ``VIRUS_GENE_ID_ALIASES``) applied
    to each bundled GTF's gene IDs, so every accession that has a legacy name
    keeps it.
``sibling_group``
    Near-identical genomes whose shared-k-mer mass EM or the equal split can
    move between them (``SIBLING_VIRUS_PAIRS``, extended with EBV type 1/2 per
    F-017). Seeded by NCBI taxid, not by name.
``risk_class``
    ``eve`` for families/genera in ``EVE_RISK_GENERA`` (empty since ANELLO-PRIOR.4),
    ``low_complexity`` for ``LOW_COMPLEXITY_RISK_GENERA``, replacing the substring
    test on virus names.
``role``
    ``decoy`` for lab-contaminant decoys (e.g. the max panel's MLV/XMRV set).
``panel``
    ``shipped`` for every accession the bundled builder emits (bundled GTF
    seqnames plus the anellovirus accession table) or that was catalogued before
    the max-panel merge (``--shipped-accessions``); ``max`` for the rest. The
    CAT-31 reconciliation guard treats only ``shipped`` rows as detection
    targets of the shipped panel.

Only empty cells are filled; a regeneration of the catalogue carries every
overlay forward (``build_virus_catalog.OVERLAY_COLUMNS``), so a human decision
is never overwritten.

Usage
-----
    python extras/seed_catalogue_overlays.py \\
        [--catalogue src/viralscan/data/virus_catalog.tsv] \\
        [--gtf-dir src/viralscan/data] [--decoys candidates.tsv]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from viralscan.constants import (  # noqa: E402
    EVE_RISK_GENERA,
    LOW_COMPLEXITY_RISK_GENERA,
    VIRUS_NAME_MAP,
)
from viralscan.virus_grouping import virus_name_for_gene  # noqa: E402

DEFAULT_CATALOGUE = REPO_ROOT / "src" / "viralscan" / "data" / "virus_catalog.tsv"
DEFAULT_GTF_DIR = REPO_ROOT / "src" / "viralscan" / "data"

#: NCBI taxid -> sibling group. Membership requires near-identity (shared k-mer
#: mass that allocation can move), which is why EBV/KSHV are *not* siblings.
SIBLING_GROUPS_BY_TAXID: dict[str, str] = {
    "10298": "HSV",  # Human alphaherpesvirus 1
    "10310": "HSV",  # Human alphaherpesvirus 2
    "32603": "HHV-6",  # Human betaherpesvirus 6A
    "32604": "HHV-6",  # Human betaherpesvirus 6B
    "10376": "HHV-4",  # Human gammaherpesvirus 4 (EBV type 1)
    "12509": "HHV-4",  # Human herpesvirus 4 type 2 (EBV type 2), F-017
}

_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def legacy_names_by_accession(gtf_dir: Path) -> dict[str, str]:
    """Majority legacy display name of each bundled GTF seqname."""
    votes: dict[str, Counter[str]] = defaultdict(Counter)
    for gtf in sorted(gtf_dir.glob("*.gtf")):
        with open(gtf, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                match = _GENE_ID_RE.search(cols[8])
                if not match:
                    continue
                gene_id = match.group(1)
                name = virus_name_for_gene(gene_id, VIRUS_NAME_MAP)
                if name != gene_id:
                    votes[cols[0]][name] += 1
    return {acc: counter.most_common(1)[0][0] for acc, counter in votes.items()}


def shipped_accessions(gtf_dir: Path, extra: Path | None) -> set[str]:
    """Base accessions the bundled builder emits, plus an optional prior list."""
    from viralscan.anellovirus import load_accession_table

    shipped = {row["accession"].strip().split(".")[0] for row in load_accession_table()}
    for gtf in gtf_dir.glob("*.gtf"):
        with open(gtf, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith("#") and "\t" in line:
                    shipped.add(line.split("\t", 1)[0].split(".")[0])
    if extra is not None:
        shipped |= {
            line.strip().split(".")[0] for line in extra.read_text().splitlines() if line.strip()
        }
    return shipped


def decoy_accessions(path: Path | None) -> set[str]:
    if path is None:
        return set()
    with open(path, newline="") as handle:
        return {
            row["accession_version"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("role") == "decoy"
        }


def seed(
    rows: list[dict[str, str]],
    legacy: dict[str, str],
    decoys: set[str],
    shipped: set[str],
) -> Counter[str]:
    filled: Counter[str] = Counter()
    for row in rows:
        acc = row["accession_version"]
        if not row.get("panel"):
            row["panel"] = "shipped" if row["accession"] in shipped else "max"
            filled[f"panel={row['panel']}"] += 1
        name = legacy.get(acc) or legacy.get(row["accession"])
        if name and not row.get("common_name"):
            row["common_name"] = name
            filled["common_name"] += 1
        group = SIBLING_GROUPS_BY_TAXID.get(row.get("taxid", ""))
        if group and not row.get("sibling_group"):
            row["sibling_group"] = group
            filled["sibling_group"] += 1
        for genera, risk in ((EVE_RISK_GENERA, "eve"), (LOW_COMPLEXITY_RISK_GENERA, "low_complexity")):
            hit = row.get("family") in genera or row.get("genus") in genera
            if hit and not row.get("risk_class"):
                row["risk_class"] = risk
                filled["risk_class"] += 1
        if acc in decoys and not row.get("role"):
            row["role"] = "decoy"
            filled["role"] += 1
    return filled


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE)
    parser.add_argument("--gtf-dir", type=Path, default=DEFAULT_GTF_DIR)
    parser.add_argument("--decoys", type=Path, help="TSV with accession_version and role columns")
    parser.add_argument(
        "--shipped-accessions",
        type=Path,
        help="Extra accessions to scope as `shipped` (e.g. the pre-merge catalogue)",
    )
    args = parser.parse_args(argv)

    with open(args.catalogue, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    missing = [
        c for c in ("taxid", "common_name", "sibling_group", "role", "panel") if c not in fieldnames
    ]
    if missing:
        parser.error(f"catalogue lacks {missing}; regenerate it with build_virus_catalog.py first")

    filled = seed(
        rows,
        legacy_names_by_accession(args.gtf_dir),
        decoy_accessions(args.decoys),
        shipped_accessions(args.gtf_dir, args.shipped_accessions),
    )
    with open(args.catalogue, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"seeded {args.catalogue}: "
        + ", ".join(f"{k}={v}" for k, v in sorted(filled.items()))
        + f" over {len(rows)} rows",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
