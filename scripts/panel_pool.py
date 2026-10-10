#!/usr/bin/env python3
"""PANEL-01 (WP2 pool): Virus-Host DB viruses with a human host, against the catalogue and candidates (gap check).

    python scripts/panel_pool.py --vhdb <virushostdb.tsv> [--candidates analysis/panel_expansion/candidates.tsv]
        [--out analysis/panel_expansion/pool_vhdb.tsv]

`virushostdb.tsv` is the table from https://www.genome.jp/virushostdb/ (not committed, 17 MB). A virus counts as
human-hosted when `host tax id` is 9606. Each row is marked by where its RefSeq/GenBank ids already are:
  panel      an id is a `shipped` catalogue row
  candidate  an id is a not-excluded row of candidates.tsv (WP1)
  excluded   an id is an excluded candidate row (twin, alias, ...)
  catalogue  an id is a catalogue row with another `panel` value (see `panel`)
  gap        no id is known to either; these are what the user reviews, none is added automatically

`tier` is a hint for review order, not a verdict: `disease` (the database names a disease), `refseq_evidence`
(host taken from RefSeq/UniProt), `literature_only` (a paper reports the infection; includes spill-over and
serology, e.g. avian paramyxoviruses, bat adenoviruses). ICTV VMR adds genus/family but no host field, so it is not
used here; `family`/`genus` come from the database's own lineage.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from viralscan.virus_catalog import load_catalogue  # noqa: E402

HUMAN_TAXID = "9606"
COLUMNS = [
    "virus_tax_id",
    "virus_name",
    "family",
    "genus",
    "accessions",
    "status",
    "panel",
    "tier",
    "disease",
    "evidence",
]
STATUS_ORDER = ["panel", "candidate", "excluded", "catalogue", "gap"]


def base_accession(accession: str) -> str:
    return accession.strip().upper().split(".")[0]


def lineage_names(lineage: str) -> tuple[str, str]:
    """(family, genus) from a taxonomy lineage: the first '*viridae' token; the last one-word '*virus' token."""
    tokens = [t.strip() for t in lineage.split(";")]
    family = next((t for t in tokens if t.endswith("viridae")), "")
    genus = next((t for t in reversed(tokens) if " " not in t and t.endswith("virus")), "")
    return family, genus


def classify(
    accessions: list[str], panel_of: dict[str, str], cand_excluded: dict[str, bool]
) -> tuple[str, str]:
    """(status, panel value) of the best-placed accession; STATUS_ORDER decides which wins."""
    found: dict[str, str] = {}
    for a in accessions:
        if panel_of.get(a) == "shipped":
            found.setdefault("panel", "shipped")
        elif a in cand_excluded:
            found.setdefault("excluded" if cand_excluded[a] else "candidate", panel_of.get(a, ""))
        elif a in panel_of:
            found.setdefault("catalogue", panel_of[a])
    status = next((s for s in STATUS_ORDER if s in found), "gap")
    return status, found.get(status, "")


def tier(row: dict) -> str:
    if row["DISEASE"]:
        return "disease"
    if re.search(r"RefSeq|UniProt", row["evidence"]):
        return "refseq_evidence"
    return "literature_only"


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1] / "analysis" / "panel_expansion"
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--vhdb", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, default=root / "candidates.tsv")
    ap.add_argument("--out", type=Path, default=root / "pool_vhdb.tsv")
    args = ap.parse_args(argv)

    panel_of = {base_accession(r["accession_version"]): r["panel"] for r in load_catalogue()}
    with args.candidates.open() as handle:
        cand_excluded = {
            base_accession(r["accession"]): bool(r["exclusion_reason"])
            for r in csv.DictReader(handle, delimiter="\t")
        }
    with args.vhdb.open() as handle:
        human = [
            r for r in csv.DictReader(handle, delimiter="\t") if r["host tax id"] == HUMAN_TAXID
        ]

    out = []
    for r in human:
        accessions = [base_accession(a) for a in re.split(r"[,\s]+", r["refseq id"]) if a]
        status, panel = classify(accessions, panel_of, cand_excluded)
        family, genus = lineage_names(r["virus lineage"])
        out.append(
            {
                "virus_tax_id": r["virus tax id"],
                "virus_name": r["virus name"],
                "family": family,
                "genus": genus,
                "accessions": ",".join(accessions),
                "status": status,
                "panel": panel,
                "tier": tier(r),
                "disease": r["DISEASE"],
                "evidence": r["evidence"],
            }
        )
    out.sort(key=lambda o: (STATUS_ORDER.index(o["status"]), o["family"], o["virus_name"].lower()))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)
    counts = {s: sum(o["status"] == s for o in out) for s in STATUS_ORDER}
    gap_tiers = {
        t: sum(o["status"] == "gap" and o["tier"] == t for o in out)
        for t in ("disease", "refseq_evidence", "literature_only")
    }
    print(f"{len(out)} human-host viruses: {counts}; gap by tier: {gap_tiers}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
