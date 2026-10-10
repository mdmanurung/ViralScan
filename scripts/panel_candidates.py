#!/usr/bin/env python3
"""PANEL-01 (WP1): candidate table for promoting human-relevant viruses into the shipped panel (read-only).

    python scripts/panel_candidates.py [--evonk analysis/panel_expansion/evonk_candidates.tsv]
        [--curated analysis/panel_expansion/human_relevant_curated.tsv]
        [--census analysis/panel_expansion/census.tsv]
        [--kmer-sharing analysis/panel_expansion/kmer_sharing.tsv] [--kmer-partners .../kmer_partners.tsv] [--panel-fasta <build>/viral.fa] [--candidate-fasta other.fasta ...]
        [--out analysis/panel_expansion/candidates.tsv]

Candidates are the catalogue rows whose `panel` is `max`, `broad` or `legacy` (listed, never shipped) plus the
accessions of `--evonk` that the catalogue does not hold. Nothing is added to the panel here: the table is what the
user reviews row by row. Output is deterministic (sorted, no timestamps).

`relevance` (why a virus may belong in a human-sample panel):
  H1             the catalogue `host` field says Homo sapiens (record-level evidence)
  H2_species     no host field, but another catalogue record of the same species has H1
  H2_name        no host field; the species name starts with "Human " (name-based, to review)
  H2_curated     a `status=accepted` curated row (--curated) names the accession or species
  H3_curated     as H2_curated, class H3: a vector/reagent virus that occurs in human samples
  proposed       a `status=proposed` curated row: suggested, not yet accepted by the user
  census_host    only NCBI's host qualifier (scripts/panel_census.py) says Homo sapiens; submitter-provided and
                 known to include phage and bacterial records, so it is a reason to look, not to promote
  unreviewed     none of the above; never added automatically

`exclusion_reason` (non-empty means: do not promote, whatever `relevance` says):
  anellovirus_max_not_promoted   Anelloviridae `max` rows are redundant at 95 % ANI (design doc 2026-09-27)
  eve_risk                       catalogue risk_class is `eve`
  refseq_alias_of_genbank_hpv:<id>  an evonk HPV RefSeq twin of GenBank <id>; <id> is a candidate in its own right
  sequence_twin_of:<id>          identical (upper-cased) sequence already in the panel under another ID
  kmer_twin_of:<id>              under 5 % of its k=31 k-mers lie outside the panel (scripts/panel_kmer_sharing.py,
                                 --kmer-sharing); <id> is its top panel partner in --kmer-partners. Same-virus
                                 duplicates, not novel genomes. Run order: this script, panel_kmer_sharing.py, this
                                 script again (panel_kmer_sharing.py keeps rescoring these rows, so the result is stable)

Sequence twins are only checked for sequences present in --panel-fasta and --candidate-fasta; `twin_checked`
says which rows could be checked, so an unchecked row is never read as "no twin".
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from viralscan.virus_catalog import load_catalogue  # noqa: E402

COLUMNS = [
    "accession",
    "source",
    "species",
    "family",
    "genus",
    "host",
    "segment",
    "panel",
    "relevance",
    "relevance_basis",
    "exclusion_reason",
    "twin_checked",
    "sequence_twin_of",
]


def base_accession(accession: str) -> str:
    return accession.strip().upper().split(".")[0]


def fasta_md5(path: Path) -> dict[str, str]:
    """base accession -> md5 of the upper-cased sequence."""
    out: dict[str, str] = {}
    name, chunks = None, []

    def flush() -> None:
        if name is not None:
            out[name] = hashlib.md5("".join(chunks).upper().encode()).hexdigest()  # noqa: S324

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                flush()
                name, chunks = base_accession(line[1:].split()[0]), []
            elif line:
                chunks.append(line)
    flush()
    return out


ANELLO_NAME = re.compile(r"torque teno|anellovir|gyrovirus|torquevirus", re.I)


KMER_TWIN_BELOW = 0.05


def load_kmer_twins(sharing_path: Path | None, partners_path: Path | None) -> dict[str, str]:
    """base accession -> top panel partner, for candidates with under KMER_TWIN_BELOW of k-mers outside the panel."""
    if sharing_path is None or not sharing_path.is_file():
        return {}
    with sharing_path.open() as handle:
        low = {
            base_accession(r["accession"])
            for r in csv.DictReader(handle, delimiter="\t")
            if float(r["frac_not_in_panel"]) < KMER_TWIN_BELOW
        }
    partner: dict[str, str] = {}
    if partners_path is not None and partners_path.is_file():
        with partners_path.open() as handle:
            for r in csv.DictReader(handle, delimiter="\t"):
                base = base_accession(r["candidate"])
                if r["partner_role"] == "panel" and base not in partner:
                    partner[base] = r["partner"]  # rows are rank-ordered
    return {base: partner.get(base, "unknown") for base in low}


def load_curated(path: Path | None) -> dict[str, dict[str, str]]:
    """key (base accession or lower-cased species) -> curated row."""
    if path is None or not path.is_file():
        return {}
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        (base_accession(r["key"]) if r["key_type"] == "accession" else r["key"].strip().lower()): r
        for r in rows
    }


def relevance(
    host: str, species: str, accession: str, human_species: set[str], curated: dict
) -> tuple[str, str]:
    base = base_accession(accession)
    row = curated.get(base) or curated.get(species.strip().lower())
    if row and row["status"] == "accepted":
        return f"{row['class']}_curated", row["basis"]
    if "homo sapiens" in host.lower():
        return "H1", "catalogue host field"
    if species.strip().lower() in human_species:
        return "H2_species", "another record of this species has a Homo sapiens host"
    if species.strip().lower().startswith("human "):
        return "H2_name", "species name starts with 'Human '"
    if row and row["status"] == "proposed":
        return "proposed", row["basis"]
    return "unreviewed", ""


def build_candidates(
    catalogue: list[dict[str, str]],
    evonk: list[dict[str, str]],
    curated: dict,
    panel_md5: dict[str, str],
    candidate_md5: dict[str, str],
    census: list[dict[str, str]] | None = None,
    kmer_twins: dict[str, str] | None = None,
) -> list[dict[str, str]]:
    human_species = {
        r["species"].strip().lower() for r in catalogue if "homo sapiens" in r["host"].lower()
    }
    in_catalogue = {base_accession(r["accession"]) for r in catalogue}
    panel_by_md5: dict[str, str] = {}
    for acc, digest in sorted(panel_md5.items()):
        panel_by_md5.setdefault(digest, acc)
    rows: list[dict[str, str]] = []
    for r in catalogue:
        if r["panel"] in ("max", "broad", "legacy"):
            rows.append({**r, "source": f"catalogue-{r['panel']}"})
    for e in evonk:
        if base_accession(e["accession"]) not in in_catalogue:
            rows.append(
                {
                    "accession": e["accession"],
                    "species": e["name"],
                    "family": "",
                    "genus": "",
                    "host": "",
                    "segment": "",
                    "panel": "",
                    "risk_class": "",
                    "group": e["group"],
                    "genbank_twin": e["genbank_twin"],
                    "source": "evonk-new",
                }
            )
    seen = {base_accession(r["accession"]) for r in rows}
    for c in census or []:
        if c["in_catalogue_panel"] or base_accession(c["accession"]) in seen:
            continue  # catalogued rows and evonk rows are already candidates
        rows.append(
            {
                "accession": c["accession"],
                "species": c["organism"],
                "family": "",
                "genus": "",
                "host": "",
                "segment": "",
                "panel": "",
                "risk_class": "",
                "source": "census-new",
            }
        )
    out: list[dict[str, str]] = []
    for r in rows:
        base = base_accession(r["accession"])
        rel, basis = relevance(r["host"], r["species"], r["accession"], human_species, curated)
        if r["source"] == "census-new" and rel == "unreviewed":
            rel, basis = (
                "census_host",
                "NCBI host qualifier is Homo sapiens (submitter-provided; verify)",
            )
            if ANELLO_NAME.search(r["species"]):
                basis += "; anellovirus: the panel total is decided separately (PLAN ANELLO-15)"
        digest = candidate_md5.get(base) or panel_md5.get(base)
        twin = panel_by_md5.get(digest, "") if digest else ""
        if twin == base:
            twin = ""
        reason = ""
        if r["family"] == "Anelloviridae" and r["panel"] == "max":
            reason = "anellovirus_max_not_promoted"
        elif r.get("risk_class") == "eve":
            reason = "eve_risk"
        elif r.get("group") == "evonk-hpv-alias":
            reason = f"refseq_alias_of_genbank_hpv:{r['genbank_twin']}"
        elif twin:
            reason = f"sequence_twin_of:{twin}"
        elif base in (kmer_twins or {}):
            reason = f"kmer_twin_of:{kmer_twins[base]}"
        out.append(
            {
                "accession": r["accession"],
                "source": r["source"],
                "species": r["species"],
                "family": r["family"],
                "genus": r.get("genus", ""),
                "host": r["host"],
                "segment": r["segment"],
                "panel": r["panel"],
                "relevance": rel,
                "relevance_basis": basis,
                "exclusion_reason": reason,
                "twin_checked": "yes" if digest else "no",
                "sequence_twin_of": twin,
            }
        )
    out.sort(key=lambda row: (row["source"], row["accession"]))
    return out


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1] / "analysis" / "panel_expansion"
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--evonk", type=Path, default=root / "evonk_candidates.tsv")
    ap.add_argument("--curated", type=Path, default=root / "human_relevant_curated.tsv")
    ap.add_argument("--census", type=Path, default=root / "census.tsv")
    ap.add_argument("--kmer-sharing", type=Path, default=root / "kmer_sharing.tsv")
    ap.add_argument("--kmer-partners", type=Path, default=root / "kmer_partners.tsv")
    ap.add_argument("--panel-fasta", type=Path)
    ap.add_argument("--candidate-fasta", type=Path, action="append", default=[])
    ap.add_argument("--out", type=Path, default=root / "candidates.tsv")
    args = ap.parse_args(argv)

    panel_md5 = fasta_md5(args.panel_fasta) if args.panel_fasta else {}
    candidate_md5: dict[str, str] = {}
    for path in args.candidate_fasta:
        candidate_md5.update(fasta_md5(path))
    evonk = list(csv.DictReader(args.evonk.open(), delimiter="\t")) if args.evonk.is_file() else []
    census = (
        list(csv.DictReader(args.census.open(), delimiter="\t")) if args.census.is_file() else []
    )
    rows = build_candidates(
        load_catalogue(),
        evonk,
        load_curated(args.curated),
        panel_md5,
        candidate_md5,
        census,
        load_kmer_twins(args.kmer_sharing, args.kmer_partners),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    keep = [r for r in rows if not r["exclusion_reason"]]
    counts: dict[str, int] = {}
    for r in keep:
        counts[r["relevance"]] = counts.get(r["relevance"], 0) + 1
    print(
        f"{len(rows)} candidates, {len(keep)} not excluded: {dict(sorted(counts.items()))}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
