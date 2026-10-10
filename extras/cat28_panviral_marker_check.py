#!/usr/bin/env python3
"""Check every ``gene_programs.tsv`` marker against the panviral annotation table (PLAN CAT-28).

Validation report only: nothing in the catalogue changes. The input is the vendored
``panviral/reference/pan_virus_annotation_plain.tsv`` (not ``pan_viral_...``; the file is
gitignored under ``extras/vendor_sources/``). Its hazards, from ``VENDOR_REPORT.md`` section 5:

* it is EC-keyed, so a gene repeats (``BWRF1`` x6, ``LMP-1`` x3, ``US33A`` x3): rows are
  de-duplicated per (virus, accession, gene) before matching;
* file line 532 carries a stray ``[`` (``[NC_001699.1``, JCV ``Jvgp6``): repaired, and
  counted, because a join on accession would drop that gene silently;
* HHV7 (``U43400.1``) and HHV8 (``MK733606.1``) use a different accession from the panel's
  RefSeq (``NC_001716.2``, ``NC_009333.1``).

Each marker's ``refseq_gene`` is matched to the panviral ``Gene`` of the same virus by exact,
case-sensitive equality (as PROG-11 requires: ``BARF1`` and ``BaRF1`` are different genes).
Status:

* ``matched``: the gene is in the panviral table under the same accession as our panel;
* ``accession_mismatch``: the gene is there, but the panviral accession is not ours;
* ``unmatched``: the virus is not in the panviral table, or the gene is not under it. A near
  miss (case-insensitive, or one half of a ``A/B`` fused name) is named in ``detail`` but is
  not promoted to ``matched``.

Usage: python extras/cat28_panviral_marker_check.py [--panviral TSV] [--out TSV]
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PANVIRAL = (
    REPO_ROOT
    / "extras"
    / "vendor_sources"
    / "panviral"
    / "reference"
    / "pan_virus_annotation_plain.tsv"
)
DEFAULT_CATALOGUE = REPO_ROOT / "src" / "viralscan" / "data" / "gene_programs.tsv"
DEFAULT_OUT = REPO_ROOT / "analysis" / "panel_expansion" / "cat28_panviral_marker_check.tsv"

PANVIRAL_COLUMNS = ("EC", "Gene", "Nuccore", "Virus")

#: Catalogue virus name -> panviral ``Virus`` label. HHV-2 and HHV-6A are not in the table
#: (it holds AAV2, HHV1/3/4/5/6B/7/8 and JCV), so their markers are ``unmatched``.
PANVIRAL_VIRUS = {
    "Epstein-Barr virus": "HHV4",
    "Human cytomegalovirus": "HHV5",
    "Human herpesvirus 1": "HHV1",
    "Human herpesvirus 6b": "HHV6B",
    "Human herpesvirus 7": "HHV7",
    "Human herpesvirus 8": "HHV8",
    "Varicella-zoster virus": "HHV3",
}

#: The accession of each virus in our bundled panel (the seqname of its GTF).
OUR_ACCESSION = {
    "Epstein-Barr virus": "NC_007605.1",
    "Human cytomegalovirus": "NC_006273.2",
    "Human herpesvirus 1": "NC_001806.2",
    "Human herpesvirus 6b": "AF157706.1",
    "Human herpesvirus 7": "NC_001716.2",
    "Human herpesvirus 8": "NC_009333.1",
    "Varicella-zoster virus": "NC_001348.1",
}

OUT_COLUMNS = ("marker", "virus", "status", "detail")


def read_panviral(path: Path) -> tuple[dict[tuple[str, str], dict[str, int]], int]:
    """``({(Virus, Nuccore): {Gene: n EC rows}}, n repaired accessions)``.

    The stray bracket is stripped from ``Nuccore`` (the repair is counted, never silent);
    a row that is still not four fields is an error, not skipped.
    """
    table: dict[tuple[str, str], dict[str, int]] = {}
    repaired = 0
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if tuple(header) != PANVIRAL_COLUMNS:
            raise ValueError(f"{path}: expected columns {PANVIRAL_COLUMNS}, got {tuple(header)}")
        for lineno, row in enumerate(reader, start=2):
            if len(row) != 4:
                raise ValueError(f"{path}:{lineno}: expected 4 fields, got {len(row)}: {row}")
            _, gene, nuccore, virus = row
            clean = nuccore.strip().strip("[]").strip()
            if clean != nuccore:
                repaired += 1
            counts = table.setdefault((virus, clean), {})
            counts[gene] = counts.get(gene, 0) + 1
    return table, repaired


def _near_miss(gene: str, genes: dict[str, int]) -> str:
    folded = [g for g in genes if g.casefold() == gene.casefold()]
    if folded:
        return f"; near miss (case differs): {', '.join(sorted(folded))}"
    fused = [g for g in genes if "/" in g and gene in g.split("/")]
    if fused:
        return f"; near miss (part of a fused name): {', '.join(sorted(fused))}"
    return ""


def check_markers(
    catalogue: list[dict[str, str]], table: dict[tuple[str, str], dict[str, int]]
) -> list[dict[str, str]]:
    """One report row per catalogue marker row, in catalogue order."""
    out: list[dict[str, str]] = []
    for row in catalogue:
        virus, gene = row["virus"], row["refseq_gene"]
        where = f"gene_id={row['gene_id_bundled']}"
        label = PANVIRAL_VIRUS.get(virus)
        entries = {acc: genes for (v, acc), genes in table.items() if v == label}
        if label is None or not entries:
            status, detail = "unmatched", f"{where}; virus is not in the panviral table"
        else:
            hit = [acc for acc, genes in entries.items() if gene in genes]
            ours = OUR_ACCESSION[virus]
            if not hit:
                near = "".join(_near_miss(gene, g) for g in entries.values())
                status, detail = "unmatched", f"{where}; {gene} is not a {label} Gene{near}"
            elif ours in hit:
                n = entries[ours][gene]
                status = "matched"
                detail = f"{where}; {label} {ours}, {n} EC row(s)"
                if virus == "Human herpesvirus 6b":
                    detail += "; this accession is the AF157706.1 pseudocontig, not a genome"
            else:
                status = "accession_mismatch"
                detail = f"{where}; panviral {label} uses {hit[0]}, our panel uses {ours}"
        out.append({"marker": gene, "virus": virus, "status": status, "detail": detail})
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--panviral", type=Path, default=DEFAULT_PANVIRAL)
    parser.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    if not args.panviral.is_file():
        sys.exit(f"ERROR: panviral table not found: {args.panviral}")
    table, repaired = read_panviral(args.panviral)
    with open(args.catalogue, newline="", encoding="utf-8") as handle:
        catalogue = list(csv.DictReader(handle, delimiter="\t"))
    rows = check_markers(catalogue, table)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter(r["status"] for r in rows)
    n_rows = sum(sum(g.values()) for g in table.values())
    print(
        f"{len(rows)} markers: "
        + ", ".join(f"{counts[s]} {s}" for s in ("matched", "unmatched", "accession_mismatch"))
        + f" ({n_rows} panviral rows, {repaired} accession(s) repaired) -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
