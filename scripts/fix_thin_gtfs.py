#!/usr/bin/env python3
"""Drop genome-spanning pseudo-transcripts and fix one malformed feature column (PLAN CAT-38, REF-13).

HAV ``NC_001489.1`` and HTLV-2 ``NC_001488.1`` carry a RefSeq ``misc_RNA`` / ``prim_transcript`` that
spans 94-100 % of the genome beside their real CDS. A transcript that covers the whole genome shares
every k-mer with every CDS transcript, which is the equivalence-class collapse measured for
anelloviruses (ANELLO-12). Those transcripts are removed and the gene row is re-spanned to the
remaining transcripts; the real CDS stay. Human astrovirus ``NC_001943.1`` writes the feature type
``Non structural gene`` in column 3 for two CDS rows; they become ``CDS``.

Idempotent. Usage: python scripts/fix_thin_gtfs.py [--data-dir DIR]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

SPANNING_TYPES = {"misc_RNA", "prim_transcript"}
MIN_SPAN_FRACTION = 0.9
# (gtf stem, sequence length): lengths of the RefSeq records NC_001489.1 and NC_001488.1
GENOMES = {
    "Hepatitis_A_virus_NC_001489": 7478,
    "Human_T-lymphotropic_virus_2_NC_001488": 8952,
}
ASTROVIRUS = "Human_astrovirus_NC_001943"


def attrs(column: str) -> dict[str, str]:
    return dict(re.findall(r'(\w+) "([^"]*)"', column))


def drop_spanning(lines: list[str], genome_length: int) -> list[str]:
    """Remove whole-genome misc_RNA/prim_transcript transcripts; re-span their genes."""
    rows = [line.split("\t") for line in lines]
    doomed: set[str] = set()
    for cols in rows:
        if len(cols) < 9 or cols[2] != "transcript":
            continue
        a = attrs(cols[8])
        spans = (int(cols[4]) - int(cols[3]) + 1) >= MIN_SPAN_FRACTION * genome_length
        if spans and a.get("gbkey") in SPANNING_TYPES:
            doomed.add(a["transcript_id"])
    kept = [c for c in rows if len(c) < 9 or attrs(c[8]).get("transcript_id") not in doomed]
    extent: dict[str, list[int]] = {}
    for cols in kept:
        if len(cols) >= 9 and cols[2] in {"transcript", "exon", "CDS"}:
            lo_hi = extent.setdefault(attrs(cols[8])["gene_id"], [10**12, 0])
            lo_hi[0], lo_hi[1] = min(lo_hi[0], int(cols[3])), max(lo_hi[1], int(cols[4]))
    out = []
    for cols in kept:
        if len(cols) >= 9 and cols[2] == "gene":
            gene = attrs(cols[8]).get("gene_id", "")
            if gene not in extent:
                continue  # every transcript of this gene was the spanning one
            cols[3], cols[4] = (str(v) for v in extent[gene])
        out.append("\t".join(cols))
    return out


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=root / "src" / "viralscan" / "data")
    args = parser.parse_args()
    for stem, length in GENOMES.items():
        path = args.data_dir / f"{stem}.gtf"
        before = path.read_text().splitlines()
        after = drop_spanning(before, length)
        path.write_text("\n".join(after) + "\n")
        print(f"{path.name}: {len(before)} -> {len(after)} rows")
    path = args.data_dir / f"{ASTROVIRUS}.gtf"
    text = path.read_text()
    fixed = text.replace("\tNon structural gene\t", "\tCDS\t")
    path.write_text(fixed)
    print(f"{path.name}: {text.count(chr(9) + 'Non structural gene' + chr(9))} feature types fixed")


if __name__ == "__main__":
    main()
