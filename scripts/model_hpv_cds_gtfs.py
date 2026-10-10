#!/usr/bin/env python3
"""Model core-gene CDS for the two thin HPV records (PLAN CAT-39).

HPV-8 ``M12737.1`` has no CDS upstream (its GTF is one whole-genome pseudo-transcript) and HPV-71
``NC_039089.1`` annotates only E1. Each missing gene is modelled from an annotated sibling
(HPV-5 ``NC_001531.1`` for HPV-8, HPV-16 ``NC_001526.4`` for HPV-71): tblastn of the sibling
protein onto the target genome, then the longest ATG ORF in the hit's frame (stop to stop,
leftmost ATG), the convention of ``model_nocds_gtfs.py``. Only unspliced core genes (E6, E7, E1,
E2, L2, L1) are modelled; spliced E1^E4 and the variable E5 are not.

A model is kept only if its protein is 70-130 % of the sibling's length and the HSP covers at
least 50 % of the sibling protein; every attempt, kept or not, is written to the QC table.

Needs ``tblastn`` and ``makeblastdb`` (BLAST+) on PATH and the cached NCBI records
(``~/.cache/viralscan/ncbi/<acc>/<acc>.fasta``).
Usage: python scripts/model_hpv_cds_gtfs.py [--data-dir DIR] [--cache DIR] [--qc OUT.tsv]
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

CORE_GENES = ("E6", "E7", "E1", "E2", "L2", "L1")
MIN_LEN_RATIO, MAX_LEN_RATIO, MIN_COVERAGE = 0.7, 1.3, 0.5
STOPS = {"TAA", "TAG", "TGA"}
CODONS = dict(
    zip(
        (a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
)


@dataclass(frozen=True)
class Target:
    accession: str  # versioned, as cached
    gtf_stem: str  # file in the data dir
    sibling: str  # versioned
    sibling_gtf_stem: str
    genes: tuple[str, ...]  # genes to model
    replace_file: bool  # True: the existing GTF has no usable gene and is rewritten


TARGETS = (
    Target(
        "M12737.1",
        "Human_papillomavirus_8_M12737",
        "NC_001531.1",
        "Human_papillomavirus_5_NC_001531",
        CORE_GENES,
        True,
    ),
    Target(
        "NC_039089.1",
        "Human_papillomavirus_71_NC_039089",
        "NC_001526.4",
        "Human_papillomavirus_16_18_NC_001526",
        ("E6", "E7", "E2", "L2", "L1"),  # E1 is annotated upstream
        False,
    ),
)


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    ).upper()


def translate(dna: str) -> str:
    return "".join(CODONS.get(dna[i : i + 3], "X") for i in range(0, len(dna) - 2, 3))


def sibling_proteins(gtf: Path, genome: str) -> dict[str, str]:
    """Protein of each unspliced core gene (single CDS/exon row per transcript), by ``gene``.

    Older RefSeq GTFs in the corpus carry ``exon`` rows only, newer ones ``CDS`` rows.
    """
    rows: dict[str, list[tuple[str, int, int, str]]] = {}
    for line in gtf.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) < 9 or cols[2] not in {"CDS", "exon"} or cols[6] != "+":
            continue
        gene, tx = (
            re.search(r'gene "([^"]+)"', cols[8]),
            re.search(r'transcript_id "([^"]+)"', cols[8]),
        )
        if gene and tx and gene.group(1) in CORE_GENES:
            rows.setdefault(tx.group(1), []).append(
                (cols[2], int(cols[3]), int(cols[4]), gene.group(1))
            )
    best: dict[str, str] = {}
    for transcript_rows in rows.values():
        kinds = {kind for kind, *_ in transcript_rows}
        use = [r for r in transcript_rows if r[0] == ("CDS" if "CDS" in kinds else "exon")]
        if len(use) != 1:
            continue  # spliced product (E1^E4, E6*): not a single ORF
        _, start, end, gene = use[0]
        protein = translate(genome[start - 1 : end]).rstrip("*")
        if "*" not in protein and len(protein) > len(best.get(gene, "")):
            best[gene] = protein
    return best


def best_hsps(queries: dict[str, str], genome_path: Path) -> dict[str, list[str]]:
    """tblastn of each protein on the genome; best-bitscore HSP row per gene."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        query = tmp_dir / "q.faa"
        query.write_text("".join(f">{g}\n{p}\n" for g, p in queries.items()))
        db = tmp_dir / "db"
        subprocess.run(
            ["makeblastdb", "-in", str(genome_path), "-dbtype", "nucl", "-out", str(db)],
            check=True,
            capture_output=True,
        )
        out = subprocess.run(
            [
                "tblastn",
                "-query",
                str(query),
                "-db",
                str(db),
                "-evalue",
                "1e-5",
                "-outfmt",
                "6 qseqid sstart send pident length qstart qend evalue bitscore",
                "-seg",
                "no",
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    best: dict[str, list[str]] = {}
    for line in out.splitlines():
        row = line.split("\t")
        if row[0] not in best or float(row[8]) > float(best[row[0]][8]):
            best[row[0]] = row
    return best


def orf_around(genome: str, hit_start: int, hit_end: int) -> tuple[int, int] | None:
    """1-based inclusive (start, end incl. stop) of the longest ATG ORF in the hit's frame.

    The ORF runs from the last in-frame stop before the hit start to the first one at or after it;
    its start is the leftmost ATG inside that stop-to-stop window. tblastn may score through a
    stop into a neighbouring frame's gene, so the hit end is not used to place the stop.
    """
    first = hit_start - 1  # 0-based, in frame with every codon below
    left = first % 3
    for i in range(first - 3, -1, -3):
        if genome[i : i + 3] in STOPS:
            left = i + 3
            break
    right = next(
        (i + 3 for i in range(first, len(genome) - 2, 3) if genome[i : i + 3] in STOPS), None
    )
    if right is None:
        return None
    for i in range(left, right, 3):
        if genome[i : i + 3] == "ATG":
            return i + 1, right
    return None


def model_gene(
    gene: str, sibling_protein: str, hsp: list[str] | None, genome: str
) -> dict[str, object]:
    """One QC row: the modelled span, or the reason it was not kept."""
    row: dict[str, object] = {"gene": gene, "status": "no_hit", "sibling_aa": len(sibling_protein)}
    if hsp is None:
        return row
    sstart, send = int(hsp[1]), int(hsp[2])
    row.update(pident=float(hsp[3]), hsp_aa=int(hsp[4]), evalue=hsp[7], bitscore=float(hsp[8]))
    if sstart > send:
        row["status"] = "minus_strand_hit"
        return row
    span = orf_around(genome, sstart, send)
    if span is None:
        row["status"] = "no_orf"
        return row
    start, end = span
    aa = (end - start + 1) // 3 - 1
    row.update(start=start, end=end, aa=aa, ratio=round(aa / len(sibling_protein), 3))
    # aligned sibling residues that fall inside the ORF (an HSP can run on past its stop)
    qstart, qend = int(hsp[5]), int(hsp[6])
    in_orf = min(qend, qstart - 1 + (end - 3 - (sstart - 1)) // 3) - qstart + 1
    coverage = max(in_orf, 0) / len(sibling_protein)
    row["coverage"] = round(coverage, 3)
    if not MIN_LEN_RATIO <= aa / len(sibling_protein) <= MAX_LEN_RATIO:
        row["status"] = "length_mismatch"
    elif coverage < MIN_COVERAGE:
        row["status"] = "low_coverage"
    else:
        row["status"] = "modelled"
    return row


def gtf_rows(acc: str, source: str, gene: str, start: int, end: int, note: str) -> list[str]:
    gid = f"{acc.split('.')[0]}_{gene}_modelled"
    base = (
        f'gene_id "{gid}"; transcript_id "{gid}_t1"; gene_name "{gene}"; gene "{gene}"; '
        f'gene_biotype "protein_coding"; note "{note}"; n_exons "1";'
    )
    return [
        f"{acc}\t{source}\t{feature}\t{start}\t{end}\t.\t+\t.\t{base}"
        + (' exon_number "1";' if feature == "exon" else "")
        for feature in ("gene", "transcript", "exon")
    ]


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=root / "src" / "viralscan" / "data")
    parser.add_argument("--cache", type=Path, default=Path.home() / ".cache" / "viralscan" / "ncbi")
    parser.add_argument(
        "--qc", type=Path, default=root / "analysis" / "panel_expansion" / "hpv_modelled_cds_qc.tsv"
    )
    parser.add_argument("--dry-run", action="store_true", help="print QC only; write no GTF")
    args = parser.parse_args()

    qc_rows: list[dict[str, object]] = []
    for target in TARGETS:
        genome_path = args.cache / target.accession / f"{target.accession}.fasta"
        genome = read_fasta(genome_path)
        sib_genome = read_fasta(args.cache / target.sibling / f"{target.sibling}.fasta")
        proteins = sibling_proteins(args.data_dir / f"{target.sibling_gtf_stem}.gtf", sib_genome)
        queries = {g: proteins[g] for g in target.genes if g in proteins}
        hsps = best_hsps(queries, genome_path)
        new_rows: list[str] = []
        for gene in target.genes:
            if gene not in queries:
                qc_rows.append(
                    {"target": target.accession, "gene": gene, "status": "sibling_gene_missing"}
                )
                continue
            row = model_gene(gene, queries[gene], hsps.get(gene), genome)
            row["target"], row["sibling"] = target.accession, target.sibling
            qc_rows.append(row)
            if row["status"] == "modelled":
                note = (
                    f"modelled CDS: GenBank has none; longest ORF {row['aa']} aa, "
                    f"{row['pident']:.0f} % identity (tblastn) to {target.sibling} {gene}"
                )
                new_rows += gtf_rows(
                    target.accession, "ViralScan", gene, int(row["start"]), int(row["end"]), note
                )
        path = args.data_dir / f"{target.gtf_stem}.gtf"
        if new_rows and not args.dry_run:
            old = [] if target.replace_file else path.read_text().splitlines()
            path.write_text("\n".join(old + new_rows) + "\n")
            print("wrote", path, f"({len(new_rows) // 3} modelled genes)")

    fields = [
        "target", "sibling", "gene", "status", "sibling_aa", "aa", "ratio", "coverage", "pident",
        "hsp_aa", "start", "end", "evalue", "bitscore",
    ]  # fmt: skip
    if not args.dry_run:
        args.qc.parent.mkdir(parents=True, exist_ok=True)
        with args.qc.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", restval="")
            writer.writeheader()
            writer.writerows(qc_rows)
    for row in qc_rows:
        print(
            {
                k: row.get(k)
                for k in (
                    "target",
                    "gene",
                    "status",
                    "start",
                    "end",
                    "aa",
                    "ratio",
                    "pident",
                    "coverage",
                )
            }
        )


if __name__ == "__main__":
    main()
