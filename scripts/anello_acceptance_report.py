#!/usr/bin/env python3
"""Score the ANDET-09e acceptance: does alignment see what kallisto misses?

Two measurements, because they answer different halves of the question.

**Per-genome assay (the head-to-head).** The planted reads carry their own
read IDs, so they can be pulled back out of the planted arm and run on their
own. Each (genome, window) set goes through `kb count` and through the branch's
own STAR command, against the same panel. With no background, every anellovirus
molecule either method reports came from a planted read, so recovery is a plain
fraction and needs no equivalence-class decoding. Read count per set is fixed,
so the fraction is what the detection floor at any planting level follows from.

**Whole-arm assay (the context).** planted, unplanted and negative run
end-to-end through `viralscan`. The negative gives the false-positive rate on
artefact classes that no tiny assay contains, and unplanted gives the real
anellovirus background of the covid library.

Pass (pre-registered, PLAN ANDET-09e):

1. alignment recovery >= kallisto recovery for every planted genome, and
   strictly greater for any genome below 95 % identity to its nearest panel
   genome;
2. 0 alignment molecules in the negative arm;
3. no planted read lost to the host filter.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from viralscan.anello_align import (  # noqa: E402
    accession_metrics,
    align_cmd,
    iter_fasta,
    parse_sam_line,
    virus_molecules,
)
from viralscan.scripts.host_filter import starsolo_barcode_args  # noqa: E402

#: Anellovirus accessions are the panel rows whose family is Anelloviridae.
FAMILY = "Anelloviridae"


def load_truth(path: Path) -> dict[str, dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return {r["read_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def split_planted(planted_dir: Path, truth: dict, out_dir: Path) -> dict[tuple, int]:
    """Write one FASTQ pair per (genome, window); return the read count of each.

    The planted reads sit at the end of the arm's FASTQs, but this filters by
    read ID rather than by position so a reordered file cannot silently shift
    which reads are treated as planted.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    handles: dict[tuple, tuple] = {}
    counts: dict[tuple, int] = defaultdict(int)
    with gzip.open(planted_dir / "R1.fastq.gz", "rt") as f1, gzip.open(
        planted_dir / "R2.fastq.gz", "rt"
    ) as f2:
        while True:
            rec1 = [f1.readline() for _ in range(4)]
            rec2 = [f2.readline() for _ in range(4)]
            if not rec1[0]:
                break
            read_id = rec1[0][1:].split()[0]
            row = truth.get(read_id)
            if row is None:
                continue
            key = (row["genome"], row["window"])
            if key not in handles:
                stem = out_dir / f"{row['genome']}_{row['window']}"
                handles[key] = (
                    gzip.open(f"{stem}_R1.fastq.gz", "wt"),
                    gzip.open(f"{stem}_R2.fastq.gz", "wt"),
                )
            h1, h2 = handles[key]
            h1.write("".join(rec1))
            h2.write("".join(rec2))
            counts[key] += 1
    for h1, h2 in handles.values():
        h1.close()
        h2.close()
    return dict(counts)


def anellovirus_genes(t2g: Path, identity: Path) -> set[str]:
    """Gene IDs of the panel's Anelloviridae records."""
    accessions = set()
    with open(identity, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("family") == FAMILY and row.get("genome_accession"):
                accessions.add(row["genome_accession"])
    genes = set()
    with open(t2g, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[4] in accessions:
                genes.add(f[1])
    return genes


def kallisto_molecules(fastq_stem: Path, index: Path, t2g: Path, genes: set[str],
                       whitelist: Path, out_dir: Path, threads: int) -> int:
    """Anellovirus molecules `kb count` reports for one planted set."""
    out_dir.mkdir(parents=True, exist_ok=True)
    # kb exits 0 with no counts when nothing pseudoaligns, and for a genome with
    # no near reference that is the expected result, not an error. Verify by
    # artifact: no counts directory means zero molecules.
    subprocess.run(
        ["kb", "count", "-i", str(index), "-g", str(t2g), "-x", "10xv3",
         "-w", str(whitelist), "-o", str(out_dir), "-t", str(threads), "--overwrite",
         f"{fastq_stem}_R1.fastq.gz", f"{fastq_stem}_R2.fastq.gz"],
        check=False, capture_output=True,
    )
    counts = out_dir / "counts_unfiltered"
    if not (counts / "cells_x_genes.mtx").is_file():
        return 0
    names = (counts / "cells_x_genes.genes.txt").read_text().split("\n")
    total = 0.0
    with open(counts / "cells_x_genes.mtx", encoding="utf-8") as fh:
        fh.readline()
        fh.readline()
        fh.readline()
        for line in fh:
            _, gene_idx, value = line.split()
            if names[int(gene_idx) - 1] in genes:
                total += float(value)
    return int(round(total))


def star_molecules(fastq_stem: Path, anello_index: Path, whitelist: Path,
                   out_dir: Path, threads: int) -> tuple[int, int, float]:
    """``(molecules, reads, median identity)`` the branch reports for one set."""
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        align_cmd("STAR", str(anello_index), f"{fastq_stem}_R2.fastq.gz",
                  f"{fastq_stem}_R1.fastq.gz",
                  starsolo_barcode_args("10xv3", str(whitelist)), f"{out_dir}/", threads),
        check=True, capture_output=True,
    )
    bam = out_dir / "Aligned.sortedByCoord.out.bam"
    if not bam.is_file():
        return 0, 0, 0.0
    subprocess.run(["samtools", "index", str(bam)], check=True)

    def alignments():
        out = subprocess.run(["samtools", "view", str(bam)], capture_output=True,
                             text=True, check=True)
        for line in out.stdout.splitlines():
            a = parse_sam_line(line)
            if a is not None:
                yield a

    lengths = {n: len(s) for n, s in iter_fasta(anello_index / "anello.fa")}
    rows = accession_metrics(alignments(), lengths)
    # Every panel contig in this index is Anelloviridae, so the whole index is
    # one group for the purpose of "did the branch see anellovirus sequence".
    molecules = virus_molecules(alignments(), {n: FAMILY for n in lengths})
    reads = len({a.qname for a in alignments()})
    idents = [float(r["median_identity"]) for r in rows if r["median_identity"] != ""]
    median = round(sorted(idents)[len(idents) // 2], 4) if idents else 0.0
    return molecules.get(FAMILY, {}).get("molecules", 0), reads, median


def planted_reads_surviving_host_filter(run_dir: Path, truth: dict) -> dict:
    """How many planted reads the STAR host filter kept (criterion 3).

    A planted read removed here never reaches either branch, so a loss would
    cap both methods equally and has to be reported rather than absorbed.
    """
    kept = 0
    filtered = list(run_dir.glob("*/host_filtered/R1.fastq.gz"))
    if not filtered:
        return {"error": f"no host_filtered/R1.fastq.gz under {run_dir}"}
    with gzip.open(filtered[0], "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0 and line[1:].split()[0] in truth:
                kept += 1
    return {"planted_total": len(truth), "planted_kept": kept,
            "planted_lost": len(truth) - kept}


def arm_summary(run_dir: Path) -> dict:
    """Anellovirus totals from one whole-arm viralscan run."""
    summaries = list(run_dir.glob("*/results/viral_summary.tsv"))
    if not summaries:
        return {"error": f"no viral_summary.tsv under {run_dir}"}
    kallisto = alignment = 0.0
    rows_seen = []
    with open(summaries[0], encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row.get("alignment_status"):
                continue  # not an anellovirus row
            kallisto += float(row["viral_molecules_total_est"] or 0)
            alignment += float(row["alignment_molecules_unique"] or 0)
            rows_seen.append(
                {
                    "virus_name": row["virus_name"],
                    "detection_source": row["detection_source"],
                    "kallisto_molecules": row["viral_molecules_total_est"],
                    "alignment_molecules": row["alignment_molecules_unique"],
                    "alignment_homopolymer_fraction": row["alignment_homopolymer_fraction"],
                }
            )
    return {
        "path": str(summaries[0]),
        "kallisto_anellovirus_molecules": int(round(kallisto)),
        "alignment_anellovirus_molecules": int(round(alignment)),
        "rows": rows_seen,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", type=Path, required=True, help="runs_anello_plant/")
    p.add_argument("--index", type=Path, required=True)
    p.add_argument("--t2g", type=Path, required=True)
    p.add_argument("--anello-index", type=Path, required=True)
    p.add_argument("--whitelist", type=Path, required=True)
    p.add_argument("--identity", type=Path, required=True,
                   help="virus_identity.tsv from any run on this panel")
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args()

    work = args.work
    truth = load_truth(work / "fastq" / "truth.tsv")
    selection = {r["genome"]: r for r in csv.DictReader(
        open(work / "fastq" / "selection.tsv", encoding="utf-8"), delimiter="\t")}
    genes = anellovirus_genes(args.t2g, args.identity)
    print(f"{len(genes)} anellovirus gene IDs in the panel", flush=True)

    print("Splitting planted reads per (genome, window) ...", flush=True)
    sets = split_planted(work / "fastq" / "planted", truth, work / "per_genome" / "fastq")

    rows = []
    for (genome, window), n_reads in sorted(sets.items()):
        stem = work / "per_genome" / "fastq" / f"{genome}_{window}"
        k = kallisto_molecules(stem, args.index, args.t2g, genes, args.whitelist,
                               work / "per_genome" / "kb" / f"{genome}_{window}", args.threads)
        s, s_reads, ident = star_molecules(stem, args.anello_index, args.whitelist,
                                           work / "per_genome" / "star" / f"{genome}_{window}",
                                           args.threads)
        meta = selection[genome]
        rows.append({
            "genome": genome,
            "genus": meta["genus"],
            "identity_band": meta["identity_band"],
            "nearest_identity": meta["nearest_identity"],
            "window": window,
            "planted_reads": n_reads,
            "kallisto_molecules": k,
            "alignment_molecules": s,
            "alignment_reads": s_reads,
            "alignment_median_identity": ident,
            "kallisto_recovery": round(k / n_reads, 4),
            "alignment_recovery": round(s / n_reads, 4),
        })
        print(f"  {genome:24s} {window:7s} planted={n_reads:5d} "
              f"kallisto={k:5d} alignment={s:5d}", flush=True)

    out_tsv = work / "acceptance_per_genome.tsv"
    with open(out_tsv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    arms = {a: arm_summary(work / "runs" / a) for a in ("planted", "unplanted", "negative")}

    # Pass criteria, evaluated per genome over the windows pooled.
    by_genome: dict[str, dict] = defaultdict(lambda: {"kallisto": 0, "alignment": 0})
    for r in rows:
        by_genome[r["genome"]]["kallisto"] += r["kallisto_molecules"]
        by_genome[r["genome"]]["alignment"] += r["alignment_molecules"]
    criterion_1 = []
    for genome, got in by_genome.items():
        band = selection[genome]["identity_band"]
        divergent = band != ">=95"
        ok = got["alignment"] > got["kallisto"] if divergent else got["alignment"] >= got["kallisto"]
        criterion_1.append({"genome": genome, "identity_band": band, **got, "pass": ok})
    host = planted_reads_surviving_host_filter(work / "runs" / "planted", truth)
    negative_alignment = arms["negative"].get("alignment_anellovirus_molecules")
    verdict = {
        "criterion_1_alignment_beats_kallisto": all(c["pass"] for c in criterion_1),
        "criterion_1_detail": criterion_1,
        "criterion_2_negative_is_zero": negative_alignment == 0,
        "criterion_2_negative_alignment_molecules": negative_alignment,
        "criterion_3_no_planted_read_lost": host.get("planted_lost") == 0,
        "criterion_3_detail": host,
        "arms": arms,
    }
    (work / "acceptance_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(json.dumps({k: v for k, v in verdict.items() if k != "arms"}, indent=2))
    print(f"\nWrote {out_tsv} and {work / 'acceptance_verdict.json'}")


if __name__ == "__main__":
    main()
