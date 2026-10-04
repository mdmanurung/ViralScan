#!/usr/bin/env python3
"""Body census over anellovirus-aligned reads (PLAN ANELLO-PRIOR.3).

A genuine 10x R2 read of a TTV mRNA 3' end is ``[complex viral body]
[untemplated poly-A]``, so the one label that does not lean on the branch's own
measures is whether the *body* — the sequence 5' of the first long homopolymer
run (``read_body``) — is itself anellovirus sequence. This script places every
body >= ``MIN_BODY_LEN`` on the anellovirus panel, then scores the read-side
measures (``is_complex_body``, ``has_tso``) against that label:

* sensitivity — does ``complex_body and not tso`` keep every body-mapping read?
* leak — how many reads it keeps have a body that maps nowhere?

The body-mapping reads, against the reads and distinct CB+UMI examined, bound
the genuine component of the call. Nothing here filters anything.

Placement is ungapped, on either strand. A body "maps" when its longest stretch
with <= 2 mismatches is >= min(25, body length) — the review's mlen >= 25,
NM <= 2, relaxed for the 20-24 nt bodies ``MIN_BODY_LEN`` admits. Seeding is
exhaustive for that criterion: a window of ``need`` bases with <= 2 mismatches
holds an exact run of ceil((need - 2) / 3), and every k-mer of that length in
the body is seeded.

Usage::

    samtools view -F 0x904 viral_reads.bam > primary.sam
    python scripts/anello_body_census.py primary.sam anello_refs.fa \\
        src/viralscan/data/anellovirus_accessions.tsv --out census.tsv
"""

import argparse
import bisect
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from viralscan.anello_align import (  # noqa: E402
    MIN_BODY_LEN,
    TSO,
    _contains,
    _revcomp,
    dinucleotide_entropy,
    has_tso,
    is_complex_body,
    parse_sam_line,
    read_body,
)

MAX_MISMATCH = 2
MIN_MATCH = 25
#: 3' 15 nt of the TSO: what survives when the read starts inside the oligo.
TSO_CORE = TSO[-15:]
#: Illumina TruSeq Read 1 primer tail (10x R1 side) — reagent, not template.
TRUSEQ_R1 = "CTACACGACGCTCTTCCGATCT"


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def longest_stretch(a, b, max_mm):
    """Longest window of equal-length ``a``/``b`` holding <= max_mm mismatches."""
    best = left = mm = 0
    for right in range(len(a)):
        mm += a[right] != b[right]
        while mm > max_mm:
            mm -= a[left] != b[left]
            left += 1
        best = max(best, right - left + 1)
    return best


class Panel:
    """Every genome, both strands, in one N-separated string."""

    def __init__(self, genomes):
        parts, self.starts, self.names, pos = [], [], [], 0
        for name, g in genomes.items():
            for strand in (g, _revcomp(g)):
                self.starts.append(pos)
                self.names.append(name)
                parts.append(strand)
                pos += len(strand) + 1
        self.seq = "N".join(parts)

    def best(self, body):
        """(longest <=2-mismatch stretch, reference) over all seeded placements."""
        need = min(MIN_MATCH, len(body))
        k = math.ceil((need - MAX_MISMATCH) / (MAX_MISMATCH + 1))
        offsets = set()
        for i in range(len(body) - k + 1):
            start = self.seq.find(body[i : i + k])
            while start != -1:
                offsets.add(start - i)
                start = self.seq.find(body[i : i + k], start + 1)
        best = (0, "")
        for off in offsets:
            lo, hi = max(0, -off), min(len(body), len(self.seq) - off)
            hit = longest_stretch(body[lo:hi], self.seq[off + lo : off + hi], MAX_MISMATCH)
            if hit > best[0]:
                best = (hit, self.names[bisect.bisect_right(self.starts, off + lo) - 1])
        return best


def reagent(seq):
    """First reagent class found in the read, else ''."""
    up = seq.upper()
    if has_tso(up):
        return "tso"
    if TSO_CORE in up or _revcomp(TSO_CORE) in up:
        return "tso_partial"
    if _contains(up, TRUSEQ_R1, 2) or _contains(up, _revcomp(TRUSEQ_R1), 2):
        return "truseq_r1"
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("sam", help="primary records (samtools view -F 0x904)")
    ap.add_argument("refs", help="anellovirus genome FASTA")
    ap.add_argument("accessions", help="anellovirus_accessions.tsv")
    ap.add_argument("--out", required=True, help="per-read TSV")
    args = ap.parse_args()

    anello = {r["accession"] for r in csv.DictReader(open(args.accessions), delimiter="\t")}

    rows = []
    for line in open(args.sam):
        aln = parse_sam_line(line)
        if aln is None or aln.secondary or aln.rname not in anello:
            continue
        seq = aln.read_seq()
        body = read_body(seq)
        cb, umi = aln.qname.split("_")[:2]
        rows.append({
            "qname": aln.qname, "cb_umi": f"{cb}_{umi}", "rname": aln.rname,
            "reverse": int(aln.reverse), "body_len": len(body),
            "body_entropy": round(dinucleotide_entropy(body), 3),
            "complex_body": int(is_complex_body(seq)), "tso": int(has_tso(seq)),
            "reagent": reagent(seq), "query_coverage": round(aln.query_coverage(), 3),
            "body_match": "", "body_match_ref": "", "body_maps": 0, "body": body,
        })

    panel = Panel(read_fasta(args.refs))
    placed = {}
    for r in rows:
        if r["body_len"] < MIN_BODY_LEN:
            continue
        body = r["body"].upper()
        if body not in placed:
            placed[body] = panel.best(body)
        r["body_match"], r["body_match_ref"] = placed[body]
        r["body_maps"] = int(r["body_match"] >= min(MIN_MATCH, r["body_len"]))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    def report(label, sel):
        print(f"{label}\t{len(sel)}\t{len({r['cb_umi'] for r in sel})}")

    kept = [r for r in rows if r["complex_body"] and not r["tso"]]
    maps = [r for r in rows if r["body_maps"]]
    print(f"unique_bodies_placed\t{len(placed)}")
    print("class\treads\tcb_umi")
    report("aligned_primary", rows)
    report(f"body_len>={MIN_BODY_LEN}", [r for r in rows if r["body_len"] >= MIN_BODY_LEN])
    report("complex_body", [r for r in rows if r["complex_body"]])
    report("tso", [r for r in rows if r["tso"]])
    report("kept(complex_and_no_tso)", kept)
    for cls in ("tso_partial", "truseq_r1", ""):
        report(f"kept|reagent={cls or 'none'}", [r for r in kept if r["reagent"] == cls])
    report("body_maps", maps)
    report("body_maps_and_kept", [r for r in maps if r["complex_body"] and not r["tso"]])
    report("body_maps_no_reagent", [r for r in maps if not r["reagent"]])
    report("kept_body_does_not_map", [r for r in kept if not r["body_maps"]])

    # Post-hoc refinement, chosen after inspecting the body_maps reads: the
    # body places over >= 90 % of its length, carries no reagent (including a
    # TSO with mismatches in its core) and is complex. It leans on
    # is_complex_body, the measure under test, so it is not the bound.
    full = [
        r for r in rows
        if r["body_match"] != "" and r["body_match"] >= 0.9 * r["body_len"]
        and not r["reagent"] and not _contains(TSO, r["body"][-15:].upper(), 3)
    ]
    report("posthoc_full_length_reagent_free", full)
    report("posthoc_full_length_reagent_free_complex", [r for r in full if r["complex_body"]])
    print(f"posthoc_max_body_len\t{max((r['body_len'] for r in full), default=0)}")


if __name__ == "__main__":
    main()
