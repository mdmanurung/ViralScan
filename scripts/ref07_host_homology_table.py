#!/usr/bin/env python3
"""REF-07: per-locus viral/host homology table, partition-split (PLAN REF-07, VAL-01 D7).

Why this is not `build_reference.measure_host_homology`. That function keeps the
single longest hit per viral record, which answers "is this genome
host-homologous at all". `VAL-03` needs the opposite: every host locus that a
viral segment resembles, because each one is a *source* of adversarial challenge
reads. Multi-copy loci are exactly the dangerous class and are the ones a
max-per-query reduction throws away.

What the levels are built from. Properties of the sequence pair -- identity,
aligned length, longest exact run, shared 31-mers, strand, host expression
context -- never the behaviour of one index. STARsolo and Viral-Track carry no
kallisto D-list, so an index-tuned band would not transfer to them.

Partition discipline (protocol `leakage_prohibitions`, VAL-01 D4). Host loci are
challenge templates, so they must not be shared across partitions. Hits are
clustered by their *viral* interval first and whole clusters are assigned, so
near-identical multi-copy loci cannot land on both sides. Cut-points for the
frozen levels must be read off the training clusters alone.

Low complexity is held apart (advisor, 2026-10-05): `low_complexity` is its own
VAL-01 factor, and simple repeats would otherwise dominate the hit counts and
confound the two. Flagged hits stay in the table and are excluded when levels
are chosen.

Usage:

    python scripts/ref07_host_homology_table.py \\
        --viral-fasta <build>/viral.fa \\
        --host-genome <ref>/fasta/genome.fa \\
        --host-gtf <ref>/genes/genes.gtf \\
        --out <dir>/host_homology_loci.tsv
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

#: Read length the identity window is evaluated over. A challenge read is a
#: read, so identity across a whole multi-kb alignment is not the quantity that
#: decides whether one read is confusable; the worst 90-nt window is.
WINDOW = 90

#: kallisto needs one exact 31-mer. An exact run shorter than this contributes
#: no shared k-mer at all, which is the regime boundary the levels must bracket.
KMER = 31

#: 3' chemistry: only the last part of a transcript is sampled, so a host locus
#: further upstream cannot produce a realistic challenge read in 10x data.
THREE_PRIME_WINDOW = 300

#: Holdout share, matching protocol partitions.holdout_fraction.
HOLDOUT_FRACTION = 0.3

#: Minimum alignment block worth recording. Below ~50 nt a hit cannot carry a
#: 90-nt read's worth of confusable sequence.
MIN_BLOCK = 50

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_ATTR = re.compile(r'(\S+) "([^"]*)"')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def dust_intervals(viral_fasta: Path) -> dict[str, list[tuple[int, int]]]:
    """Low-complexity intervals per viral record, from ``dustmasker``."""
    dustmasker = shutil.which("dustmasker")
    if dustmasker is None:
        raise RuntimeError("dustmasker not found; it ships in the full ViralScan env")
    proc = subprocess.run(
        [dustmasker, "-in", str(viral_fasta), "-outfmt", "acclist"],
        check=True,
        capture_output=True,
        text=True,
    )
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for line in proc.stdout.splitlines():
        if not line.startswith(">"):
            continue
        name, start, end = line[1:].split("\t")
        out[name.split()[0]].append((int(start), int(end) + 1))
    return out


def overlap(intervals: list[tuple[int, int]], start: int, end: int) -> int:
    """Bases of ``[start, end)`` covered by sorted, merged ``intervals``."""
    if not intervals:
        return 0
    starts = [s for s, _ in intervals]
    i = max(0, bisect.bisect_right(starts, start) - 1)
    total = 0
    for s, e in intervals[i:]:
        if s >= end:
            break
        total += max(0, min(e, end) - max(s, start))
    return total


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def host_annotation(gtf: Path) -> tuple[dict, dict]:
    """``(exons, three_prime)`` merged interval lists per contig.

    ``three_prime`` is the terminal ``THREE_PRIME_WINDOW`` of each transcript's
    genomic span on its own strand. Taking the genomic span rather than walking
    the exon chain overstates the window on spliced 3' UTRs; it is an
    approximation, and it is recorded as such in the table header.
    """
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    spans: dict[str, tuple[str, int, int, str]] = {}
    with gtf.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            contig, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            exons[contig].append((start, end))
            attrs = dict(_ATTR.findall(f[8]))
            tid = attrs.get("transcript_id")
            if tid is None:
                continue
            if tid in spans:
                _, s0, e0, _ = spans[tid]
                spans[tid] = (contig, min(s0, start), max(e0, end), strand)
            else:
                spans[tid] = (contig, start, end, strand)
    three: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for contig, start, end, strand in spans.values():
        if strand == "-":
            three[contig].append((start, min(end, start + THREE_PRIME_WINDOW)))
        else:
            three[contig].append((max(start, end - THREE_PRIME_WINDOW), end))
    return (
        {c: merge(v) for c, v in exons.items()},
        {c: merge(v) for c, v in three.items()},
    )


def exact_runs(cigar: str) -> list[int]:
    """Lengths of exact-match runs, from an ``=``/``X`` CIGAR.

    minimap2 ``--eqx`` emits ``=`` for matches and ``X`` for mismatches, so a run
    of ``=`` is an exact stretch. Without ``--eqx`` every aligned base is ``M``
    and the longest run would be the whole block, which is wrong.
    """
    return [int(n) for n, op in _CIGAR.findall(cigar) if op == "="]


def window_identity(cigar: str, window: int = WINDOW) -> float:
    """Best identity over any ``window`` of query+reference aligned positions.

    A single read is the unit of confusion, so this reports the most confusable
    window rather than the whole-alignment average.
    """
    ops: list[tuple[int, bool]] = []
    for n, op in _CIGAR.findall(cigar):
        if op in "=XMID":
            ops.append((int(n), op in "=M"))
    column = []
    for n, is_match in ops:
        column.extend([is_match] * min(n, window * 4))
    if len(column) <= window:
        return sum(column) / len(column) if column else 0.0
    best = run = sum(column[:window])
    for i in range(window, len(column)):
        run += column[i] - column[i - window]
        best = max(best, run)
    return best / window


def run_minimap2(viral: Path, genome: Path, threads: int, extra: list[str]) -> str:
    """Align with the **viral panel as the reference** and the genome as query.

    The obvious direction (genome as reference) was tried first and returned
    zero alignments for all 2,343 genomes: an ``asm`` preset wants long colinear
    high-identity blocks, while viral/host homology is short diverged patches,
    and indexing 3.1 Gbp costs 78 s and 13 GB before a single hit. Indexing the
    11 MB panel instead costs 0.5 GB, runs ~70x faster, and finds the real
    signal (HHV-6 telomeric repeats, HSV-2 patches). Query and target are
    therefore swapped relative to the usual convention: PAF query = host contig,
    PAF target = viral record.
    """
    minimap2 = shutil.which("minimap2")
    if minimap2 is None:
        raise RuntimeError("minimap2 not found; it ships in the full ViralScan env")
    cmd = [
        minimap2, "-c", "--eqx", "-k", "15", "-w", "10",
        "-N", "200", "-p", "0.01", "-s", "40", "-m", "20",
        "-t", str(threads), *extra, str(viral), str(genome),
    ]
    print("  " + " ".join(cmd), file=sys.stderr)
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return proc.stdout


def partition_of(cluster_key: str, seed: int) -> str:
    """Whole-cluster assignment, deterministic from the key and the split seed.

    Hashing rather than apportioning: the cluster count is not known until the
    alignment runs, and the protocol's largest-remainder algorithm applies to
    biological samples, not to challenge templates.
    """
    digest = hashlib.sha256(f"{seed}:{cluster_key}".encode()).hexdigest()
    return "holdout" if int(digest[:8], 16) / 0xFFFFFFFF < HOLDOUT_FRACTION else "training"


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--viral-fasta", type=Path, required=True)
    p.add_argument("--host-genome", type=Path, required=True)
    p.add_argument("--host-gtf", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--paf", type=Path, help="reuse an existing PAF instead of aligning")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--split-seed", type=int, default=20260727001,
                   help="protocol seeds.split; hashed with the cluster key")
    p.add_argument("--minimap2-arg", action="append", default=[],
                   help="extra minimap2 argument, repeatable")
    args = p.parse_args()

    print("Low-complexity mask ...", file=sys.stderr)
    dust = dust_intervals(args.viral_fasta)
    dust = {k: merge(v) for k, v in dust.items()}

    print("Host annotation ...", file=sys.stderr)
    exons, three_prime = host_annotation(args.host_gtf)

    if args.paf and args.paf.exists():
        print(f"Reusing {args.paf}", file=sys.stderr)
        paf_text = args.paf.read_text()
    else:
        print("minimap2 ...", file=sys.stderr)
        paf_text = run_minimap2(
            args.viral_fasta, args.host_genome, args.threads, args.minimap2_arg
        )
        if args.paf:
            args.paf.write_text(paf_text)

    rows = []
    for line in paf_text.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        block = int(f[10])
        if block < MIN_BLOCK:
            continue
        tags = dict(t.split(":", 2)[::2] for t in f[12:] if t.count(":") >= 2)
        cigar = tags.get("cg", "")
        runs = exact_runs(cigar)
        # Swapped orientation (see run_minimap2): query = host, target = viral.
        h_contig, h_start, h_end = f[0], int(f[2]), int(f[3])
        v_record, v_start, v_end = f[5], int(f[7]), int(f[8])
        host_len = h_end - h_start
        rows.append({
            "viral_record": v_record,
            "viral_start": v_start,
            "viral_end": v_end,
            "viral_length": int(f[6]),
            "host_contig": h_contig,
            "host_start": h_start,
            "host_end": h_end,
            "strand": f[4],
            "aligned_bases": block,
            "identity_block": round(int(f[9]) / block, 4),
            "identity_best_90nt_window": round(window_identity(cigar), 4),
            "longest_exact_run": max(runs) if runs else 0,
            "shared_31mers": sum(max(0, r - KMER + 1) for r in runs),
            "low_complexity_fraction": round(
                overlap(dust.get(v_record, []), v_start, v_end) / max(1, v_end - v_start), 4
            ),
            "host_exon_fraction": round(
                overlap(exons.get(h_contig, []), h_start, h_end) / max(1, host_len), 4
            ),
            "host_3p_fraction": round(
                overlap(three_prime.get(h_contig, []), h_start, h_end) / max(1, host_len), 4
            ),
            "mapq": int(f[11]),
            "alignment_type": tags.get("tp", ""),
        })

    # Cluster so that a template family cannot straddle the split. Two rows join
    # when they share a viral segment OR a host locus. Viral-side merging alone
    # is not enough: MPXV and Cowpox both hit chr9:109765432-109765517 from
    # different records, and HHV-6A/6B/7 all hit the same human telomeres, so a
    # per-record clustering would hash one host locus into both partitions and
    # breach D4.
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    def group_keys(prefix: str, key_of, start_of, end_of) -> list[str]:
        """One interval-group key per row, merging overlaps within each key."""
        buckets: dict[str, list[int]] = defaultdict(list)
        for i, r in enumerate(rows):
            buckets[key_of(r)].append(i)
        out = [""] * len(rows)
        for key, idx in buckets.items():
            merged = merge([(start_of(rows[i]), end_of(rows[i])) for i in idx])
            for i in idx:
                for n, (s, e) in enumerate(merged):
                    if start_of(rows[i]) < e and s < end_of(rows[i]):
                        out[i] = f"{prefix}:{key}:{n}"
                        break
        return out

    viral_keys = group_keys(
        "v", lambda r: r["viral_record"], lambda r: r["viral_start"], lambda r: r["viral_end"]
    )
    host_keys = group_keys(
        "h", lambda r: r["host_contig"], lambda r: r["host_start"], lambda r: r["host_end"]
    )
    for vk, hk in zip(viral_keys, host_keys):
        union(vk, hk)
    for i, r in enumerate(rows):
        key = find(viral_keys[i])
        r["cluster_id"] = key
        r["partition"] = partition_of(key, args.split_seed)

    clusters = {r["cluster_id"] for r in rows}
    header = {
        "generated": "REF-07 host-homology locus table",
        "viral_fasta": str(args.viral_fasta),
        "viral_fasta_sha256": sha256(args.viral_fasta),
        "host_genome": str(args.host_genome),
        "host_genome_sha256": sha256(args.host_genome),
        "host_gtf": str(args.host_gtf),
        "host_gtf_sha256": sha256(args.host_gtf),
        "split_seed": args.split_seed,
        "holdout_fraction": HOLDOUT_FRACTION,
        "window_nt": WINDOW,
        "kmer": KMER,
        "min_block": MIN_BLOCK,
        "three_prime_window_nt": THREE_PRIME_WINDOW,
        "caveats": [
            "host_3p_fraction uses the transcript's genomic span, which overstates "
            "the window on spliced 3' UTRs",
            "levels must be chosen from partition == training rows only",
            "low_complexity_fraction > 0 rows belong to the low_complexity factor, "
            "not to host_virus_homology",
        ],
        "n_alignments": len(rows),
        "n_clusters": len(clusters),
    }
    if not rows:
        raise SystemExit("no alignment survived the filters; nothing to write")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        for line in json.dumps(header, indent=2).splitlines():
            fh.write(f"# {line}\n")
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(
            rows, key=lambda r: (r["viral_record"], r["viral_start"], r["host_contig"])
        ))
    print(json.dumps(header, indent=2))


if __name__ == "__main__":
    main()
