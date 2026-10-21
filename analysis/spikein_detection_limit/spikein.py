"""Spike real viral reads into a null 10x library to measure the detection limit.

EXPLORATORY (user override of discovery Q21 / LOCK-01 A2, 2026-10-21). Standard library only; calls ``samtools``.

``bg``     keep the first N read pairs of a background library (once; the spikes are layered on it).
``spike``  draw viral R2 reads from a STARsolo BAM, stratified by flank-audit class at the source's real
           proportions, give each its own random UMI in one of K background barcodes, and write
           ``spike_L<level>_s<seed>_R{1,2}.fastq.gz`` plus a truth table. Spikes are separate files:
           ``cat`` them after the background (concatenated gzip is valid gzip) or pass them as a
           comma-separated list to STAR / ``kb count``.

Source R2 is restored to its original orientation (BAM SEQ is reverse-complemented for flag 16), trimmed to
the background R2 length keeping the 5' end, and reverse-complemented only when the source and background
cDNA read orientations differ. One spiked read = one spiked UMI (conservative: real data carry PCR duplicates).
"""

import argparse
import collections
import gzip
import itertools
import json
import random
import re
import subprocess
import sys
from pathlib import Path

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
FLANK = 1000  # same window as scripts/flank_audit.py
CLASSES = ("sense_in_gene", "flank3_unannotated", "antisense_in_gene", "intergenic")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def load_genes(gtf: str, contig: str) -> list[tuple[int, int, str]]:
    genes = []
    with open(gtf) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene" and f[0] == contig:
                genes.append((int(f[3]), int(f[4]), f[6]))
    return genes


def classify(flag: int, pos: int, cigar: str, genes, orient: str) -> str:
    """Read class; the same rule as ``flank_audit.py`` (orient: RNA strand is opposite the read when 'reverse')."""
    span = sum(int(a) for a, b in CIG.findall(cigar) if b in "MDN=X")
    s, e = pos, pos + span - 1
    rna = ("+" if flag & 16 else "-") if orient == "reverse" else ("-" if flag & 16 else "+")
    sense = anti = flank3 = False
    for gs, ge, st in genes:
        if gs <= e and s <= ge:
            if st == rna:
                sense = True
            else:
                anti = True
        elif st == rna and (
            (st == "+" and ge < s <= ge + FLANK) or (st == "-" and gs - FLANK <= e < gs)
        ):
            flank3 = True
    return (
        "sense_in_gene"
        if sense
        else "flank3_unannotated"
        if flank3
        else "antisense_in_gene"
        if anti
        else "intergenic"
    )


def allocate(counts: dict[str, int], n: int) -> dict[str, int]:
    """Largest-remainder split of n over the classes in proportion to counts."""
    tot = sum(counts.values())
    raw = {k: n * v / tot for k, v in counts.items()}
    out = {k: int(x) for k, x in raw.items()}
    for k in sorted(raw, key=lambda k: -(raw[k] - out[k]))[: n - sum(out.values())]:
        out[k] += 1
    return out


def source_pool(
    bam: str, contig: str, genes, orient: str, r2_len: int, seed: int, target: int = 20000
):
    """~target primary unique reads from *contig* (``samtools view -s``), original R2 orientation restored."""
    total = int(
        next(
            ln.split("\t")[2]
            for ln in subprocess.run(
                ["samtools", "idxstats", bam], capture_output=True, text=True, check=True
            ).stdout.splitlines()
            if ln.split("\t")[0] == contig
        )
    )
    frac = min(1.0, target / max(total, 1))
    p = subprocess.Popen(
        ["samtools", "view", "-F", "0x904", "-s", f"{seed + frac:.6f}", bam, contig],
        stdout=subprocess.PIPE,
        text=True,
    )
    pool = []
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        flag, seq, qual = int(f[1]), f[9], f[10]
        if "NH:i:1" not in f[11:] or len(seq) < r2_len:
            continue
        if flag & 16:
            seq, qual = revcomp(seq), qual[::-1]
        pool.append(
            {
                "seq": seq,
                "qual": qual,
                "cls": classify(flag, int(f[3]), f[5], genes, orient),
                "pos": int(f[3]),
                "flag": flag,
            }
        )
    p.wait()
    return pool


def top_barcodes(r1: str, cb_len: int, k: int, n: int = 1_000_000) -> list[str]:
    c = collections.Counter()
    with gzip.open(r1, "rt") as fh:
        for seq in itertools.islice(itertools.islice(fh, 1, None, 4), n):
            c[seq[:cb_len]] += 1
    return [b for b, _ in c.most_common(k)]


def make_spike(read: dict, name: str, cb: str, umi: str, r1_len: int, r2_len: int, rc: bool):
    """(R1 record, R2 record) for one spiked read; R1 = barcode + UMI, padded with A to the background R1 length."""
    seq, qual = read["seq"][:r2_len], read["qual"][:r2_len]
    if rc:
        seq, qual = revcomp(seq), qual[::-1]
    r1 = (cb + umi).ljust(r1_len, "A")
    return f"@{name}\n{r1}\n+\n{'F' * len(r1)}\n", f"@{name}\n{seq}\n+\n{qual}\n"


def cmd_bg(a):
    for src, dst in ((a.r1, a.out_prefix + "_R1.fastq.gz"), (a.r2, a.out_prefix + "_R2.fastq.gz")):
        with gzip.open(src, "rt") as fi, gzip.open(dst, "wt", compresslevel=3) as fo:
            fo.writelines(itertools.islice(fi, 4 * a.n_pairs))
    print(
        json.dumps(
            {
                "n_pairs": a.n_pairs,
                "out_prefix": a.out_prefix,
                "note": "first N reads, not a random subsample",
            }
        )
    )


def cmd_spike(a):
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    genes = load_genes(a.gtf, a.contig)
    cbs = top_barcodes(a.bg_r1, a.cb_len, a.cells)
    rc = a.source_orient != a.bg_orient
    for seed in a.seeds:
        pool = source_pool(a.bam, a.contig, genes, a.source_orient, a.r2_len, seed)
        by_cls = {c: [r for r in pool if r["cls"] == c] for c in CLASSES}
        mix = {c: len(v) for c, v in by_cls.items()}
        print(f"seed {seed}: pool {len(pool)} class mix {mix}", file=sys.stderr)
        for level in a.levels:
            rng = random.Random(f"{seed}-{level}")
            used, rows = set(), []
            i = 0
            r1_path = out / f"spike_L{level}_s{seed}_R1.fastq.gz"
            r2_path = out / f"spike_L{level}_s{seed}_R2.fastq.gz"
            with gzip.open(r1_path, "wt") as r1o, gzip.open(r2_path, "wt") as r2o:
                for cls, k in allocate(mix, level).items():
                    for read in rng.sample(by_cls[cls], k):
                        while True:
                            umi = "".join(rng.choices("ACGT", k=a.umi_len))
                            if (cbs[i % len(cbs)], umi) not in used:
                                break
                        cb = cbs[i % len(cbs)]
                        used.add((cb, umi))
                        name = f"SPIKE:{cls}:{i}"
                        rec1, rec2 = make_spike(read, name, cb, umi, a.r1_len, a.r2_len, rc)
                        r1o.write(rec1)
                        r2o.write(rec2)
                        rows.append((name, cb, umi, cls, a.contig, read["pos"], read["flag"]))
                        i += 1
            with open(out / f"spike_L{level}_s{seed}.truth.tsv", "w") as t:
                t.write("name\tcb\tumi\tclass\tcontig\tsrc_pos\tsrc_flag\n")
                t.writelines("\t".join(map(str, r)) + "\n" for r in rows)
            (out / f"spike_L{level}_s{seed}.json").write_text(
                json.dumps(
                    {**vars(a), "seed": seed, "level": level, "rc": rc, "pool_class_mix": mix},
                    indent=1,
                )
            )


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("bg")
    b.add_argument("--r1", required=True)
    b.add_argument("--r2", required=True)
    b.add_argument("--n-pairs", type=int, required=True)
    b.add_argument("--out-prefix", required=True)
    s = sub.add_parser("spike")
    s.add_argument("--bam", required=True)
    s.add_argument("--gtf", required=True)
    s.add_argument("--bg-r1", required=True, help="background R1 (for barcode choice)")
    s.add_argument("--out-dir", required=True)
    s.add_argument("--contig", default="VIRUS|AF157706.1")
    s.add_argument("--levels", type=lambda x: [int(v) for v in x.split(",")], required=True)
    s.add_argument("--seeds", type=lambda x: [int(v) for v in x.split(",")], default=[1])
    s.add_argument("--cells", type=int, default=5)
    s.add_argument("--cb-len", type=int, required=True)
    s.add_argument("--umi-len", type=int, required=True)
    s.add_argument("--r1-len", type=int, required=True)
    s.add_argument("--r2-len", type=int, required=True)
    s.add_argument("--source-orient", choices=("forward", "reverse"), default="reverse")
    s.add_argument("--bg-orient", choices=("forward", "reverse"), required=True)
    a = p.parse_args(argv)
    {"bg": cmd_bg, "spike": cmd_spike}[a.cmd](a)


if __name__ == "__main__":
    main()
