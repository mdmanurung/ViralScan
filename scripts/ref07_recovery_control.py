#!/usr/bin/env python3
"""REF-07 sensitivity control: can the aligner see diverged non-repeat homology?

F-023 claims that viral/host homology in the panel is almost entirely
low-complexity. That claim is only as strong as the aligner's ability to find
the *other* kind -- a diverged, non-repeat, exonic match -- and the 8/8
genome-level control cannot show that, for two reasons: it matches at genome
level (HHV-6 passes on 1,439 telomeric hits whatever happens at `p23`), and
every genome in it produced low-complexity hits in the first place.

So: take real non-dust human exons, mutate them to known identities, write them
as a fake viral panel, and run the identical minimap2 command. Recovery per
identity tier is the sensitivity curve. The 85-90 % band matters most: that is
where a comparator such as STARsolo would still align while kallisto would not,
and it is the band F-023 says is empty.

Usage:

    python scripts/ref07_recovery_control.py \\
        --host-genome <ref>/fasta/genome.fa --host-gtf <ref>/genes/genes.gtf \\
        --out-dir <dir>
"""

from __future__ import annotations

import argparse
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

#: Identity tiers to plant. 0.95 should be trivially recovered; 0.85 is the
#: hard edge; 0.80 is below where a 90-nt read carries an exact 31-mer at all.
TIERS = (0.95, 0.90, 0.85, 0.80)

#: Exons per tier. Enough to read a recovery rate, small enough to stay cheap.
PER_TIER = 20

#: Minimum exon length: a decoy must be able to carry a 90-nt read plus margin.
MIN_EXON = 300

#: Chromosomes the decoys are drawn from, so the query side stays small.
CHROMS = ("chr20", "chr21", "chr22")

BASES = "ACGT"
_ATTR = re.compile(r'gene_id "([^"]*)"')


def read_chroms(genome: Path, wanted: set[str]) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with genome.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in wanted:
                    seqs[name] = []
            elif name in seqs:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def exons(gtf: Path, wanted: set[str]) -> list[tuple[str, int, int, str]]:
    out = []
    with gtf.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "exon" or f[0] not in wanted:
                continue
            start, end = int(f[3]) - 1, int(f[4])
            if end - start < MIN_EXON:
                continue
            gene = _ATTR.search(f[8])
            out.append((f[0], start, end, gene.group(1) if gene else "?"))
    return out


def dust(path: Path) -> set[str]:
    """Record ids carrying any dustmasker interval."""
    dustmasker = shutil.which("dustmasker")
    if dustmasker is None:
        raise RuntimeError("dustmasker not found; it ships in the full ViralScan env")
    proc = subprocess.run(
        [dustmasker, "-in", str(path), "-outfmt", "acclist"],
        check=True, capture_output=True, text=True,
    )
    return {ln[1:].split("\t")[0].split()[0] for ln in proc.stdout.splitlines()
            if ln.startswith(">")}


def mutate(seq: str, identity: float, rng: random.Random) -> str:
    out = list(seq)
    for i in rng.sample(range(len(out)), int(round(len(out) * (1 - identity)))):
        out[i] = rng.choice([b for b in BASES if b != out[i]])
    return "".join(out)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--host-genome", type=Path, required=True)
    p.add_argument("--host-gtf", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--seed", type=int, default=20261005)
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args()

    rng = random.Random(args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    wanted = set(CHROMS)

    print("Reading chromosomes and exons ...", file=sys.stderr)
    seqs = read_chroms(args.host_genome, wanted)
    pool = [e for e in exons(args.host_gtf, wanted) if e[0] in seqs]
    print(f"  {len(seqs)} chromosomes, {len(pool)} exons >= {MIN_EXON} nt", file=sys.stderr)

    # Draw candidates, then drop any the masker touches: a decoy must be
    # non-repeat, or the control would measure the same thing as the finding.
    rng.shuffle(pool)
    cand = args.out_dir / "candidates.fa"
    picked: list[tuple[str, str]] = []
    with cand.open("w", encoding="utf-8") as fh:
        for contig, start, end, gene in pool[: PER_TIER * len(TIERS) * 4]:
            s = seqs[contig][start:end]
            if "N" in s:
                continue
            name = f"{contig}_{start}_{end}_{gene}"
            picked.append((name, s))
            fh.write(f">{name}\n{s}\n")
    dusty = dust(cand)
    clean = [(n, s) for n, s in picked if n not in dusty]
    print(f"  {len(clean)} of {len(picked)} candidate exons carry no dust", file=sys.stderr)
    if len(clean) < PER_TIER * len(TIERS):
        raise SystemExit(f"only {len(clean)} clean exons; need {PER_TIER * len(TIERS)}")

    decoys = args.out_dir / "decoy_panel.fa"
    truth: list[tuple[str, float, str]] = []
    with decoys.open("w", encoding="utf-8") as fh:
        it = iter(clean)
        for tier in TIERS:
            for _ in range(PER_TIER):
                name, seq = next(it)
                rec = f"decoy_{tier:.2f}_{name}"
                fh.write(f">{rec}\n{mutate(seq, tier, rng)}\n")
                truth.append((rec, tier, name))

    query = args.out_dir / "query.fa"
    with query.open("w", encoding="utf-8") as fh:
        for contig in CHROMS:
            fh.write(f">{contig}\n{seqs[contig]}\n")

    minimap2 = shutil.which("minimap2")
    cmd = [minimap2, "-c", "--eqx", "-k", "15", "-w", "10",
           "-N", "200", "-p", "0.01", "-s", "40", "-m", "20",
           "-t", str(args.threads), str(decoys), str(query)]
    print("  " + " ".join(cmd), file=sys.stderr)
    paf = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    (args.out_dir / "recovery.paf").write_text(paf)

    hit: dict[str, int] = {}
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12 or int(f[10]) < 50:
            continue
        hit[f[5]] = max(hit.get(f[5], 0), int(f[10]))

    print(f"\n== recovery of non-repeat decoys ({PER_TIER} per tier, "
          f"blocks >= 50 nt, same minimap2 settings as the REF-07 table)")
    for tier in TIERS:
        recs = [r for r, t, _ in truth if t == tier]
        found = [r for r in recs if r in hit]
        blocks = sorted(hit[r] for r in found)
        med = blocks[len(blocks) // 2] if blocks else 0
        print(f"   identity {tier:.2f}: recovered {len(found):2d}/{len(recs)} "
              f"({100*len(found)/len(recs):5.1f}%)  median block {med} nt")
    print("\nRead this as the sensitivity floor of the REF-07 measurement: a tier "
          "recovered poorly here cannot be claimed absent in the panel.")


if __name__ == "__main__":
    main()
