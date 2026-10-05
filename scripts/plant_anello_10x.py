#!/usr/bin/env python3
"""Plant held-out anellovirus reads into a real 10x library (PLAN ANDET-09e).

Builds the three FASTQ sets the acceptance run compares:

* ``planted/``   — background + reads simulated from held-out genomes;
* ``unplanted/`` — the identical background, nothing added;
* ``negative/``  — synthetic host + artefact reads, no viral sequence at all.

Why held-out genomes. The panel's anellovirus set was clustered at 95 % ANI
(F-013), so a divergent strain shares few exact 31-mers with its nearest
reference and kallisto cannot see it. None of the 829 Modha et al. 2025
genomes (TPA BK068993-BK069821) is in the panel — verified by accession here —
so they measure exactly the regime the alignment branch exists for.

Why planted minus unplanted. The covid background is real human blood and may
carry genuine anellovirus, so an absolute count in the planted arm is not
recovery. The same background in both arms makes the difference the planted
signal.

Why a synthetic negative. A nonzero count in a real library is not proof of an
artefact (it could be real commensal anellovirus). The negative contains no
viral sequence, so anything the branch reports there is a false positive. Its
artefact classes are the measured ones: poly-G no-signal reads (F-019), poly-A
tail capture (F-021), CAG repeats and TSO+poly-A joins (F-010). They are
planted far above their real rate, which makes it a stress test.

Sequences and CDS coordinates both come from ``Modha_genomes_annotated.gbk``.
The repository's ``Modha_contigs.fas`` is **ORF1 only** -- all 829 of its
records match the metadata's ``ORF1_len`` exactly -- so planting from it would
sample only the hypervariable ORF1 and leave the windows undefined.

Usage (see PLAN ANDET-09e for the pass criterion):

    python scripts/plant_anello_10x.py --out-dir <dir> \\
        --gbk <dir>/heldout/Modha_genomes_annotated.gbk \\
        --metadata <dir>/heldout/Modha_allsequence_data.csv \\
        --paf <dir>/heldout/modha_vs_panel.paf \\
        --background-r1 <run>/host_filtered/R1.fastq.gz \\
        --background-r2 <run>/host_filtered/R2.fastq.gz \\
        --whitelist <ref>/cellranger_whitelist.txt \\
        --host-cdna <build>/cdna.fa --host-genome <ref>/genome.fa
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import re
import sys
from pathlib import Path

# ── Simulation constants ─────────────────────────────────────────────────────
#: 10x v3: R1 is 16 nt barcode + 12 nt UMI; the x213 library's R2 is 90 nt.
CB_LEN, UMI_LEN, READ_LEN = 16, 12, 90

#: Per-base substitution rate for planted viral reads. F-010's held-out plant
#: used 1 %, and this plant is read against the same finding.
ERROR_RATE = 0.01

#: Host reads carry a lower rate: they stand in for ordinary sequencing error,
#: not for divergence.
HOST_ERROR_RATE = 0.001

#: Molecules planted per genome per window. One read per molecule, so a
#: molecule is one distinct (CB, UMI) and recovery is countable directly.
LEVELS = (10, 100, 1000)

#: 5' window: the 500 nt before ORF1, as in F-010. 3' window: the 300 nt after
#: the last CDS ends, where a 3' 10x library's reads sit (the anellovirus polyA
#: site is downstream of the ORFs). Genomes are circular, so both wrap.
WINDOW_5P_LEN, WINDOW_3P_LEN = 500, 300

#: 10x template-switch oligo; TSO+poly-A joins are a measured artefact (F-010).
TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"

#: Shortest held-out genome worth planting. These are assembled genomes and
#: 609 of 829 are partial, so a strict whole-genome floor would leave the
#: under-sampled genera with nothing to pick.
MIN_LENGTH = 1500

#: A homopolymer this long marks the artefact classes the negative plants.
ARTEFACT_RUN = 40

#: Barcodes the planted molecules are spread over, so a plant is not one cell.
N_CELLS = 2000

BASES = "ACGT"


# ── Inputs ───────────────────────────────────────────────────────────────────
#: Ensembl transcript IDs. The negative must contain no viral sequence, and a
#: pipeline's own `cdna.fa` is the COMBINED host+viral cDNA that kb ref built:
#: the cat42b one holds 471,944 records of which 6,175 are viral, anelloviruses
#: among them. Sampling "host" reads from it puts real virus in the negative
#: control, which is how the first ANDET-09e negative arm reported ~16,000
#: anellovirus molecules.
HOST_RECORD_PREFIX = "ENST"


def read_fasta(path: Path, prefix: str = "") -> dict[str, str]:
    """Records of *path*, keeping only those whose ID starts with *prefix*."""
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks)
    if prefix:
        kept = {k: v for k, v in seqs.items() if k.startswith(prefix)}
        if not kept:
            raise ValueError(
                f"{path} holds no record whose ID starts with {prefix!r}. The "
                "negative control must be built from host sequence only."
            )
        if len(kept) != len(seqs):
            print(f"  {path.name}: kept {len(kept)} {prefix}* of {len(seqs)} records "
                  f"({len(seqs) - len(kept)} non-host excluded)")
        return kept
    return seqs


_CDS_RE = re.compile(r"^     CDS             (\d+)\.\.(\d+)\s*$")
_SEQ_RE = re.compile(r"[^acgtnACGTN]")


def read_gbk(gbk: Path) -> dict[str, dict]:
    """``accession -> {seq, orf1_start, last_cds_end}`` from the annotated GenBank.

    ORF1 is taken as the longest CDS: every product in this file is labelled
    "hypothetical protein", but ORF1 is several times the length of the other
    ORFs in every anellovirus genome. Joined and complemented CDS lines are
    skipped, and a genome left with no simple CDS is dropped.
    """
    out: dict[str, dict] = {}
    name: str | None = None
    spans: list[tuple[int, int]] = []
    chunks: list[str] = []
    in_origin = False

    def flush() -> None:
        if name and spans and chunks:
            orf1 = max(spans, key=lambda s: s[1] - s[0])
            out[name] = {
                "seq": "".join(chunks).upper(),
                "orf1_start": orf1[0],
                "last_cds_end": max(e for _, e in spans),
            }

    with open(gbk, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                flush()
                name, spans, chunks, in_origin = line.split()[1], [], [], False
            elif line.startswith("ORIGIN"):
                in_origin = True
            elif line.startswith("//"):
                in_origin = False
            elif in_origin:
                chunks.append(_SEQ_RE.sub("", line))
            else:
                m = _CDS_RE.match(line.rstrip("\n"))
                if m:
                    spans.append((int(m.group(1)) - 1, int(m.group(2))))
    flush()
    return out


def nearest_identity(paf: Path) -> dict[str, tuple[float, str]]:
    """``contig -> (identity, panel genome)`` for the best panel alignment.

    Identity is minimap2's matching bases over alignment block length. The best
    hit maximises identity times query coverage, so a short high-identity patch
    does not outrank a genome-length alignment.
    """
    best: dict[str, tuple[float, str, float]] = {}
    with open(paf, encoding="utf-8") as fh:
        for line in fh:
            f = line.split("\t")
            query, qlen, target = f[0], int(f[1]), f[5]
            block = int(f[10])
            identity, coverage = int(f[9]) / block, block / qlen
            if query not in best or identity * coverage > best[query][2]:
                best[query] = (identity, target, identity * coverage)
    return {q: (v[0], v[1]) for q, v in best.items()}


def band(identity: float | None) -> str:
    """Identity band; kallisto needs an exact 31-mer, so 0.95**31 ~ 0.20 and
    0.85**31 ~ 0.007 bracket where pseudoalignment stops working."""
    if identity is None:
        return "no_alignment"
    for cut, label in ((0.95, ">=95"), (0.90, "90-95"), (0.85, "85-90")):
        if identity >= cut:
            return label
    return "<85"


def select_genomes(genomes, metadata, completeness, paf, rng) -> list[dict]:
    """One genome per (genus, identity band) cell, over the cells that exist.

    Spanning both axes is the point: the identity axis measures the divergence
    where the branch should beat kallisto, and the genus axis covers the four
    genera the panel barely represents (Gyro-, Mem-, Samek-, Hetorquevirus;
    PLAN ANELLO-PRIOR, CAT-20).
    """
    meta = {r["label"]: r for r in csv.DictReader(open(metadata, encoding="utf-8"))}
    complete = {
        r["ContigID"]: r["PredictionResult"]
        for r in csv.DictReader(open(completeness, encoding="utf-8"))
    }
    nearest = nearest_identity(paf)
    wanted = [
        ("alpha", ">=95"),
        ("alpha", "85-90"),
        ("beta", "90-95"),
        ("beta", "<85"),
        ("gamma", "no_alignment"),
        ("samek", "<85"),
        ("mem", "<85"),
        ("he", "<85"),
    ]
    pools: dict[tuple[str, str], list[str]] = {}
    for label, rec in genomes.items():
        row = meta.get(label)
        if row is None or row["dataset"] != "this study":
            continue
        if len(rec["seq"]) < MIN_LENGTH:
            continue
        identity = nearest.get(label, (None, "-"))[0]
        pools.setdefault((row["Clade"], band(identity)), []).append(label)
    chosen = []
    for genus, want in wanted:
        pool = sorted(pools.get((genus, want), []))
        if not pool:
            print(f"  WARNING: no held-out genome for {genus}/{want}", file=sys.stderr)
            continue
        # A full-length genome gives both windows their intended meaning; fall
        # back to a partial contig only when the cell has nothing else.
        full = [c for c in pool if complete.get(c) == "Full-length"]
        label = rng.choice(full or pool)
        identity, target = nearest.get(label, (None, "-"))
        chosen.append(
            {
                "genome": label,
                "genus": genus,
                "identity_band": want,
                "nearest_identity": round(identity, 4) if identity else "",
                "nearest_panel_genome": target,
                "completeness": complete.get(label, "unknown"),
                "length": len(genomes[label]["seq"]),
                "orf1_start": genomes[label]["orf1_start"],
                "last_cds_end": genomes[label]["last_cds_end"],
            }
        )
    return chosen


# ── Read simulation ──────────────────────────────────────────────────────────
def circular_slice(seq: str, start: int, length: int) -> str:
    """*length* bases from *start*, wrapping; anellovirus genomes are circular."""
    start %= len(seq)
    if start + length <= len(seq):
        return seq[start : start + length]
    return seq[start:] + seq[: length - (len(seq) - start)]


def cleavage_site(row: dict) -> int:
    """Where the transcript is cleaved and the untemplated poly-A tail begins.

    The CDS-only gene models carry no polyA annotation, so this takes the far
    edge of the 3' window — the furthest downstream point the model knows about.
    It is an estimate, and the only thing that rests on it is how long a planted
    3'-end read's templated body is, which ``truth.tsv`` records per read.
    """
    return row["last_cds_end"] + WINDOW_3P_LEN


def polya_read(seq: str, start: int, cleavage: int, rng: random.Random) -> tuple[str, int]:
    """A genuine 10x 3'-end read: ``[templated body][untemplated poly-A]``.

    Planting a clean genome slice here would make any poly-A or homopolymer gate
    look free, because no planted read would carry the tail such a gate removes
    (Biomni review, 2026-10-04). A read starting *start* bases in runs out of
    transcript at *cleavage* and is filled with A's from there; one starting far
    enough upstream never reaches the site and gets no tail at all. Returns the
    read and its templated body length, so recovery can be scored against it.
    """
    body_len = max(0, min(READ_LEN, cleavage - start))
    body = circular_slice(seq, start, body_len) if body_len else ""
    return body + "A" * (READ_LEN - body_len), body_len


def window_range(row: dict, window: str) -> tuple[int, int]:
    """``(first start, span)`` of read start positions for *window*."""
    length, orf1, last_end = row["length"], row["orf1_start"], row["last_cds_end"]
    if window == "5p":
        return orf1 - WINDOW_5P_LEN, WINDOW_5P_LEN
    if window == "3p":
        return last_end, WINDOW_3P_LEN
    return 0, length


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    out = list(seq)
    for i, base in enumerate(out):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in BASES if b != base])
    return "".join(out)


def umi(rng: random.Random) -> str:
    """A UMI STARsolo will keep: it discards homopolymer UMIs, and such a read
    then carries CB:Z:-/UB:Z:- and counts as no molecule at all."""
    while True:
        candidate = "".join(rng.choice(BASES) for _ in range(UMI_LEN))
        if len(set(candidate)) > 1 and not re.search(r"(.)\1{5,}", candidate):
            return candidate


def fastq(handle, name: str, seq: str) -> None:
    handle.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def plant_reads(selection, genomes, barcodes, rng, r1, r2, truth) -> int:
    """Write planted pairs and their ground-truth rows; return the read count."""
    n = 0
    for row in selection:
        seq = genomes[row["genome"]]["seq"]
        cleavage = cleavage_site(row)
        for window in ("5p", "3p", "uniform"):
            first, span = window_range(row, window)
            for level in LEVELS:
                for _ in range(level):
                    start = first + rng.randrange(span)
                    # Only the 3' window models the transcript end, so only it
                    # carries an untemplated tail; 5p and uniform reads are
                    # internal fragments and are fully templated.
                    if window == "3p":
                        read, body_len = polya_read(seq, start, cleavage, rng)
                    else:
                        read, body_len = circular_slice(seq, start, READ_LEN), READ_LEN
                    read = mutate(read, ERROR_RATE, rng)
                    cb, ub = rng.choice(barcodes), umi(rng)
                    name = f"plant_{n}"
                    fastq(r1, f"{name} 1", cb + ub)
                    fastq(r2, f"{name} 2", read)
                    truth.writerow(
                        [name, row["genome"], row["genus"], row["identity_band"],
                         window, level, cb, ub, start % len(seq),
                         body_len, READ_LEN - body_len]
                    )
                    n += 1
    return n


def artefact_read(kind: str, rng: random.Random) -> str:
    """One read of a measured artefact class (F-010, F-019, F-021)."""
    if kind == "polyG":  # two-colour no-signal reads
        return "G" * READ_LEN
    if kind == "polyA":  # poly-A tail capture
        pad = rng.randrange(READ_LEN - ARTEFACT_RUN)
        return "".join(rng.choice(BASES) for _ in range(pad)) + "A" * (READ_LEN - pad)
    if kind == "CAG":  # trinucleotide repeat
        return ("CAG" * READ_LEN)[:READ_LEN]
    tail = READ_LEN - len(TSO)  # TSO joined to a poly-A stretch
    return TSO + "A" * tail


def write_negative(out_dir: Path, cdna, genome_fa, barcodes, n_pairs, rng) -> dict:
    """Synthetic negative: host sequence plus artefact reads, no viral sequence."""
    # Deliberately artefact-rich: 20 % of reads are artefacts, far above any
    # real library, so a branch that calls them cannot pass by rarity.
    counts = {"host_cdna": int(n_pairs * 0.6), "host_genome": int(n_pairs * 0.2)}
    per_artefact = (n_pairs - sum(counts.values())) // 4
    for kind in ("polyG", "polyA", "CAG", "TSO_polyA"):
        counts[kind] = per_artefact
    cdna_names = [k for k, v in cdna.items() if len(v) >= READ_LEN]
    genome_names = [k for k, v in genome_fa.items() if len(v) >= READ_LEN]
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(out_dir / "R1.fastq.gz", "wt") as r1, gzip.open(
        out_dir / "R2.fastq.gz", "wt"
    ) as r2:
        for kind, count in counts.items():
            for _ in range(count):
                if kind == "host_cdna":
                    seq = cdna[rng.choice(cdna_names)]
                    start = rng.randrange(len(seq) - READ_LEN + 1)
                    read = mutate(seq[start : start + READ_LEN], HOST_ERROR_RATE, rng)
                elif kind == "host_genome":
                    seq = genome_fa[rng.choice(genome_names)]
                    start = rng.randrange(len(seq) - READ_LEN + 1)
                    read = mutate(seq[start : start + READ_LEN].upper(), HOST_ERROR_RATE, rng)
                    if "N" * 10 in read:  # unplaced/telomeric N runs are not sequence
                        continue
                else:
                    read = artefact_read(kind, rng)
                fastq(r1, f"neg_{n} 1", rng.choice(barcodes) + umi(rng))
                fastq(r2, f"neg_{n} 2", read)
                n += 1
    counts["written"] = n
    return counts


def copy_background(src_r1, src_r2, n_pairs, handles) -> int:
    """Stream the first *n_pairs* of the background into every handle pair.

    The first N is a contiguous block of the real library rather than a random
    sample; it is the same block in both arms, which is what the planted-minus-
    unplanted difference needs.
    """
    written = 0
    with gzip.open(src_r1, "rt") as f1, gzip.open(src_r2, "rt") as f2:
        while written < n_pairs:
            rec1 = [f1.readline() for _ in range(4)]
            rec2 = [f2.readline() for _ in range(4)]
            if not rec1[0] or not rec2[0]:
                break
            for r1, r2 in handles:
                r1.write("".join(rec1))
                r2.write("".join(rec2))
            written += 1
    return written


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--gbk", type=Path, required=True)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--completeness", type=Path, required=True)
    p.add_argument("--paf", type=Path, required=True)
    p.add_argument("--background-r1", type=Path, required=True)
    p.add_argument("--background-r2", type=Path, required=True)
    p.add_argument("--whitelist", type=Path, required=True)
    p.add_argument("--host-cdna", type=Path, required=True)
    p.add_argument("--host-genome", type=Path, required=True)
    p.add_argument("--n-background", type=int, default=5_000_000)
    p.add_argument("--seed", type=int, default=20261003)
    p.add_argument("--negative-only", action="store_true",
                   help="Rebuild only the synthetic negative, leaving the planted and "
                        "unplanted arms (and truth.tsv) untouched.")
    args = p.parse_args()

    rng = random.Random(args.seed)
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    if args.negative_only:
        barcodes = [ln.strip() for ln in open(args.whitelist, encoding="utf-8") if ln.strip()]
        counts = write_negative(
            out / "negative",
            read_fasta(args.host_cdna, HOST_RECORD_PREFIX),
            read_fasta(args.host_genome),
            rng.sample(barcodes, N_CELLS),
            args.n_background,
            rng,
        )
        print(json.dumps(counts, indent=2))
        return

    print("Reading held-out genomes and annotation ...")
    genomes = read_gbk(args.gbk)
    print(f"  {len(genomes)} annotated genomes")
    selection = select_genomes(genomes, args.metadata, args.completeness, args.paf, rng)
    with open(out / "selection.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(selection[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(selection)
    for row in selection:
        print(
            f"  {row['genus']:6s} {row['identity_band']:12s} "
            f"{row['completeness']:11s} {row['length']:5d}  {row['genome']}"
        )

    barcodes = [ln.strip() for ln in open(args.whitelist, encoding="utf-8") if ln.strip()]
    barcodes = rng.sample(barcodes, N_CELLS)

    planted, unplanted = out / "planted", out / "unplanted"
    planted.mkdir(exist_ok=True)
    unplanted.mkdir(exist_ok=True)
    print(f"Copying {args.n_background:,} background pairs into both arms ...")
    with gzip.open(planted / "R1.fastq.gz", "wt") as p1, gzip.open(
        planted / "R2.fastq.gz", "wt"
    ) as p2, gzip.open(unplanted / "R1.fastq.gz", "wt") as u1, gzip.open(
        unplanted / "R2.fastq.gz", "wt"
    ) as u2, open(out / "truth.tsv", "w", encoding="utf-8", newline="") as tf:
        n_background = copy_background(
            args.background_r1, args.background_r2, args.n_background, [(p1, p2), (u1, u2)]
        )
        truth = csv.writer(tf, delimiter="\t")
        truth.writerow(
            ["read_id", "genome", "genus", "identity_band", "window", "level", "cb", "umi",
             "pos", "body_len", "tail_len"]
        )
        print("Planting held-out reads into the planted arm ...")
        n_planted = plant_reads(selection, genomes, barcodes, rng, p1, p2, truth)

    print(f"Writing the synthetic negative ({args.n_background:,} pairs) ...")
    negative_counts = write_negative(
        out / "negative",
        read_fasta(args.host_cdna, HOST_RECORD_PREFIX),
        read_fasta(args.host_genome),
        barcodes,
        args.n_background,
        rng,
    )

    manifest = {
        "seed": args.seed,
        "read_length": READ_LEN,
        "error_rate": ERROR_RATE,
        "levels": list(LEVELS),
        "windows": ["5p", "3p", "uniform"],
        "polya_tail_on_3p": True,
        "n_background_pairs": n_background,
        "n_planted_reads": n_planted,
        "n_genomes": len(selection),
        "negative_composition": negative_counts,
        "background_r1": str(args.background_r1),
        "background_r2": str(args.background_r2),
    }
    (out / "plant_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
