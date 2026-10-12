"""Tiny real-style reference inputs for the truth-panel generator's AUD-18 options.

Written into ``tmp_path`` by the tests (nothing is committed): a multi-line, mixed-case FASTA of
protocol-named viruses, a gene-span TSV (overlapping genes on both strands) and a barcode on-list.
Four viruses carry two complete protocol pairs (HSV-1/HSV-2, EBV/KSHV); HHV-6A has no partner in the
FASTA, so it is an absent (negative) virus.
"""

from __future__ import annotations

import random
from pathlib import Path

VIRUSES = ("EBV", "HHV-6A", "HSV-1", "HSV-2", "KSHV")
PLANTED = ("EBV", "HSV-1", "HSV-2", "KSHV")
LENGTH = 1500
#: 1-based inclusive; G1/G2 overlap on "+", G3 is on "-" and overlaps G2 on the opposite strand
SPANS = (("G1", 1, 700, "+"), ("G2", 601, 1300, "+"), ("G3", 1101, 1500, "-"))


def write_real_inputs(out: Path, n_barcodes: int = 600) -> tuple[Path, Path, Path]:
    rng = random.Random(20261011)
    out.mkdir(parents=True, exist_ok=True)
    fasta, spans, onlist = out / "viral.fa", out / "spans.tsv", out / "onlist.txt"
    with open(fasta, "w") as handle:
        for index, virus in enumerate(VIRUSES):
            seq = "".join(rng.choice("ACGT") for _ in range(LENGTH))
            seq = seq.lower() if index == 0 else seq  # soft-masked input is upper-cased
            handle.write(f">{virus} synthetic stand-in for a RefSeq record\n")
            handle.writelines(seq[i : i + 60] + "\n" for i in range(0, LENGTH, 60))
    lines = ["virus_id\tgene_id\tstart\tend\tstrand"]
    lines += [f"{v}\t{v}_{g}\t{a}\t{b}\t{sign}" for v in VIRUSES for g, a, b, sign in SPANS]
    spans.write_text("\n".join(lines) + "\n")
    barcodes: set[str] = set()
    while len(barcodes) < n_barcodes:
        barcodes.add("".join(rng.choice("ACGT") for _ in range(16)))
    onlist.write_text("\n".join(sorted(barcodes)) + "\n")
    return fasta, spans, onlist


def write_background(
    out: Path, n_cells: int = 20, n_empty: int = 8, n_stragglers: int = 30, r1_len: int = 28
) -> tuple[Path, Path]:
    """A tiny stand-in for a real 10x v3 library (D1): FASTQ pair with a clear cell / empty split.

    ``n_cells`` barcodes carry 40 reads, ``n_empty`` carry 3 and the stragglers 1, so the count rank
    puts the cells first. Qualities are not all ``I``, so a test can tell they were kept.
    """
    import gzip

    rng = random.Random(20261012)
    out.mkdir(parents=True, exist_ok=True)
    barcodes: list[str] = []
    while len(barcodes) < n_cells + n_empty + n_stragglers:
        barcode = "".join(rng.choice("ACGT") for _ in range(16))
        if barcode not in barcodes:
            barcodes.append(barcode)
    depths = [40] * n_cells + [3] * n_empty + [1] * n_stragglers
    paths = out / "bg_R1.fastq.gz", out / "bg_R2.fastq.gz"
    with gzip.open(paths[0], "wt") as f1, gzip.open(paths[1], "wt") as f2:
        n = 0
        for barcode, depth in zip(barcodes, depths):
            for _ in range(depth):
                umi = "".join(rng.choice("ACGT") for _ in range(r1_len - 16))
                f1.write(f"@bg{n} 1:N:0\n{barcode}{umi}\n+\n{'F' * (r1_len - 4)}:,#\n")
                f2.write(f"@bg{n} 2:N:0\n{_rand_read(rng)}\n+\n{_rand_qual(rng)}\n")
                n += 1
    return paths


def _rand_read(rng: random.Random, n: int = 91) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _rand_qual(rng: random.Random, n: int = 91) -> str:
    return "".join(rng.choice("F:,#") for _ in range(n))
