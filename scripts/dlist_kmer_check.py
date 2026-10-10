#!/usr/bin/env python3
"""DLIST-01 / I0: are the viral k-mers that the false-positive reads share with the target virus in a D-list?

    python scripts/dlist_kmer_check.py <viral_reads.fasta> <virus.fa> <dlist.fa|cdna.fa[.gz]> [...]

Prints the k-mers (k=31, both strands) shared by the reads and the virus FASTA, and for each D-list how many
of them occur in it. A k-mer absent from the D-list stays in the index, so D-listing cannot remove that call.
"""
import gzip
import sys

K = 31


def records(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, parts = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts).upper()


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def kmers(seq):
    return {seq[i:i + K] for i in range(len(seq) - K + 1)}


reads_fa, virus_fa, *dlists = sys.argv[1:]
virus = set()
for _, s in records(virus_fa):
    virus |= kmers(s) | kmers(rc(s))
shared = set()
for _, s in records(reads_fa):
    shared |= kmers(s) & virus
print(f"{len(shared)} distinct {K}-mers shared by the reads and the virus FASTA")
for k in sorted(shared)[:5]:
    print("  ", k)
for path in dlists:
    found = set()
    for _, s in records(path):
        found |= {k for k in shared if k in s or rc(k) in s}
    print(f"{path}: {len(found)}/{len(shared)} shared k-mers present")
