#!/usr/bin/env python3
"""READS-01 smoke test on a real run: python scripts/check_reads_api.py <run_dir> <virus> <out_dir>."""
import sys

from viralscan.reads import extract_virus_reads

run_dir, virus, out_dir = sys.argv[1:4]
vr = extract_virus_reads(run_dir, [virus], out_dir, cores=4)[virus]
seqs = vr.sequences()
print(f"{virus}: {len(seqs)} reads in {vr.fasta}; lineage rows {len(vr.lineage())}")
