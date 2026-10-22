"""Score a spike-in run against its truth table (exploratory; see spikein.py).

``star``        fate of each spiked read in a STARsolo host+virus BAM, per class: any alignment on any truth contig
                (headline), primary on a virus / another virus / host, or absent (unmapped).
``hostfilter``  how many spiked reads survive the real STAR host filter (``host_filtered/R2.fastq.gz``).
``kb``          viral UMIs in the spike barcodes, spiked run minus the unspiked run on the same background.
"""

import argparse
import collections
import csv
import gzip
import itertools
import json
import subprocess
import sys
from pathlib import Path


def truth(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def by_class(rows, hit):
    n = collections.Counter(r["class"] for r in rows)
    k = collections.Counter(r["class"] for r in rows if r["name"] in hit)
    return n, k


def truth_contigs(contig, rows, truth_path):
    """Truth contig set: --contig (comma list), else the truth table's ``contigs`` column, else the sibling meta json."""
    if contig:
        return set(contig.split(","))
    if rows and rows[0].get("contigs"):
        return set(rows[0]["contigs"].split(","))
    meta = Path(str(truth_path).removesuffix(".truth.tsv") + ".json")
    if meta.exists() and json.loads(meta.read_text()).get("contig"):
        return set(json.loads(meta.read_text())["contig"].split(","))
    sys.exit(
        "score: no truth contig (pass --contig or use a truth table / meta json from spikein.py)"
    )


def star_fates(lines, contigs):
    """(primary-fate per spiked read, names with ANY alignment, secondary included, on ANY of *contigs*)."""
    contigs = set(contigs.split(",")) if isinstance(contigs, str) else set(contigs)
    fate, on_contig = {}, set()
    for line in lines:
        if not line.startswith("SPIKE"):
            continue
        f = line.split("\t", 4)
        if f[2] in contigs:
            on_contig.add(f[0])
        if int(f[1]) & 0x100:
            continue
        fate[f[0]] = (
            "source"
            if f[2] in contigs
            else "other_virus"
            if f[2].startswith("VIRUS|")
            else "host"
            if f[2].startswith("HOST|")
            else "other"
        )
    return fate, on_contig


def cmd_star(a):
    rows = truth(a.truth)
    p = subprocess.Popen(["samtools", "view", a.bam], stdout=subprocess.PIPE, text=True)
    fate, on_contig = star_fates(p.stdout, truth_contigs(a.contig, rows, a.truth))
    p.wait()
    # n_source_contig_any is the headline column (pre-set rule: >= 1 STAR read on any truth contig);
    # n_virus = primary on any VIRUS| contig (source or other).
    print("class\tn_spiked\tn_source_contig_any\tn_virus\tn_other_virus\tn_host\tn_absent")
    for c in sorted({r["class"] for r in rows}) + ["ALL"]:
        names = [r["name"] for r in rows if c == "ALL" or r["class"] == c]
        cnt = collections.Counter(fate.get(n, "absent") for n in names)
        print(
            f"{c}\t{len(names)}\t{sum(n in on_contig for n in names)}\t"
            f"{cnt['source'] + cnt['other_virus']}\t{cnt['other_virus']}\t{cnt['host']}\t{cnt['absent']}"
        )


def cmd_hostfilter(a):
    rows = truth(a.truth)
    with gzip.open(a.filtered_r2, "rt") as fh:
        kept = {
            ln[1:].split()[0] for ln in itertools.islice(fh, 0, None, 4) if ln.startswith("@SPIKE")
        }
    n, k = by_class(rows, kept)
    print("class\tn_spiked\tn_survive_host_filter")
    for c in sorted(n):
        print(f"{c}\t{n[c]}\t{k[c]}")
    print(f"ALL\t{sum(n.values())}\t{sum(k.values())}")


def viral_umis_in(counts_dir, cbs):
    import scipy.io

    m = scipy.io.mmread(counts_dir + "/cells_x_genes.mtx").tocsr()
    genes = Path(counts_dir, "cells_x_genes.genes.txt").read_text().split()
    bcs = Path(counts_dir, "cells_x_genes.barcodes.txt").read_text().split()
    rows = [i for i, b in enumerate(bcs) if b in cbs]
    vi = [i for i, g in enumerate(genes) if not g.startswith("ENSG")]
    return int(m[rows][:, vi].sum()) if rows and vi else 0


def cmd_kb(a):
    cbs = {r["cb"] for r in truth(a.truth)}
    n = len(truth(a.truth))
    s, u = viral_umis_in(a.spiked_dir, cbs), viral_umis_in(a.unspiked_dir, cbs)
    print(
        f"n_spiked\tviral_umi_spiked\tviral_umi_unspiked\trecovered\nspiked={n}\t{s}\t{u}\t{s - u}"
    )


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("star")
    s.add_argument("--bam", required=True)
    s.add_argument("--truth", required=True)
    s.add_argument(
        "--contig",
        default=None,
        help="truth contig or comma-separated set; default: read from the truth table / meta json",
    )
    h = sub.add_parser("hostfilter")
    h.add_argument("--filtered-r2", required=True)
    h.add_argument("--truth", required=True)
    k = sub.add_parser("kb")
    k.add_argument("--spiked-dir", required=True, help="kb output dir's counts_unfiltered/")
    k.add_argument("--unspiked-dir", required=True)
    k.add_argument("--truth", required=True)
    a = p.parse_args(argv)
    {"star": cmd_star, "hostfilter": cmd_hostfilter, "kb": cmd_kb}[a.cmd](a)


if __name__ == "__main__":
    main()
