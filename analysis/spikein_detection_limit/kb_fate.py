"""Per-spike fate in a kb run: pseudoaligned or not, and how many genes its equivalence class spans.

Reads ``output.unfiltered.bus`` (via ``bustools text``), ``matrix.ec`` and ``transcripts.txt`` of a ``kb count``
run, plus the transcript-to-gene map, and the spike truth table. The gene-level matrix drops UMIs whose class
spans several genes, so ``multi_gene`` is the mechanism behind a spike that pseudoaligns but is not counted.
"""

import argparse
import collections
import csv
import subprocess
from pathlib import Path


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--kb-dir", required=True)
    p.add_argument("--t2g", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--bustools", default="bustools")
    a = p.parse_args(argv)
    kb = Path(a.kb_dir)
    with open(a.truth) as fh:
        truth = {(r["cb"], r["umi"]): r["class"] for r in csv.DictReader(fh, delimiter="\t")}
    cbs = {cb for cb, _ in truth}
    tx = (kb / "transcripts.txt").read_text().split()
    t2g = {}
    with open(a.t2g) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            t2g[f[0]] = f[1]
    ec_genes = {}
    with open(kb / "matrix.ec") as fh:
        for ln in fh:
            ec, ts = ln.split()
            ec_genes[int(ec)] = {t2g.get(tx[int(t)], "?") for t in ts.split(",")}
    bus = kb / "output.unfiltered.bus"
    fate = {}
    proc = subprocess.Popen([a.bustools, "text", "-p", str(bus)], stdout=subprocess.PIPE, text=True)
    for ln in proc.stdout:
        cb, umi, ec, _ = ln.split()
        if (cb, umi) in truth:
            fate[(cb, umi)] = ec_genes[int(ec)]
    proc.wait()
    out = collections.defaultdict(collections.Counter)
    for key, cls in truth.items():
        genes = fate.get(key)
        if genes is None:
            k = "not_pseudoaligned"
        else:
            viral = {g for g in genes if not g.startswith("ENSG")}
            k = (
                "host_only"
                if not viral
                else "viral_single_gene"
                if len(genes) == 1
                else "multi_gene"
            )
        out[cls][k] += 1
        out["ALL"][k] += 1
    keys = ("viral_single_gene", "multi_gene", "host_only", "not_pseudoaligned")
    print("class\tn_spiked\t" + "\t".join(keys))
    for cls in sorted(out):
        print(f"{cls}\t{sum(out[cls].values())}\t" + "\t".join(str(out[cls][k]) for k in keys))


if __name__ == "__main__":
    main()
