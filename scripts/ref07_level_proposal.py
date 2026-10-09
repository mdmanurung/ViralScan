#!/usr/bin/env python3
"""REF-07 follow-up: sensitivity control check and homology level proposal.

Two jobs, both read-only over the REF-07 table:

1. **Control check.** REF-06 showed spurious unique viral calls disappearing
   once GRCh38 entered the D-list. Every one of those genomes must appear in
   this table, because a host-derived read is exactly what a host-homologous
   locus produces. A miss means the alignment dropped a real locus -- most
   likely a multi-copy one, which minimap2's minimizer filtering discards
   first, and which is the most dangerous class to miss.

2. **Level proposal.** Cut-points are read off the **training** clusters only
   (protocol: the homology characterization is training-only), with
   low-complexity hits held out because `low_complexity` is its own VAL-01
   factor. Three levels including `none`, matching the protocol's 63-strata
   illustration: each extra level costs 84 samples x 2 twins.

Nothing here writes to protocol.yaml. The output is a DEF-00 input.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

#: Genomes whose spurious unique calls vanished under the genome D-list
#: (REF-06, 2026-10-05). Gene -> genome accession, resolved from combined.gtf.
CONTROLS = {
    "HHV-6 p23": "NC_001664.4",
    "MPXV gp132/gp028": "NC_003310.1",
    "MOCV gp001": "NC_001731.1",
    "HPV9 gp4": "NC_001596.1",
    "HPV77 L2": "Y15175.1",
    "EBV EBNA-2": "NC_007605.1",
    "CPXV159": "NC_003663.2",
    "HHV-1 gp00p39/p61": "NC_001806.2",
}

#: A hit below this cannot carry a confusable 90-nt read even in principle.
READ_LEN = 90


def load(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        rows = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", type=Path, required=True)
    args = p.parse_args()
    rows = load(args.table)
    for r in rows:
        for k in ("identity_best_90nt_window", "low_complexity_fraction",
                  "host_exon_fraction", "host_3p_fraction", "identity_block"):
            r[k] = float(r[k])
        for k in ("aligned_bases", "longest_exact_run", "shared_31mers"):
            r[k] = int(r[k])

    print(f"== table: {len(rows)} alignments, "
          f"{len({r['viral_record'] for r in rows})} viral genomes, "
          f"{len({r['cluster_id'] for r in rows})} clusters")
    part = {}
    for r in rows:
        part[r["partition"]] = part.get(r["partition"], 0) + 1
    print(f"   partitions: {part}")

    print("\n== control check (REF-06 genomes whose calls the D-list removed)")
    missing = []
    for label, acc in CONTROLS.items():
        hits = [r for r in rows if r["viral_record"] == acc]
        if not hits:
            missing.append(label)
            print(f"   MISS  {label:20s} {acc}")
            continue
        best = max(hits, key=lambda r: r["shared_31mers"])
        print(f"   ok    {label:20s} {acc}  hits={len(hits):4d}  "
              f"best: {best['aligned_bases']}nt id90={best['identity_best_90nt_window']:.3f} "
              f"exact={best['longest_exact_run']} 31mers={best['shared_31mers']} "
              f"exon={best['host_exon_fraction']:.2f}")
    print(f"   -> {len(CONTROLS) - len(missing)}/{len(CONTROLS)} present"
          + (f"; MISSING {missing}" if missing else ""))

    # Levels are proposed from training rows only, low-complexity excluded.
    pool = [r for r in rows
            if r["partition"] == "training"
            and r["low_complexity_fraction"] < 0.2
            and r["aligned_bases"] >= READ_LEN]
    print(f"\n== level pool: {len(pool)} training, non-low-complexity, >={READ_LEN}nt hits")
    if not pool:
        print("   empty: cannot propose levels")
        return

    exact = sorted(r["longest_exact_run"] for r in pool)
    ident = sorted(r["identity_best_90nt_window"] for r in pool)

    def pct(xs, q):
        return xs[min(len(xs) - 1, int(q * len(xs)))]

    print(f"   longest_exact_run  p50={pct(exact,.5)} p90={pct(exact,.9)} "
          f"p99={pct(exact,.99)} max={exact[-1]}")
    print(f"   id90               p50={pct(ident,.5):.3f} p90={pct(ident,.9):.3f} "
          f"p99={pct(ident,.99):.3f} max={ident[-1]:.3f}")

    # The regime boundary is the exact 31-mer: a locus with a >=31 exact run can
    # be pseudoaligned, one without it cannot, whatever its overall identity.
    with_kmer = [r for r in pool if r["shared_31mers"] > 0]
    without = [r for r in pool if r["shared_31mers"] == 0]
    print(f"\n   shared 31-mer present: {len(with_kmer)}  absent: {len(without)}")
    expressed = [r for r in with_kmer if r["host_exon_fraction"] > 0.5]
    three_p = [r for r in expressed if r["host_3p_fraction"] > 0]
    print(f"   of those with a 31-mer: {len(expressed)} are >50% exonic, "
          f"{len(three_p)} touch a transcript 3' end")

    print("\n== proposed levels (DEF-00 input, NOT frozen here)")
    print("   none      : no planted host-homologous read")
    print("   subkmer   : longest exact run < 31 -- confusable by alignment, "
          "invisible to pseudoalignment")
    print("   exactkmer : longest exact run >= 31 and >50% exonic -- the regime "
          "that produced the F-005 artefact")
    print(f"   candidate loci: subkmer={len(without)}  exactkmer={len(expressed)}")


if __name__ == "__main__":
    main()
