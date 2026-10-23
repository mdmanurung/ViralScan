"""Real reproduction of STAR 2.7.11b's exit-139 Solo crash and the star_run fallback (PLAN STAR-02 / STAR-ACCEPT-01).

Random reads that map to neither host nor virus, with a chemistry that has no on-list (Drop-seq), make STAR run
``--soloCBwhitelist None`` with zero mapped reads. star_hostvirus must recover (dummy whitelist retry, zero mapped
accepted), leave ``star139_fallback.txt``, and still write an (empty) per-virus table.
"""

import argparse
import gzip
import random
from pathlib import Path

from viralscan import star_hostvirus as sh


def write_random_pairs(
    out: Path, n: int, cb_umi: int, r2_len: int, seed: int = 7
) -> tuple[Path, Path]:
    rng = random.Random(seed)
    r1, r2 = out / "R1.fastq.gz", out / "R2.fastq.gz"
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        for i in range(n):
            f1.write(f"@r{i}\n{''.join(rng.choices('ACGT', k=cb_umi))}\n+\n{'F' * cb_umi}\n")
            f2.write(f"@r{i}\n{''.join(rng.choices('ACGT', k=r2_len))}\n+\n{'F' * r2_len}\n")
    return r1, r2


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("out", type=Path)
    p.add_argument("--index", type=Path, required=True)
    p.add_argument("--technology", default="dropseq")
    p.add_argument("--n", type=int, default=2000)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    cb_len, umi_len = sh.cb_umi_geometry(a.technology)
    r1, r2 = write_random_pairs(a.out, a.n, cb_len + umi_len, 60)
    work = a.out / "star"
    bam = sh.align_host_unmapped(r1, r2, a.index, work, a.technology, threads=4, prune=True)
    marker = work / "star139_fallback.txt"
    print("pruned bam:", bam)
    print("fallback marker present:", marker.exists())
    print("summary:", (work / "star_support_by_virus.tsv").read_text().strip().splitlines())
    return 0 if bam.exists() else 1


if __name__ == "__main__":
    raise SystemExit(main())
