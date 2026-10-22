"""Spike-in detection-limit grid: ALL paths and parameters (CONFIG block) plus the conditions table (PLAN SPIKE-02).

``python grid.py [--out conditions.tsv]`` writes one row per condition: 5 levels x 3 seeds x 2 layouts per
background (60 spiked in total for two backgrounds) plus N_NULL_REPS unspiked replicates per background (66 rows).
Rows are arm-independent: the old index, the flank index and the STAR arm are all scored on the same condition.

A condition's background subsample replicate ``bg_rep`` is its spike seed (spiked) or its replicate number (null), so
spike seed s and null replicate s share one background and the null floor can be subtracted. Condition ids contain no
underscore: ``viralscan run`` derives the sample id from the R1 name up to the first underscore.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

# ----------------------------------------------------------------------------------------------- CONFIG
M = Path("/exports/archive/hg-funcgenom-research/mdmanurung")
W = M / "viralscan_work"
X = W / "trackB_explore_2026-10"
S = W / "star_hostvirus"
ROOT = Path(
    os.environ.get("SPIKE_GRID_ROOT", W / "spikein_grid_2026-10")
)  # all outputs; overridable
ENV_BIN = M / "conda/envs/viralscan_test_full/bin"  # kb, snakemake, STAR, samtools, Rscript, python
SRC = Path(__file__).resolve().parents[2] / "src"  # this worktree's viralscan package
HOST_INDEX = M / "ViralScan/ebv_latest_ref_2026-09-27/star_host_2.7.11b"
STAR_HV_INDEX = W / "star_confirm_accept_2026-10/index"  # holds star_hostvirus_index.json
ARMS = {  # kb arms; the STAR arm is the --star-confirm output of the FIRST arm (it depends on the host-filtered reads only)
    "old": dict(
        idx=M / "viral_ref_cat42d/build/panel.idx",
        t2g=M / "viral_ref_cat42d/build/panel.t2g",
        gtf=M / "viral_ref_cat42d/build/viral_panel.gtf",
        star_confirm=True,
    ),
    "flank": dict(
        idx=X / "ref_flank_v2/panel_flank.idx",
        t2g=X / "ref_flank_v2/panel_flank.t2g",
        gtf=X / "ref_flank_v2/viral_panel_flank.gtf",
        star_confirm=False,
    ),
}
SOURCES = {  # spike sources: STARsolo host+virus BAM, annotation, truth contig set (first = class reference)
    "HHV6B": dict(
        bam=S / "GSM6416022/Aligned.sortedByCoord.out.bam",
        gtf=S / "ref/hostvirus.gtf",
        contigs="VIRUS|AF157706.1",
        orient="reverse",  # 5' library
        max_nh=1,
    ),
    "EBV": dict(  # 72 percent of reads are NH=2 on type 1 + type 2 (decisions.md, SPIKE-02 / I1)
        bam=S / "GSM4796271/Aligned.sortedByCoord.out.bam",
        gtf=S / "ref/hostvirus.gtf",
        contigs="VIRUS|NC_007605.1,VIRUS|NC_009334.1",
        orient="forward",  # 3' library
        max_nh=2,
    ),
}
BACKGROUNDS = {
    "D03": dict(
        lib="SRR29320037",
        r1=X / "inputs/SRR29320037/SRR29320037_R1.fastq.gz",
        r2=X / "inputs/SRR29320037/SRR29320037_R2.fastq.gz",
        technology="10xv2",
        strand="reverse",
        orient="reverse",
        cb_len=16,
        umi_len=10,
        r1_len=26,
        r2_len=90,
        source="HHV6B",
    ),
    "R01": dict(
        lib="R01_SOC4",
        r1=X / "inputs/R01_SOC4/R01_SOC4_R1.fastq.gz",
        r2=X / "inputs/R01_SOC4/R01_SOC4_R2.fastq.gz",
        technology="10xv3",
        strand="forward",
        orient="forward",
        cb_len=16,
        umi_len=12,
        r1_len=28,
        r2_len=90,
        source="EBV",
    ),
}
# on-list per technology: resolved at run time with viralscan.chemistry.onlist_path(technology) (needs kb-python)
LEVELS = (3, 10, 30, 100, 300)  # spiked UMIs
SEEDS = (1, 2, 3)
LAYOUTS = {"c005": 5, "c100": 100}  # spiked cells
N_NULL_REPS = 3
BG_TARGET_PAIRS = (
    45_000_000  # seeded random subsample of each background (frac = target / counted total)
)
N_CALLED = 100  # frozen called-cell list = top-N on-list barcodes of the background replicate
CORES = 8
KB_MIN_MOLECULES = 3  # call rule: >= 3 kb molecules AND >= 1 STAR read on the same contig
STAR_MIN_READS = 1
# -------------------------------------------------------------------------------------------------------

COLUMNS = ["cond_id", "background", "source", "level", "seed", "layout", "cells", "bg_rep"]


def build_conditions() -> list[dict]:
    rows = []
    for bg, b in BACKGROUNDS.items():
        for level in LEVELS:
            for seed in SEEDS:
                for layout, cells in LAYOUTS.items():
                    rows.append(
                        dict(
                            cond_id=f"{bg}-L{level:03d}-s{seed}-{layout}",
                            background=bg,
                            source=b["source"],
                            level=level,
                            seed=seed,
                            layout=layout,
                            cells=cells,
                            bg_rep=seed,
                        )
                    )
        for rep in range(1, N_NULL_REPS + 1):
            rows.append(
                dict(
                    cond_id=f"{bg}-NULL-r{rep}",
                    background=bg,
                    source="",
                    level=0,
                    seed=0,
                    layout="",
                    cells=0,
                    bg_rep=rep,
                )
            )
    return rows


def load_conditions(path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--out", default=str(ROOT / "conditions.tsv"))
    a = p.parse_args(argv)
    rows = build_conditions()
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    (ROOT / "slurm").mkdir(parents=True, exist_ok=True)  # sbatch -o directory must exist
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(
        f"{len(rows)} conditions -> {a.out}; sbatch --array=1-{len(rows)}%6 grid.sbatch",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
