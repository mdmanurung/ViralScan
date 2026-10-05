"""COST-01: per-sample run cost of the combined and STAR two-step arms.

Reads `sacct` accounting for the 2026-09-28/29 full-depth runs on three real
libraries and normalises allocated and used core-hours to 50 M read pairs (the
VAL-01 D5 sample: 2,000 cells x 25 k reads). Read counts are kb `n_processed`
(combined) and STAR `Number of input reads` (two-step), both from the run dirs.

    python scripts/cost01_per_sample_cost.py <outdir>
"""

import csv
import subprocess
import sys
from pathlib import Path

PAIRS_PER_SAMPLE = 50_000_000
# job id -> (arm, dataset, chemistry, read pairs)
JOBS = {
    "25652878_0": ("combined", "SRR12682296", "10xv2", 127_045_580),
    "25652878_1": ("combined", "SRR20710641", "10xv3", 115_458_803),
    "25652878_2": ("combined", "SRR8315713", "dropseq", 99_862_579),
    "25666447_0": ("combined_corrected", "SRR12682296", "10xv2", 127_045_580),
    "25666447_2": ("combined_corrected", "SRR8315713", "dropseq", 99_862_579),
    "25652882_0": ("two_step", "SRR12682296", "10xv2", 127_045_580),
    "25652882_2": ("two_step", "SRR8315713", "dropseq", 99_862_579),
}


def _seconds(text: str) -> float:
    days, _, clock = text.rpartition("-")
    parts = [float(p) for p in clock.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    return (int(days) if days else 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def main() -> int:
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    ids = sorted({j.split("_")[0] for j in JOBS})
    raw = subprocess.run(
        [
            "sacct",
            "-j",
            ",".join(ids),
            "--parsable2",
            "--noheader",
            "-o",
            "JobID,AllocCPUS,Elapsed,CPUTimeRAW,TotalCPU,MaxRSS,State",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    (out / "sacct_snapshot.txt").write_text(raw)

    batch = {}
    for line in raw.splitlines():
        job, cpus, elapsed, cputime, totalcpu, maxrss, state = line.split("|")
        if job.endswith(".batch") and job[: -len(".batch")] in JOBS:
            batch[job[: -len(".batch")]] = (
                int(cpus),
                elapsed,
                int(cputime),
                totalcpu,
                maxrss,
                state,
            )
    missing = set(JOBS) - set(batch)
    if missing:
        sys.exit(f"sacct has no batch step for {sorted(missing)}")

    rows = []
    for job, (arm, srr, chem, pairs) in JOBS.items():
        cpus, elapsed, cputime, totalcpu, maxrss, state = batch[job]
        assert state == "COMPLETED", (job, state)
        scale = PAIRS_PER_SAMPLE / pairs
        rows.append(
            {
                "job": job,
                "arm": arm,
                "dataset": srr,
                "chemistry": chem,
                "read_pairs": pairs,
                "alloc_cpus": cpus,
                "elapsed": elapsed,
                "alloc_core_h": round(cputime / 3600, 3),
                "used_cpu_h": round(_seconds(totalcpu) / 3600, 3),
                "max_rss_gb": round(int(maxrss.rstrip("K")) / 1024**2, 1),
                "alloc_core_h_per_50M": round(cputime / 3600 * scale, 3),
                "used_cpu_h_per_50M": round(_seconds(totalcpu) / 3600 * scale, 3),
            }
        )
    with open(out / "per_sample_cost.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(str(v) for v in r.values()))
    return 0


if __name__ == "__main__":
    assert _seconds("01:17:36") == 4656 and _seconds("48:37.098") == 2917.098
    assert _seconds("1-00:00:01") == 86401
    sys.exit(main())
