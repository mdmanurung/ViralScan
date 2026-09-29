#!/usr/bin/env bash
# EXPLORATORY (grill R2.6/R2.10, 2026-09-29): strand test for 5' libraries.
# kb count on the first 4M read pairs under --strand forward/reverse/unstranded;
# compares pseudoalignment rate and host/viral UMI totals. kallisto's 10x default
# is forward; 5' R2 reads run antisense, so reverse/unstranded should win if the
# default loses reads.
#SBATCH -J vs_strand
#SBATCH --array=0-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH -o logs/vs_strand_%A_%a.log
#SBATCH -e logs/vs_strand_%A_%a.err
set -euo pipefail
P=/exports/para-lipg-hpc/mdmanurung/ViralScan
C=$P/covid_viralscan; V1=$P/ebv_latest_ref_2026-09-27/full_panel
case $SLURM_ARRAY_TASK_ID in
  0) NAME=SRR20710641; ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
     R1=$P/benchmark_inputs/reference_strategy/SRR20710641/SRR20710641_1.fastq.gz
     R2=$P/benchmark_inputs/reference_strategy/SRR20710641/SRR20710641_2.fastq.gz
     IDX=$V1/index.idx; T2G=$V1/t2g.txt; TECH=10xv2; WL=() ;;
  1|2) S=(x x LUM-SJ-x213-g LUM-SJ-x216-g); NAME=${S[$((SLURM_ARRAY_TASK_ID+1))]}
     ENV=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
     R1=$C/data/$NAME/${NAME}_merged_R1.fastq.gz; R2=$C/data/$NAME/${NAME}_merged_R2.fastq.gz
     IDX=$C/viralscan_ref/index.idx; T2G=$C/viralscan_ref/t2g.txt; TECH=10xv3
     WL=(-w "$C/viralscan_ref/cellranger_whitelist.txt") ;;
esac
export PATH="$ENV/bin:$PATH"
OUT=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/strand_test/$NAME
mkdir -p "$OUT"; cd "$OUT"
set +o pipefail  # zcat gets SIGPIPE when head stops reading
zcat "$R1" | head -n 16000000 | gzip -1 > sub_R1.fastq.gz
zcat "$R2" | head -n 16000000 | gzip -1 > sub_R2.fastq.gz
set -o pipefail
[ "$(zcat sub_R1.fastq.gz | wc -l)" -eq 16000000 ] || { echo "subsample short" >&2; exit 1; }
for st in forward reverse unstranded; do
  kb count -i "$IDX" -g "$T2G" -x "$TECH" "${WL[@]}" --strand "$st" -t 8 -o "kb_$st" \
    sub_R1.fastq.gz sub_R2.fastq.gz
done
python - "$T2G" <<'PY'
import json, sys, csv
from pathlib import Path
import numpy as np, scipy.io
rows = []
for st in ("forward", "reverse", "unstranded"):
    d = Path(f"kb_{st}")
    info = json.load(open(d / "run_info.json"))
    m = scipy.io.mmread(d / "counts_unfiltered" / "cells_x_genes.mtx").tocsc()
    genes = [l.strip() for l in open(d / "counts_unfiltered" / "cells_x_genes.genes.txt")]
    tot = np.asarray(m.sum(axis=0)).ravel()
    host = sum(t for g, t in zip(genes, tot) if g.startswith("ENSG"))
    viral = {g: t for g, t in zip(genes, tot) if not g.startswith("ENSG") and t > 0}
    top = sorted(viral.items(), key=lambda kv: -kv[1])[:5]
    rows.append({"strand": st, "n_processed": info.get("n_processed"),
                 "p_pseudoaligned": info.get("p_pseudoaligned"),
                 "host_umis": int(host), "viral_umis": int(sum(viral.values())),
                 "top_viral_genes": ";".join(f"{g}={int(t)}" for g, t in top)})
with open("strand_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(open("strand_summary.tsv").read())
PY
