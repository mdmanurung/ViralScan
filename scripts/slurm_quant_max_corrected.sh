#!/usr/bin/env bash
# EXPLORATORY (grill R2.10 = D, 2026-09-29): F-017 redo with barcode correction on
# (SW-13/SW-20 fixed). Dilution regression for the max viral panel: rerun the combined (host+virus
# index) strategy on the two datasets with a known answer and compare against
# the 2,215-genome panel runs in ebv_latest_ref_2026-09-27/runs/combined/.
# Pass criteria (plan): EBV-1 within +-5% of 906,202; EBV-2 not a separate call
# >1% of EBV-1; HSV-1 within +-5% of 32,404; HHV-2 still possible_em_bleed;
# anellovirus ~0.
#
# Submit (after the max-panel kb-ref job, e.g. jobid 12345):
#   sbatch --dependency=afterok:12345 scripts/slurm_quant_max_regression.sh
#
#SBATCH -J vs_max_corr
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH -o logs/vs_max_corr_%A_%a.log
#SBATCH -e logs/vs_max_corr_%A_%a.err

set -euo pipefail

REPO=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
REF=/exports/para-lipg-hpc/mdmanurung/ViralScan/viral_panel_max_2026-09-28
INPUTS=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/reference_strategy

export PATH="$ENV/bin:$PATH"
# emptyDrops (cell calling, decided 2026-10-01, PLAN SW-23) needs Rscript with DropletUtils.
export PATH="$PATH:/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

NAMES=(ebv hsv1)
SRRS=(SRR12682296 SRR8315713)
TECHS=(10xv2 dropseq)

i=$SLURM_ARRAY_TASK_ID
NAME=${NAMES[$i]}
SRR=${SRRS[$i]}
TECH=${TECHS[$i]}

OUT=$REF/runs/combined_max_corrected/$SRR
mkdir -p "$OUT"

echo "=== dataset: $NAME ($SRR, $TECH) — combined strategy, max panel ==="
git -C "$REPO" rev-parse HEAD
git -C "$REPO" status --short

viralscan \
  -o "$OUT" \
  -s1 "$INPUTS/$SRR/${SRR}_1.fastq.gz" -s2 "$INPUTS/$SRR/${SRR}_2.fastq.gz" \
  -i "$REF/index_max.idx" -t "$REF/t2g_max.txt" \
  -gtf "$REF/viral_max.gtf" \
  -x "$TECH" \
  -c 8 \
  --cell-calling emptydrops \
  --anellovirus-gene-ids \
  --yes \
  --verbose

echo "Done: $NAME. Results: $OUT/results/viral_summary.tsv"
