#!/usr/bin/env bash
# Strategy 2 (STARsolo host-filter, then virus-only index) quantification,
# 3 datasets, as a SLURM job array. Depends on scripts/slurm_kbref_full_panel.sh
# having already built full_panel/virus_index.idx + virus_t2g.txt.
#
# v3 only supports --host-filter starsolo (kallisto host-filter removed).
#
# Submit (after the kb-ref job, e.g. jobid 12345):
#   sbatch --dependency=afterok:12345 scripts/slurm_quant_twostep_array.sh
#
#SBATCH -J vs_twostep
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH -o logs/vs_twostep_%A_%a.log
#SBATCH -e logs/vs_twostep_%A_%a.err

set -euo pipefail

REPO=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
WORK=/exports/para-lipg-hpc/mdmanurung/ViralScan/ebv_latest_ref_2026-09-27
REF=$WORK/full_panel
INPUTS=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/reference_strategy

export PATH="$ENV/bin:$PATH"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

HOST_STAR_DIR=/exports/para-lipg-hpc/mdmanurung/ViralScan/ebv_latest_ref_2026-09-27/star_host_2.7.11b

NAMES=(ebv hhv6b hsv1)
SRRS=(SRR12682296 SRR20710641 SRR8315713)
TECHS=(10xv2 10xv3 dropseq)

i=$SLURM_ARRAY_TASK_ID
NAME=${NAMES[$i]}
SRR=${SRRS[$i]}
TECH=${TECHS[$i]}

OUT=$WORK/runs/two_step/$SRR
mkdir -p "$OUT"

R1=$INPUTS/$SRR/${SRR}_1.fastq.gz
R2=$INPUTS/$SRR/${SRR}_2.fastq.gz

echo "=== dataset: $NAME ($SRR, $TECH) — two-step (STARsolo host-filter + virus-only index) ==="
git -C "$REPO" rev-parse HEAD
git -C "$REPO" status --short
echo "STAR version:"
STAR --version

viralscan \
  -o "$OUT" \
  -s1 "$R1" -s2 "$R2" \
  -i "$REF/virus_index.idx" -t "$REF/virus_t2g.txt" \
  -gtf "$REF/viral_final.gtf" \
  -x "$TECH" \
  -c 16 \
  --host-filter starsolo --host-index "$HOST_STAR_DIR" \
  --cell-calling knee \
  --anellovirus-gene-ids \
  --yes \
  --verbose

echo "Done: $NAME. Results: $OUT/results/viral_summary.tsv"
