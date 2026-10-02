#!/usr/bin/env bash
# SW-13/SW-20 impact check: the combined-index runs for EBV and HSV-1 repeated
# on the same v1 index with barcode correction on (runs/combined/ are the
# uncorrected `-w None` baseline). Derived from slurm_quant_combined_array.sh.
# Strategy 1 (combined host+virus index) quantification, 3 datasets, as a
# SLURM job array. Depends on scripts/slurm_kbref_full_panel.sh having
# already built full_panel/index.idx + t2g.txt.
#
# Submit (after the kb-ref job, e.g. jobid 12345):
#   sbatch --dependency=afterok:12345 scripts/slurm_quant_combined_array.sh
#
#SBATCH -J vs_corrected
#SBATCH --array=0,2
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH -o logs/vs_corrected_%A_%a.log
#SBATCH -e logs/vs_corrected_%A_%a.err

set -euo pipefail

REPO=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
WORK=/exports/para-lipg-hpc/mdmanurung/ViralScan/ebv_latest_ref_2026-09-27
REF=$WORK/full_panel
INPUTS=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/reference_strategy

export PATH="$ENV/bin:$PATH"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

# dataset table: name, SRR, technology
NAMES=(ebv hhv6b hsv1)
SRRS=(SRR12682296 SRR20710641 SRR8315713)
TECHS=(10xv2 10xv3 dropseq)

i=$SLURM_ARRAY_TASK_ID
NAME=${NAMES[$i]}
SRR=${SRRS[$i]}
TECH=${TECHS[$i]}

OUT=$WORK/runs/combined_corrected/$SRR
mkdir -p "$OUT"

R1=$INPUTS/$SRR/${SRR}_1.fastq.gz
R2=$INPUTS/$SRR/${SRR}_2.fastq.gz

echo "=== dataset: $NAME ($SRR, $TECH) — combined strategy ==="
git -C "$REPO" rev-parse HEAD
git -C "$REPO" status --short

viralscan \
  -o "$OUT" \
  -s1 "$R1" -s2 "$R2" \
  -i "$REF/index.idx" -t "$REF/t2g.txt" \
  -gtf "$REF/viral_final.gtf" \
  -x "$TECH" \
  -c 8 \
  --cell-calling knee \
  --anellovirus-gene-ids \
  --yes \
  --verbose

echo "Done: $NAME. Results: $OUT/results/viral_summary.tsv"
