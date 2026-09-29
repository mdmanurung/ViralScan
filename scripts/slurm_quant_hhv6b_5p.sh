#!/usr/bin/env bash
# EXPLORATORY (grill R2.10 = D, 2026-09-29): HHV-6B SRR20710641 is 10x 5' v1/v2
# (TSO at R1 base 27 -> 16 bp barcode + 10 bp UMI), but runs/combined/ ran it as
# -x 10xv3. Re-run the combined strategy on the same v1 index as -x 10xv2.
# The strand is kallisto's 10x default (forward); the strand test
# (slurm_strand_test.sh) measures whether that loses 5' reads.
#SBATCH -J vs_hhv6b_5p
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH -o logs/vs_hhv6b_5p_%j.log
#SBATCH -e logs/vs_hhv6b_5p_%j.err
set -euo pipefail
REPO=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
WORK=/exports/para-lipg-hpc/mdmanurung/ViralScan/ebv_latest_ref_2026-09-27
REF=$WORK/full_panel
IN=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/reference_strategy/SRR20710641
export PATH="$ENV/bin:$PATH"; export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
OUT=$WORK/runs/combined_corrected_5p/SRR20710641; mkdir -p "$OUT"
git -C "$REPO" rev-parse HEAD
viralscan -o "$OUT" -s1 "$IN/SRR20710641_1.fastq.gz" -s2 "$IN/SRR20710641_2.fastq.gz" \
  -i "$REF/index.idx" -t "$REF/t2g.txt" -gtf "$REF/viral_final.gtf" \
  -x 10xv2 -c 8 --cell-calling knee --anellovirus-gene-ids --yes --verbose
