#!/usr/bin/env bash
# EXPLORATORY (PLAN DEF-02 sub-item, F-020): rerun the 5' positives with an explicit
# kallisto strand. Each library keeps its baseline index, so strand is the only
# change. Every run writes to a fresh -o: resume refuses a strand mismatch. Sets no
# default; the default stays unset until DEF-01 (unstranded re-admits poly-G).
#
# CODE must point at a detached worktree pinned at a recorded SHA (PLAN coordination
# rule 5), with the gitignored src/viralscan/data/*.gtf copied in:
#   sbatch --export=ALL,CODE=/path/to/pinned/worktree scripts/slurm_quant_5p_strand_rerun.sh
# Tasks: 0-1 covid x213 {reverse,unstranded}; 2-3 covid x216; 4-5 HHV-6B SRR20710641.
#SBATCH -J vs_5p_strand
#SBATCH --array=0-5
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=2-00:00:00
#SBATCH -o logs/vs_5p_strand_%A_%a.log
#SBATCH -e logs/vs_5p_strand_%A_%a.err
set -euo pipefail
: "${CODE:?set CODE to the pinned worktree}"
R_BIN=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin
STRANDS=(reverse unstranded)
i=$SLURM_ARRAY_TASK_ID
STRAND=${STRANDS[$((i % 2))]}
DATASET=$((i / 2))
export PYTHONPATH="$CODE/src"
echo "code_sha $(git -C "$CODE" rev-parse HEAD) strand $STRAND task $i"

if [ "$DATASET" -lt 2 ]; then
  # covid: index built with kallisto 0.51.1 (test_viralscan env), CellRanger whitelist.
  ENV=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
  C=/exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan
  R=$C/viralscan_ref
  SAMPLES=(LUM-SJ-x213-g LUM-SJ-x216-g); s=${SAMPLES[$DATASET]}
  OUT=$C/results_v3_strand_$STRAND/$s
  ARGS=(-s1 "$C/data/$s/${s}_merged_R1.fastq.gz" -s2 "$C/data/$s/${s}_merged_R2.fastq.gz"
        -i "$R/index.idx" -t "$R/t2g.txt"
        -gtf "$R/viral/viral_whole_genome.gtf,/exports/para-lipg-hpc/mdmanurung/ViralScan/references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf"
        -w "$R/cellranger_whitelist.txt" -x 10xv3 -c 16)
else
  # HHV-6B: 10x 5' v1/v2 on the v1 full panel (as slurm_quant_hhv6b_5p.sh).
  ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
  WORK=/exports/para-lipg-hpc/mdmanurung/ViralScan/ebv_latest_ref_2026-09-27
  REF=$WORK/full_panel
  IN=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/reference_strategy/SRR20710641
  OUT=$WORK/runs/combined_corrected_5p_$STRAND/SRR20710641
  ARGS=(-s1 "$IN/SRR20710641_1.fastq.gz" -s2 "$IN/SRR20710641_2.fastq.gz"
        -i "$REF/index.idx" -t "$REF/t2g.txt" -gtf "$REF/viral_final.gtf"
        -x 10xv2 -c 16 --anellovirus-gene-ids)
fi
export PATH="$ENV/bin:$PATH:$R_BIN"
mkdir -p "$OUT"
git -C "$CODE" rev-parse HEAD > "$OUT/code_sha.txt"
python -m viralscan.menu -o "$OUT" "${ARGS[@]}" --strand "$STRAND" \
  --cell-calling emptydrops --yes --verbose
echo "Done: task $i -> $OUT"
