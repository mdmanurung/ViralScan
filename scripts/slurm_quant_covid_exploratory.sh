#!/usr/bin/env bash
# EXPLORATORY (grill R2.10 = D, 2026-09-29): covid x213/x216 re-run with current
# code on the published index/t2g/GTF/CellRanger whitelist. R1 is 28 bp with no
# TSO (16 bp barcode + 12 bp UMI), so -x 10xv3 geometry is correct; correction
# was already on (user whitelist). Strand stays kallisto's forward default until
# --strand exists; no read filter yet (R2.0), so the F-019 poly-G signal is
# expected to reproduce. Uses test_viralscan (kallisto 0.51.1) because the index
# was built with it.
#SBATCH -J vs_covid_expl
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=2-00:00:00
#SBATCH -o logs/vs_covid_expl_%A_%a.log
#SBATCH -e logs/vs_covid_expl_%A_%a.err
set -euo pipefail
REPO=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
ENV=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
C=/exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan
R=$C/viralscan_ref
export PATH="$ENV/bin:$PATH"; export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
S=(LUM-SJ-x213-g LUM-SJ-x216-g); s=${S[$SLURM_ARRAY_TASK_ID]}
OUT=$C/results_v3_exploratory/$s; mkdir -p "$OUT"
git -C "$REPO" rev-parse HEAD
python -m viralscan.menu -o "$OUT" \
  -s1 "$C/data/$s/${s}_merged_R1.fastq.gz" -s2 "$C/data/$s/${s}_merged_R2.fastq.gz" \
  -i "$R/index.idx" -t "$R/t2g.txt" \
  -gtf "$R/viral/viral_whole_genome.gtf,/exports/para-lipg-hpc/mdmanurung/ViralScan/references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf" \
  -w "$R/cellranger_whitelist.txt" -x 10xv3 -c 16 --cell-calling knee --yes --verbose
