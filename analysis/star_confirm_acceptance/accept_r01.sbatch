#!/usr/bin/env bash
# 3' v3 / forward-strand acceptance of the star_confirm rule on R01_SOC4 (RAW reads, 172 M pairs; the pilot STAR run on the
# same reads, no whitelist and unstranded, is the manual reference). Exercises the 10xv3 on-list, 12-nt UMI, --strand forward.
#SBATCH -J star_acc_r01
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -t 06:00:00
#SBATCH -p highmem
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/star_confirm_accept_2026-10/slurm_%x_%j.log
set -euo pipefail
M=/exports/archive/hg-funcgenom-research/mdmanurung
W=$M/viralscan_work; X=$W/trackB_explore_2026-10
export PATH=$M/conda/envs/viralscan_test_full/bin:$PATH PYTHONPATH=$M/vs_work/star-confirm/src
date; STAR --version
python $M/vs_work/star-confirm/analysis/star_confirm_acceptance/accept_ebv.py $W/star_confirm_accept_2026-10 --name R01_SOC4 \
  --technology 10xv3 --strand forward --r1 $X/inputs/R01_SOC4/R01_SOC4_R1.fastq.gz --r2 $X/inputs/R01_SOC4/R01_SOC4_R2.fastq.gz \
  --manual-bam $X/runs/R01_SOC4/star/Aligned.sortedByCoord.out.bam
date; echo STAR_ACCEPT_R01_DONE
