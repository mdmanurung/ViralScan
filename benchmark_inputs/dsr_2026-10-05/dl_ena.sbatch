#!/usr/bin/env bash
#SBATCH -J dsr_ena
#SBATCH --array=0-13%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH -o /exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/dsr_2026-10-05/logs/%x_%A_%a.log
# DSR plan row-2 downloads: ENA FASTQ files with the ENA md5 checked after download.
set -euo pipefail
D=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/dsr_2026-10-05
IFS=$'\t' read -r DS RUN URL MD5 < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $D/ena_files.tsv)
O=$D/$DS/$RUN; mkdir -p $O; F=$O/$(basename $URL)
if [ -s $F ] && [ "$(md5sum < $F | cut -d' ' -f1)" = "$MD5" ]; then echo "already ok: $F"; exit 0; fi
for i in 1 2 3; do
  curl -fsS --retry 5 --retry-delay 30 -C - -o $F "https://$URL" && break || sleep 60
done
GOT=$(md5sum < $F | cut -d' ' -f1)
[ "$GOT" = "$MD5" ] || { echo "MD5 MISMATCH $F: $GOT != $MD5"; exit 1; }
echo "$MD5  $(basename $F)" >> $O/md5.txt
echo "ok: $F $(zcat $F | sed -n 2p | tr -d '\n' | wc -c)bp"
