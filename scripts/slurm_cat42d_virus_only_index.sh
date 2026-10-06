#!/usr/bin/env bash
# Virus-only cat42d index for the STAR two-step arm (host reads removed by STARsolo
# first). Same viral.fa (homopolymer + low-complexity masked), same GTF (DSR-01:
# combined.gtf minus host_cdna rows, == cat42 viral_panel.gtf) and same GRCh38 genome
# D-list as panel.idx; only the host cDNA is left out. Same kb/kallisto 0.50.1 as the
# cat42d panel build.
#
#   sbatch scripts/slurm_cat42d_virus_only_index.sh
#
#SBATCH -J cat42d_virus_idx
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viral_ref_cat42d/logs/%x_%j.log
set -euo pipefail
C=/exports/archive/hg-funcgenom-research/mdmanurung/viral_ref_cat42d
B=$C/build
O=$C/build_virus_only
export PATH=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_env_pinned/bin:$PATH
mkdir -p "$O"
kb ref -i "$O/virus_panel.idx" -g "$O/virus_panel.t2g" -f1 "$O/virus_cdna.fa" --overwrite \
  -t 8 --d-list "$C/dlist.fa" "$B/viral.fa" "$B/viral_panel.gtf"
test -s "$O/virus_panel.idx" && test -s "$O/virus_panel.t2g"
sha256sum "$B/viral.fa" "$B/viral_panel.gtf" "$C/dlist.fa" "$O/virus_panel.t2g" | tee "$O/inputs.sha256"
echo BUILD_OK
