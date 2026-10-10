#!/usr/bin/env bash
# DLIST-01 / I1: cat42d virus-only index with a genome + cDNA D-list (the Pachter macaque notebook
# 7_virus_host_captured_dlist_cdna_dna recipe). Same viral.fa, GTF, kb/kallisto pin and flags as
# slurm_cat42d_virus_only_index.sh; only the D-list differs: cat42d/dlist.fa (GRCh38 genome) + the
# GRCh38 cDNA, which adds spliced exon-junction k-mers the genome lacks. New directory; the frozen
# cat42d reference is untouched.
#
#   sbatch scripts/slurm_cat42d_virus_only_dlist2_index.sh
#
#SBATCH -J cat42d_dlist2_idx
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viral_ref_cat42d_dlist2/logs/%x_%j.log
set -euo pipefail
M=/exports/archive/hg-funcgenom-research/mdmanurung
C=$M/viral_ref_cat42d
B=$C/build
O=$M/viral_ref_cat42d_dlist2/build_virus_only
export PATH=$M/viralscan_env_pinned/bin:$PATH
mkdir -p "$O"
DL=$M/viral_ref_cat42d_dlist2/dlist_cdna_dna.fa
if [ ! -s "$DL" ]; then
  { cat "$C/dlist.fa"; zcat "$B/host/Homo_sapiens.GRCh38.cdna.all.fa.gz"; } > "$DL.tmp" && mv "$DL.tmp" "$DL"
fi
sha256sum "$DL" > "$M/viral_ref_cat42d_dlist2/dlist_cdna_dna.sha256"
kb ref -i "$O/virus_panel.idx" -g "$O/virus_panel.t2g" -f1 "$O/virus_cdna.fa" --overwrite \
  -t 8 --d-list "$DL" "$B/viral.fa" "$B/viral_panel.gtf"
test -s "$O/virus_panel.idx" && test -s "$O/virus_panel.t2g"
sha256sum "$B/viral.fa" "$B/viral_panel.gtf" "$DL" "$O/virus_panel.t2g" | tee "$O/inputs.sha256"
echo BUILD_OK
