#!/usr/bin/env bash
# Corrected v2 build: reverts the mistaken HHV-6B swap. AF157706.1 ("Human
# herpesvirus 6B strain Z29, complete genome" -- confirmed via the source
# Nature paper s41586-023-06704-2, which uses this exact accession as its
# HHV-6B reference) is restored; NC_000898.1 (mistakenly substituted in,
# based on a misreading of PLAN.md's CAT-28 note, which was actually about a
# different accounting issue in a virus-catalog completeness audit, not a
# defect in AF157706.1 itself) is removed. Written to new v2 filenames so the
# currently-running EBV/HSV-1 two-step jobs (which read index.idx/
# virus_index.idx built from the v1 panel) are not disturbed mid-run --
# those results are unaffected by this HHV-6B-only correction anyway.
#
# Submit:
#   sbatch scripts/slurm_kbref_full_panel_v2.sh
#
#SBATCH -J kbref_full_panel_v2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH -o logs/kbref_full_panel_v2_%j.log
#SBATCH -e logs/kbref_full_panel_v2_%j.err

set -euo pipefail

REPO=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
WORK=/exports/para-lipg-hpc/mdmanurung/ViralScan/ebv_latest_ref_2026-09-27
REF=$WORK/full_panel

export PATH="$ENV/bin:$PATH"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

echo "=== HEAD at build time ==="
git -C "$REPO" rev-parse HEAD
git -C "$REPO" status --short

echo "=== Combined host+virus index v2 (AF157706.1 restored) ==="
kb ref -i "$REF/index_v2.idx" -g "$REF/t2g_v2.txt" -f1 "$REF/cdna_v2.fa" \
  --tmp "$WORK/tmp_kbref_combined_v2" \
  "$REF/combined_v2.fa" "$REF/combined_v2.gtf"

echo "=== Virus-only index v2 ==="
kb ref -i "$REF/virus_index_v2.idx" -g "$REF/virus_t2g_v2.txt" -f1 "$REF/virus_cdna_v2.fa" \
  --tmp "$WORK/tmp_kbref_virusonly_v2" \
  "$REF/viral_final_v2.fa" "$REF/viral_final_v2.gtf"

echo "=== Spot checks ==="
t2g_lines=$(wc -l < "$REF/t2g_v2.txt")
host_lines=$(grep -c "^ENST" "$REF/t2g_v2.txt" || true)
echo "combined t2g_v2: $t2g_lines total, $host_lines host ENST*"
af157706=$(grep -c "AF157706" "$REF/t2g_v2.txt" || true)
nc000898=$(grep -c "NC_000898" "$REF/t2g_v2.txt" || true)
echo "AF157706 (HHV-6B, restored) entries: $af157706   NC_000898 (should be 0 now): $nc000898"

echo "Done."
