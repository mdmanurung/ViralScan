#!/usr/bin/env bash
# One-time kb ref build for the 2,215-genome full-panel reference (194 bundled
# + 2,021 clareaulab anellovirus, CAT-05 duplicate dropped, CAT-17 masked,
# HHV-6B pseudocontig AF157706.1 swapped for the real NC_000898.1 genome).
#
# Builds BOTH indices used downstream:
#   - full_panel/index.idx + t2g.txt        (combined host+virus, strategy 1)
#   - full_panel/virus_index.idx + virus_t2g.txt  (virus-only, strategy 2)
#
# The quant array jobs (slurm_quant_combined_array.sh,
# slurm_quant_twostep_array.sh) depend on this job via --dependency=afterok.
#
# Submit:
#   sbatch scripts/slurm_kbref_full_panel.sh
#
#SBATCH -J kbref_full_panel
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH -o logs/kbref_full_panel_%j.log
#SBATCH -e logs/kbref_full_panel_%j.err

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

echo "=== Combined host+virus index (2,215 viral + ~465k host transcripts) ==="
kb ref -i "$REF/index.idx" -g "$REF/t2g.txt" -f1 "$REF/cdna.fa" \
  --tmp "$WORK/tmp_kbref_combined_full" \
  "$REF/combined.fa" "$REF/combined.gtf"

echo "=== Virus-only index (2,215 viral transcripts) ==="
kb ref -i "$REF/virus_index.idx" -g "$REF/virus_t2g.txt" -f1 "$REF/virus_cdna.fa" \
  --tmp "$WORK/tmp_kbref_virusonly_full" \
  "$REF/viral_final.fa" "$REF/viral_final.gtf"

echo "=== Spot checks ==="
t2g_lines=$(wc -l < "$REF/t2g.txt")
host_lines=$(grep -c "^ENST" "$REF/t2g.txt" || true)
echo "combined t2g: $t2g_lines total, $host_lines host ENST*"
epstein=$(grep -c "EPSTEIN" "$REF/t2g.txt" || true)
hhv6b=$(grep -c "NC_000898" "$REF/t2g.txt" || true)
echo "EPSTEIN (EBV) entries: $epstein   NC_000898 (HHV-6B) entries: $hhv6b"

echo "Done."
