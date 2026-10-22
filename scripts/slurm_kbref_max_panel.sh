#!/usr/bin/env bash
# kallisto indices for the max viral panel (viral_panel_max_2026-09-28/):
# baseline v2 (2,215) + every deduplicated, host-filtered, hard-masked addition
# from 01_collect.py .. 04_assemble.py in that directory. Same two builds and
# same implicit D-list as scripts/slurm_kbref_full_panel_v2.sh, so results are
# comparable with the v2 panel.
#
# Submit:
#   sbatch scripts/slurm_kbref_max_panel.sh
#
#SBATCH -J kbref_max_panel
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH -o logs/kbref_max_panel_%j.log
#SBATCH -e logs/kbref_max_panel_%j.err

SITE_ROOT=${VIRALSCAN_SITE_ROOT:?set VIRALSCAN_SITE_ROOT to the site work root}
SITE_SCRATCH=${VIRALSCAN_SITE_SCRATCH:?set VIRALSCAN_SITE_SCRATCH to the site scratch root}

set -euo pipefail

REPO=$SITE_ROOT/ViralScan
ENV=$SITE_ROOT/conda/envs/viralscan_bench
REF=$SITE_SCRATCH/ViralScan/viral_panel_max_2026-09-28

export PATH="$ENV/bin:$PATH"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

echo "=== HEAD at build time ==="
git -C "$REPO" rev-parse HEAD
git -C "$REPO" status --short
kb --version 2>&1 | head -1 || true

echo "=== Combined host+virus index (max panel) ==="
kb ref -i "$REF/index_max.idx" -g "$REF/t2g_max.txt" -f1 "$REF/cdna_max.fa" \
  --tmp "$REF/tmp_kbref_combined_max" \
  "$REF/combined_max.fa" "$REF/combined_max.gtf"

echo "=== Virus-only index (max panel) ==="
kb ref -i "$REF/virus_index_max.idx" -g "$REF/virus_t2g_max.txt" -f1 "$REF/virus_cdna_max.fa" \
  --tmp "$REF/tmp_kbref_virusonly_max" \
  "$REF/viral_max.fa" "$REF/viral_max.gtf"

echo "=== Spot checks ==="
t2g_lines=$(wc -l < "$REF/t2g_max.txt")
host_lines=$(grep -c "^ENST" "$REF/t2g_max.txt" || true)
virus_only=$(wc -l < "$REF/virus_t2g_max.txt")
echo "combined t2g_max: $t2g_lines total, $host_lines host ENST* (v2: 465769), viral $((t2g_lines - host_lines)); virus-only t2g: $virus_only"
for acc in AF157706 NC_000898 AB303562 NC_009334 NC_001802 NC_045512; do
  printf '%-10s %s\n' "$acc" "$(grep -c "$acc" "$REF/virus_t2g_max.txt" || true)"
done

echo "Done."
