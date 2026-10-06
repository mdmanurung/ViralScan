#!/usr/bin/env bash
# Grid-cost check vs Q10d: time one `rerun-multimap em-cell` on the COST-01
# full-depth combined run of SRR12682296 (10xv2, 127,045,580 pairs; F-024).
# Same 8 cores as the COST-01 combined jobs so allocated core-h compare 1:1.
# 64G OOM-killed multimap after ~20 min (job 25702291, MaxRSS 63.6 GB): em-cell
# keeps a dense n_genes theta per barcode (46,238 x 8 B x up to 793,308 ≈ 293 GB).
#
#   sbatch scripts/slurm_rerun_multimap_em_cell_cost.sh
#
#SBATCH -J vs_rerun_emcell
#SBATCH --cpus-per-task=8
#SBATCH --mem=480G
#SBATCH --time=16:00:00
#SBATCH -o logs/vs_rerun_emcell_%j.log
#SBATCH -e logs/vs_rerun_emcell_%j.err

set -euo pipefail

REPO=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_test_full
SRC=/exports/para-lipg-hpc/mdmanurung/ViralScan/ebv_latest_ref_2026-09-27/runs/combined/SRR12682296
OUT=/exports/para-lipg-hpc/mdmanurung/ViralScan/cost01/rerun_em_cell/SRR12682296

export PATH="$ENV/bin:$PATH"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

git -C "$REPO" rev-parse HEAD
git -C "$REPO" status --short --untracked-files=no
mkdir -p "$(dirname "$OUT")"

date +%s > "$(dirname "$OUT")/SRR12682296.start"
python -c "from viralscan.menu import main; main()" rerun-multimap \
  --run-dir "$SRC" -o "$OUT" \
  --multimap-method em-cell \
  -c 8 \
  --verbose
date +%s > "$(dirname "$OUT")/SRR12682296.end"

# kb_count must not have rerun: the copied BUS file keeps the source mtime.
stat -c '%y %n' "$SRC/SRR12682296/kb-python/output.bus" "$OUT/SRR12682296/kb-python/output.bus"
ls -la --time-style=full-iso "$OUT/SRR12682296/log/"
