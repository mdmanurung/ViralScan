#!/usr/bin/env bash
#SBATCH --job-name=vs3-ebv-mm
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=analysis/v3_ebv_baseline/slurm-%j.out
#SBATCH --error=analysis/v3_ebv_baseline/slurm-%j.err
set -euo pipefail

REPO="${VIRALSCAN_REPO:-$PWD}"
SAMPLE_DIR="${VIRALSCAN_EBV_SAMPLE_DIR:-$REPO/benchmark_runs/reference_strategy_2026-06-28_fresh12b/runs/ebv__viralscan__combined/SRR12682296}"
T2G="${VIRALSCAN_EBV_T2G:-}"
OUT_DIR="$REPO/analysis/v3_ebv_baseline"

if [[ -z "$T2G" ]]; then
  echo "Set VIRALSCAN_EBV_T2G to the checksum-pinned transcript-to-gene map." >&2
  exit 2
fi

mkdir -p "$OUT_DIR"
cd "$REPO"
export PYTHONPATH="$REPO/src"
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/viralscan-numba-cache"

/usr/bin/time -v python3 scripts/benchmark_v3_multimap.py \
  --sample-dir "$SAMPLE_DIR" \
  --t2g "$T2G" \
  --method host-conservative \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --output "$OUT_DIR/host_conservative.json"
