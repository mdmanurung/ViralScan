#!/usr/bin/env bash
# Run the multimap profiling analysis.
# Usage: bash analysis/multimap_profiling/scripts/run.sh
# Must be run from the repo root: /exports/archive/hg-funcgenom-research/mdmanurung/ViralScan

set -euo pipefail

PYTHON=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench/bin/python
REPO=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan

echo "[run.sh] Working directory: $(pwd)"
echo "[run.sh] Python: $PYTHON"
echo "[run.sh] Python version: $($PYTHON --version)"

PYTHONPATH="$REPO/src" \
    "$PYTHON" \
    analysis/multimap_profiling/scripts/profile_multimap.py \
    "$PYTHON" \
    2>&1 | tee analysis/multimap_profiling/outputs/profile_run.log

echo "[run.sh] done."
