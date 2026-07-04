#!/usr/bin/env bash
# run.sh — execute harmonize_2x2.py with the viralscan_bench conda env
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../../.."; pwd)"
PYTHON=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench/bin/python

cd "$REPO"
PYTHONPATH=src "$PYTHON" analysis/reference_strategy_benchmark/scripts/harmonize_2x2.py "$@" 2>&1 | tee analysis/reference_strategy_benchmark/outputs/harmonize_2x2.log
