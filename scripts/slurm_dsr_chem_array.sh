#!/usr/bin/env bash
# DSR round 1, step 1: one `viralscan check-chemistry` per manifest row. -x is NOT passed, so the
# call is unbiased. The JSON (chem/<sample>.json) is the only source of -x / --strand for the
# rerun arrays; a nonzero exit (unresolved or conflicting) is recorded, not hidden.
#
#   sbatch --array=0-18 scripts/slurm_dsr_chem_array.sh
# ROUND_DIR / CODE_DIR / REF_BUILD default to round 1. Round 2 reuses round 1's frozen chem JSONs
# and does not run this step; if it is ever rerun, pass `-o $ROUND_DIR/logs/%x_%A_%a.log`.
#
#SBATCH -J dsr_chem
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round1/logs/%x_%A_%a.log
set -uo pipefail
# >>> dsr-env (tests/test_dsr_round_env.py evaluates this block alone)
M=/exports/archive/hg-funcgenom-research/mdmanurung
R=${ROUND_DIR:-$M/viralscan_work/dsr_round1}
CODE=${CODE_DIR:-$M/vs_pinned/bbf1821}
REF=${REF_BUILD:-$M/viral_ref_cat42d/build}
# <<< dsr-env
export PATH="$M/conda/envs/viralscan_test_full/bin:$PATH" PYTHONPATH="$CODE/src" MPLBACKEND=Agg
[ -z "$(git -C "$CODE" status --porcelain)" ] || { echo "dirty pinned tree" >&2; exit 2; }
IFS=$'\t' read -r DS SAMPLE GRP R1 R2 WL PAIRS HINT < <(tail -n +2 "$R/manifest.tsv" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
[ "$WL" != - ] || WL=""   # manifest uses - for an empty field (tab IFS would collapse it)
echo "code $(git -C "$CODE" describe --always --dirty)  $DS/$SAMPLE"
ARGS=(-s1 "$R1" -s2 "$R2" -i "$REF/panel.idx" -t "$REF/panel.t2g" -c 4 --json "$R/chem/$SAMPLE.json")
[ -n "$WL" ] && ARGS+=(-w "$WL")
python -m viralscan.menu check-chemistry "${ARGS[@]}" 2>&1 | tee "$R/chem/$SAMPLE.txt"
echo "exit=${PIPESTATUS[0]}" | tee -a "$R/chem/$SAMPLE.txt"
