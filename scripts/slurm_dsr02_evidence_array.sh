#!/usr/bin/env bash
# DSR-02: `viralscan evidence` for one row of calls.tsv (scripts/dsr02_enumerate_calls.py output).
#   python scripts/dsr02_enumerate_calls.py <round> > <round>/evidence/calls.tsv
#   CALLS=<round>/evidence/calls.tsv sbatch --array=0-$(( $(wc -l < $CALLS) - 2 ))%8 scripts/slurm_dsr02_evidence_array.sh
# Output: <round>/evidence/<dataset>__<sample>__<arm>__<virus slug>/ (read by scripts/dsr02_verdicts.py).
#SBATCH -J dsr02_ev
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr02/logs/%x_%A_%a.log
set -euo pipefail
: "${CALLS:?calls.tsv}"
OFFSET=${OFFSET:-0}   # MaxArraySize is 125: rows >=125 go in a second array with OFFSET=125
M=/exports/archive/hg-funcgenom-research/mdmanurung
CODE=${CODE:-$M/vs_pinned/ef0c0c8}
HOST=${HOST_FA:-/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa}
PANEL_FA=${PANEL_FA:-$M/viral_ref_cat42d/build/viral.fa}
export PATH="$M/conda/envs/viralscan_test_full/bin:$PATH" PYTHONPATH="$CODE/src" MPLBACKEND=Agg
IFS=$'\t' read -r DS ARM SAMPLE RUN_DIR VIRUS _ < <(tail -n +2 "$CALLS" | sed -n "$((SLURM_ARRAY_TASK_ID + OFFSET + 1))p")
SLUG=$(printf '%s' "$VIRUS" | tr -cs 'A-Za-z0-9' '_' | sed 's/^_//; s/_$//')
OUT=$(dirname "$CALLS")/${DS}__${SAMPLE}__${ARM}__${SLUG}
[ -s "$OUT/interpretation_flags.tsv" ] && { echo "already done: $OUT"; exit 0; }
echo "code_sha $(git -C "$CODE" rev-parse HEAD) $DS $SAMPLE $ARM virus=$VIRUS"
python -m viralscan.menu evidence --run-dir "$RUN_DIR" -o "$OUT" --virus "$VIRUS" \
  --viral-fasta "$PANEL_FA" --host-fasta "$HOST" --blast --cell-tags --cores "${SLURM_CPUS_PER_TASK:-8}" --verbose
test -s "$OUT/interpretation_flags.tsv"
rm -rf "$OUT/competitive_host_target.fasta" "$OUT"/blast/competitive_db*   # ~GB each; the TSVs/BAMs are kept
echo "DONE $OUT"
