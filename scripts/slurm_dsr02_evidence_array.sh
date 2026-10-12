#!/usr/bin/env bash
# DSR-02: read-validate one row of calls.tsv (scripts/dsr02_enumerate_calls.py output) with
# scripts/dsr02_validate_call.py (`viralscan evidence` + per-read verdicts + repeat class).
#   python scripts/dsr02_enumerate_calls.py <round> > <round>/evidence/calls.tsv
#   CALLS=<round>/evidence/calls.tsv sbatch --array=0-$(( $(wc -l < $CALLS) - 2 ))%8 scripts/slurm_dsr02_evidence_array.sh
# Output: <round>/evidence/<dataset>__<sample>__<arm>__<virus slug>/ (read by scripts/dsr02_verdicts.py).
# Dataset-specific panels: PANEL_FOR="hpv16=/path/a.fa ebv=/path/b.fa" (default PANEL_FA for the rest).
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
PANEL_ARGS=()
for kv in ${PANEL_FOR:-}; do PANEL_ARGS+=(--panel-for "$kv"); done
echo "code_sha $(git -C "$CODE" rev-parse HEAD) row=$((SLURM_ARRAY_TASK_ID + OFFSET))"
python "$CODE/scripts/dsr02_validate_call.py" --calls "$CALLS" --row $((SLURM_ARRAY_TASK_ID + OFFSET)) \
  --panel-fasta "$PANEL_FA" ${PANEL_ARGS[@]+"${PANEL_ARGS[@]}"} --host-fasta "$HOST" --cores "${SLURM_CPUS_PER_TASK:-8}"
