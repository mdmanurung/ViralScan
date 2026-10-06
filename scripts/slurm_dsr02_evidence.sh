#!/usr/bin/env bash
# DSR-02: read-level check of one (run, virus) call with the native `viralscan evidence`
# (exact CB+UMI read replay, competitive alignment vs GRCh38 + panel, alignment QC, BLAST,
# interpretation flags). Replaces the hand-edited run.sh copies used for F-025 / F-026.
#
#   RUN_DIR=<sample run dir> VIRUS=<exact label or gene id> PANEL_FA=<viral fasta of that run's index> \
#   OUT=<evidence dir> sbatch scripts/slurm_dsr02_evidence.sh
#
#SBATCH -J dsr02_evidence
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr02/logs/%x_%j.log
set -euo pipefail
: "${RUN_DIR:?}" "${VIRUS:?}" "${PANEL_FA:?}" "${OUT:?}"
CODE=${CODE:-/exports/archive/hg-funcgenom-research/mdmanurung/vs_pinned/4346dc8}
ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_test_full
HOST=${HOST_FA:-/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa}
export PATH="$ENV/bin:$PATH" PYTHONPATH="$CODE/src" MPLBACKEND=Agg
echo "code_sha $(git -C "$CODE" rev-parse HEAD) virus=$VIRUS run=$RUN_DIR"
python -m viralscan.menu evidence --run-dir "$RUN_DIR" -o "$OUT" --virus "$VIRUS" \
  --viral-fasta "$PANEL_FA" --host-fasta "$HOST" --blast --cell-tags --cores 8 --verbose
ls "$OUT"
