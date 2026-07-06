#!/usr/bin/env bash
# Publication checklist — submit all PR-22 validation jobs in one command.
#
# Chains the three SLURM submissions needed to close the empirical gap
# for publication (P22.4, P22.5, P22.6):
#
#   Step 1  sbatch --array=0-2  slurm_full_depth_validation.sh
#           Three datasets: HHV-6/CAR-T (SRR20710641), EBV/LCL (SRR12682296),
#           HSV-1/fibroblast (SRR8315713).
#
#   Step 2  afterok dependency — emit BENCHMARK_COMPARISON_full_depth.tsv once
#           all three array tasks succeed.
#
#   Step 3  slurm_starsolo_ebv_comparison.sh — independent of the array; runs
#           in parallel.  Produces comparison_starsolo_vs_viralscan.tsv.
#
# Usage:
#   bash scripts/publication_checklist.sh          # submit everything
#   bash scripts/publication_checklist.sh --dry-run  # print sbatch commands only
#
# After all jobs finish, transcribe the two TSVs into BENCHMARK_COMPARISON.md
# and docs/manuscript_draft.md (Tables 3.2-3.3 + §3.4), then flip P22.4-P22.7
# checkboxes in PLAN.md in the same commit (CLAUDE.md contract).

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS="$REPO/scripts"
DRY_RUN=false

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
        --help|-h)
            sed -n '2,/^set -euo/p' "${BASH_SOURCE[0]}" | grep '^#' | sed 's/^# \?//'
            exit 0
            ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# ── preflight ─────────────────────────────────────────────────────────────────
if ! command -v sbatch &>/dev/null; then
    echo "ERROR: sbatch not on PATH — are you on the cluster login node?" >&2
    exit 1
fi

echo "=== ViralScan publication checklist (PR 22) ==="
echo "Repo: $REPO"
echo ""

# ── Step 1: full-depth SLURM array (P22.4) ───────────────────────────────────
FD_SCRIPT="$SCRIPTS/slurm_full_depth_validation.sh"
echo "Step 1: submit full-depth validation array (3 datasets)"
echo "  sbatch --parsable --array=0-2 $FD_SCRIPT"

if [[ $DRY_RUN == false ]]; then
    FD_JOB=$(sbatch --parsable --array=0-2 "$FD_SCRIPT")
    echo "  → array job ID: $FD_JOB"
else
    FD_JOB="DRY_RUN_FD"
    echo "  [dry-run] would submit; using placeholder job ID $FD_JOB"
fi
echo ""

# ── Step 2: emit benchmark TSV after array succeeds (P22.4 summary) ──────────
# Uses afterok so this runs only if ALL three array tasks exit 0.
# The --summarize flag now has an early dispatch at the top of the script
# (no kb/snakemake required; safe on login node or minimal compute node).
SUMMARIZE_CMD="bash $FD_SCRIPT --summarize"
echo "Step 2: emit BENCHMARK_COMPARISON_full_depth.tsv (afterok:$FD_JOB)"
echo "  sbatch --dependency=afterok:$FD_JOB --job-name=vs_summarize"
echo "         --cpus-per-task=1 --mem=1G --time=00:10:00"
echo "         --wrap \"$SUMMARIZE_CMD\""

if [[ $DRY_RUN == false ]]; then
    SUM_JOB=$(sbatch \
        --parsable \
        --dependency="afterok:$FD_JOB" \
        --job-name=vs_summarize \
        --cpus-per-task=1 \
        --mem=1G \
        --time=00:10:00 \
        --wrap "$SUMMARIZE_CMD")
    echo "  → summarize job ID: $SUM_JOB"
else
    SUM_JOB="DRY_RUN_SUM"
    echo "  [dry-run] would submit; using placeholder job ID $SUM_JOB"
fi
echo ""

# ── Step 3: STARsolo EBV comparison (P22.6) — independent ───────────────────
STAR_SCRIPT="$SCRIPTS/slurm_starsolo_ebv_comparison.sh"
echo "Step 3: submit STARsolo EBV comparison (independent of array)"
echo "  sbatch $STAR_SCRIPT"

if [[ $DRY_RUN == false ]]; then
    STAR_JOB=$(sbatch --parsable "$STAR_SCRIPT")
    echo "  → STARsolo job ID: $STAR_JOB"
else
    STAR_JOB="DRY_RUN_STAR"
    echo "  [dry-run] would submit; using placeholder job ID $STAR_JOB"
fi
echo ""

# ── summary ───────────────────────────────────────────────────────────────────
VS_ROOT=/exports/para-lipg-hpc/mdmanurung/viralscan_showcase
STAR_WORKDIR="$REPO/starsolo_p22_6"

echo "=== Jobs submitted ==="
echo "  Full-depth array  : $FD_JOB   (squeue -j $FD_JOB)"
echo "  Summarize (afterok): $SUM_JOB"
echo "  STARsolo           : $STAR_JOB  (squeue -j $STAR_JOB)"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  tail -f logs/vs_val_${FD_JOB}_*.out  # full-depth array logs"
echo "  tail -f logs/ebv_starsolo_${STAR_JOB}.log  # STARsolo log"
echo ""
echo "When all jobs complete, review:"
echo "  $VS_ROOT/BENCHMARK_COMPARISON_full_depth.tsv  (written by job $SUM_JOB)"
echo "  $STAR_WORKDIR/comparison_starsolo_vs_viralscan.tsv"
echo ""
echo "Then transcribe both TSVs into:"
echo "  BENCHMARK_COMPARISON.md      (update per-study blocks + Summary Table)"
echo "  docs/manuscript_draft.md     (Tables 3.2-3.3 + §3.4)"
echo "  PLAN.md                      (flip P22.4/P22.5/P22.6/P22.7 → [x]; same commit)"
