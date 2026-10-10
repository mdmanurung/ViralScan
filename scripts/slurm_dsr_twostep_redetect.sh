#!/usr/bin/env bash
# DSR-15: rerun detection only (no STAR, no kb) on a finished twostep run, with cells called on
# the STARsolo host matrix. Writes runs/<dataset>/twostep_v2/<sample>; the source run is never touched.
# Detection-only needs `rerun-triggers: mtime` (profile in $R/redetect_cwd/<ds>-<sample>): the source ran on an older
# commit, so snakemake's default "code changed" trigger reran STARsolo + kb and the kb `mv` collided
# with the copied kb-python/counts_unfiltered (25715231, 25715437).
#   sbatch --array=<manifest rows> -J dsr_two_redet --cpus-per-task=4 --mem=16G --time=4:00:00 \
#       scripts/slurm_dsr_twostep_redetect.sh
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round1/logs/%x_%A_%a.log
set -euo pipefail
M=/exports/archive/hg-funcgenom-research/mdmanurung
R=$M/viralscan_work/dsr_round1
CODE=${CODE:?pinned worktree with the DSR-15 patch}
export PATH="$M/conda/envs/viralscan_bench/bin:$M/conda/envs/R4_51/bin:$PATH" PYTHONPATH="$CODE/src" MPLBACKEND=Agg
SHA=$(git -C "$CODE" describe --always --dirty)
case "$SHA" in *-dirty) echo "dirty pinned tree: $SHA" >&2; exit 2;; esac
IFS=$'\t' read -r DS SAMPLE _ < <(tail -n +2 "$R/manifest.tsv" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
SRC=$R/runs/$DS/twostep/$SAMPLE
OUT=$R/runs/$DS/twostep_v2/$SAMPLE
# the source may have failed in detection (emptyDrops n/a); rerun-multimap needs only multimap.done
[ -f "$SRC/$SAMPLE/log/multimap.done" ] || { echo "no multimap checkpoint: $SRC" >&2; exit 3; }
# the copy inherits the source's results, so only the run marker proves a v2 finished
if [ -s "$OUT/run_complete.json" ]; then echo "already done: $OUT"; exit 0; fi
rm -rf "$OUT" "$OUT.meta"; mkdir -p "$OUT.meta" "$(dirname "$OUT")"
mkdir -p "$R/redetect_cwd/$DS-$SAMPLE/profiles/default"
printf 'rerun-triggers:\n  - mtime\n' > "$R/redetect_cwd/$DS-$SAMPLE/profiles/default/config.v9+.yaml"
cd "$R/redetect_cwd/$DS-$SAMPLE"
echo "$SHA" > "$OUT.meta/code_sha.txt"
python -m viralscan.menu rerun-multimap --run-dir "$SRC" -o "$OUT" --multimap-method host-conservative --cores "${SLURM_CPUS_PER_TASK:-4}" --verbose
test -s "$OUT/$SAMPLE/results/viral_summary.tsv"
echo "DONE $DS/$SAMPLE twostep_v2 $SHA"
