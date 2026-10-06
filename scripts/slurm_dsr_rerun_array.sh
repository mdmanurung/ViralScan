#!/usr/bin/env bash
# DSR round 1, step 2: one ViralScan run per (manifest row, arm). Arms differ ONLY by flags:
#   combined_off      -i panel.idx                      (no filter flag)
#   combined_artefact -i panel.idx  --read-filter artefact
#   twostep           -i virus_panel.idx  --host-filter starsolo --host-index <STAR GRCh38-2024-A>
# Everything else (-x, --strand, -w, --cell-calling, reference, code) is the same for every arm
# and comes from the frozen manifest + chem/<sample>.json. Resources are set on the sbatch line.
#
#   ARM=combined_off sbatch --array=<rows> --cpus-per-task=4 --mem=32G --time=12:00:00 \
#       -J dsr_off_S scripts/slurm_dsr_rerun_array.sh
#
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round1/logs/%x_%A_%a.log
set -euo pipefail
: "${ARM:?combined_off|combined_artefact|twostep}"
M=/exports/archive/hg-funcgenom-research/mdmanurung
R=$M/viralscan_work/dsr_round1
CODE=$M/vs_pinned/bbf1821
B=$M/viral_ref_cat42d/build
STAR_IDX=$M/ViralScan/ebv_latest_ref_2026-09-27/star_host_2.7.11b
export PATH="$M/conda/envs/viralscan_bench/bin:$M/conda/envs/R4_51/bin:$PATH" PYTHONPATH="$CODE/src" MPLBACKEND=Agg
SHA=$(git -C "$CODE" describe --always --dirty)
case "$SHA" in *-dirty) echo "dirty pinned tree: $SHA" >&2; exit 2;; esac
IFS=$'\t' read -r DS SAMPLE GRP R1 R2 WL PAIRS HINT < <(tail -n +2 "$R/manifest.tsv" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
[ "$WL" != - ] || WL=""   # manifest uses - for an empty field (tab IFS would collapse it)
CHEM=$R/chem/$SAMPLE.json
[ -s "$CHEM" ] || { echo "no chemistry JSON for $SAMPLE" >&2; exit 3; }
read -r TECH STRAND < <(python - "$CHEM" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
print(d["chemistry"]["chemistry"] or "NONE", (d.get("strand_pilot") or {}).get("choice") or "NONE")
PY
)
[ "$TECH" != NONE ] || { echo "chemistry unresolved for $SAMPLE" >&2; exit 3; }
case "$ARM" in
  combined_off)      IDX=(-i "$B/panel.idx" -t "$B/panel.t2g"); EXTRA=();;
  combined_artefact) IDX=(-i "$B/panel.idx" -t "$B/panel.t2g"); EXTRA=(--read-filter artefact);;
  twostep)           IDX=(-i "$M/viral_ref_cat42d/build_virus_only/virus_panel.idx" -t "$M/viral_ref_cat42d/build_virus_only/virus_panel.t2g")
                     EXTRA=(--host-filter starsolo --host-index "$STAR_IDX");;
  *) echo "unknown ARM $ARM" >&2; exit 2;;
esac
OUT=$R/runs/$DS/$ARM/$SAMPLE
if [ -s "$OUT/$SAMPLE/results/viral_summary.tsv" ]; then echo "already done: $OUT"; exit 0; fi
# run metadata goes beside $OUT: viralscan refuses a non-empty -o directory
META=$OUT.meta
rm -rf "$OUT" "$META"; mkdir -p "$META"
echo "$SHA" > "$META/code_sha.txt"; cp "$CHEM" "$META/check_chemistry.json"
ARGS=(-o "$OUT" -s1 "$R1" -s2 "$R2" "${IDX[@]}" -gtf "$B/viral_panel.gtf" -x "$TECH" -c "${SLURM_CPUS_PER_TASK:-4}"
      --cell-calling emptydrops --anellovirus-gene-ids --yes --verbose)
[ "$STRAND" = NONE ] || ARGS+=(--strand "$STRAND")
[ -z "$WL" ] || ARGS+=(-w "$WL")
[ -z "${FORCE_TECH:-}" ] || ARGS+=(--force-technology)   # GSM4682311 only: documented 10xv2, check-chemistry unresolved (user, 2026-10-06)
echo "code $SHA  $DS/$SAMPLE  arm=$ARM  -x $TECH --strand $STRAND  ${WL:+-w $WL}"
/usr/bin/time -v -o "$META/time.txt" python -m viralscan.menu "${ARGS[@]}" "${EXTRA[@]}"
test -s "$OUT/$SAMPLE/results/viral_summary.tsv"
echo "DONE $DS/$SAMPLE $ARM"
