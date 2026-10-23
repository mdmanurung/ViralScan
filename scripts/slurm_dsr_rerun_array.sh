#!/usr/bin/env bash
# DSR rounds 1-2, step 2: one ViralScan run per (manifest row, arm). Arms differ ONLY by flags:
#   combined_off      -i panel.idx                      (no filter flag)
#   combined_artefact -i panel.idx  --read-filter artefact
#   twostep           -i virus_panel.idx  --host-filter starsolo --host-index <STAR GRCh38-2024-A>
# DEF-05 two-step cell calling needs no extra flag: with --host-filter starsolo and --cell-calling
# emptydrops, detection.py calls cells on the STARsolo host matrix (solo_raw_dir -> host_called_cells).
# Round 2 therefore needs no twostep_v2 redetect pass.
# Everything else (-x, --strand, -w, --cell-calling, reference, code) is the same for every arm
# and comes from the frozen manifest + chem/<sample>.json. Resources are set on the sbatch line.
# ROUND_DIR / CODE_DIR / REF_BUILD default to round 1 (dsr_round1, vs_pinned/bbf1821, cat42d build).
# Round 2: export ROUND_DIR CODE_DIR first and pass `-o $ROUND_DIR/logs/%x_%A_%a.log` to sbatch
# (the #SBATCH -o below is round 1's and cannot read env). Outside round 1 every finished run is
# appended to $ROUND_DIR/explore_ledger.tsv (scripts/explore_ledger.py); EXPLORE_LEDGER overrides it.
#
#   ARM=combined_off sbatch --array=<rows> --cpus-per-task=4 --mem=32G --time=12:00:00 \
#       -J dsr_off_S scripts/slurm_dsr_rerun_array.sh
#
#SBATCH -o slurm_%x_%j.log
set -euo pipefail
: "${ARM:?combined_off|combined_artefact|twostep}"
# >>> dsr-env (tests/test_dsr_round_env.py evaluates this block alone)
M=${VIRALSCAN_SITE_ROOT:?set VIRALSCAN_SITE_ROOT to the site work root}
ROUND1=$M/viralscan_work/dsr_round1
R=${ROUND_DIR:-$ROUND1}
CODE=${CODE_DIR:-$M/vs_pinned/bbf1821}
B=${REF_BUILD:-$M/viral_ref_cat42d/build}
VIRUS_ONLY=${VIRUS_ONLY_BUILD:-$(dirname "$B")/build_virus_only}
STAR_IDX=${STAR_IDX:-$M/ViralScan/ebv_latest_ref_2026-09-27/star_host_2.7.11b}
# <<< dsr-env
REF_ID=${REF_ID:-$(basename "$(dirname "$B")")}
LEDGER=${EXPLORE_LEDGER-}
[ -n "$LEDGER" ] || [ "$R" = "$ROUND1" ] || LEDGER=$R/explore_ledger.tsv
LEDGER_PY=${EXPLORE_LEDGER_PY:-$CODE/scripts/explore_ledger.py}
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
  twostep)           IDX=(-i "$VIRUS_ONLY/virus_panel.idx" -t "$VIRUS_ONLY/virus_panel.t2g")
                     EXTRA=(--host-filter starsolo --host-index "$STAR_IDX");;
  *) echo "unknown ARM $ARM" >&2; exit 2;;
esac
# Track B guard: never run a deny-listed accession (checked against the manifest row, not just the ledger)
if [ -n "$LEDGER" ]; then
  DENY=$(dirname "$LEDGER_PY")/../config/track_b_denylist.txt
  [ -f "$DENY" ] || { echo "no deny-list at $DENY" >&2; exit 2; }
  if printf '%s\n' "$DS $SAMPLE $R1 $R2" | grep -iqFf <(sed 's/#.*//;/^[[:space:]]*$/d' "$DENY"); then
    echo "deny-listed accession in manifest row: $DS/$SAMPLE" >&2; exit 4
  fi
fi
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
[ -z "$LEDGER" ] || python "$LEDGER_PY" --ledger "$LEDGER" append --dataset "$DS" --sample "$SAMPLE" \
  --arm "$ARM" --code-sha "$SHA" --reference-id "$REF_ID"
echo "DONE $DS/$SAMPLE $ARM"
