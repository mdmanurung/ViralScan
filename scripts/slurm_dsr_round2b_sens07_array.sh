#!/usr/bin/env bash
# DSR-16 round 2b / SENS-07: per-dataset planted-read sensitivity. One array task per round-2 manifest
# row (array index i = data row i of dsr_round2/manifest.tsv, 0-based). Per row:
#   1. sens07 `plant`   background FASTQ pair -> planted/ and unplanted/ twin (exact genomes only)
#   2. viralscan        planted/ and unplanted/, with round 2's `combined_off` settings (no filter flag)
#   3. sens07 `score`   planted run minus twin -> recovery_by_condition.tsv, family_summary.tsv
# Every viralscan run is appended to $ROUND_DIR/explore_ledger.tsv (label `explore`). Exploratory,
# unscored, one library per dataset (n = 1); not VAL-01.
#
# Reuses round 2's manifest.tsv and chem/<sample>.json read-only; writes only under ROUND_DIR
# (runs/ status/ logs/ explore_ledger.tsv). Per-row choices live in a config TSV (below).
#
# Setup after the user approves (the user or the assistant, once):
#   R=$VIRALSCAN_SITE_ROOT/viralscan_work/dsr_round2b
#   mkdir -p $R/logs $R/status
#   cp config/dsr_round2b_sens07.template.tsv $R/config.tsv        # fill in; freeze before submitting
#   cp docs/plans/dsr_round2b_sens07_GUARD.md $R/GUARD.md
# Submit (from the repo dir; resources are the only thing an agent may change later; throttle %2):
#   ROUND_DIR=$R CODE_DIR=$VIRALSCAN_SITE_ROOT/vs_pinned/35de249 \
#     sbatch -o $R/logs/%x_%A_%a.log -J dsr2b_sens07 --array=0-16%2 --cpus-per-task=8 --mem=64G \
#       --time=48:00:00 scripts/slurm_dsr_round2b_sens07_array.sh
# Rows whose config `run` is `no` exit 0 at once. Dry run (prints commands, creates nothing):
#   DRY_RUN=1 bash scripts/slurm_dsr_round2b_sens07_array.sh 0 2       # rows 0 and 2
#
# Round 2c (PLAN SENS-07, user Q63/Q64) is this same script with SELECTION_SEED_MODE=per-family, a
# round-2c ROUND_DIR/CODE_DIR/CONFIG_TSV and config/dsr_round2c_sens07.template.tsv (a `trim_r1`
# column); see docs/plans/dsr_round2c_sens07_GUARD.md. Round 2b is unchanged: with neither, the
# `plant` arguments are exactly those of round 2b.
#
# Config TSV (CONFIG_TSV, default $ROUND_DIR/config.tsv): `-` means empty, as in manifest.tsv.
#   dataset  sample  run  expected_family  n_background  [trim_r1]  note   (columns are read by header
#   name; `trim_r1` is optional and defaults to `no`, which is the round-2b 6-column layout)
#   run              yes|no   (no = skip the row, recorded n/a; the note says why)
#   expected_family  family label(s) of the dataset's own virus, comma-separated, added to PLANT_FAMILIES;
#                    `-` = only PLANT_FAMILIES. Labels are those of virus_identity.tsv (EBV is
#                    `Herpesviridae`, the other human herpesviruses `Orthoherpesviridae`).
#   trim_r1          yes|no  yes = pass --trim-r1-to-geometry to `plant`: R1 longer than barcode + UMI is
#                    trimmed to that length in the planted AND the unplanted twin (a declared read change,
#                    recorded in plant_manifest.json). R1 shorter than barcode + UMI is still refused.
#   n_background     first N background read pairs to copy and plant into (`-` = N_BACKGROUND_DEFAULT).
#                    The plant step holds every (CB, UMI) of those pairs in memory (about 10 GB per 50M).
#
# Parameters (environment; defaults = the user's pending Q48/Q49 recommendations):
#   ROUND_DIR  CODE_DIR  REF_BUILD  ROUND2_DIR  CONFIG_TSV  EXACT_FASTA  FAMILY_IDENTITY
#   PLANT_FAMILIES="Orthoherpesviridae Herpesviridae Papillomaviridae Polyomaviridae"
#   PER_FAMILY=2   SEED=20261012   N_BACKGROUND_DEFAULT=50000000
#   PLANT_EXTRA=   extra `plant` arguments, e.g. "--levels 10 100 --cells-per-condition 10" or
#                  "--heldout-fasta H.fa --paf H.paf" (Q48 says exact only: leave empty)
#   SELECTION_SEED_MODE=global   global (round 2b: one RNG over all families) | per-family (round 2c, Q64:
#                  each family's genomes are drawn from sha256(SEED, family), the same in every dataset)
#   FORCE_TECH=1   pass --force-technology (GSM4682311 only)   DRY_RUN=1  print, do not run
#SBATCH -J dsr2b_sens07
#SBATCH -o slurm_%x_%j.log
set -euo pipefail
# >>> dsr-env (tests/test_dsr_round2b_sens07.py evaluates this block alone)
M=${VIRALSCAN_SITE_ROOT:?set VIRALSCAN_SITE_ROOT to the site work root}
ROUND2=${ROUND2_DIR:-$M/viralscan_work/dsr_round2}
R=${ROUND_DIR:-$M/viralscan_work/dsr_round2b}
CODE=${CODE_DIR:-$M/vs_pinned/35de249}
B=${REF_BUILD:-$M/viral_ref_cat42d/build}
CONFIG=${CONFIG_TSV:-$R/config.tsv}
EXACT=${EXACT_FASTA:-$B/viral.fa}
FAMILY_IDENTITY=${FAMILY_IDENTITY:-$ROUND2/runs/hsv1/combined_off/SRR8315713/SRR8315713/results/virus_identity.tsv}
# <<< dsr-env
PLANT_FAMILIES=${PLANT_FAMILIES:-Orthoherpesviridae Herpesviridae Papillomaviridae Polyomaviridae}
PER_FAMILY=${PER_FAMILY:-2}
SEED=${SEED:-20261012}
NBG_DEFAULT=${N_BACKGROUND_DEFAULT:-50000000}
SEL_MODE=${SELECTION_SEED_MODE:-global}
DRY=${DRY_RUN-}
REF_ID=${REF_ID:-$(basename "$(dirname "$B")")}
LEDGER=${EXPLORE_LEDGER:-$R/explore_ledger.tsv}
LEDGER_PY=$CODE/scripts/explore_ledger.py
SENS=$CODE/scripts/sens07_planted_sensitivity.py
export PATH="$M/conda/envs/viralscan_bench/bin:$M/conda/envs/R4_51/bin:$PATH" PYTHONPATH="$CODE/src" MPLBACKEND=Agg
if [ -n "$DRY" ]; then SHA=$(git -C "$CODE" describe --always --dirty 2>/dev/null || echo unknown); else
  SHA=$(git -C "$CODE" describe --always --dirty)
  case "$SHA" in *-dirty) echo "dirty pinned tree: $SHA" >&2; exit 2;; esac
fi

run() { printf '+'; printf ' %q' "$@"; echo; [ -n "$DRY" ] || "$@"; }   # echo the command, run it unless dry

one_row() {
  local IDX=$1 DS SAMPLE GRP R1 R2 WL PAIRS HINT
  IFS=$'\t' read -r DS SAMPLE GRP R1 R2 WL PAIRS HINT < <(tail -n +2 "$ROUND2/manifest.tsv" | sed -n "$((IDX + 1))p")
  [ -n "${SAMPLE:-}" ] || { echo "no manifest row $IDX" >&2; return 2; }
  [ "$WL" != - ] || WL=""
  local RUNFLAG EXP NBG TRIM NOTE
  # columns by header name (round 2b's 6-column config has no trim_r1: it reads as `no`)
  IFS=$'\t' read -r RUNFLAG EXP NBG TRIM NOTE < <(awk -F'\t' -v d="$DS" -v s="$SAMPLE" '
    NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    $1 == d && $2 == s { printf "%s\t%s\t%s\t%s\t%s\n", $c["run"], $c["expected_family"], $c["n_background"], ("trim_r1" in c ? $c["trim_r1"] : "no"), $c["note"] }' "$CONFIG") || true
  [ -n "${RUNFLAG:-}" ] || { echo "no config row for $DS/$SAMPLE in $CONFIG" >&2; return 3; }
  [ "$RUNFLAG" = yes ] || { echo "== row $IDX $DS/$SAMPLE skipped (config run=$RUNFLAG): $NOTE"; return 0; }
  [ "$EXP" != - ] || EXP=""
  [ "${NBG:--}" != - ] || NBG=$NBG_DEFAULT
  local XF=() F
  # the round-2c flags are added only when asked for, after PLANT_EXTRA, so round 2b's command line is unchanged
  case "$SEL_MODE" in global) ;; per-family) XF+=(--selection-seed-mode per-family);; *) echo "SELECTION_SEED_MODE must be global or per-family, not $SEL_MODE" >&2; return 2;; esac
  case "$TRIM" in yes) XF+=(--trim-r1-to-geometry);; no) ;; *) echo "config trim_r1 must be yes or no, not $TRIM ($DS/$SAMPLE)" >&2; return 3;; esac
  # round 2b may still be running: round-2c settings must never write into its directory (rm -rf below)
  case "$R" in */dsr_round2b | */dsr_round2b/) [ ${#XF[@]} -eq 0 ] || { echo "round-2c settings (${XF[*]}) with ROUND_DIR=$R: set ROUND_DIR to dsr_round2c" >&2; return 2; };; esac
  for F in ${XF[@]+"${XF[@]}"}; do
    case "$F" in --*) [ ! -f "$SENS" ] || grep -q -e "$F" "$SENS" || { echo "$SENS has no $F: pin CODE_DIR to a commit that does" >&2; return 2; };; esac
  done
  echo "== row $IDX $DS/$SAMPLE expected_family=${EXP:-none} n_background=$NBG per_family=$PER_FAMILY seed=$SEED selection=$SEL_MODE trim_r1=$TRIM code=$SHA ${DRY:+(DRY RUN)}"
  local CHEM=$ROUND2/chem/$SAMPLE.json
  [ -s "$CHEM" ] || { echo "no chemistry JSON for $SAMPLE" >&2; return 3; }
  local TECH STRAND
  read -r TECH STRAND < <(python - "$CHEM" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
print(d["chemistry"]["chemistry"] or "NONE", (d.get("strand_pilot") or {}).get("choice") or "NONE")
PY
)
  [ "$TECH" != NONE ] || { echo "chemistry unresolved for $SAMPLE" >&2; return 3; }
  # plant R2 sense to the transcript, or antisense for a --strand reverse library (5' chemistry)
  local LIBSTRAND=forward; [ "$STRAND" != reverse ] || LIBSTRAND=reverse
  [ -n "$DRY" ] || [ -s "$B/cdna.fa" ] || { echo "no $B/cdna.fa (sense transcripts to plant from)" >&2; return 3; }
  # Track B guard, as in slurm_dsr_rerun_array.sh
  local DENY=$CODE/config/track_b_denylist.txt
  [ -f "$DENY" ] || { echo "no deny-list at $DENY" >&2; return 2; }
  if printf '%s\n' "$DS $SAMPLE $R1 $R2" | grep -iqFf <(sed 's/#.*//;/^[[:space:]]*$/d' "$DENY"); then
    echo "deny-listed accession in manifest row: $DS/$SAMPLE" >&2; return 4
  fi
  local OUT=$R/runs/$DS/sens07/$SAMPLE META=$R/runs/$DS/sens07/$SAMPLE.meta
  if [ -s "$OUT/score/family_summary.tsv" ]; then echo "already done: $OUT"; return 0; fi

  # families.tsv (accession<TAB>family): panel genomes of the planted families, from the reference's identity table
  local FAMS="$PLANT_FAMILIES ${EXP//,/ }" FAMTSV
  # shellcheck disable=SC2086
  FAMTSV=$(python - "$FAMILY_IDENTITY" "$EXACT" $FAMS <<'PY'
import csv, sys
ident, fasta, want = sys.argv[1], sys.argv[2], set(sys.argv[3:])
have = {ln[1:].split()[0] for ln in open(fasta) if ln.startswith(">")}
rows = {}
for r in csv.DictReader(open(ident), delimiter="\t"):
    if r["viral"] == "true" and r["family"] in want and r["genome_accession"] in have:
        rows[r["genome_accession"]] = r["family"]
missing = want - set(rows.values())
if missing:
    sys.exit(f"family not in {ident} / the exact FASTA: {sorted(missing)}")
for acc in sorted(rows):
    print(f"{acc}\t{rows[acc]}")
PY
)
  echo "families: $(printf '%s\n' "$FAMTSV" | cut -f2 | sort | uniq -c | awk '{printf "%s %s ", $1, $2}')"
  echo "would create (dry run: nothing is): $OUT/{plant,run_planted,run_unplanted,score} and $META"
  if [ -z "$DRY" ]; then
    rm -rf "$OUT" "$META"; mkdir -p "$META" "$OUT"
    printf '%s\n' "$FAMTSV" > "$META/families.tsv"
    echo "$SHA" > "$META/code_sha.txt"; cp "$CHEM" "$META/check_chemistry.json"
    sha256sum "$FAMILY_IDENTITY" "$EXACT" "$CONFIG" > "$META/input_sha256.txt"
  fi
  local PX=(); read -ra PX <<< "${PLANT_EXTRA-}"
  run /usr/bin/time -v -o "$META/time_plant.txt" python "$SENS" plant --out-dir "$OUT/plant" --technology "$TECH" \
    --background-r1 "$R1" --background-r2 "$R2" --exact-fasta "$EXACT" --families "$META/families.tsv" \
    --transcripts-fasta "$B/cdna.fa" --library-strand "$LIBSTRAND" \
    --per-family "$PER_FAMILY" --seed "$SEED" --n-background "$NBG" ${PX[@]+"${PX[@]}"} ${XF[@]+"${XF[@]}"}
  # viralscan names the sample dir from the file name before the first "_"
  local ARM
  for ARM in planted unplanted; do
    local D=$OUT/plant/$ARM
    run ln -s R1.fastq.gz "$D/${SAMPLE}_R1.fastq.gz"; run ln -s R2.fastq.gz "$D/${SAMPLE}_R2.fastq.gz"
    # the same arguments as round 2 `combined_off` (slurm_dsr_rerun_array.sh)
    local ARGS=(-o "$OUT/run_$ARM" -s1 "$D/${SAMPLE}_R1.fastq.gz" -s2 "$D/${SAMPLE}_R2.fastq.gz"
                -i "$B/panel.idx" -t "$B/panel.t2g" -gtf "$B/viral_panel.gtf" -x "$TECH"
                -c "${SLURM_CPUS_PER_TASK:-4}" --cell-calling emptydrops --anellovirus-gene-ids --yes --verbose)
    [ "$STRAND" = NONE ] || ARGS+=(--strand "$STRAND")
    [ -z "$WL" ] || ARGS+=(-w "$WL")
    [ -z "${FORCE_TECH:-}" ] || ARGS+=(--force-technology)
    run /usr/bin/time -v -o "$META/time_$ARM.txt" python -m viralscan.menu "${ARGS[@]}"
    [ -n "$DRY" ] || test -s "$OUT/run_$ARM/$SAMPLE/results/viral_summary.tsv"
    run python "$LEDGER_PY" --ledger "$LEDGER" append --dataset "$DS" --sample "$SAMPLE" \
      --arm "combined_off_$ARM" --code-sha "$SHA" --reference-id "$REF_ID"
  done
  local RUN=$OUT/run_planted/$SAMPLE TWIN=$OUT/run_unplanted/$SAMPLE
  run python "$SENS" score --out-dir "$OUT/score" --truth "$OUT/plant/truth.tsv" \
    --counts-h5ad "$RUN/kb-python/counts_unfiltered/adata_multimap.h5ad" \
    --twin-h5ad "$TWIN/kb-python/counts_unfiltered/adata_multimap.h5ad" \
    --layer X --identity "$RUN/results/virus_identity.tsv"
  echo "${DRY:+DRY RUN: }DONE $DS/$SAMPLE sens07"
}

if [ "$#" -gt 0 ]; then ROWS=("$@"); else ROWS=("${SLURM_ARRAY_TASK_ID:?array index or row arguments}"); fi
for i in "${ROWS[@]}"; do one_row "$i"; done
