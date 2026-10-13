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
#   R=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round2b
#   mkdir -p $R/logs $R/status
#   cp config/dsr_round2b_sens07.template.tsv $R/config.tsv        # fill in; freeze before submitting
#   cp docs/plans/dsr_round2b_sens07_GUARD.md $R/GUARD.md
# Submit (from the repo dir; resources are the only thing an agent may change later; throttle %2):
#   ROUND_DIR=$R CODE_DIR=/exports/archive/hg-funcgenom-research/mdmanurung/vs_pinned/35de249 \
#     sbatch -o $R/logs/%x_%A_%a.log -J dsr2b_sens07 --array=0-16%2 --cpus-per-task=8 --mem=64G \
#       --time=48:00:00 scripts/slurm_dsr_round2b_sens07_array.sh
# Rows whose config `run` is `no` exit 0 at once. Dry run (prints commands, creates nothing):
#   DRY_RUN=1 bash scripts/slurm_dsr_round2b_sens07_array.sh 0 2       # rows 0 and 2
#
# Config TSV (CONFIG_TSV, default $ROUND_DIR/config.tsv): `-` means empty, as in manifest.tsv.
#   dataset  sample  run  expected_family  n_background  note
#   run              yes|no   (no = skip the row, recorded n/a; the note says why)
#   expected_family  family label(s) of the dataset's own virus, comma-separated, added to PLANT_FAMILIES;
#                    `-` = only PLANT_FAMILIES. Labels are those of virus_identity.tsv (EBV is
#                    `Herpesviridae`, the other human herpesviruses `Orthoherpesviridae`).
#   n_background     first N background read pairs to copy and plant into (`-` = N_BACKGROUND_DEFAULT).
#                    The plant step holds every (CB, UMI) of those pairs in memory (about 10 GB per 50M).
#
# Parameters (environment; defaults = the user's pending Q48/Q49 recommendations):
#   ROUND_DIR  CODE_DIR  REF_BUILD  ROUND2_DIR  CONFIG_TSV  EXACT_FASTA  FAMILY_IDENTITY
#   PLANT_FAMILIES="Orthoherpesviridae Herpesviridae Papillomaviridae Polyomaviridae"
#   PER_FAMILY=2   SEED=20261012   N_BACKGROUND_DEFAULT=50000000
#   PLANT_EXTRA=   extra `plant` arguments, e.g. "--levels 10 100 --cells-per-condition 10" or
#                  "--heldout-fasta H.fa --paf H.paf" (Q48 says exact only: leave empty)
#   FORCE_TECH=1   pass --force-technology (GSM4682311 only)   DRY_RUN=1  print, do not run
#SBATCH -J dsr2b_sens07
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round2b/logs/%x_%A_%a.log
set -euo pipefail
# >>> dsr-env (tests/test_dsr_round2b_sens07.py evaluates this block alone)
M=/exports/archive/hg-funcgenom-research/mdmanurung
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
  local RUNFLAG EXP NBG NOTE
  IFS=$'\t' read -r _ _ RUNFLAG EXP NBG NOTE < <(awk -F'\t' -v d="$DS" -v s="$SAMPLE" 'NR > 1 && $1 == d && $2 == s' "$CONFIG") || true
  [ -n "${RUNFLAG:-}" ] || { echo "no config row for $DS/$SAMPLE in $CONFIG" >&2; return 3; }
  [ "$RUNFLAG" = yes ] || { echo "== row $IDX $DS/$SAMPLE skipped (config run=$RUNFLAG): $NOTE"; return 0; }
  [ "$EXP" != - ] || EXP=""
  [ "${NBG:--}" != - ] || NBG=$NBG_DEFAULT
  echo "== row $IDX $DS/$SAMPLE expected_family=${EXP:-none} n_background=$NBG per_family=$PER_FAMILY seed=$SEED code=$SHA ${DRY:+(DRY RUN)}"
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
  echo "families: $(printf '%s\n' "$FAMTSV" | cut -f2 | sort | uniq -c | tr '\n' ' ')"
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
    --per-family "$PER_FAMILY" --seed "$SEED" --n-background "$NBG" ${PX[@]+"${PX[@]}"}
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
