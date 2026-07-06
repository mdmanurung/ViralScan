#!/usr/bin/env bash
# PR 23 — bulk panel index build, resumed from kb ref with all covid Stage-2 lessons baked in.
#
# Background: build_bundled_panel_ref.py (job 25138575) completed Steps 1-6 — it fetched
# host cDNA + 194 curated + 2022 anello and wrote combined.fa (1.4 GB, 465,769 host ENST +
# 2216 viral) — but died at Step 7 (`kb ref`) because `kb` was not on PATH. Two further
# latent bugs would also have bitten (both hit covid Stage 2):
#   (1) combined.gtf carried Ensembl *chromosomal* seqnames (1, 2, X) that don't match the
#       ENST cDNA FASTA headers → kb ref hangs forever at "Splitting genome".
#   (2) `kallisto index` aborts on any repeated FASTA target name.
#
# This job uses a pre-generated cDNA-level GTF (combined_cdna.gtf, seqname = ENST/viral
# accession, validated locally) and adds a keep-first dedup fallback + artifact verification.
# It does NOT re-fetch (combined.fa/combined_cdna.gtf already on disk).
#
#SBATCH -J build_panel_kbref
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH -o /exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref/build_panel_kbref_%j.log
#SBATCH -e /exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref/build_panel_kbref_%j.err

set -euo pipefail

VS_CONDA_ENV=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
export PATH="$VS_CONDA_ENV/bin:$PATH"          # <- the fix the failed run was missing
KALLISTO=$VS_CONDA_ENV/lib/python3.12/site-packages/kb_python/bins/linux/kallisto/kallisto

OUT=/exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref/ref
COMBINED_FA=$OUT/combined.fa
CDNA_GTF=$OUT/combined_cdna.gtf
INDEX=$OUT/panel.idx
T2G=$OUT/panel.t2g
CDNA=$OUT/cdna.fa

echo "[$(date)] Panel kb ref on $(hostname)"
for f in "$COMBINED_FA" "$CDNA_GTF"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
command -v kb >/dev/null || { echo "ERROR: kb not on PATH" >&2; exit 1; }

# Fail-fast guard against the covid-#1 hang: the cDNA GTF must have NO chromosomal seqnames.
CHR_ROWS=$(awk '{print $1}' "$CDNA_GTF" | grep -cE '^[0-9]+$' || true)
[[ "$CHR_ROWS" -eq 0 ]] || { echo "ERROR: $CDNA_GTF has $CHR_ROWS chromosomal seqname rows — would hang kb ref." >&2; exit 1; }
echo "  GTF seqname guard OK (0 chromosomal rows)."

# ── kb ref (extraction + index). May fail at the index step on duplicate names. ──
rm -f "$INDEX"
set +e
kb ref -i "$INDEX" -g "$T2G" -f1 "$CDNA" "$COMBINED_FA" "$CDNA_GTF"
KB_RC=$?
set -e
echo "[$(date)] kb ref returned $KB_RC."

# ── Verify by ARTIFACT, not exit code. If the index is missing/empty but cdna.fa exists,
#    the failure is almost certainly duplicate FASTA names → dedup + index directly. ──
if [[ ! -s "$INDEX" ]]; then
    echo "[$(date)] panel.idx missing/empty — applying keep-first dedup fallback (covid-#2)."
    [[ -s "$CDNA" ]] || { echo "ERROR: cdna.fa also missing — kb ref failed earlier (check .err)." >&2; exit 1; }
    DUPES=$(grep '^>' "$CDNA" | awk '{print $1}' | sort | uniq -d | wc -l)
    echo "  duplicate target names in cdna.fa: $DUPES"
    awk '/^>/{keep=!seen[$1]++} keep' "$CDNA"        > "$CDNA.dedup"        && mv -f "$CDNA.dedup" "$CDNA"
    awk '!seen[$1]++'                 "$T2G"         > "$T2G.dedup"         && mv -f "$T2G.dedup" "$T2G"
    awk '/^>/{keep=!seen[$1]++} keep' "$COMBINED_FA" > "$COMBINED_FA.dedup" && mv -f "$COMBINED_FA.dedup" "$COMBINED_FA"
    echo "[$(date)] rebuilding index on deduped inputs ..."
    "$KALLISTO" index -i "$INDEX" -k 31 -t 8 -d "$COMBINED_FA" "$CDNA"
fi

[[ -s "$INDEX" ]] || { echo "ERROR: panel.idx still empty after fallback — aborting." >&2; exit 1; }

echo ""
echo "[$(date)] Verify by artifact:"
echo "  panel.idx size: $(du -sh "$INDEX" | cut -f1)"
"$KALLISTO" inspect "$INDEX" 2>&1 | sed 's/^/    /' || true
echo "  t2g entries:  $(wc -l < "$T2G")"
echo "  ENST rows:    $(grep -c ENST  "$T2G" || echo 0)   (host; expect >=100k)"
echo "  _gene rows:   $(grep -c _gene "$T2G" || echo 0)   (anello; expect >0)"

echo ""
echo "[$(date)] Panel build complete."
echo "  Index: $INDEX"
echo "  T2G:   $T2G"
echo "  Next: P23.op3 format probe (kb count -x BULK on one GSE128078 sample)."
