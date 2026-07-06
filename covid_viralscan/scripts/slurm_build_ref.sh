#!/usr/bin/env bash
# Stage 2 — Build combined kallisto reference for SARS-CoV-2 covid_viralscan analysis
#
# Produces a combined index covering:
#   - Human transcriptome (Homo_sapiens, Ensembl current release)
#   - SARS-CoV-2 (NC_045512.2, downloaded from NCBI)
#   - Full Serratus + expanded anellovirus viral panel (spliced in from disk)
#
# Steps:
#   1. viralscan build-ref --no-kb-ref --no-anellovirus
#      → emits combined.fa (host cDNA + SARS-CoV-2 genome)
#        and combined.gtf (host GTF + SARS-CoV-2 whole-genome GTF)
#   2. Append Serratus+anellovirus panel to combined.{fa,gtf}
#   3. kb ref → index.idx, t2g.txt, cdna.fa
#
# Uses:
#   - test_viralscan conda env (provides kb + viralscan package)
#   - PYTHONPATH override to run from repo (for --no-kb-ref flag)
#
# Submit:
#   sbatch covid_viralscan/scripts/slurm_build_ref.sh
#
#SBATCH -J covid_vs_buildref
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH -o covid_viralscan/logs/build_ref_%j.log
#SBATCH -e covid_viralscan/logs/build_ref_%j.err

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO=/exports/para-lipg-hpc/mdmanurung/ViralScan
VS_CONDA_ENV=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan

export PATH="$VS_CONDA_ENV/bin:$PATH"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

OUTDIR=$REPO/covid_viralscan/viralscan_ref
PANEL_FA=$REPO/references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.fa
PANEL_GTF=$REPO/references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf

# Ensembl download cache — keep outside repo to avoid re-downloading on reruns
CACHE_DIR=/exports/para-lipg-hpc/mdmanurung/viralscan_cache

mkdir -p "$OUTDIR" "$CACHE_DIR" "$REPO/covid_viralscan/logs"

# ── Verify prerequisites ───────────────────────────────────────────────────────
for f in "$PANEL_FA" "$PANEL_GTF"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
command -v kb   || { echo "ERROR: kb not on PATH"; exit 1; }
command -v python || { echo "ERROR: python not on PATH"; exit 1; }

echo "[$(date)] Prerequisites OK."
echo "  Serratus panel FA:  $PANEL_FA  ($(wc -l < "$PANEL_FA") lines)"
echo "  Serratus panel GTF: $PANEL_GTF  ($(wc -l < "$PANEL_GTF") lines)"
echo "  Output dir:         $OUTDIR"

# ── Step 1: build-ref (no indexing, no anellovirus — panel covers it) ─────────
COMBINED_FA=$OUTDIR/combined.fa
COMBINED_GTF=$OUTDIR/combined.gtf

if [[ -f "$COMBINED_FA" && -f "$COMBINED_GTF" ]]; then
    echo "[$(date)] SKIP step 1 — combined.fa/gtf already exist (rerun protection)"
else
    echo "[$(date)] Step 1: viralscan build-ref (NC_045512.2, --no-kb-ref, --no-anellovirus) ..."
    python "$REPO/src/viralscan/menu.py" build-ref \
        --host human \
        --virus-accessions NC_045512.2 \
        --no-anellovirus \
        --no-kb-ref \
        --output "$OUTDIR" \
        --cache-dir "$CACHE_DIR" \
        --ncbi-email mikhael.manurung@gmail.com
    echo "[$(date)] Step 1 complete. combined.fa: $(du -sh "$COMBINED_FA" | cut -f1)"
fi

# ── Step 2: splice in Serratus+anellovirus panel ──────────────────────────────
# Check if panel was already appended (idempotency: look for a Serratus accession)
if grep -q "SERRATUS\|ANELLO\|HUM_HERP6B" "$COMBINED_FA" 2>/dev/null; then
    echo "[$(date)] SKIP step 2 — Serratus panel already in combined.fa"
else
    echo "[$(date)] Step 2: appending Serratus+anellovirus panel ..."
    cat "$PANEL_FA"  >> "$COMBINED_FA"
    cat "$PANEL_GTF" >> "$COMBINED_GTF"
    echo "[$(date)] Step 2 complete. combined.fa now: $(du -sh "$COMBINED_FA" | cut -f1)"
    echo "           combined.gtf now: $(du -sh "$COMBINED_GTF" | cut -f1)"
fi

# ── Step 2.5: sanitise combined GTF for kb ref ───────────────────────────────
# Two issues from the Serratus panel:
#   1. Non-unique transcript IDs (unassigned_transcript_N shared across many
#      different viral genes) — prefix with seqname to make globally unique.
#      Pattern is exact so this is idempotent on repeated runs.
#   2. Blank lines and GFF3 '###' flush directives — ngs_tools (used by kb ref)
#      cannot parse these and raises GtfEntryError.
echo "[$(date)] Step 2.5: sanitising combined.gtf (dedup transcript IDs, strip blank/### lines) ..."
awk '
  /^$/ || /^###/ { next }
  /transcript_id "unassigned_transcript/ {
    gsub(/transcript_id "unassigned_transcript/, "transcript_id \"" $1 "-unassigned_transcript")
  }
  { print }
' "$COMBINED_GTF" > "${COMBINED_GTF}.clean"
mv "${COMBINED_GTF}.clean" "$COMBINED_GTF"
echo "[$(date)] Step 2.5 complete. Lines: $(wc -l < "$COMBINED_GTF")"

# ── Step 2.6: dedup combined.fa by sequence name (keep-first) ─────────────────
# The Serratus panel contains at least one accession (NC_002076.2, Torque teno
# virus 1) TWICE — it is in both the anellovirus and Serratus sets. Duplicate
# FASTA headers make `kallisto index` abort with "repeated name in FASTA file",
# and ngs_tools extracts each overlapping gene model once per duplicate copy.
# Removing byte-identical duplicates is a k-mer no-op. Idempotent (dedup of a
# deduped file is a fixed point).
DUPES=$(grep '^>' "$COMBINED_FA" | awk '{print $1}' | sort | uniq -d | wc -l)
if [[ "$DUPES" -gt 0 ]]; then
    echo "[$(date)] Step 2.6: removing $DUPES duplicate accession(s) from combined.fa ..."
    awk '/^>/{keep=!seen[$1]++} keep' "$COMBINED_FA" > "${COMBINED_FA}.dedup"
    mv "${COMBINED_FA}.dedup" "$COMBINED_FA"
    echo "[$(date)] Step 2.6 complete. combined.fa now: $(grep -c '^>' "$COMBINED_FA") sequences"
else
    echo "[$(date)] Step 2.6: no duplicate accessions in combined.fa — skipping"
fi

# ── Step 3: kb ref ────────────────────────────────────────────────────────────
INDEX=$OUTDIR/index.idx
T2G=$OUTDIR/t2g.txt
CDNA=$OUTDIR/cdna.fa

if [[ -f "$INDEX" && -f "$T2G" ]]; then
    echo "[$(date)] SKIP step 3 — index.idx and t2g.txt already exist"
else
    echo "[$(date)] Step 3: kb ref (threads: $SLURM_CPUS_PER_TASK) ..."
    kb ref \
        -i "$INDEX" \
        -g "$T2G" \
        -f1 "$CDNA" \
        "$COMBINED_FA" \
        "$COMBINED_GTF"
    echo "[$(date)] Step 3 complete."
    echo "  index.idx: $(du -sh "$INDEX" | cut -f1)"
    echo "  t2g.txt:   $(wc -l < "$T2G") transcript-to-gene entries"
    echo "  cdna.fa:   $(du -sh "$CDNA" | cut -f1)"
fi

echo "[$(date)] Reference build complete."
echo "  Index:  $INDEX"
echo "  T2G:    $T2G"
echo "  Use viral GTF for -gtf: $PANEL_GTF"
