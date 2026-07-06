#!/usr/bin/env bash
# Stage 2 v2 — Repair kb ref using cDNA-level GTF
#
# Background: job 25137178 hung at "Splitting genome" for 3h+ with 0 bytes output.
# Root cause: combined.gtf had chromosomal seqnames (1, 2, X, …) that don't match
# the cDNA FASTA headers (ENST transcript IDs), so ngs_tools could not find any
# human sequences.
#
# This job:
#   A. Generates combined_cdna.gtf from FASTA headers (seqname = ENST transcript ID)
#   B. Runs kb ref with combined.fa + combined_cdna.gtf
#
# Inputs (completed by job 25137178 Steps 1+2):
#   combined.fa     — host cDNA + SARS-CoV-2 + Serratus/anellovirus panel (1.4 GB)
#   viral_whole_genome.gtf — SARS-CoV-2 gene models
#   viral_genome.gtf       — Serratus+anellovirus panel gene models
#
# Outputs:
#   combined_cdna.gtf   — new compatible GTF
#   index.idx           — kallisto index
#   t2g.txt             — transcript-to-gene table
#   cdna.fa             — extracted cDNA sequences
#
#SBATCH -J covid_vs_buildref_v2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH -o covid_viralscan/logs/build_ref_v2_%j.log
#SBATCH -e covid_viralscan/logs/build_ref_v2_%j.err

set -euo pipefail

REPO=/exports/para-lipg-hpc/mdmanurung/ViralScan
VS_CONDA_ENV=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
export PATH="$VS_CONDA_ENV/bin:$PATH"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"

OUTDIR=$REPO/covid_viralscan/viralscan_ref
VIRAL_GTF=$OUTDIR/viral/viral_whole_genome.gtf
PANEL_GTF=$REPO/references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf
COMBINED_FA=$OUTDIR/combined.fa
COMBINED_CDNA_GTF=$OUTDIR/combined_cdna.gtf
INDEX=$OUTDIR/index.idx
T2G=$OUTDIR/t2g.txt
CDNA=$OUTDIR/cdna.fa

mkdir -p "$REPO/covid_viralscan/logs"

echo "[$(date)] Stage 2 v2 starting on $(hostname)"

# ── Prerequisite checks ───────────────────────────────────────────────────────
for f in "$COMBINED_FA" "$VIRAL_GTF" "$PANEL_GTF"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
command -v kb >/dev/null || { echo "ERROR: kb not on PATH" >&2; exit 1; }
command -v python >/dev/null || { echo "ERROR: python not on PATH" >&2; exit 1; }

echo "  combined.fa:   $(du -sh "$COMBINED_FA" | cut -f1)"
echo "  viral GTF:     $VIRAL_GTF  ($(wc -l < "$VIRAL_GTF") lines)"
echo "  panel GTF:     $PANEL_GTF  ($(wc -l < "$PANEL_GTF") lines)"
echo ""

# ── Step A: Generate cDNA-level GTF from FASTA headers ───────────────────────
if [[ -f "$COMBINED_CDNA_GTF" ]]; then
    echo "[$(date)] SKIP step A — combined_cdna.gtf already exists ($(wc -l < "$COMBINED_CDNA_GTF") lines)"
else
    echo "[$(date)] Step A: generating combined_cdna.gtf ..."
    python "$REPO/covid_viralscan/scripts/gen_combined_cdna_gtf.py" \
        "$COMBINED_FA" \
        "$VIRAL_GTF" \
        "$PANEL_GTF" \
        "$COMBINED_CDNA_GTF"
    echo "[$(date)] Step A complete."
    echo "  combined_cdna.gtf: $(wc -l < "$COMBINED_CDNA_GTF") lines, $(du -sh "$COMBINED_CDNA_GTF" | cut -f1)"
fi

# Verify the cDNA GTF has ENST entries (sanity check)
ENST_LINES=$(grep -c "^ENST" "$COMBINED_CDNA_GTF" || true)
if [[ "$ENST_LINES" -eq 0 ]]; then
    echo "ERROR: combined_cdna.gtf has no ENST entries — GTF generation failed." >&2
    exit 1
fi
echo "  ENST lines in new GTF: $ENST_LINES"
NC_LINES=$(grep -c "^NC_" "$COMBINED_CDNA_GTF" || true)
echo "  NC_ (viral) lines:     $NC_LINES"

# ── Step B: kb ref ────────────────────────────────────────────────────────────
if [[ -f "$INDEX" && -f "$T2G" && -f "$CDNA" ]]; then
    echo "[$(date)] SKIP step B — index.idx, t2g.txt, cdna.fa already exist"
else
    echo "[$(date)] Step B: kb ref (threads: ${SLURM_CPUS_PER_TASK:-16}) ..."
    kb ref \
        -i "$INDEX" \
        -g "$T2G" \
        -f1 "$CDNA" \
        "$COMBINED_FA" \
        "$COMBINED_CDNA_GTF"
    echo "[$(date)] Step B complete."
    echo "  index.idx: $(du -sh "$INDEX"    | cut -f1)"
    echo "  t2g.txt:   $(wc -l  < "$T2G")  transcript-to-gene entries"
    echo "  cdna.fa:   $(du -sh "$CDNA"    | cut -f1)"
fi

# ── Post-build verification ──────────────────────────────────────────────────
echo ""
echo "[$(date)] Post-build sanity checks ..."
COV2=$(grep -c "NC_045512" "$T2G"  || echo 0)
ENST_T2G=$(grep -c "ENST"      "$T2G"  || echo 0)
GENE_T2G=$(grep -c "_gene"     "$T2G"  || echo 0)
SARSP=$(grep -c "sarsp"        "$T2G"  || echo 0)
echo "  NC_045512 rows: $COV2    (SARS-CoV-2; expect ≥1)"
echo "  ENST rows:      $ENST_T2G  (human transcripts; expect ≥100k)"
echo "  _gene rows:     $GENE_T2G  (anellovirus gene entries; expect ≥1000)"
echo "  sarsp rows:     $SARSP     (SARS-CoV-1 locus tags; expect ≥1)"

[[ "$COV2"    -ge 1      ]] || echo "WARNING: no SARS-CoV-2 in t2g"
[[ "$ENST_T2G" -ge 100000 ]] || echo "WARNING: very few ENST entries in t2g"
[[ "$SARSP"   -ge 1      ]] || echo "WARNING: no sarsp entries in t2g"

echo ""
echo "[$(date)] Stage 2 v2 complete."
echo "  Index:  $INDEX"
echo "  T2G:    $T2G"
echo "  Proceed to Stage 3: sbatch --array=0-1 covid_viralscan/scripts/slurm_viralscan_quant.sh"
