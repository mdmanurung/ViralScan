#!/usr/bin/env bash
# P22.10 — STARsolo EBV comparison WITH 10x v2 whitelist (matched-barcode run)
#
# Reuses the GRCh38+EBV genome index and FASTQs built by slurm_starsolo_ebv_comparison.sh
# (P22.6). Adds --soloCBwhitelist so barcode correction matches the paper's CellRanger space.
#
# Submit:
#   sbatch scripts/slurm_starsolo_ebv_wl.sh
#
# After completion, find the raw matrix (for matched-barcode intersection) at:
#   $WORKDIR/starsolo_ebv_wl/Solo.out/GeneFull/raw/
#
#SBATCH -J ebv_starsolo_wl_p22.10
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH -o logs/ebv_starsolo_wl_%j.log
#SBATCH -e logs/ebv_starsolo_wl_%j.err

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO=/exports/para-lipg-hpc/mdmanurung/ViralScan

# Reuse P22.6 genome index and FASTQs (already built/downloaded)
GENOME_DIR=$REPO/starsolo_p22_6/genome_GRCh38_EBV
FASTQ_DIR=$REPO/starsolo_p22_6/fastq
SRR=SRR12682296
R1=$FASTQ_DIR/${SRR}_1.fastq.gz
R2=$FASTQ_DIR/${SRR}_2.fastq.gz

# New output directory for the whitelist run (preserve P22.6 no-whitelist output)
WORKDIR=$REPO/starsolo_p22_6b
STARSOLO_OUT=$WORKDIR/starsolo_ebv_wl

# 10x v2 whitelist (737,280 barcodes; extracted from ngs_tools in Step 0)
WHITELIST=$REPO/ref/10x_version2_whitelist.txt

STAR=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/starsolo/bin/STAR

mkdir -p "$WORKDIR" "$REPO/logs"

# ── Verify prerequisites ───────────────────────────────────────────────────────
if [[ ! -f "$GENOME_DIR/SA" ]]; then
    echo "ERROR: Genome index not found at $GENOME_DIR" >&2
    echo "  Run slurm_starsolo_ebv_comparison.sh first to build the P22.6 genome." >&2
    exit 1
fi
if [[ ! -f "$R2" ]]; then
    echo "ERROR: FASTQ not found: $R2" >&2
    exit 1
fi
if [[ ! -f "$WHITELIST" ]]; then
    echo "ERROR: Whitelist not found: $WHITELIST" >&2
    echo "  Run: zcat <ngs_tools>/chemistry/whitelists/10x_version2_whitelist.txt.gz > $WHITELIST" >&2
    exit 1
fi
echo "[$(date)] Prerequisites verified."
echo "  Genome: $GENOME_DIR"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Whitelist: $WHITELIST ($(wc -l < "$WHITELIST") barcodes)"

# ── STARsolo (10xv2, GeneFull, WITH 10x v2 whitelist) ─────────────────────────
echo "[$(date)] STEP 1: Run STARsolo with 10x v2 whitelist"

RAW_DIR="$STARSOLO_OUT/Solo.out/GeneFull/raw"

if [[ -d "$RAW_DIR" ]]; then
    echo "  [SKIP] STARsolo output already exists at $STARSOLO_OUT"
else
    echo "  Running STARsolo (10xv2, GeneFull, soloCBwhitelist) ..."
    "$STAR" \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 --soloCBlen 16 \
        --soloUMIstart 17 --soloUMIlen 10 \
        --soloCBwhitelist "$WHITELIST" \
        --soloFeatures GeneFull \
        --soloCellFilter CellRanger2.2 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$R2" "$R1" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --runThreadN 8 \
        --outFileNamePrefix "$STARSOLO_OUT/"

    echo "  [DONE] STARsolo complete."
fi

# ── Quick summary ──────────────────────────────────────────────────────────────
FILTERED_DIR="$STARSOLO_OUT/Solo.out/GeneFull/filtered"
if [[ -f "$FILTERED_DIR/barcodes.tsv" ]] || [[ -f "$FILTERED_DIR/barcodes.tsv.gz" ]]; then
    n_filtered=$(zcat "$FILTERED_DIR/barcodes.tsv.gz" 2>/dev/null | wc -l || wc -l < "$FILTERED_DIR/barcodes.tsv")
    echo "[$(date)] Filtered cells (CellRanger2.2): $n_filtered"
fi
RAW_BC=$(zcat "$RAW_DIR/barcodes.tsv.gz" 2>/dev/null | wc -l || wc -l < "$RAW_DIR/barcodes.tsv" || echo "?")
echo "[$(date)] Raw barcodes (whitelist-corrected): $RAW_BC"
echo "[$(date)] P22.10 STARsolo-wl run complete."
echo "  Raw matrix for matched-barcode intersection: $RAW_DIR"
