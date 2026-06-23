#!/usr/bin/env bash
# P22.6 — STARsolo comparison on EBV dataset SRR12682296 (GSE158275)
#
# Builds a combined GRCh38 + EBV STAR genome, runs STARsolo in 10xv2 mode,
# then calls compare_starsolo_viralscan.py to compare infected-cell rates
# against the ViralScan kallisto output.
#
# Submit:
#   sbatch scripts/slurm_starsolo_ebv_comparison.sh
#
# After the job completes, find results at:
#   $WORKDIR/comparison_starsolo_vs_viralscan.tsv
#
#SBATCH -J ebv_starsolo_p22.6
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH -o logs/ebv_starsolo_%j.log
#SBATCH -e logs/ebv_starsolo_%j.err

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO=/exports/para-lipg-hpc/mdmanurung/ViralScan
WORKDIR=$REPO/starsolo_p22_6

STAR=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/starsolo/bin/STAR
PYTHON=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/scale_py/bin/python

GRCH38_FA=/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa
GRCH38_GTF=/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf
EBV_FA=/exports/archive/hg-funcgenom-research/evonk/old/intern/fasta_viruses/Serratus/fasta_split/Epstein-Barr_virus_NC_007605.fasta
EBV_GTF=/exports/archive/hg-funcgenom-research/evonk/old/intern/fasta_viruses/Serratus/split_gtf/Epstein_Barr_virus_NC_007605.gtf

SRR=SRR12682296
GENOME_DIR=$WORKDIR/genome_GRCh38_EBV
COMBINED_FA=$WORKDIR/GRCh38_EBV_combined.fa
COMBINED_GTF=$WORKDIR/GRCh38_EBV_combined.gtf
FASTQ_DIR=$WORKDIR/fastq
STARSOLO_OUT=$WORKDIR/starsolo_ebv

mkdir -p "$WORKDIR" "$FASTQ_DIR" "$GENOME_DIR" "$REPO/logs"

# ── STEP 1: Build combined GRCh38 + EBV STAR genome index ─────────────────────
echo "[$(date)] STEP 1: Build combined genome (skip if SA file exists)"

if [[ -f "$GENOME_DIR/SA" ]]; then
    echo "  [SKIP] Genome already built at $GENOME_DIR"
else
    echo "  Concatenating GRCh38 + EBV FASTA..."
    cat "$GRCH38_FA" "$EBV_FA" > "$COMBINED_FA"

    echo "  Concatenating GRCh38 + EBV GTF..."
    cat "$GRCH38_GTF" "$EBV_GTF" > "$COMBINED_GTF"

    echo "  Building STAR genome index (8 threads, 64 GB) ..."
    "$STAR" \
        --runMode genomeGenerate \
        --genomeDir "$GENOME_DIR" \
        --genomeFastaFiles "$COMBINED_FA" \
        --sjdbGTFfile "$COMBINED_GTF" \
        --genomeSAindexNbases 14 \
        --runThreadN 8

    echo "  [DONE] Genome built: $GENOME_DIR"
fi

# ── STEP 2: Download SRR12682296 FASTQs ───────────────────────────────────────
echo "[$(date)] STEP 2: Download $SRR FASTQs (skip if present)"

R1=$FASTQ_DIR/${SRR}_1.fastq.gz
R2=$FASTQ_DIR/${SRR}_2.fastq.gz

if [[ -f "$R2" ]]; then
    echo "  [SKIP] FASTQs already present."
else
    if command -v fasterq-dump >/dev/null 2>&1; then
        echo "  Downloading via fasterq-dump..."
        fasterq-dump "$SRR" \
            --outdir "$FASTQ_DIR" \
            --threads 8 \
            --temp "$FASTQ_DIR/tmp_sra"
        gzip "$FASTQ_DIR/${SRR}_1.fastq" "$FASTQ_DIR/${SRR}_2.fastq"
    elif command -v prefetch >/dev/null 2>&1; then
        echo "  Trying prefetch + fasterq-dump..."
        prefetch "$SRR" -O "$FASTQ_DIR"
        fasterq-dump "$FASTQ_DIR/$SRR" --outdir "$FASTQ_DIR" --threads 8
        gzip "$FASTQ_DIR/${SRR}_1.fastq" "$FASTQ_DIR/${SRR}_2.fastq"
    else
        echo "ERROR: neither fasterq-dump nor prefetch found. Install sra-tools." >&2
        exit 1
    fi
    echo "  [DONE] FASTQs downloaded."
fi

# ── STEP 3: Run STARsolo (10xv2 chemistry) ────────────────────────────────────
# 10xv2 layout: R1 = 26 bp (16 bp CB + 10 bp UMI), R2 = cDNA
# STARsolo readFilesIn expects: <cDNA_R2> <barcode_R1>
echo "[$(date)] STEP 3: Run STARsolo (skip if filtered matrix exists)"

FILTERED_DIR="$STARSOLO_OUT/Solo.out/GeneFull/filtered"

if [[ -d "$FILTERED_DIR" ]]; then
    echo "  [SKIP] STARsolo output already exists at $STARSOLO_OUT"
else
    echo "  Running STARsolo (10xv2, GeneFull) ..."
    "$STAR" \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 --soloCBlen 16 \
        --soloUMIstart 17 --soloUMIlen 10 \
        --soloCBwhitelist None \
        --soloFeatures GeneFull \
        --soloCellFilter CellRanger2 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$R2" "$R1" \
        --readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --runThreadN 8 \
        --outFileNamePrefix "$STARSOLO_OUT/"

    echo "  [DONE] STARsolo complete."
fi

# ── STEP 4: Compare STARsolo vs ViralScan ─────────────────────────────────────
echo "[$(date)] STEP 4: Run comparison script"

# Find ViralScan summary from the 1M-read dry-run if present
VS_SUMMARY=""
VS_CANDIDATE=$(ls "$REPO"/results*/${SRR}*/viral_detection/viral_summary_*.tsv 2>/dev/null | head -1 || true)
if [[ -n "$VS_CANDIDATE" ]]; then
    VS_SUMMARY="--viralscan-summary $VS_CANDIDATE"
    echo "  Using ViralScan summary: $VS_CANDIDATE"
else
    echo "  No ViralScan summary found; comparison will show STARsolo numbers only."
fi

"$PYTHON" "$REPO/scripts/compare_starsolo_viralscan.py" \
    --starsolo-dir "$STARSOLO_OUT/Solo.out" \
    --out "$WORKDIR/comparison_starsolo_vs_viralscan.tsv" \
    $VS_SUMMARY

echo "[$(date)] P22.6 complete."
echo "  Comparison table: $WORKDIR/comparison_starsolo_vs_viralscan.tsv"
echo "  STARsolo logs:    $STARSOLO_OUT/Log.final.out"
