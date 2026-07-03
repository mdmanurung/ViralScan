#!/usr/bin/env bash
# STARsolo combined host+viral run on the covid_viralscan samples.
#
# Implements the "CellRanger-with-combined-references" idea (STARsolo is the
# open-source CellRanger equivalent; CellRanger isn't on this cluster). Builds a
# combined GRCh38 + viral-panel + SARS-CoV-2 STAR genome and quantifies host+viral
# in one pass with proper cell-calling (EmptyDrops_CR). Complements the kb/ViralScan
# run: STARsolo counts unique reads only (loses ViralScan's multimap-recovered viral
# signal) but gives a correct barcode space + built-in cell calls — a cross-check and
# a source of an external called-cell list for `viralscan cell_calling=external`.
#
# Two-phase (submit index build first, then the run array with an afterok dependency):
#   sbatch covid_viralscan/scripts/slurm_starsolo_covid.sh build      # STEP 1: index
#   sbatch --array=0-1 --dependency=afterok:<JOB> \
#          covid_viralscan/scripts/slurm_starsolo_covid.sh            # STEP 2: STARsolo
#
#SBATCH -J covid_starsolo
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=12:00:00
#SBATCH -o covid_viralscan/logs/starsolo_%x_%A_%a.log
#SBATCH -e covid_viralscan/logs/starsolo_%x_%A_%a.err

set -euo pipefail

REPO=/exports/para-lipg-hpc/mdmanurung/ViralScan
STAR=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/starsolo/bin/STAR

# ── inputs ────────────────────────────────────────────────────────────────────
GRCH38_FA=/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa
GRCH38_GTF=/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf
PANEL_FA=$REPO/references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.fa
PANEL_GTF=$REPO/references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf
SARS_FA=$REPO/covid_viralscan/viralscan_ref/viral/reference.fasta
SARS_GTF=$REPO/covid_viralscan/viralscan_ref/viral/viral_whole_genome.gtf

# GEM-X-5' barcode whitelist (CellRanger-derived; the bundled 10x lists don't match) + geometry
WHITELIST=$REPO/covid_viralscan/viralscan_ref/cellranger_whitelist.txt
FASTQ_DATA=$REPO/covid_viralscan/data

OUT=/exports/para-lipg-hpc/mdmanurung/covid_viralscan_starsolo
GENOME_DIR=$OUT/genome_GRCh38_viral
COMBINED_FA=$OUT/GRCh38_viral_combined.fa
COMBINED_GTF=$OUT/GRCh38_viral_combined.gtf
mkdir -p "$OUT" "$REPO/covid_viralscan/logs"

# ── STEP 1: build combined index ("build" arg) ───────────────────────────────
if [[ "${1:-}" == "build" ]]; then
    if [[ -f "$GENOME_DIR/SAindex" ]]; then
        echo "[$(date)] SKIP build — $GENOME_DIR/SAindex exists"; exit 0
    fi
    for f in "$GRCH38_FA" "$GRCH38_GTF" "$PANEL_FA" "$PANEL_GTF" "$SARS_FA" "$SARS_GTF"; do
        [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
    done
    echo "[$(date)] Concatenating GRCh38 + viral panel + SARS-CoV-2 ..."
    cat "$GRCH38_FA" "$PANEL_FA" "$SARS_FA" > "$COMBINED_FA"
    cat "$GRCH38_GTF" "$PANEL_GTF" "$SARS_GTF" > "$COMBINED_GTF"
    mkdir -p "$GENOME_DIR"
    echo "[$(date)] STAR genomeGenerate (16 threads) ..."
    "$STAR" --runMode genomeGenerate --runThreadN 16 \
        --genomeDir "$GENOME_DIR" \
        --genomeFastaFiles "$COMBINED_FA" \
        --sjdbGTFfile "$COMBINED_GTF" \
        --sjdbOverhang 89 \
        --genomeSAsparseD 3
    echo "[$(date)] Index build complete: $GENOME_DIR"
    exit 0
fi

# ── STEP 2: STARsolo per sample (array 0-1) ──────────────────────────────────
SAMPLES=("LUM-SJ-x213-g" "LUM-SJ-x216-g")
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID:-0}]}
[[ -f "$GENOME_DIR/SAindex" ]] || { echo "ERROR: index not built — run '$0 build' first" >&2; exit 1; }
[[ -f "$WHITELIST" ]] || { echo "ERROR: missing whitelist $WHITELIST" >&2; exit 1; }

R1=$FASTQ_DATA/$SAMPLE/${SAMPLE}_merged_R1.fastq.gz
R2=$FASTQ_DATA/$SAMPLE/${SAMPLE}_merged_R2.fastq.gz
for f in "$R1" "$R2"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done

SOUT=$OUT/$SAMPLE/
mkdir -p "$SOUT"
echo "[$(date)] STARsolo $SAMPLE (GeneFull, EmptyDrops_CR, GEM-X 16+12) ..."
# 10x geometry: R2 = cDNA (genomic read), R1 = CB(16)+UMI(12). soloBarcodeReadLength 0
# disables the R1-length check (GEM-X R1 can differ from 26/28).
"$STAR" --runMode alignReads --runThreadN 16 \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$R2" "$R1" \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloFeatures GeneFull \
    --soloCellFilter EmptyDrops_CR \
    --soloStrand Forward \
    --outSAMtype None \
    --outFileNamePrefix "$SOUT"
echo "[$(date)] STARsolo $SAMPLE done. Filtered cells: $SOUT/Solo.out/GeneFull/filtered/"
echo "  Use barcodes.tsv there as the external called-cell list for viralscan cell_calling=external."
