#!/usr/bin/env bash
# Build a genome-discriminated ViralScan panel reference (F-005 fix).
#
# The standard cDNA-only host reference cannot suppress reads from GRCh38
# intronic / intergenic regions that share k-mers with anellovirus or other
# viral sequences, producing a ~90% spurious "viral" signal in whole-blood
# RNA-seq (finding F-005, 2026-07-06).
#
# This script adds the GRCh38 primary-assembly genome FASTA as a kallisto
# D-list.  k-mers shared between the genome and any viral sequence are masked
# in the index, so reads whose k-mers appear anywhere in GRCh38 (including
# introns and intergenic regions) are NOT counted as viral.
#
# Estimated resource requirements: 64 GB RAM, 8 CPUs, ~8 h wall time.
# Output: /exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref_genomic/ref/
#
# Usage (from repo root):
#   sbatch scripts/build_genome_panel_ref.sh

#SBATCH --job-name=vs_panel_genomic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --chdir=/exports/para-lipg-hpc/mdmanurung/ViralScan
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref_genomic/build_%j.out
#SBATCH --error=/exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref_genomic/build_%j.err

set -euo pipefail

# ── paths ─────────────────────────────────────────────────────────────────────
REPO=/exports/para-lipg-hpc/mdmanurung/ViralScan

# GRCh38 2024-A primary assembly from the CellRanger reference package.
# Contains all chromosomal sequences including introns and intergenic regions.
GENOME_FA=/exports/para-lipg-hpc/mdmanurung/malaria_bcells/data/refgenome/refdata-gex-GRCh38-2024-A/fasta/genome.fa

OUTDIR=/exports/para-lipg-hpc/mdmanurung/viralscan_panel_ref_genomic/ref

# ── conda env ─────────────────────────────────────────────────────────────────
CONDA_INIT=/share/software/tools/miniconda/3.10/23.3.1/etc/profile.d/conda.sh
CONDA_ENV=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
if [[ -f $CONDA_INIT ]]; then
    # shellcheck source=/dev/null
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV"
fi
command -v kb || { echo "ERROR: kb not on PATH — activate test_viralscan conda env"; exit 1; }

# ── preflight ─────────────────────────────────────────────────────────────────
[[ -f $GENOME_FA ]] || { echo "ERROR: genome FASTA not found: $GENOME_FA"; exit 1; }
mkdir -p "$OUTDIR"

echo "[$(date)] Building genome-discriminated panel reference"
echo "  Genome D-list : $GENOME_FA ($(du -sh "$GENOME_FA" | cut -f1))"
echo "  Output dir    : $OUTDIR"

# ── build ─────────────────────────────────────────────────────────────────────
PYTHONPATH="$REPO/src" NCBI_EMAIL="${NCBI_EMAIL:-mikhael.manurung@gmail.com}" \
    python "$REPO/scripts/build_bundled_panel_ref.py" \
    --out "$OUTDIR" \
    --genome-dlist "$GENOME_FA"

echo "[$(date)] Build complete."
echo ""
echo "To use this index for the B5 bulk scan:"
echo "  REFDIR=$OUTDIR sbatch --array=0-5 scripts/bulk_viral_scan.sh   # pilot first"
echo "  REFDIR=$OUTDIR sbatch --array=0-98 scripts/bulk_viral_scan.sh  # full cohort"
