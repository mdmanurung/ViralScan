#!/usr/bin/env bash
# Build STAR genomeGenerate directories for the reference-strategy benchmark.
#
# Submit examples:
#   sbatch --parsable --export=ALL,REF_KIND=human scripts/slurm_build_starsolo_reference_strategy_refs.sh
#   sbatch --parsable --export=ALL,REF_KIND=all_virus scripts/slurm_build_starsolo_reference_strategy_refs.sh
#   sbatch --parsable --export=ALL,REF_KIND=combined scripts/slurm_build_starsolo_reference_strategy_refs.sh

#SBATCH --job-name=build_ref_strategy_star
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_runs/reference_strategy_2026-06-27_provenance_gate/logs/%x_%j.out
#SBATCH --error=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_runs/reference_strategy_2026-06-27_provenance_gate/logs/%x_%j.err

set -euo pipefail

ROOT=/exports/para-lipg-hpc/mdmanurung/ViralScan
STAR_BIN=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/starsolo/bin/STAR
THREADS=${SLURM_CPUS_PER_TASK:-8}
REF_KIND=${REF_KIND:?set REF_KIND to all_virus or combined}

mkdir -p "$ROOT/benchmark_runs/reference_strategy_2026-06-27_provenance_gate/logs"

case "$REF_KIND" in
  human)
    GENOME_DIR=$ROOT/references/starsolo/human_GRCh38_2024A
    FASTA=/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa
    GTF=/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf
    SA_NBASES=14
    ;;
  all_virus)
    GENOME_DIR=$ROOT/references/starsolo/all_virus_serratus_plus_anellovirus
    FASTA=$GENOME_DIR/viral_genome.fa
    GTF=$GENOME_DIR/viral_genome.gtf
    SA_NBASES=10
    ;;
  combined)
    GENOME_DIR=$ROOT/references/starsolo/combined_GRCh38_2024A_serratus_plus_anellovirus
    FASTA=$GENOME_DIR/combined_genome.fa
    GTF=$GENOME_DIR/combined.gtf
    SA_NBASES=14
    ;;
  *)
    echo "Unsupported REF_KIND=$REF_KIND" >&2
    exit 2
    ;;
esac

mkdir -p "$GENOME_DIR"
test -x "$STAR_BIN"
test -s "$FASTA"
test -s "$GTF"

"$STAR_BIN" \
  --runMode genomeGenerate \
  --genomeDir "$GENOME_DIR" \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 99 \
  --genomeSAindexNbases "$SA_NBASES" \
  --runThreadN "$THREADS"
