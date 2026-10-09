#!/usr/bin/env bash
# Rebuild the pure-human (GRCh38 2024-A) STAR genome index under the
# viralscan_bench env's STAR 2.7.11b. The existing index at
# evonk/.../refdata-gex-GRCh38-2024-A/star was built with STAR 2.7.1a, which
# STAR 2.7.11b refuses to load ("Genome version: 2.7.1a is INCOMPATIBLE with
# running STAR version: 2.7.11b"). Same source FASTA/GTF, same build
# parameters (read from the original genomeParameters.txt), just regenerated
# under the current STAR binary. Written to our own scratch dir since we
# don't own the evonk reference directory.
#
# Submit:
#   sbatch scripts/slurm_star_host_index_rebuild.sh
#
#SBATCH -J star_host_rebuild
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH -o logs/star_host_rebuild_%j.log
#SBATCH -e logs/star_host_rebuild_%j.err

set -euo pipefail

ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
export PATH="$ENV/bin:$PATH"

HR=/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A
OUT=/exports/para-lipg-hpc/mdmanurung/ViralScan/ebv_latest_ref_2026-09-27/star_host_2.7.11b
mkdir -p "$OUT"

echo "STAR version:"
STAR --version

STAR --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir "$OUT" \
  --genomeFastaFiles "$HR/fasta/genome.fa" \
  --sjdbGTFfile "$HR/genes/genes.gtf" \
  --genomeSAindexNbases 14 \
  --genomeChrBinNbits 18 \
  --genomeSAsparseD 3 \
  --sjdbOverhang 100 \
  --limitGenomeGenerateRAM 60000000000

echo "Done. Index at $OUT"
ls -la "$OUT"
