#!/bin/bash
#SBATCH --job-name=evidence_rerun
#SBATCH --partition=all
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan/covid_viralscan/logs/evidence_rerun_%j.log
#SBATCH --error=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan/covid_viralscan/logs/evidence_rerun_%j.err

# Reproducibility fix for T5: regenerate viral_reads.bam + coverage.tsv from
# the committed pipeline, closing the gap left by the undocumented hf_align
# job 25181135.
#
# Context: the original `viralscan evidence` run (job 25180994) crashed at
# `samtools sort` because the viral reference (combined.fa) contained duplicate
# NC_002076.2 headers. viral_reads.fasta was already extracted before the crash.
# viral_genome.dedup.fa was created with the dedup fix (1 copy of NC_002076.2).
# This script re-aligns from the existing viral_reads.fasta to the deduped
# reference, producing BAM + coverage.tsv via the same code path as
# `viralscan evidence --viral-fasta`.
#
# Usage: sbatch scripts/slurm_evidence_rerun.sh

set -euo pipefail

# ── Conda ─────────────────────────────────────────────────────────────────────
__conda_setup="$('/share/software/tools/miniconda/3.10/23.3.1/bin/conda' 'shell.bash' 'hook' 2>/dev/null)"
eval "$__conda_setup"
conda activate evtools

# ── Paths ─────────────────────────────────────────────────────────────────────
REPODIR=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
VIRAL_FA=/exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/viralscan_ref/viral_genome.dedup.fa
THREADS=8

echo "=== Evidence rerun $(date) ==="
echo "Viral reference: ${VIRAL_FA}"
echo "Sequences in reference: $(grep -c '^>' ${VIRAL_FA})"
echo "NC_002076.2 copies: $(grep -c 'NC_002076.2' ${VIRAL_FA})"

for SAMPLE in x213-g x216-g; do
    EVDIR=${REPODIR}/covid_viralscan/results_hostfilter/LUM-SJ-${SAMPLE}/evidence
    READS_FA=${EVDIR}/viral_reads.fasta
    OUT_BAM=${EVDIR}/viral_reads.bam
    COV_TSV=${EVDIR}/coverage.tsv

    if [[ ! -s "${READS_FA}" ]]; then
        echo "  ${SAMPLE}: viral_reads.fasta missing or empty — skipping"
        continue
    fi

    NREADS=$(grep -c '^>' "${READS_FA}")
    echo ""
    echo "=== ${SAMPLE}: ${NREADS} reads in ${READS_FA} ==="

    # Align to deduped viral reference (minimap2 short-read mode)
    echo "  Aligning to ${VIRAL_FA} ..."
    minimap2 -ax sr -t ${THREADS} \
        "${VIRAL_FA}" "${READS_FA}" 2>/dev/null \
        | samtools sort -@ 4 -o "${OUT_BAM}" -
    samtools index "${OUT_BAM}"
    echo "  BAM written: ${OUT_BAM}"

    # Coverage table (replicates viralscan.evidence.coverage_table())
    echo "  Computing coverage ..."
    PYTHONPATH="${REPODIR}/src" python3 - <<PYEOF
import csv, sys
sys.path.insert(0, "${REPODIR}/src")
from viralscan.evidence import coverage_table
cov = coverage_table("${OUT_BAM}")
if not cov:
    print("  WARNING: 0 reads aligned", file=sys.stderr)
    sys.exit(0)
with open("${COV_TSV}", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(cov[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(cov)
print(f"  Wrote {len(cov)} covered references to ${COV_TSV}")
for r in cov[:5]:
    print(f"  {r.get('rname')}: reads={r.get('numreads')} "
          f"coverage={r.get('coverage')}% meandepth={r.get('meandepth')}")
PYEOF

done

echo ""
echo "=== Evidence rerun complete $(date) ==="
echo "Outputs:"
for SAMPLE in x213-g x216-g; do
    EVDIR=${REPODIR}/covid_viralscan/results_hostfilter/LUM-SJ-${SAMPLE}/evidence
    echo "  ${SAMPLE}: ${EVDIR}/viral_reads.bam + coverage.tsv"
done
