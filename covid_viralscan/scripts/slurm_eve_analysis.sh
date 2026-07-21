#!/bin/bash
#SBATCH --job-name=eve_analysis
#SBATCH --partition=all
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=eve_analysis_%j.log
#SBATCH --error=eve_analysis_%j.err

# EVE (Endogenous Viral Element) characterisation analysis
# Three-phase analysis to identify human genomic loci with anellovirus homology:
#   Phase A: reads → GRCh38 (minimap2) — which human loci do artifact reads come from?
#   Phase B: covered viral seqs → NT (BLAST, human-only) — known EVEs / gene annotation?
#   Phase C: full viral panel → GRCh38 (minimap2 asm20) — genome-wide EVE screen

set -euo pipefail

# ── Conda ────────────────────────────────────────────────────────────────────
CONDA_EXE=${CONDA_EXE:-conda}
EVE_CONDA_ENV=${EVE_CONDA_ENV:-evtools}
if ! __conda_setup="$("${CONDA_EXE}" shell.bash hook 2>/dev/null)"; then
    echo "ERROR: set CONDA_EXE to a usable conda executable" >&2
    exit 2
fi
eval "$__conda_setup"
conda activate "${EVE_CONDA_ENV}"

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPODIR=${REPODIR:-$(cd -- "${SCRIPT_DIR}/../.." && pwd)}
OUTDIR=${OUTDIR:-${REPODIR}/covid_viralscan/results_hostfilter/eve_analysis}
VIRAL_FA=${VIRAL_FA:-${REPODIR}/covid_viralscan/viralscan_ref/viral_genome.dedup.fa}
BAM_213=${BAM_213:-${REPODIR}/covid_viralscan/results_hostfilter/LUM-SJ-x213-g/evidence/viral_reads.bam}
BAM_216=${BAM_216:-${REPODIR}/covid_viralscan/results_hostfilter/LUM-SJ-x216-g/evidence/viral_reads.bam}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-8}}

: "${GRCh38_FA:?Set GRCh38_FA to the GRCh38 FASTA path}"
: "${GENES_GTF:?Set GENES_GTF to the GRCh38 genes.gtf path}"
: "${NT_DB:?Set NT_DB to the local BLAST nt database prefix}"

for tool in samtools minimap2 blastn python3; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: required tool not found on PATH: ${tool}" >&2
        exit 2
    fi
done

for input_file in "${VIRAL_FA}" "${GRCh38_FA}" "${GENES_GTF}" "${BAM_213}" "${BAM_216}"; do
    if [[ ! -s "${input_file}" ]]; then
        echo "ERROR: required input file missing or empty: ${input_file}" >&2
        exit 2
    fi
done

# High-depth accessions from both samples (meandepth > 5x post-host-filter)
# Includes anellovirus genera, EMCV (NC_001479.1), HCV (NC_004102.1), OC43 (NC_002645.1)
KEY_ACCS=(
    NC_001479.1   # EMCV — 841/1621x depth, 156 bases, identical in both samples
    MW455365.1    # Samektorquevirus — 32/8x depth, 63/59 bases
    MW455373.1    # Samektorquevirus — 48/15x depth, 57/51 bases
    MW455378.1    # Samektorquevirus — 40/9x depth, 60/53 bases
    MZ286238.1    # Gammatorquevirus — 70/16x depth, 80/45 bases
    KP343825.1    # Gammatorquevirus — 19/7x depth, 88/90 bases (3.41% breadth)
    MW455439.1    # Alphatorquevirus — 10x depth, 42/40 bases
    MN771265.1    # Betatorquevirus — 5/5x depth, 59/49 bases
    NC_004102.1   # HCV — 76/50x depth, 112/118 bases (unexpected)
    NC_002645.1   # OC43 coronavirus — 2.5/0.7x depth, 45/46 bases
    AB303557.1    # Anellovirus — 1.2x depth, 44/37 bases
)

mkdir -p "${REPODIR}/covid_viralscan/logs" "${OUTDIR}/"{reads,depth,blast,minimap2_panel,annotation}

echo "=== EVE Analysis $(date) ==="
echo "Samples: x213-g (${BAM_213}) x216-g (${BAM_216})"

# ── Phase A: Extract reads per accession → map to GRCh38 ─────────────────────
echo ""
echo "=== Phase A: reads → GRCh38 mapping ==="

for ACC in "${KEY_ACCS[@]}"; do
    echo "  Extracting reads for ${ACC}..."

    # Extract reads from both BAMs. Keep stderr (do not discard to /dev/null) and
    # surface a visible WARN on real failure (e.g. missing index) instead of a
    # silent `|| true`, which would be indistinguishable from "region absent".
    for pair in "x213:${BAM_213}" "x216:${BAM_216}"; do
        s="${pair%%:*}"; bam="${pair#*:}"
        errlog="${OUTDIR}/reads/${ACC}_${s}.samtools.err"
        if ! samtools view -u "${bam}" "${ACC}" 2>"${errlog}" \
                | samtools fasta - > "${OUTDIR}/reads/${ACC}_${s}.fa" 2>>"${errlog}"; then
            echo "    WARN: read extraction failed for ${ACC} (${s}); see ${errlog}" >&2
        fi
    done

    for SAMPLE in x213 x216; do
        READS=${OUTDIR}/reads/${ACC}_${SAMPLE}.fa
        if [[ ! -s "${READS}" ]]; then
            echo "    ${SAMPLE}: no reads extracted — skipping"
            continue
        fi
        NREADS=$(grep -c "^>" "${READS}" || true)
        echo "    ${SAMPLE}: ${NREADS} reads → mapping to GRCh38"

        minimap2 -ax sr -t "${THREADS}" \
            "${GRCh38_FA}" "${READS}" 2>/dev/null \
            | samtools sort -@ "${THREADS}" \
            | samtools view -F 4 -b \
            > "${OUTDIR}/reads/${ACC}_${SAMPLE}_grch38.bam"

        samtools index "${OUTDIR}/reads/${ACC}_${SAMPLE}_grch38.bam"

        # Summary: which chromosomes / regions are hit?
        samtools idxstats "${OUTDIR}/reads/${ACC}_${SAMPLE}_grch38.bam" \
            | awk '$3>0 {print}' \
            > "${OUTDIR}/reads/${ACC}_${SAMPLE}_grch38.idxstats.txt"

        # Per-base coverage on the human genome at these positions
        samtools depth -a "${OUTDIR}/reads/${ACC}_${SAMPLE}_grch38.bam" \
            | awk '$3>0' \
            > "${OUTDIR}/reads/${ACC}_${SAMPLE}_grch38.depth.txt"
    done
done
echo "Phase A complete."

# ── Phase B: Extract covered viral positions → BLAST vs NT (human only) ──────
echo ""
echo "=== Phase B: covered viral seqs → NT BLAST (human only) ==="

# Get per-base depth across the whole viral BAM for both samples
echo "  Computing per-base depth..."
samtools depth -a "${BAM_213}" > "${OUTDIR}/depth/per_base_213.tsv"
samtools depth -a "${BAM_216}" > "${OUTDIR}/depth/per_base_216.tsv"
echo "  Depth files written."

# Ensure viral FASTA is indexed
if [[ ! -f "${VIRAL_FA}.fai" ]]; then
    samtools faidx "${VIRAL_FA}"
fi

# Python: extract covered regions (depth ≥ 10 in either sample), write query FASTA
DEPTH_213="${OUTDIR}/depth/per_base_213.tsv"
DEPTH_216="${OUTDIR}/depth/per_base_216.tsv"
OUTFA="${OUTDIR}/blast/covered_viral_regions.fa"
export VIRAL_FA DEPTH_213 DEPTH_216 OUTFA
python3 - <<'PYEOF'
import collections
import os
import sys

VIRAL_FA = os.environ["VIRAL_FA"]
DEPTH_213 = os.environ["DEPTH_213"]
DEPTH_216 = os.environ["DEPTH_216"]
OUTFA = os.environ["OUTFA"]
DEPTH_THRESH = 10

# Accumulate covered positions per accession
covered = collections.defaultdict(set)
for fpath in [DEPTH_213, DEPTH_216]:
    with open(fpath) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            acc, pos, depth = parts[0], int(parts[1]), int(parts[2])
            if depth >= DEPTH_THRESH:
                covered[acc].add(pos)

# For each accession, merge positions into contiguous regions (±5 bp tolerance)
def merge_positions(positions, slop=5):
    if not positions:
        return []
    spos = sorted(positions)
    regions = []
    start = end = spos[0]
    for p in spos[1:]:
        if p <= end + slop:
            end = p
        else:
            regions.append((start, end))
            start = end = p
    regions.append((start, end))
    return regions

# Load viral FASTA
seqs = {}
with open(VIRAL_FA) as fh:
    cur_id = None
    buf = []
    for line in fh:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if cur_id:
                seqs[cur_id] = ''.join(buf)
            cur_id = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if cur_id:
        seqs[cur_id] = ''.join(buf)

# Write query FASTA for BLAST
written = 0
with open(OUTFA, 'w') as out:
    for acc, positions in sorted(covered.items()):
        if acc not in seqs:
            print(f"  WARNING: {acc} not in viral FASTA — skipping", file=sys.stderr)
            continue
        seq = seqs[acc]
        regions = merge_positions(positions)
        for start, end in regions:
            region_len = end - start + 1
            if region_len < 20:   # skip tiny fragments
                continue
            # samtools depth positions are 1-based; Python slices are 0-based/end-exclusive.
            s0 = max(0, start - 1 - 20)
            e0 = min(len(seq), end + 20)
            subseq = seq[s0:e0]
            header = f">{acc}:{s0 + 1}-{e0} len={len(subseq)} depth_pos={start}-{end} n_pos={len(positions)}"
            out.write(f"{header}\n{subseq}\n")
            written += 1

print(f"  Wrote {written} covered regions from {len(covered)} accessions to {OUTFA}")
PYEOF

echo "  Running BLAST against NT (human only, taxid 9606)..."
if [[ -s "${OUTDIR}/blast/covered_viral_regions.fa" ]]; then
    blastn \
        -db "${NT_DB}" \
        -query "${OUTDIR}/blast/covered_viral_regions.fa" \
        -taxids 9606 \
        -outfmt "6 qseqid sseqid stitle pident length qlen qstart qend sstart send evalue bitscore" \
        -perc_identity 80 \
        -evalue 1e-5 \
        -max_target_seqs 10 \
        -num_threads "${THREADS}" \
        -out "${OUTDIR}/blast/covered_vs_nt_human.tsv" \
        > "${OUTDIR}/blast/blastn.stdout.log" \
        2> "${OUTDIR}/blast/blastn.stderr.log"
else
    : > "${OUTDIR}/blast/covered_vs_nt_human.tsv"
    echo "  No covered viral regions passed the BLAST query threshold; skipping BLAST."
fi

NHITS=$(wc -l < "${OUTDIR}/blast/covered_vs_nt_human.tsv")
echo "  Phase B complete: ${NHITS} BLAST hits."

# ── Phase C: Full viral panel → GRCh38 (minimap2 asm20) ──────────────────────
echo ""
echo "=== Phase C: full viral panel → GRCh38 (EVE screen) ==="

# asm20 mode: up to ~20% divergence, appropriate for ancient EVE integrations
minimap2 -x asm20 -t "${THREADS}" \
    --cs \
    -o "${OUTDIR}/minimap2_panel/panel_vs_grch38.paf" \
    "${GRCh38_FA}" \
    "${VIRAL_FA}" \
    2>"${OUTDIR}/minimap2_panel/minimap2.log"

NHITS_C=$(wc -l < "${OUTDIR}/minimap2_panel/panel_vs_grch38.paf")
echo "  Phase C complete: ${NHITS_C} minimap2 alignments."

# Quick summary: which viral accessions have GRCh38 hits?
awk '$11>=1000 && $12>=10 {print $1}' "${OUTDIR}/minimap2_panel/panel_vs_grch38.paf" \
    | sort | uniq -c | sort -rn \
    > "${OUTDIR}/minimap2_panel/accessions_with_grch38_hits.txt"
N_PANEL_HITS=$(wc -l < "${OUTDIR}/minimap2_panel/accessions_with_grch38_hits.txt")
echo "  Accessions with >=1000 bp, MAPQ>=10 hits: ${N_PANEL_HITS}"

# ── Phase D: Annotate results with GRCh38 gene GTF ───────────────────────────
echo ""
echo "=== Phase D: annotation ==="

python3 "${REPODIR}/covid_viralscan/scripts/annotate_eve.py" \
    --phase-a-dir "${OUTDIR}/reads" \
    --phase-b-blast "${OUTDIR}/blast/covered_vs_nt_human.tsv" \
    --phase-c-paf "${OUTDIR}/minimap2_panel/panel_vs_grch38.paf" \
    --gtf "${GENES_GTF}" \
    --outdir "${OUTDIR}/annotation" \
    --key-accs "${KEY_ACCS[*]}"

echo ""
echo "=== EVE analysis complete $(date) ==="
echo "Output: ${OUTDIR}"
echo ""
echo "Key files:"
echo "  Phase A: ${OUTDIR}/reads/{ACC}_{sample}_grch38.depth.txt"
echo "  Phase B: ${OUTDIR}/blast/covered_vs_nt_human.tsv"
echo "  Phase C: ${OUTDIR}/minimap2_panel/panel_vs_grch38.paf"
echo "  Summary: ${OUTDIR}/annotation/eve_summary_report.txt"
