#!/usr/bin/env bash
# Stage 2 v2 — step B resume: kallisto index on deduplicated inputs.
#
# Background: job 25138052 extracted cdna.fa + t2g.txt correctly (2 h), then
# `kallisto index` failed with "repeated name in FASTA file" because the source
# panel viral_genome.fa contains the accession NC_002076.2 (Torque teno virus 1)
# TWICE (byte-identical). ngs_tools extracted its 4 gene/transcript records twice
# → 4 duplicate names in cdna.fa, which kallisto rejects.
#
# Fix applied before this job: keep-first dedup of cdna.fa (470472→470468),
# t2g.txt (→470468), and the D-list combined.fa (468083→468082). All removed
# records are byte-identical, so this equals a clean build from a deduped panel.
#
# This job just resumes the pipeline at the kallisto index step (the literal
# command from kb ref's own traceback) and swaps the deduped files into place.
#
#SBATCH -J covid_vs_kalidx
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH -o covid_viralscan/logs/kallisto_index_%j.log
#SBATCH -e covid_viralscan/logs/kallisto_index_%j.err

set -euo pipefail

REPO=/exports/para-lipg-hpc/mdmanurung/ViralScan
VS_CONDA_ENV=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
KALLISTO=$VS_CONDA_ENV/lib/python3.12/site-packages/kb_python/bins/linux/kallisto/kallisto
export PATH="$VS_CONDA_ENV/bin:$PATH"

OUTDIR=$REPO/covid_viralscan/viralscan_ref
INDEX=$OUTDIR/index.idx

echo "[$(date)] kallisto index resume on $(hostname)"
[[ -f "$OUTDIR/cdna.dedup.fa"     ]] || { echo "ERROR: missing cdna.dedup.fa" >&2; exit 1; }
[[ -f "$OUTDIR/combined.dedup.fa" ]] || { echo "ERROR: missing combined.dedup.fa" >&2; exit 1; }
[[ -f "$OUTDIR/t2g.dedup.txt"     ]] || { echo "ERROR: missing t2g.dedup.txt" >&2; exit 1; }
[[ -x "$KALLISTO" ]] || { echo "ERROR: kallisto binary not found at $KALLISTO" >&2; exit 1; }

rm -f "$INDEX"   # clear the 0-byte artifact from the failed run

echo "[$(date)] building index (k=31, t=8, D-list=combined.dedup.fa) ..."
"$KALLISTO" index \
    -i "$INDEX" \
    -k 31 \
    -t 8 \
    -d "$OUTDIR/combined.dedup.fa" \
    "$OUTDIR/cdna.dedup.fa"
echo "[$(date)] kallisto index finished (exit $?)."

# Verify by ARTIFACT, not exit code (previous run logged success with a 0-byte index)
if [[ ! -s "$INDEX" ]]; then
    echo "ERROR: index.idx is empty after build — aborting, files NOT swapped." >&2
    exit 1
fi
echo "  index.idx size: $(du -sh "$INDEX" | cut -f1)"
echo "[$(date)] kallisto inspect:"
"$KALLISTO" inspect "$INDEX" 2>&1 | sed 's/^/    /' || true

# Promote deduped files to canonical names for Stage 3
mv -f "$OUTDIR/cdna.dedup.fa"     "$OUTDIR/cdna.fa"
mv -f "$OUTDIR/combined.dedup.fa" "$OUTDIR/combined.fa"
mv -f "$OUTDIR/t2g.dedup.txt"     "$OUTDIR/t2g.txt"
echo "[$(date)] deduped files promoted to canonical names."

echo ""
echo "[$(date)] Post-build sanity checks on t2g.txt ..."
echo "  NC_045512 rows: $(grep -c NC_045512 "$OUTDIR/t2g.txt" || echo 0)  (SARS-CoV-2; expect >=1)"
echo "  ENST rows:      $(grep -c ENST      "$OUTDIR/t2g.txt" || echo 0)  (human; expect >=100k)"
echo "  _gene rows:     $(grep -c _gene     "$OUTDIR/t2g.txt" || echo 0)  (anello; expect >=1000)"
echo "  sarsp rows:     $(grep -c sarsp     "$OUTDIR/t2g.txt" || echo 0)  (SARS-CoV-1; expect >=1)"
echo "  NC_002076.2_tx1 rows: $(grep -cP '^NC_002076.2_tx1\t' "$OUTDIR/t2g.txt" || echo 0)  (dedup check; expect 1)"

echo ""
echo "[$(date)] Stage 2 v2 (index resume) complete."
echo "  Proceed to Stage 3: sbatch --array=0-1 covid_viralscan/scripts/slurm_viralscan_quant.sh"
