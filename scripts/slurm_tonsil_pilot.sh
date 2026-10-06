#!/usr/bin/env bash
# TONSIL-02 pilot: SFL tonsil x223 (10x 5' v3, R1 28 bp, R2 90 bp, 1.6 B pairs) on
# the first 4 M pairs. Settles the inputs of the full runs before spending them:
#   1. `viralscan check-whitelist` against cellranger's raw barcode universe;
#   2. native `viralscan` (combined, cat42d, emptydrops) under forward / reverse /
#      unstranded, with per-arm wall time and peak RSS (/usr/bin/time -v).
# Same flag set as viral_ref_cat42d/covid_x213.sbatch (same lab, same chemistry).
# Code: detached worktree at the round's pinned commit, with the ignored GTFs copied in.
#
#   sbatch scripts/slurm_tonsil_pilot.sh
#
#SBATCH -J tonsil_pilot
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/sfl_tonsil/logs/%x_%j.log
set -euo pipefail
CODE=/exports/archive/hg-funcgenom-research/mdmanurung/vs_pinned/4346dc8
ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_bench
R_BIN=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin
REF=/exports/archive/hg-funcgenom-research/mdmanurung/viral_ref_cat42d/build
W=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/sfl_tonsil
WL=$W/ref/cellranger_raw_whitelist.txt
FQ=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan/benchmark_runs/sfl_tonsil_screen_2026-09-26/2025-3178-LUM-SJ-x223-x226/raw_fastq/LUM-SJ-x223_fastq
export PATH="$ENV/bin:$PATH:$R_BIN" PYTHONPATH="$CODE/src" MPLBACKEND=Agg
echo "code_sha $(git -C "$CODE" rev-parse HEAD)"

P=$W/pilot; mkdir -p "$P"; cd "$P"
if [ ! -s x223_R2.fastq.gz ]; then
  set +o pipefail  # zcat gets SIGPIPE when head stops reading
  zcat "$FQ/LUM-SJ-x223_S4_L004_R1_001.fastq.gz" | head -n 16000000 | gzip -1 > x223_R1.fastq.gz
  zcat "$FQ/LUM-SJ-x223_S4_L004_R2_001.fastq.gz" | head -n 16000000 | gzip -1 > x223_R2.fastq.gz
  set -o pipefail
fi
[ "$(zcat x223_R2.fastq.gz | wc -l)" -eq 16000000 ] || { echo "subsample short" >&2; exit 1; }

python -m viralscan.menu check-whitelist -s1 x223_R1.fastq.gz -w "$WL" -x 10xv3 \
  --n-sample 1000000 2>&1 | tee check_whitelist.txt

for st in forward reverse unstranded; do
  OUT=$P/$st
  [ -s "$OUT/x223/results/viral_summary.tsv" ] && continue
  rm -rf "$OUT"
  /usr/bin/time -v -o "time_$st.txt" python -m viralscan.menu -o "$OUT" \
    -s1 x223_R1.fastq.gz -s2 x223_R2.fastq.gz \
    -i "$REF/panel.idx" -t "$REF/panel.t2g" -gtf "$REF/viral_panel.gtf" \
    -w "$WL" -x 10xv3 -c 4 --strand "$st" \
    --cell-calling emptydrops --anellovirus-gene-ids --yes --verbose
done
grep -H "Elapsed\|Maximum resident" time_*.txt
echo PILOT_OK
