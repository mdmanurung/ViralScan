#!/usr/bin/env bash
# Diagnostic + salvage: reproduce the bustools sort→correct→sort→count chain on the
# existing output.bus for covid quant task 0 (LUM-SJ-x213-g).
#
# Why: the Snakefile kb_count rule captures `kb count ... 2>&1` into a shell variable
# that is discarded, so kb's real failure never reached disk — only the downstream
# `mv counts_unfiltered/` error did. kb produced output.bus (2.47 GB) + inspect.json but
# NO counts_unfiltered/, so a bustools step after `inspect` failed silently. This script
# re-runs that chain with stderr → file to expose the error, and (if it succeeds) salvages
# the ~40-min pseudoalignment by producing counts_unfiltered/ directly.
#
#SBATCH -J covid_diag_count
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH -o covid_viralscan/logs/diag_count_%j.log
#SBATCH -e covid_viralscan/logs/diag_count_%j.err

set -uo pipefail   # NOT -e: we want to run every step and see which one fails

VS=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
BT=$VS/lib/python3.12/site-packages/kb_python/bins/linux/bustools/bustools
WL=$VS/lib/python3.12/site-packages/ngs_tools/chemistry/whitelists/10x_version3_whitelist.txt.gz
T2G=/exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/viralscan_ref/t2g.txt

KB=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan/covid_viralscan/results/LUM-SJ-x213-g/kb-python
TMP=$KB/diag_tmp
mkdir -p "$TMP"

echo "[$(date)] Diagnostic bustools chain on $(hostname)"
echo "  bus: $(du -sh "$KB/output.bus" | cut -f1)"
echo "  free on archive: $(df -h "$KB" | tail -1 | awk '{print $4}')"

echo "[$(date)] STEP 1: bustools sort ..."
$BT sort -t 8 -T "$TMP/s1" -m 48G -o "$TMP/sorted.bus" "$KB/output.bus"
echo "  rc=$?  → $(ls -la "$TMP/sorted.bus" 2>/dev/null | awk '{print $5}') bytes"

echo "[$(date)] STEP 2: bustools correct (whitelist) ..."
$BT correct -w "$WL" -o "$TMP/corrected.bus" "$TMP/sorted.bus"
echo "  rc=$?  → $(ls -la "$TMP/corrected.bus" 2>/dev/null | awk '{print $5}') bytes"

echo "[$(date)] STEP 3: bustools sort (post-correct) ..."
$BT sort -t 8 -T "$TMP/s2" -m 48G -o "$TMP/corrected.sorted.bus" "$TMP/corrected.bus"
echo "  rc=$?"

echo "[$(date)] STEP 4: bustools count (--genecounts) ..."
mkdir -p "$TMP/counts_unfiltered"
$BT count -o "$TMP/counts_unfiltered/cells_x_genes" -g "$T2G" \
    -e "$KB/matrix.ec" -t "$KB/transcripts.txt" --genecounts "$TMP/corrected.sorted.bus"
echo "  rc=$?"

echo "[$(date)] RESULT:"
ls -la "$TMP/counts_unfiltered/" 2>/dev/null || echo "  counts_unfiltered NOT produced"
echo "[$(date)] Done."
