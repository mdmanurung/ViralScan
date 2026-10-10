#!/usr/bin/env bash
# DSR round 1: concatenate multi-lane / multi-run FASTQs per GSM into dsr_round1/inputs/.
# Lanes are NextSeq splits of one library (GSE190558) or replicate runs of one GSM (GSE154900);
# they must be one input per sample so UMIs are not double counted. Order is by SRR on both mates.
# Concatenated gzip members are valid gzip. Map files by read length, never by suffix:
#   GSE190558: _1 = I1 8 bp (dropped), _2 = R1 28 bp, _3 = cDNA 55 bp
#   GSE154900: _1 / _2 as downloaded (28/89 for SRR..51/52; 101/101 for SRR..53/54, flagged)
#
#   sbatch scripts/slurm_dsr_concat_inputs.sh      # array over the 6 GSMs
#
#SBATCH -J dsr_concat
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH -o /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round1/logs/%x_%A_%a.log
set -euo pipefail
D=/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan/benchmark_inputs/dsr_2026-10-05
O=/exports/archive/hg-funcgenom-research/mdmanurung/viralscan_work/dsr_round1/inputs
# name | study | mate-suffix for R1 | mate-suffix for R2 | SRRs (lane order)
ROWS=(
"GSM5725698|GSE190558|2|3|SRR17180386 SRR17180387 SRR17180388 SRR17180389"
"GSM5725697|GSE190558|2|3|SRR17180390 SRR17180391 SRR17180392 SRR17180393"
"GSM5725696|GSE190558|2|3|SRR17180394 SRR17180395 SRR17180396 SRR17180397"
"GSM5725695|GSE190558|2|3|SRR17180398 SRR17180399 SRR17180400 SRR17180401"
"GSM4682311|GSE154900|1|2|SRR12287751 SRR12287752"
"GSM4682312|GSE154900|1|2|SRR12287753 SRR12287754"
)
IFS='|' read -r NAME STUDY M1 M2 SRRS <<< "${ROWS[$SLURM_ARRAY_TASK_ID]}"
mkdir -p "$O/$NAME"
for pair in "R1:$M1" "R2:$M2"; do
  mate=${pair%%:*}; suf=${pair##*:}
  out="$O/$NAME/${NAME}_${mate}.fastq.gz"
  [ -s "$out" ] && continue
  files=(); for s in $SRRS; do files+=("$D/$STUDY/$s/${s}_${suf}.fastq.gz"); done
  cat "${files[@]}" > "$out.tmp"
  gzip -t "$out.tmp"
  mv "$out.tmp" "$out"
done
for m in R1 R2; do echo "$NAME $m $(zcat "$O/$NAME/${NAME}_$m.fastq.gz" | awk 'END{print NR/4}') reads"; done
