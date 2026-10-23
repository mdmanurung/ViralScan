#!/usr/bin/env bash
# Spike-in SMOKE test (exploratory; PLAN SPIKE-01): 1 background (D03 SRR29320037, 5' 10xv2, reverse), levels 10 and 100,
# seed 1, 5 cells. Checks the spike machinery recovers in kb (old + flank), STARsolo host+virus, and the real STAR host filter.
# Order: spike-only arms first (minutes, so a failure is cheap), then the full-background runs.
# Background = FIRST 5M pairs (plumbing only; the full grid uses ~45M pairs via `bg --frac --seed` and comma-separated
# inputs, not cat'd copies).
# usage: sbatch smoke.sbatch   (no job is submitted by the commit that adds this file)
#SBATCH -J spike_smoke
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -t 08:00:00
#SBATCH -p highmem
#SBATCH -o slurm_%x_%j.log  # relative to the submit directory
set -euo pipefail
M=${VIRALSCAN_SITE_ROOT:?set VIRALSCAN_SITE_ROOT to the site work root}
W=$M/viralscan_work; X=$W/trackB_explore_2026-10; S=$W/star_hostvirus
SRC=$M/vs_work/spikein; A=$SRC/analysis/spikein_detection_limit
O=$W/spikein_smoke_2026-10; LIB=SRR29320037; HOST=$M/ViralScan/ebv_latest_ref_2026-09-27/star_host_2.7.11b
CONDA=$M/conda/envs/viralscan_test_full/bin; KB=$M/viralscan_env_pinned/bin
BUSTOOLS=$M/viralscan_env_pinned/lib/python3.8/site-packages/kb_python/bins/linux/bustools/bustools
export PATH=$CONDA:$PATH PYTHONPATH=$SRC/src
mkdir -p $O/runs $O/scores; cd $O
# 10xv2 barcode on-list (kb-python's copy); spike barcodes must be on it or STARsolo/kb drop them
ONLIST=$(python -c 'from viralscan.chemistry import onlist_path; print(onlist_path("10xv2"))')
date; STAR --version

[ -s bg_R2.fastq.gz ] || python -I $A/spikein.py bg --r1 $X/inputs/$LIB/${LIB}_R1.fastq.gz --r2 $X/inputs/$LIB/${LIB}_R2.fastq.gz --n-pairs 5000000 --out-prefix $O/bg
python -I $A/spikein.py spike --bam $S/GSM6416022/Aligned.sortedByCoord.out.bam --gtf $S/ref/hostvirus.gtf --bg-r1 $O/bg_R1.fastq.gz \
  --whitelist $ONLIST --out-dir $O/spikes --levels 10,100 --seeds 1 --cells 5 --cb-len 16 --umi-len 10 --r1-len 26 --r2-len 90 --bg-orient reverse

# --- spike-only runs (first, minutes): (a) the real STAR host filter (no viral competitor), (b) host+virus index.
# Both get a plain-text whitelist of the spike barcodes: STAR 2.7.11b exits 139 in Solo with --soloCBwhitelist None
# when zero reads map, and HHV-6B reads vs a host-only index may well map nothing (ANDET-09-STAR139).
# Compare placements with the full runs below.
for TAG in L10_s1 L100_s1; do
  D=$O/runs/spikeonly_$TAG; mkdir -p $D/hf
  tail -n +2 spikes/spike_$TAG.truth.tsv | cut -f2 | sort -u > $D/spike_cbs.txt
  python - <<PYEOF
from pathlib import Path
from viralscan.scripts.host_filter import _starsolo_filter
od = Path("$D/hf")
_starsolo_filter("$O/spikes/spike_${TAG}_R1.fastq.gz", "$O/spikes/spike_${TAG}_R2.fastq.gz", "$HOST", "10xv2", "$D/spike_cbs.txt",
                 od, str(od / "R1.fastq.gz"), str(od / "R2.fastq.gz"), 8)
PYEOF
  python -I $A/score.py hostfilter --filtered-r2 $D/hf/R2.fastq.gz --truth spikes/spike_$TAG.truth.tsv > scores/hostfilter_$TAG.tsv
  # accounting (informational, never fails the job): input = mapped + unmapped; survivors = reads in the filtered FASTQ
  LF=$(find $D/hf -name Log.final.out | head -1 || true)
  if [ -n "$LF" ]; then
    echo "hostfilter $TAG accounting:" \
      "input=$(awk -F'\t' '/Number of input reads/{print $2}' $LF)" \
      "uniquely_mapped=$(awk -F'\t' '/Uniquely mapped reads number/{print $2}' $LF)" \
      "multi_mapped=$(awk -F'\t' '/Number of reads mapped to multiple loci/{print $2}' $LF)" \
      "survivors=$(( $(zcat $D/hf/R2.fastq.gz | wc -l) / 4 ))"
  else
    echo "hostfilter $TAG accounting: Log.final.out not found under $D/hf (STAR log location unverified)"
  fi
  python -m viralscan.star_hostvirus align --r1 spikes/spike_${TAG}_R1.fastq.gz --r2 spikes/spike_${TAG}_R2.fastq.gz --index $S/index --out-dir $D/star --technology 10xv2 --strand reverse --whitelist $D/spike_cbs.txt --threads 8
  python -I $A/score.py star --bam $D/star/Aligned.sortedByCoord.out.bam --truth spikes/spike_$TAG.truth.tsv > scores/star_spikeonly_$TAG.tsv
  rm -rf $D/star/_STARtmp $D/hf/star_tmp
done

# --- full-background runs: base (unspiked), L10, L100. Concatenated gzip is valid gzip; delete the copies when scored.
for TAG in base L10_s1 L100_s1; do
  D=$O/runs/$TAG; mkdir -p $D
  if [ $TAG = base ]; then cat bg_R1.fastq.gz > $D/in_R1.fastq.gz; cat bg_R2.fastq.gz > $D/in_R2.fastq.gz
  else cat bg_R1.fastq.gz spikes/spike_${TAG}_R1.fastq.gz > $D/in_R1.fastq.gz; cat bg_R2.fastq.gz spikes/spike_${TAG}_R2.fastq.gz > $D/in_R2.fastq.gz; fi
  python -m viralscan.star_hostvirus align --r1 $D/in_R1.fastq.gz --r2 $D/in_R2.fastq.gz --index $S/index --out-dir $D/star --technology 10xv2 --strand reverse --threads 8
  BAM=$D/star/Aligned.sortedByCoord.out.bam
  samtools idxstats $BAM > $D/star/idxstats.tsv
  [ $TAG = base ] || python -I $A/score.py star --bam $BAM --truth spikes/spike_$TAG.truth.tsv > scores/star_$TAG.tsv
  # keep VIRUS| contigs + SPIKE reads before deleting the full BAM, so the unspiked base run keeps its false-positive floor
  awk -F'\t' '$1 ~ /^VIRUS\|/ {print $1}' $D/star/idxstats.tsv > $D/star/virus_contigs.txt
  samtools view -b -o $D/star/virus_only.bam $BAM $(cat $D/star/virus_contigs.txt)
  samtools view -h $BAM | awk '/^@/ || $1 ~ /^SPIKE/' | samtools view -b -o $D/star/spike_reads.bam -
  rm -rf $BAM* $D/star/_STARtmp   # disk is at 93% (cap 95%)
  for IDX in old flank; do
    if [ $IDX = old ]; then I=$M/viral_ref_cat42d/build/panel.idx; T=$M/viral_ref_cat42d/build/panel.t2g
    else I=$X/ref_flank_v2/panel_flank.idx; T=$X/ref_flank_v2/panel_flank.t2g; fi
    PATH=$KB:$PATH kb count -i $I -g $T -x 10XV2 --parity single --strand unstranded -t 8 --overwrite --tmp $D/tmp_$IDX -o $D/kb_$IDX $D/in_R1.fastq.gz $D/in_R2.fastq.gz
    rm -rf $D/tmp_$IDX
    [ $TAG = base ] || python -I $A/kb_fate.py --kb-dir $D/kb_$IDX --t2g $T --truth spikes/spike_$TAG.truth.tsv --bustools $BUSTOOLS > scores/kbfate_${IDX}_$TAG.tsv
    [ $TAG = base ] || python -I $A/score.py kb --spiked-dir $D/kb_$IDX/counts_unfiltered --unspiked-dir runs/base/kb_$IDX/counts_unfiltered --truth spikes/spike_$TAG.truth.tsv > scores/kb_${IDX}_$TAG.tsv
  done
  rm -f $D/in_R1.fastq.gz $D/in_R2.fastq.gz
done
date; echo SPIKE_SMOKE_DONE
