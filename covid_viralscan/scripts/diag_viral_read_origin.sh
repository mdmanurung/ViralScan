#!/usr/bin/env bash
# Decisive test: are ViralScan's anellovirus reads GENUINELY viral, or host homology?
# STARsolo combined-ref gave 0 viral UMI, but that's a GTF artifact (anello genes lack
# exon records → GeneFull can't count them). This test sidesteps counting: align an x213
# R2 subsample to the combined GRCh38+viral STAR genome and inspect where reads that hit
# viral contigs go — uniquely to viral (NH==1 → real) or multimapping to host too
# (NH>1 → homology / possible ViralScan cDNA-only-host over-call).
#
#SBATCH -J covid_readorigin
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH -o covid_viralscan/logs/readorigin_%j.log
#SBATCH -e covid_viralscan/logs/readorigin_%j.err

set -euo pipefail
STAR=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/starsolo/bin/STAR
SAMTOOLS=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/starsolo/bin/samtools
SS=/exports/para-lipg-hpc/mdmanurung/covid_viralscan_starsolo
GENOME=$SS/genome_GRCh38_viral
R2=/exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/data/LUM-SJ-x213-g/LUM-SJ-x213-g_merged_R2.fastq.gz
WORK=$SS/readorigin
mkdir -p "$WORK"

echo "[$(date)] subsampling 5M R2 reads ..."
zcat "$R2" | head -n 20000000 | gzip > "$WORK/sub_R2.fq.gz"   # 5M reads

echo "[$(date)] aligning subsample to combined GRCh38+viral (multimappers reported) ..."
"$STAR" --runMode alignReads --runThreadN 12 \
    --genomeDir "$GENOME" \
    --readFilesIn "$WORK/sub_R2.fq.gz" --readFilesCommand zcat \
    --outFilterMultimapNmax 50 \
    --outSAMtype BAM Unsorted \
    --outSAMattributes NH HI AS nM \
    --outFileNamePrefix "$WORK/sub_"

BAM=$WORK/sub_Aligned.out.bam
echo "[$(date)] analysing viral-contig read origins ..."
# viral contig set
awk '{print $1}' "$GENOME/chrNameLength.txt" | grep -vE '^(chr)?([0-9]+|X|Y|MT?)$' | grep -vE '^(KI|GL|chrUn|chrM)' > "$WORK/viral_contigs.txt"

# primary alignments only (exclude secondary 0x100, supplementary 0x800, unmapped 0x4)
echo "  total primary-aligned reads: $($SAMTOOLS view -c -F 0x904 "$BAM")"
# reads whose primary alignment is on a viral contig
$SAMTOOLS view -F 0x904 "$BAM" | awk -v vf="$WORK/viral_contigs.txt" '
  BEGIN{while((getline c < vf)>0) v[c]=1}
  { rname=$3; nh=1; for(i=12;i<=NF;i++) if($i ~ /^NH:i:/){split($i,a,":"); nh=a[3]}
    if(rname in v){ vt++; if(nh==1) vu++; else vm++ } }
  END{ printf("  viral-primary reads: %d\n", vt);
       printf("    NH==1 (UNIQUELY viral -> real):        %d (%.1f%%)\n", vu, vt?100*vu/vt:0);
       printf("    NH>1  (multimaps host too -> homology): %d (%.1f%%)\n", vm, vt?100*vm/vt:0) }'
echo "[$(date)] done."
