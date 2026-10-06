#!/usr/bin/env bash
#SBATCH -J dsr_sra
#SBATCH --array=0-17%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH -o /exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/dsr_2026-10-05/logs/%x_%A_%a.log
# DSR plan row-2 downloads (docs/plans/2026-10-05-dataset-run-plan.md): SRA runs whose
# barcode read is only in SRA (ENA serves R2 only). Output _1 = I1, _2 = R1 (CB+UMI), _3 = R2.
set -euo pipefail
module load bioinformatics/tools/ncbi/sra/3.0.10
D=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/dsr_2026-10-05
read -r DS SRR < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $D/sra_runs.tsv)
O=$D/$DS/$SRR; mkdir -p $O $D/sra
[ -e $O/DONE ] && { echo "already done: $SRR"; exit 0; }
prefetch --max-size 200G -O $D/sra $SRR
vdb-validate $D/sra/$SRR/$SRR.sra
fasterq-dump --include-technical --split-files -e 8 -t $O/tmp -O $O $D/sra/$SRR/$SRR.sra
rm -rf $O/tmp
for f in $O/*.fastq; do echo "$(basename $f) $(sed -n 2p $f | tr -d '\n' | wc -c)bp"; done | tee $O/read_lengths.txt
pigz -p 8 $O/*.fastq
(cd $O && md5sum *.fastq.gz > md5.txt)
rm -f $D/sra/$SRR/$SRR.sra
touch $O/DONE
