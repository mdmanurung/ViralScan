#!/usr/bin/env bash
# DEF-02 (d): the --strand auto pilot on a 3' library, for the measured 3' ratio.
# EBV LCL SRR12682296 (10xv2 3', detected 2026-10-03), cat42b host+viral index,
# first 1M pairs, three explicit strands; prints the strand_inference block.
#SBATCH -J vs_strand3p
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH -o logs/vs_strand3p_%j.log
#SBATCH -e logs/vs_strand3p_%j.err
set -euo pipefail
ENV=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
B=/exports/archive/hg-funcgenom-research/mdmanurung/viral_ref_cat42b/build
I=/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/reference_strategy/SRR12682296
export PATH=$ENV/bin:$PATH PYTHONPATH=$PWD/src
python - <<PY
import json
from viralscan import strand
rates = strand.run_pilot("$I/SRR12682296_1.fastq.gz", "$I/SRR12682296_2.fastq.gz",
                         "$B/panel.idx", "$B/panel.t2g", "10xv2", None, 8)
print(json.dumps(strand.inference_block(rates), indent=2))
PY
