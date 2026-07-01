# covid_viralscan — Operational Runbook

Broad single-cell viral survey of two 10x 5′ v3 GEX libraries against a combined
human + SARS-CoV-2 + full Serratus/anellovirus reference (~2,313 viral sequences).

Tick each checkbox as it completes. All paths are absolute for cluster use.

---

## Stage 0 — Reference table (no actions)

| Item | Value |
|---|---|
| **Aim** | Broad viral-reactivation survey; SARS-CoV-2 (NC_045512.2) primary; SARS-CoV-1 (NC_004718.3) negative control |
| **Samples** | `LUM-SJ-x213-g` (batch 1), `LUM-SJ-x216-g` (batch 2) |
| **Chemistry** | 10x 5′ v3 GEX (CB 16 bp, UMI 12 bp) |
| **REPO** | `/exports/para-lipg-hpc/mdmanurung/ViralScan` |
| **Conda env** | `/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan` |
| **FASTQ root** | `/exports/para-lipg-hpc/mdmanurung/covid_viralscan_fastqs/2025-3089-LUM-SJ-x213-x218-gtm/raw_fastq` |
| **Ref dir** | `covid_viralscan/viralscan_ref/` → `index.idx`, `t2g.txt`, `cdna.fa` |
| **SARS-CoV-2 GTF** | `covid_viralscan/viralscan_ref/viral/viral_whole_genome.gtf` |
| **Serratus+anello GTF** | `references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf` |
| **CellRanger cells** | `/exports/para-lipg-hpc/Youvika/20250605_scRNAseq_YS/20250814_tino_scRNAseq_batch2_YS/data_raw/202502341a_count_v2/outs` |

**Gene-id naming notes (important for interpreting results):**
- SARS-CoV-2: `virus_name = "NC_045512.2_gene1"` in `viral_summary.tsv` (whole-genome
  single-gene GTF; no VIRUS_NAME_MAP key — falls through to raw gene_id).
- SARS-CoV-1: `virus_name` starts with `"sarsp"` (Serratus RefSeq locus tags from
  NC_004718.3 in `viral_genome.gtf`).
- Anellovirus genera (Alphatorquevirus etc.): resolved by `merged_name_map()` via the
  clareaulab accession table.

---

## Stage 2 — Reference build

### History
- **Job 25137178** (2026-07-01): Ran Steps 1+2 successfully (combined.fa 1.4 GB ready).
  Step 3 (`kb ref`) **hung for 3h29m with 0 bytes output** — cancelled.
  **Root cause**: `combined.gtf` had chromosomal seqnames (`1`, `2`, `X`, …) that don't
  match the cDNA FASTA headers (ENST transcript IDs). ngs_tools's "Splitting genome"
  step scanned 1.4 GB of FASTA trying to find chromosomal sequences that don't exist.
  **Fix**: `covid_viralscan/scripts/gen_combined_cdna_gtf.py` generates a cDNA-level GTF
  (seqname = ENST transcript ID, coords = 1 to transcript length). Viral GTFs are appended
  with the same Step-2.5 dedup sanitization.

- **Job 25138039** (2026-07-01): Re-run using `slurm_build_ref_v2.sh` + the cDNA GTF fix.
  Steps 1+2 skipped (combined.fa already exists). Steps A+B in progress.

Build script (repair): `covid_viralscan/scripts/slurm_build_ref_v2.sh`  
GTF generator:         `covid_viralscan/scripts/gen_combined_cdna_gtf.py`

- [ ] **2.1 — Confirm `kb ref` completed without errors**

  ```bash
  # Check log for completion marker
  tail -20 /exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/logs/build_ref_v2_25138039.log
  # Expect: "Stage 2 v2 complete." (includes post-build sanity check output)
  
  # Check error file for Python tracebacks
  cat /exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/logs/build_ref_v2_25138039.err
  ```

- [ ] **2.2 — Verify artifacts exist and are non-trivial**

  ```bash
  ls -lh /exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/viralscan_ref/{index.idx,t2g.txt,cdna.fa}
  # Expect: index.idx ≥ 1 GB; t2g.txt ≥ 100k lines; cdna.fa ≥ 500 MB
  ```

- [ ] **2.3 — Sanity-check t2g.txt**

  ```bash
  REF=/exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/viralscan_ref
  grep -c "NC_045512" $REF/t2g.txt     # SARS-CoV-2 must be present (≥ 1)
  grep -c "ENST"      $REF/t2g.txt     # Human transcripts (≥ 100k)
  grep -c "_gene"     $REF/t2g.txt     # Anellovirus geneN entries (≥ 1000)
  grep -c "sarsp"     $REF/t2g.txt     # SARS-CoV-1 Serratus locus tags (≥ 1)
  ```

  If any check returns 0 → re-check the build log and GTF (see ngs_tools fix notes above).

---

## Stage 3 — Quant (2-sample SLURM array)

Quant script: `covid_viralscan/scripts/slurm_viralscan_quant.sh`  
Resources: 16 CPU / 256 GB / 48 h per sample (deep 5′ libraries)

- [ ] **3.1 — Submit quant array** *(only after Stage 2 artifacts verified)*

  ```bash
  cd /exports/para-lipg-hpc/mdmanurung/ViralScan
  sbatch --array=0-1 covid_viralscan/scripts/slurm_viralscan_quant.sh
  # Records job ID (call it JQ): squeue -u mdmanurung | grep covid_vs_quant
  ```

- [ ] **3.2 — Monitor until both tasks complete**

  ```bash
  # Watch log tails (replace JOBARRAY with the actual job array ID from 3.1)
  tail -30 covid_viralscan/logs/quant_JOBARRAY_0.log   # LUM-SJ-x213-g
  tail -30 covid_viralscan/logs/quant_JOBARRAY_1.log   # LUM-SJ-x216-g
  # Expect: "viralscan quant complete." then "Stage 3 done."
  
  # Check for failures
  cat covid_viralscan/logs/quant_JOBARRAY_0.err
  cat covid_viralscan/logs/quant_JOBARRAY_1.err
  ```

- [ ] **3.3 — Check pseudoalignment rates** *(record both values)*

  ```bash
  RESULTS=/exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/results
  for SAMPLE in LUM-SJ-x213-g LUM-SJ-x216-g; do
      echo "=== $SAMPLE ==="
      python -c "import json; d=json.load(open('$RESULTS/$SAMPLE/kb-python/run_info.json')); \
                 print('p_pseudoaligned:', d['p_pseudoaligned'], '%')"
  done
  # Typical for host-dominant samples: 20–60% (human transcriptome captures most reads)
  ```

- [ ] **3.4 — Verify output TSVs exist**

  ```bash
  for SAMPLE in LUM-SJ-x213-g LUM-SJ-x216-g; do
      ls -lh $RESULTS/$SAMPLE/results/{viral_summary.tsv,per_cell_viral.tsv}
  done
  ```

---

## Stage 4 — Analysis & interpretation

Helper script: `covid_viralscan/scripts/summarize_survey.py`  
Output: `covid_viralscan/results/SURVEY_SUMMARY.md`

- [ ] **4.1 — Run the survey summary helper**

  ```bash
  cd /exports/para-lipg-hpc/mdmanurung/ViralScan
  export PYTHONPATH=$PWD/src
  CONDA_PY=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan/bin/python
  CR_OUTS=/exports/para-lipg-hpc/Youvika/20250605_scRNAseq_YS/20250814_tino_scRNAseq_batch2_YS/data_raw/202502341a_count_v2/outs

  $CONDA_PY covid_viralscan/scripts/summarize_survey.py \
      --results-dir covid_viralscan/results \
      --samples LUM-SJ-x213-g LUM-SJ-x216-g \
      --cellranger-outs "$CR_OUTS" \
      --output covid_viralscan/results/SURVEY_SUMMARY.md
  ```

  The script outputs:
  - **Ranked viral panel**: all viruses sorted by total UMI across both samples
  - **Specificity control**: NC_045512.2 (SARS-CoV-2) vs NC_004718.3 / sarsp\* (SARS-CoV-1)
  - **Cell-barcode overlap**: how many viral+ barcodes are in CellRanger's called-cell whitelist

- [ ] **4.2 — Inspect ranked panel (broad survey readout)**

  ```bash
  # Quick table preview — rank all viruses by combined UMI
  python -c "
  import pandas as pd
  dfs = []
  for s in ['LUM-SJ-x213-g', 'LUM-SJ-x216-g']:
      df = pd.read_csv(f'covid_viralscan/results/{s}/results/viral_summary.tsv', sep='\t')
      df['sample'] = s
      dfs.append(df)
  combined = pd.concat(dfs).groupby('virus_name')['total_umi'].sum().sort_values(ascending=False)
  print(combined.head(30).to_string())
  "
  ```

- [ ] **4.3 — Check SARS-CoV-2 vs SARS-CoV-1 specificity**

  ```bash
  # SARS-CoV-2 signal (expected: detectable if any COVID-positive cells)
  grep "NC_045512" covid_viralscan/results/LUM-SJ-x213-g/results/viral_summary.tsv
  grep "NC_045512" covid_viralscan/results/LUM-SJ-x216-g/results/viral_summary.tsv

  # SARS-CoV-1 signal (expected: ~0 — this is the negative specificity control)
  grep "sarsp" covid_viralscan/results/LUM-SJ-x213-g/results/viral_summary.tsv
  grep "sarsp" covid_viralscan/results/LUM-SJ-x216-g/results/viral_summary.tsv
  ```

  If SARS-CoV-1 > 0 UMI: inspect the reads (possible cross-mapping; apply `viralscan evidence`).

- [ ] **4.4 — Review CellRanger barcode overlap in SURVEY_SUMMARY.md**

  ```bash
  cat covid_viralscan/results/SURVEY_SUMMARY.md
  # Key metric: "% viral+ barcodes that are CellRanger-called cells"
  # High % = confident cell-level signal; low % = ambient RNA or empty droplets
  ```

- [ ] **4.5 — Commit results**

  ```bash
  git add covid_viralscan/results/SURVEY_SUMMARY.md
  git add covid_viralscan/scripts/summarize_survey.py
  git commit -m "feat(covid): broad viral survey results — SURVEY_SUMMARY.md (closes covid Stage 4)"
  ```

---

## Troubleshooting notes

| Symptom | Likely cause | Fix |
|---|---|---|
| `kb ref` produces empty `cdna.fa` | ngs_tools GTF parse error (tracebacks in `.err`) | Check Step 2.5 ran; inspect `.err` for `GtfEntryError` |
| `t2g.txt` missing ENST entries | GTF seqname/FASTA header mismatch | Check `cdna.fa` headers match GTF seqnames |
| `kb count` exits but no `.bus` file | Memory OOM (256 GB requested; reduce if node limit lower) | Check `sacct -j <JOB> --format=MaxRSS` |
| Stage 3 fails at prereq check | `index.idx` or GTF not found | Re-run Stage 2 checks (2.1–2.3) |
| SARS-CoV-2 = 0 UMI | Not necessarily wrong (uninfected cells/PBMC) | Cross-ref with clinical metadata; check `p_pseudoaligned` |
