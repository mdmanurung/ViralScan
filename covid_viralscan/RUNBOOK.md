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

- **Job 25138052** (2026-07-01→02): Re-run using `slurm_build_ref_v2.sh` + the cDNA GTF fix.
  Steps A (GTF gen) + cDNA extraction **succeeded** (`cdna.fa` 1.2 GB, `t2g.txt` 470,472
  entries, all sanity greps passed). But the final `kallisto index` **failed**:
  `Error: repeated name in FASTA file cdna.fa`.
  **Root cause**: the source panel `viral_genome.fa` contains accession `NC_002076.2`
  (Torque teno virus 1) **twice** (byte-identical — it is in both the anellovirus and
  Serratus sets). ngs_tools extracted its 4 gene/transcript models once per copy →
  4 duplicate names (`TTVgp1/2/3`, `NC_002076.2_tx1`) in `cdna.fa`, which kallisto rejects.
  **Fix**: keep-first dedup of `cdna.fa` (470,472→470,468), `t2g.txt` (→470,468), and the
  D-list `combined.fa` (468,083→468,082), then run `kallisto index` directly on the
  deduped inputs (equals a clean build from a deduped panel since the removed records are
  byte-identical). A dedup guard (Step 2.6) was also added to `slurm_build_ref.sh` so a
  from-scratch rebuild handles any duplicate accession automatically.

- **Job 25138556** (2026-07-02): `slurm_kallisto_index.sh` — resumes step B on the deduped
  inputs and promotes the deduped files to canonical names. Verifies by artifact
  (`index.idx` non-zero + `kallisto inspect`), not exit code.

Build script (repair):   `covid_viralscan/scripts/slurm_build_ref_v2.sh`  
Index resume:            `covid_viralscan/scripts/slurm_kallisto_index.sh`  
GTF generator:           `covid_viralscan/scripts/gen_combined_cdna_gtf.py`

- [x] **2.1 — Confirm `kb ref` completed without errors** — job 25138556 finished exit 0;
  `kallisto inspect` reports 470,468 targets / 85.1M k-mers / 621,183 D-list k-mers.

  ```bash
  # Check log for completion marker (index resume job)
  tail -25 /exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/logs/kallisto_index_25138556.log
  # Expect: "Stage 2 v2 (index resume) complete." + kallisto inspect output

  # Check error file for tracebacks
  cat /exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/logs/kallisto_index_25138556.err
  ```

- [x] **2.2 — Verify artifacts exist and are non-trivial** — `index.idx` 499 MB,
  `t2g.txt` 470,468 lines, `cdna.fa` 1.2 GB.

  ```bash
  ls -lh /exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/viralscan_ref/{index.idx,t2g.txt,cdna.fa}
  # Expect: index.idx ≥ 1 GB; t2g.txt ≥ 100k lines; cdna.fa ≥ 500 MB
  ```

- [x] **2.3 — Sanity-check t2g.txt** — NC_045512=1, ENST=465,769, _gene=2,054, sarsp=1,
  NC_002076.2_tx1=1 (dedup confirmed). All checks pass.

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

### History
- **Job 25138573/25138586** (2026-07-02): failed. First (25138573) died on a `zcat|head`
  SIGPIPE in the read-length sanity line (fixed: command substitution + `|| true`).
  Re-run (25138586) then ran `kb count` (1.2 B reads, 77.4 M pseudoaligned = **6.4 %**) but
  produced **no `counts_unfiltered/`** and failed at `mv counts_unfiltered/`. Root cause
  (found via `diag_bustools_count.sh`): `bustools correct` cannot read a **gzipped**
  whitelist — it read the compressed bytes as barcodes and died `on-list file malformed`.
  kb count swallowed the error (Snakefile captured `kb count 2>&1` into a discarded var).
  **Fix**: Snakefile `kb_count` rule now decompresses a `*.gz` whitelist, tees kb output to
  `<output>kb_count.log`, and verifies `counts_unfiltered/` exists. Re-run: **job 25139346**.
- **Job 25139346** (2026-07-02): whitelist fix let counting complete, BUT the result was
  **invalid** — the matrix was ~all empty droplets (163k barcodes, median 1 UMI, max 5,311)
  vs CellRanger's 28,922 real cells. Root cause: the bundled **10x v3 whitelist is WRONG** for
  this chemistry — only 0.4 % of raw R1 barcodes (and 0.5 % of CellRanger cells) are in it, so
  `bustools correct` dropped 96.5 % of reads. Raw R1 barcodes match no bundled ngs_tools list
  (best v4/GEM-X 4.7 %) but match CellRanger's cells 56.5 % → this is a GEM-X-5′-like chemistry.
- **Job 25140008** (2026-07-02): re-run with the **correct whitelist** — CellRanger's raw
  barcode universe (2,974,869 barcodes) extracted from the matched `cellranger-multi` run:
  ```bash
  RAW=<cellranger>/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz
  zcat "$RAW" | sed 's/-[0-9]*$//' | sort -u > covid_viralscan/viralscan_ref/cellranger_whitelist.txt
  ```
  Raw R1 barcodes match this whitelist **68.1 %** (vs 0.4 % for v3). Same chemistry whitelist
  applies to both samples. `slurm_viralscan_quant.sh` `WHITELIST` now points at it.
  ⚠️ Until 25140008 verifies (real cells with sane per-barcode UMI, viral+ ∩ CellRanger cells),
  treat all per-cell covid numbers as provisional.

Quant script: `covid_viralscan/scripts/slurm_viralscan_quant.sh`  
Resources: 16 CPU / 256 GB / 48 h per sample (deep 5′ libraries)

- [x] **3.1 — Submit quant array** *(only after Stage 2 artifacts verified)* —
  submitted 2026-07-02 as **job 25138573** (array 0–1).

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
