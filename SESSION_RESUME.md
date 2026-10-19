# Session resume — 2026-09-26

Branch `codex/viralscan-v3` @ `886cca0`. Working tree **dirty: 46 modified, 18 untracked.
Nothing has been committed or pushed.** Read `CLAUDE.md` first, then `PLAN.md`.

Env interpreter for every command below:

```
PY=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan/bin/python
```

---

## 1. DONE and verified

### 1.1 Layer 2 — gene programmes (complete, 1124 tests pass)

Shipped: `src/viralscan/gene_programs.py`, `src/viralscan/scripts/gene_programs.py`,
`src/viralscan/data/gene_programs.tsv` (79 rows / 9 viruses), generator
`extras/build_gene_programs.py`, optional Snakemake rule, `viralscan rerun-programs`,
`--gene-programs` / `--programme-min-breadth`, `tests/test_gene_programs.py` (51 tests).

Registered in `pyproject.toml` package-data, `MANIFEST.in`,
`config/public_ship_scope.json` (wheel + sdist), plus `sensitivity.py`.

Measured on the real EBV LCL run (`SRR12682296`), unique vs allocated layer:

| layer | latent | productive | mixed |
|---|---:|---:|---:|
| `counts_unique_viral` | **2,240** | 102 | 445 |
| multimap-allocated | 1,277 | 5 | — |

**0 cells** invert latent-on-unique → productive-on-allocated. The gain is **+75% latent
calls** (1,263 cells fall to `indeterminate` on the allocated layer) — lost sensitivity,
NOT a directional correction. This correction is already applied to the module docstring,
`report.html.j2`, `docs/output_reference.md`, `docs/faq.md`, `PLAN.md` and
`TestEbvLclRegression`. Do not reintroduce the older "selected layer cannot call latent" claim.

Overlap groups in EBV: `g5={EBNA-1,EBNA-2,EBNA-LP}`, `g50={BNLF2a,BNLF2b,LMP-1}`,
`g39={BTRF1,BcLF1}`. 14 groups, largest 6. Computed by exonic-interval union-find in
`extras/build_gene_programs.py` — **never hand-assign these**; the sweep and the bounding-box
variants were both wrong and only generating them found it.

Catalogue scope: `complete` = EBV, CMV, HHV-6A, HHV-7. `partial` (`latent` unreachable) =
HSV-1, HSV-2, HHV-6B, VZV, KSHV.

### 1.2 Prior P0–P3 work (also done)

`sensitivity.py` + `results/sensitivity.tsv`, positive-control reporting, anellovirus gene-ID
visibility, viral-name aliases, pinned STAR parameters, `pct_infected_comparable`,
benchmark claim retraction. Documented in `.living/` and `PLAN.md`.

### 1.3 Reference audit — findings established this session (NEW, all measured)

Measured on the existing built reference
`/exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/viralscan_ref/`:

| | count | gene structure |
|---|---:|---|
| anellovirus placeholders `<acc>_gene1` | **2,043** | 1 pseudo-gene each |
| real-gene accessions, **≥2 genes** | **33** | genuine — all herpesviruses |
| real-gene accessions, exactly 1 gene | 1,561 | 1 gene each |
| **HPV** | **4** of >200 types | generic `HpV16gp1…gp8` |

- **Only 33 of 1,637 viral genomes are quantifiable at gene level.** This is why the
  catalogue is exactly the 9 herpesviruses — not a coincidence, the only family the
  reference can resolve.
- **COVID artifact:** 4,647 viral gene IDs, all single-gene-per-accession. Total viral
  signal 1,169,272 UMI, of which **1,167,103 (99.8%) landed in one bucket**
  `MW455439.1_gene1`, reported as "Alphatorquevirus". The per-virus breakdown
  (Alphatorque 1.17M vs Betatorque 1,241) is a cross-mapping sink, not a measurement.
- **ROOT CAUSE FOUND for the anellovirus problem:** NCBI *does* carry real CDS for
  anelloviruses. Live fetch of `NC_002076.2` GenBank shows `gene:3 CDS:3 intron:3
  regulatory:3`, a contiguous `589..2901` ORF1/Rep and two spliced CDS
  (`join(353..711,2564..3077)`, `join(353..711,2374..2875)`), products VP1/VP2/VP3.
  `src/viralscan/scripts/ncbi_fetch.py` has TWO GTF paths — `_genbank_to_gtf` (real CDS,
  handles spliced `join()`) and `_whole_genome_gtf_from_fasta` (placeholder). **The
  anelloviruses took the FASTA path.** This is a code-path bug, not missing data.
  Cache at `~/.cache/viralscan/ncbi/<acc>/` holds only `.fasta` + placeholder `.gtf`;
  the GenBank flatfile was never retained.

### 1.4 SFL tonsil CITE-seq — data characterised, FASTQ extracted

Path `/exports/archive/hg-funcgenom-research/mdmanurung/_misc/sfl_tonsil_citeseq`.

24 donors, hash `s1`–`s24`, cell-hash multiplexed. 18y+ non-active tonsillitis (`s1`–`s18`);
2y–9y tonsillar hypertrophy (`s19`–`s24`). 12F/12M, AFR 14 / EUR 10.
55,239 barcodes → **36,114 singlets**; souporcell recovers **24 genotype clusters**;
ambient RNA **13.29%**.

**The count matrix cannot answer a viral question — verified, not assumed:**

| | |
|---|---|
| viral gene IDs in the reference GTFs | 2,691 |
| features in `features.tsv.gz` | 38,767 — **99.6% Ensembl `ENSG*`** |
| **intersection** | **0** |

Cellranger counted against a human-only reference, so viral reads were never quantified.
Any viral call from that matrix would be a *manufactured* negative. This is the same
failure class as the TTV work, in absolute form.

Chemistry: 10x 3′, R2 = 90 bp (20 bp CB+UMI + ~70 bp insert).

**Expected-value calculation (governs the whole exercise):**
tonsil B cells latently infected with EBV ≈ 1–10 per 10⁶; 36,114 singlets → 14–21k B cells
⇒ **E[EBV⁺ cells] ≈ 0.02–0.15**. A zero result is uninterpretable — it is what the
arithmetic predicts regardless of infection status. This is `SENS-06` (still open).

**Extraction COMPLETE, clean exit (empty `extract.log`), 202 GB** in
`benchmark_runs/sfl_tonsil_screen_2026-09-26/2025-3178-LUM-SJ-x223-x226/raw_fastq/`:

| file | size |
|---|---:|
| `LUM-SJ-x223_fastq/LUM-SJ-x223_S4_L004_R1_001.fastq.gz` | 31.3 GB |
| `LUM-SJ-x223_fastq/LUM-SJ-x223_S4_L004_R2_001.fastq.gz` | 70.5 GB |
| `LUM-SJ-x225_fastq/LUM-SJ-x225_S6_L006_R1_001.fastq.gz` | 23.6 GB |
| `LUM-SJ-x225_fastq/LUM-SJ-x225_S6_L006_R2_001.fastq.gz` | 54.2 GB |

**Note: only 2 of the 4 samples exist** (x223, x225), one lane each. The tarball is named
`x223-x226` but contains no x224/x226. Do not assume 4 samples.

---

## 2. NOT DONE

### 2.1 Option A — HPV reference with named genes (NOT STARTED)

Subagent task was written and issued; **both subagent calls were cancelled before any work
ran.** No files created.

Needed: fetch ~14 high-risk types (16/18/31/33/35/39/45/51/52/56/58/59/66/69) via the
`_genbank_to_gtf` path; determine whether NCBI `/product=` names are semantic (`E6`,`E7`) or
generic; if generic, map to canonical gene identity **by genomic coordinate** and cite the
coordinate source. Currently 4 genotypes with `gpN` names — oncogene E6/E7 activity cannot be
distinguished from passive L1. **Constraint: must NOT edit `src/viralscan/scripts/ncbi_fetch.py`
— that file is owned by Option B.**

### 2.2 Option B — anellovirus real gene structure (NOT STARTED)

Same status: **nothing ran.** Deliverable is to route the 2,042-accession panel through
`_genbank_to_gtf` instead of `_whole_genome_gtf_from_fasta`, retain the GenBank flatfile in
the cache, handle **circular genome topology** (anelloviruses are circular ssDNA — a wrong
wrap-around coordinate silently corrupts every gene call), emit real gene IDs matching repo
convention, and decide whether genogroup can be assigned. **Owns `ncbi_fetch.py`.**
Central regression test: representative anellovirus accessions now yield ≥1 real multi-exon
gene, not `<acc>_gene1`.

### 2.3 Tonsil screen (NOT RUN)

`kb count` has **never been run** on the SFL data. Reusable index is
`.../covid_viralscan/viralscan_ref/index.idx` (522 MB, 470,468 t2g pairs, 2,656 real viral
gene symbols, 2,043 placeholders) with `t2g.txt` and `combined.gtf` alongside.
`kb` and `snakemake` exist in the test_viralscan env but are **not on `PATH`** by default.

Intended design (agreed): run pooled on x223+x225; **trust gene-level herpesvirus/HPV calls;
treat anelloviruses as a BINARY technical control only — do not report per-torquevirus
abundance;** record the EBV expected-value calculation rather than a result.

### 2.4 COVID claim retraction (AWAITING USER AGREEMENT — do not edit published results yet)

`covid_viralscan/results/LUM-SJ-x213-g/results/viral_summary.tsv` publishes
`Alphatorquevirus = 1,167,103` with no caveat. On the evidence in §1.3 it should be marked
`not_applicable` for quantitative interpretation, matching the 8 excluded viruses. I asked
for agreement and did not get it. **Blocked pending explicit approval.**

### 2.5 Commit-time only

`analysis/v3_artifact_inventory.tsv` needs rows for the new files. `governance_utils` resolves
paths from the git tree at the recorded SHA, so they cannot validate uncommitted — which is
why `tests/test_artifact_inventory.py::test_repository_artifact_inventory_is_valid` fails.
It fails identically at clean `HEAD`; treat as pre-existing.

---

## 3. Blockers / constraints

- `G0` — PyPA `build` frontend missing; real wheel/sdist membership only partially verifiable.
- `REF-11` — Zenodo DOI `10.5281/zenodo.20112332` unregistered; blocks clean-install
  validation of the bundled reference path.
- `SENS-06` — nothing plants a control, so `informative_negative` cannot become `true`.
  This is why the tonsil EBV negative would be meaningless.
- `HOST-05` — the v3 STARsolo path's viral-signal cost is unmeasured.
- `PROG-08` / `PROG-09` — 5 of 9 viruses are `partial`; scope of excluded viruses.
- `VAL-01…10` — truth panel; FPR/FNR formally uncharacterised.
- Two subagent tasks were cancelled. Re-dispatch, or do the work inline.

---

## 4. Resume commands

```bash
cd /exports/archive/hg-funcgenom-research/mdmanurung/ViralScan
PY=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan/bin/python
export PATH=/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan/bin:$PATH

# state
git status --porcelain | head -60
git log --oneline -1

# tests / lint  (expect 1 known failure: test_artifact_inventory)
PYTHONPATH=src $PY -m pytest tests/ -q
ruff check $(git status --porcelain | awk '{print $2}' | grep -E '\.py$')

# SFL FASTQ (extraction already done, nothing to redo)
ls -la benchmark_runs/sfl_tonsil_screen_2026-09-26/2025-3178-LUM-SJ-x223-x226/raw_fastq/*/

# reusable index
REF=/exports/para-lipg-hpc/mdmanurung/ViralScan/covid_viralscan/viralscan_ref
ls -la $REF/index.idx $REF/t2g.txt $REF/combined.gtf
```

### Suggested order

1. Re-dispatch Option A and Option B as **parallel** subagents (A must not touch
   `ncbi_fetch.py`; B owns it). Both were briefed and ready to run.
2. While they work, run the pooled `kb count` screen on x223 + x225.
3. After B lands, re-examine whether the covid anellovirus call is an artifact — that is the
   decisive test for §2.4.
4. Get user approval for the §2.4 retraction.
5. Update `PLAN.md` in the same change as any implementation (mandatory per CLAUDE.md).

---

## 5. Standing lessons — do not regress these

- Every major defect in this repo was a **confident, clean-looking result for the wrong
  reason**. Nothing crashed and nothing warned.
- Depth was never the limit. Poisson floor ~3 molecules vs 5–20M quantified. **Reference
  capture** is the limit (substitution-only k-mer capture heuristic, 90 bp k=31: 0.25 @10% divergence, 0.063 @15%, 0.013 @20%, 0.0003 @30%).
- Where a count comes from, write it in the same breath as the count. A per-gene aggregate
  over overlapping ORFs is uninformative even when every individual number is correct.
- Generate reference annotations; do not hand-assign. Both wrong overlap-group algorithms
  (monotonic sweep, bounding boxes) were only caught by generating.
- `latent`/`mixed` are unreachable when `latency_observable_in_rna=false`. DNA-level latency
  (HIV cccDNA, integrated HPV, HBV) is **invisible to scRNA-seq, not latent**.
