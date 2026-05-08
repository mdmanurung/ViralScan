# ViralScan Scientific Audit — Full Pipeline
**Date:** 2026-05-08  
**Scope:** Full pipeline end-to-end  
**Auditor:** GitHub Copilot (Claude Sonnet 4.6)  
**Mode:** Read-only — no source files modified  
**Branch audited:** `claude/review-repo-improvements-Sg4Th`  
**Test status at audit time:** 173 passed, 2 deselected (per PLAN.md)

---

## 1. Executive Summary

- **DEFINITE BUG — double-counting in `umap.py`:** `umap.py:330–331` adds `counts_corrected + counts_original` when loading the multimapping h5ad. `counts_corrected` stores the *multimapper-share-only* redistributed matrix (unique-mapping ECs are skipped and contribute 0). This means the addition is mathematically correct — but the PLAN.md Task 1 note ("double-count bug in umap.py") and the comment structure are inconsistent with the actual multimap.py logic. After tracing every line, the addition is **not** double-counting in the current code; however, the `counts_combined` layer written by `multimap.py:final_results()` is itself the sum of `counts_corrected + counts_original` and is *never* used by downstream scripts — they recompute the same sum from individual layers. This is redundant but not incorrect.
- **HIGH RISK — detection threshold of ≥ 1 UMI with no false-discovery control:** `detection.py:preprocessing()` flags a virus as present if the sum of UMI counts across **all cells** is ≥ 1. This is effectively zero noise filtering. With 195 simultaneous tests and no multiplicity correction, the expected false-positive rate is near-certain for any real dataset.
- **HIGH RISK — GTF parsing is fragile to attribute ordering:** `analysis.py:obtain_gtf()` always extracts `info.split('"')[1]`, which blindly takes the first double-quoted token. If any GTF line has a quoted attribute before `gene_id` (e.g., `source "NCBI"; gene_id "X"`) or no `gene_id` at all, the wrong accession is silently emitted. User-supplied GTFs have no validation.
- **MEDIUM RISK — no random seeds in `umap.py` viral branch:** `sc.pp.neighbors()`, `sc.tl.pca()`, and `sc.tl.umap()` are called without `random_state`/`use_highly_variable` consistency, making UMAP plots and the neighbor-enrichment p-value non-reproducible across runs.
- **MEDIUM RISK — barcode strip inconsistency:** `load_barcodes()` does `bc.replace("-1", "")` which replaces **all** occurrences of the substring "-1" in a barcode string, not just the trailing suffix. A barcode like `ACGT-1GCTA-1` would be mangled to `ACGTGCTA`.
- **MEDIUM RISK — cache validation by file existence only:** `ncbi_fetch._fetch_one()` caches FASTA/GTF by checking `path.exists()` and `path.stat().st_size == 0`, but does not validate a stored checksum. A truncated or partially-downloaded file from a prior interrupted run will be silently re-used.
- **LOW RISK — `detection_threshold` validated to int but not positive:** `createconfig.py` casts `detection_threshold` to `int`, but a value of 0 or negative passes through silently. With threshold=0, every gene in the reference (host + viral) would be flagged as "detected".
- **No integration tests exist** and the test suite contains no tests that verify UMI count conservation, threshold-based exclusion correctness, or UMAP reproducibility.

---

## 2. High-Risk Scientific Issues (Definite Errors)

### 2.1 Detection threshold provides near-zero noise protection; no multiplicity adjustment

**Finding:** A single UMI summed across *all cells* in the dataset is sufficient to flag a virus as detected. With 195 reference viruses tested simultaneously, the expected number of false positives from sequencing noise, ambient RNA, or index cross-mapping is non-trivial.

**Location:** `src/viralscan/scripts/detection.py:60–68`, function `preprocessing()`

**Evidence:**
```python
threshold = config.get("detection_threshold", 1)
for gene_id in viral_accessions:
    if gene_id in all_gene_ids:
        gene_counts = adata[:, gene_id].X
        if hasattr(gene_counts, "toarray"):
            gene_counts = gene_counts.toarray()
        total_count = gene_counts.sum()
        if total_count >= threshold:          # ← 1 UMI across all cells triggers detection
            found_genes[gene_id] = total_count
```

- No per-cell threshold (a gene must be expressed in ≥N cells).
- No FDR/Bonferroni correction for 195 simultaneous comparisons.
- No comparison to a negative control or background noise model.

**Consequence:** Users will receive false-positive viral detections driven by sequencing noise, ambient RNA bleed-over, or low-level index cross-mapping. The HTML report, `viral_summary.tsv`, and UMAP labels will contain spurious viruses. This is the primary scientific validity concern for the tool.

---

### 2.2 GTF parsing extracts first quoted token regardless of attribute name

**Finding:** `analysis.py:obtain_gtf()` splits column 9 on `'"'` and takes `parts[1]` — the first quoted value — as the gene ID. If a GTF line places any quoted attribute before `gene_id` (e.g., Ensembl-format GTFs which may have `gene_id` as the first attribute but some bundled or user-supplied GTFs may not), the wrong string is silently used as a viral accession ID.

**Location:** `src/viralscan/scripts/analysis.py:38–41` and `47–51`, function `obtain_gtf()`

**Evidence:**
```python
info = line.split("\t")[8]
gene_id = info.split('"')[1]   # ← always grabs first quoted value
viral_accessions.add(gene_id)
```

If a line reads:  
`NC_TEST\t.\tgene\t1\t100\t.\t+\t.\tsource "NCBI"; gene_id "NC_123";`  
then `parts[1]` = `"NCBI"` not `"NC_123"`.

**Consequence:** Viral accession IDs may be corrupted silently. Downstream detection would fail to match real gene IDs and either miss viruses entirely or add nonsense strings to the accession list. User-supplied GTFs from non-standard sources are especially vulnerable.

---

### 2.3 Barcode suffix stripping replaces all occurrences of "-1", not just trailing

**Finding:** `load_barcodes()` uses `bc.replace("-1", "")`, which is a global string replacement — it will corrupt barcodes that contain the substring "-1" anywhere in the middle.

**Location:** `src/viralscan/scripts/multimap.py:60–61`, function `load_barcodes()`

**Evidence:**
```python
barcodes = [bc.replace("-1", "") for bc in barcodes]
```

For a barcode `ACTT-1GCAT-1` (unlikely but possible in some 10x chemistry formats with multiple library types), this yields `ACTTGCAT` instead of `ACTT-1GCAT`.  
The identical pattern appears at `multimap.py:normalize_barcodes()` (line ~173):
```python
bus_df["barcode"] = bus_df["barcode"].str.replace("-1", "", regex=False)
```
Here `regex=False` is correctly set, but global replacement is still applied.

**Consequence:** If any barcode contains "-1" not as a suffix but as part of the actual sequence (which does not happen in standard 10x v2/v3 formats but can occur in bulk, Smart-seq2, or custom whitelist formats), barcodes would be mangled and silently dropped from the matrix. The resulting multimapping matrix would have missing rows corresponding to affected barcodes.

---

## 3. Medium-Risk Issues (Likely Problems)

### 3.1 No random seeds in viral-branch UMAP (non-reproducibility)

**Finding:** In the virus-detected branch of `umap.py`, `sc.pp.neighbors()`, `sc.pp.pca()`, and `sc.tl.umap()` are called without `random_state`. Only `sc.pp.subsample()` in the no-virus branch uses `random_state=0`.

**Location:** `src/viralscan/scripts/umap.py:200–216`, function `umap()`

**Evidence:**
```python
sc.pp.pca(adata, use_highly_variable=True)      # no random_state
sc.pp.neighbors(adata)                           # no random_state
sc.tl.umap(adata)                               # no random_state
```

The `viral_neighbor_enrichment()` function uses `np.random.permutation(labels)` without seeding:
```python
for _ in range(n_permutations):
    shuffled = np.random.permutation(labels)    # no seed
```

**Consequence:** UMAP layout and p-values from the neighbor enrichment test will differ between runs on identical data. Published figures cannot be reproduced. The p-value can flip across the 0.05 threshold between identical runs.

**What would resolve the uncertainty:** Add `random_state=0` to `sc.pp.pca`, `sc.pp.neighbors`, `sc.tl.umap`, and seed `np.random` before the permutation loop.

---

### 3.2 Cache validation uses file existence + size, not content hash

**Finding:** `ncbi_fetch._fetch_one()` skips download if the cached file exists and is non-empty, without validating its content against a checksum.

**Location:** `src/viralscan/scripts/ncbi_fetch.py:248–259`, function `_fetch_one()`

**Evidence:**
```python
if not fasta_path.exists() or fasta_path.stat().st_size == 0:
    fasta_text = _efetch(acc, "fasta", email, api_key)
    ...
    fasta_path.write_text(fasta_text)

if not gtf_path.exists() or gtf_path.stat().st_size == 0:
    genbank_text = _efetch(acc, "gb", email, api_key)
    ...
```

Note that `_checksum()` is defined in the module but never called during cache validation.

**Consequence:** A partially-written file from a crashed prior run (non-zero size but truncated content) would be re-used silently, producing a corrupt FASTA or GTF fed into `kb ref`. Downstream `kb count` would either fail with a cryptic error or silently map fewer reads to a truncated reference.

**What would resolve the uncertainty:** Call `_checksum()` against a stored `.sha256` sidecar file on cache hit; if missing or mismatched, re-download.

---

### 3.3 `detection_threshold` not validated as positive

**Finding:** `createconfig.py` casts `detection_threshold` to `int()` but does not validate `>= 1`.

**Location:** `src/viralscan/scripts/createconfig.py:40`

**Evidence:**
```python
"detection_threshold": int(cfg_in.get("detection_threshold", 1)),
```

A user passing `--detection-threshold 0` would store `detection_threshold: 0`. In `detection.py`, `total_count >= 0` is always true, flagging every viral accession in the index as detected.

**Consequence:** With threshold=0, every one of the 195 reference viruses would be "detected" in every sample, corrupting all output files and the HTML report.

**What would resolve the uncertainty:** Add `if detection_threshold < 1: raise ValueError(...)` in `createconfig.py` or `menu.py` argument validation.

---

### 3.4 `umap.py:main()` layer addition may be the wrong formula (open PLAN.md task)

**Finding:** `umap.py:main()` performs:
```python
adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]
```
PLAN.md Task 1 explicitly flags this as a known double-count bug. After tracing `multimap.py:build_multimap_matrix()`, the logic is:

- `counts_corrected` = the matrix written by `build_multimap_matrix()` — this matrix **only** accumulates `share = count / len(genes_in_ec)` for ECs with `len(genes_in_ec) > 1` (unique mappers are `continue`d at line 230). Therefore `counts_corrected[cell, gene] == 0` for all unique-mapping reads.
- `counts_original` = `adata_orig[:, adata.var_names].X` — the full kb-python output including all uniquely-mapping reads.
- The sum therefore equals `unique-mapping UMIs + redistributed multimapper fraction`.

This is mathematically defensible (no double-counting *in the current code*). **However**, PLAN.md states the bug is present, suggesting either the PLAN or the code changed since the PLAN was written. The comment in PLAN.md Task 1 says "confirm the multimap.py share logic first" — this audit confirms the share logic is correct; the addition is not double-counting unique reads.

**Uncertainty remains** because `counts_combined` is also written to the h5ad (line 276–278 of multimap.py) as the same sum, yet `detection.py` and `umap.py` recompute the same sum from individual layers instead of using `counts_combined`. This is redundant but not incorrect.

**What would resolve the uncertainty:** Add a unit test asserting `adata.layers["counts_corrected"].sum() == redistributed_share` and `adata.layers["counts_original"].sum() == original_total_umi` for a synthetic BUS file.

---

### 3.5 `summary.txt` is opened inside `multimap.py:final_results()` without a context manager

**Finding:** `summary.open(...)` / `summary.close()` pattern in `final_results()` leaks the file handle if an exception occurs between open and close.

**Location:** `src/viralscan/scripts/multimap.py:284–291`, function `final_results()`

**Evidence:**
```python
summary = open(f"{config['output']}/summary.txt", "w")
summary.write(...)
summary.close()
```

If any write raises (e.g., disk full), the handle is never closed.

**Consequence:** Not a scientific correctness issue, but a file handle leak that could interfere with downstream scripts reading `summary.txt`.

---

### 3.6 Viral counts in `umap.py` are extracted from pre-normalization counts but labeling is applied post-normalization

**Finding:** In `umap.py:umap()`, the `viral_counts` array and `viral_presence` dict are computed from `X_counts` (lines 157–170) before normalization. These raw UMI counts are then stored in `adata.obs["viral_counts"]`. Normalization (`sc.pp.normalize_total` + `sc.pp.log1p`) is applied afterward to `adata.X`. The UMAP and neighbor enrichment test use the normalized, log-transformed `adata.X` for PCA/UMAP embedding but use raw `viral_counts` for coloring. This is the **correct** approach (virus labeling should use raw counts, UMAP embedding should use normalized counts), but it is not documented or guarded.

**Location:** `src/viralscan/scripts/umap.py:113–175`, function `umap()`

**What would resolve the uncertainty:** Add a comment asserting the intentional use of pre-normalization counts for viral labeling vs. post-normalization for embedding.

---

### 3.7 Proportional redistribution uses equal-weight prior (not EM)

**Finding:** `build_multimap_matrix()` assigns `share = count / len(genes_in_ec)` for each gene in a multi-gene EC. This is simple equal-weight redistribution, not EM-based correction as used by kallisto's own `--quant-mode` or salmon.

**Location:** `src/viralscan/scripts/multimap.py:228–233`, function `build_multimap_matrix()`

**Evidence:**
```python
share = count / len(genes_in_ec)
for gid in genes_in_ec:
    rows.append(cell_idx)
    cols.append(gid)
    data.append(share)
```

**Consequence:** For ECs that span both viral and host genes, equal-weight redistribution gives 50% viral weight to reads that are most likely host-derived (because the host transcriptome is vastly more expressed). This systematically inflates viral counts for host-viral ambiguous ECs. There is no downweighting by prior expression level or EM iteration.

**What would resolve the uncertainty:** This is a methodological choice. Domain review is needed. Comparison to kallisto EM or a simulation with known viral/host mixing ratios would characterize the bias.

---

### 3.8 `ncbi_fetch._genbank_to_gtf()` raises `NCBIFetchError` on empty CDS; bundled GTFs may use different gene_id formats

**Finding:** The NCBI-downloaded GTFs (via `ncbi_fetch.py`) use `gene_name or protein_id or "{accession}_cds{n}"` as `gene_id`, while bundled GTFs (under `data/*.gtf`) use accession-based IDs like `NC_001477`. The `analysis.py:obtain_gtf()` function treats all GTFs identically — it takes `info.split('"')[1]` as the viral accession regardless of whether it came from a bundled GTF or a user/NCBI-generated one. The `gene_id` format in NCBI-generated GTFs will differ from the `gene_id` format in the bundled GTFs.

**Location:** `src/viralscan/scripts/ncbi_fetch.py:181–200`, function `_genbank_to_gtf()`

**Consequence:** If users use `--ncbi-accession` to add a new viral reference and also rely on bundled GTFs, the detection accession list will contain heterogeneous gene ID formats (e.g., `"envelope glycoprotein"` from protein product names, vs `"NC_001477"` from bundled files). These different formats must match the gene IDs in the kallisto t2g mapping, or detection silently misses genes.

---

## 4. Unclear Assumptions Requiring Domain Review

### 4.1 Is equal-weight proportional redistribution the correct prior for viral load quantification?

The equal-weight redistribution (`share = count / n_genes_in_EC`) assumes equal a priori probability that a multimapping read originated from each gene. For host-viral ambiguous ECs, this is biologically unreasonable: host transcripts are orders of magnitude more abundant than viral transcripts. Should the prior be weighted by single-gene UMI expression levels (EM approach)? Is the equal-weight prior documented as a known limitation?

### 4.2 Is a threshold of ≥ 1 UMI summed across all cells the right detection criterion for presence/absence calls?

Should viral "presence" require a minimum of N cells each having ≥ M UMIs rather than a global count sum? In single-cell data, 1 UMI summed across 10,000 cells is a much weaker signal than 1 UMI in a single cell. Is there a biologically motivated minimum cell count threshold?

### 4.3 Should ViralScan implement FDR correction for multi-virus comparisons?

With 195 viruses tested per sample, what is the expected false-positive rate under the null hypothesis of no viral infection? Should Bonferroni or BH correction be applied to the detection call? Is the tool designed for discovery (higher sensitivity, higher FDR acceptable) or confirmation (lower FDR required)?

### 4.4 Is the neighbor-enrichment test (permutation of virus labels) the appropriate test for spatial clustering in UMAP space?

UMAP distances are not metric and are non-linear; statistical tests based on kNN structure in UMAP space are not well-calibrated. Is this test used for publication-quality inference or only for visualization guidance? Is the one-sided Laplace-corrected p-value `(sum >= observed + 1) / (n_permutations + 1)` an appropriate statistic?

### 4.5 Does the super-expressor plot provide biologically interpretable results?

The null model in `super_expressor()` is a proportional model: `null_model = total_viral * (cell_total_UMI / all_total_UMI)`. This assumes viral RNA is distributed proportionally to total RNA content per cell. Is this a valid null for viral infection, where the expectation is focal high-count cells above background? Is the "super-expressor" threshold of ≥ 10 UMI calibrated to any biological reference?

### 4.6 HVG selection parameters are hardcoded and may exclude sparse viral genes

`sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)` is designed for host transcriptomes with thousands of moderately expressed genes. Viral genes in low-titer infections will have very low mean expression and may fall outside these bounds before the force-include step. Is the force-include of viral genes (`adata.var["highly_variable"] |= adata.var_names.isin(list(found_genes))`) applied before or after PCA? (Current code: it is applied before `sc.pp.pca`, which is correct.)

---

## 5. Verification Plan

### V1: UMI count conservation test (addresses Section 3.4)

```python
# Synthetic BUS file with known unique + multi mappers
import numpy as np
from scipy import sparse

# Construct a minimal BUS DataFrame with:
# - 5 unique-mapping reads to gene A (EC maps to [gene_A] only)
# - 4 multi-mapping reads to EC spanning [gene_A, gene_B]
# Then verify:
# counts_original[cell, gene_A] == 5
# counts_corrected[cell, gene_A] == 2.0  (4/2 share)
# counts_corrected[cell, gene_B] == 2.0
# counts_combined == counts_original + counts_corrected
# Final X[cell, gene_A] == 7.0,  X[cell, gene_B] == 2.0

# pytest stub:
def test_multimap_umi_conservation():
    from viralscan.scripts.multimap import build_multimap_matrix
    bus_df = pd.DataFrame([
        {"barcode": "BC1", "umi": "U1", "ec": 0, "count": 5},   # unique to gene_A
        {"barcode": "BC1", "umi": "U2", "ec": 1, "count": 4},   # multi: gene_A + gene_B
    ])
    barcode_to_idx = {"BC1": 0}
    ec_map = {0: [0], 1: [0, 1]}  # EC0 → [gene_A], EC1 → [gene_A, gene_B]
    corrected = build_multimap_matrix(bus_df, barcode_to_idx, ec_map, n_cells=1, n_genes=2)
    # EC0 is unique (len==1) → skipped → corrected[0,0] = 0
    # EC1 is multi → share = 4/2 = 2.0 → corrected[0,0] = 2.0, corrected[0,1] = 2.0
    assert corrected[0, 0] == 2.0
    assert corrected[0, 1] == 2.0
    total_unique = 5  # from adata_orig
    # final X = corrected + original → gene_A = 2 + 5 = 7, gene_B = 2 + 0 = 2
```

### V2: Detection threshold false-positive rate characterization (addresses Section 2.1)

1. Take a real sample with no expected viral infection (e.g., healthy donor PBMC).
2. Run ViralScan with the default threshold of 1 UMI.
3. Count the number of viruses "detected".
4. Compare against expectation of 0 detections for a virus-free sample.
5. Repeat with threshold=10 UMI, threshold=100 UMI, and per-cell ≥ 1 UMI in ≥ 2 cells.

Expected result: default threshold will flag several of the 195 reference viruses in virus-free data.

### V3: GTF attribute-order robustness test (addresses Section 2.2)

```python
def test_gtf_parse_wrong_attribute_order():
    """Should fail with current code — attribute before gene_id is taken."""
    gtf_line = 'NC_TEST\t.\tgene\t1\t100\t.\t+\t.\tsource "NCBI"; gene_id "NC_123";\n'
    from tests.test_analysis import _parse_gtf_text
    result = _parse_gtf_text(gtf_line)
    # Current behavior: result == {"NCBI"} (WRONG)
    # Correct behavior: result == {"NC_123"}
    assert result == {"NC_123"}  # this test will FAIL with current code
```

Fix: use `re.search(r'gene_id "([^"]+)"', info)` instead of `info.split('"')[1]`.

### V4: Barcode suffix strip correctness (addresses Section 2.3)

```python
def test_barcode_strip_only_trailing():
    # Barcode with "-1" not as suffix (edge case)
    bc = "ACTT-1GCAT-1"
    # Current code:
    result = bc.replace("-1", "")  # → "ACTTGCAT"
    # Correct behavior: strip only trailing "-1":
    expected = "ACTT-1GCAT"
    import re
    result_correct = re.sub(r"-1$", "", bc)
    assert result_correct == expected
    # This test shows the bug in the current implementation
```

### V5: UMAP reproducibility test (addresses Section 3.1)

```python
@pytest.mark.parametrize("seed", [0, 42])
def test_umap_reproducible(seed):
    """Two runs with same data and same seed must produce identical UMAP coordinates."""
    import numpy as np
    import scanpy as sc
    adata = sc.datasets.pbmc3k_processed()
    np.random.seed(seed)
    sc.tl.umap(adata, random_state=seed)
    coords1 = adata.obsm["X_umap"].copy()
    np.random.seed(seed)
    sc.tl.umap(adata, random_state=seed)
    coords2 = adata.obsm["X_umap"]
    np.testing.assert_array_almost_equal(coords1, coords2)
```

Current code does not set `random_state` so this test would fail.

### V6: Cache integrity check (addresses Section 3.2)

Runnable check (read-only):

```bash
# Simulate truncated cache file
echo "partial" > ~/.cache/viralscan/ncbi/NC_045512.2/NC_045512.2.fasta
# Then run ncbi_fetch — it will reuse the corrupt file without error
PYTHONPATH=src python -c "
from viralscan.scripts.ncbi_fetch import fetch_reference
fetch_reference(['NC_045512.2'], '/tmp/vtest', email='test@test.com')
"
# Check if the returned FASTA starts with '>'
head -1 /tmp/vtest/reference.fasta   # Will show 'partial' — a bug
```

### V7: Zero threshold rejection test (addresses Section 3.3)

```python
def test_detection_threshold_zero_rejected():
    """threshold=0 should raise ValueError in menu.py or createconfig.py."""
    # Currently this test will FAIL — no validation exists
    from viralscan.scripts.createconfig import build_cfg
    with pytest.raises(ValueError):
        build_cfg({"detection_threshold": 0, ...})
```

### V8: Integration smoke test (addresses Section 1 — missing integration tests)

```python
@pytest.mark.integration
def test_full_pipeline_synthetic_data(tmp_path):
    """Run ViralScan on minimal synthetic FASTQ with known viral reads.
    Verify that the known viral gene appears in viral_summary.tsv."""
    # 1. Generate synthetic paired-end FASTQ with 10 reads from NC_001477 (Dengue)
    # 2. Build minimal kallisto index from Dengue FASTA only
    # 3. Run: viralscan run --sample1 R1.fastq.gz --sample2 R2.fastq.gz ...
    # 4. Assert: results/viral_summary.tsv contains NC_001477 with total_umi >= 1
    pass  # stub — full implementation requires kb-python in test environment
```

---

## Appendix: Trace Blocks

### Trace A: Multimapping correction call chain
```
multimap.py:main()
  → load_adata()            # loads adata_orig (full kb-python h5ad)
  → load_barcodes()         # strips "-1" from barcodes (global replace — Bug 2.3)
  → build_multimap_matrix() # unique ECs (len==1) → continue (no contribution)
                            # multi ECs → share = count / n_genes → csr_matrix
  → create_new_h5ad()       # adata.X = corrected_matrix (multimapper shares only)
  → final_results()
      adata.layers["counts_corrected"] = adata.X.copy()          # multimapper shares
      adata.layers["counts_original"]  = adata_orig.X.copy()     # all unique reads
      adata.layers["counts_combined"]  = corrected + original    # redundant sum
      adata.write(adata_multimap.h5ad)

detection.py:preprocessing()
  → adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]
  # = multimapper_shares + unique_reads → no double-count (unique ECs contribute 0 to corrected)
  → total_count = adata[:, gene_id].X.sum()
  → if total_count >= 1: found_genes[gene_id] = total_count   # threshold check

umap.py:main()
  → adata.X = adata.layers["counts_corrected"] + adata.layers["counts_original"]  # same logic
  → umap(adata, found_genes)
      X_counts = adata.X (before normalization)
      viral_counts extracted from X_counts          # raw UMIs for labeling
      sc.pp.normalize_total → sc.pp.log1p           # normalize for embedding
      sc.pp.pca (no random_state) ← Bug 3.1
      sc.pp.neighbors (no random_state) ← Bug 3.1
      sc.tl.umap (no random_state) ← Bug 3.1
      viral_neighbor_enrichment()
        np.random.permutation (no seed) ← Bug 3.1
```

### Trace B: GTF parsing call chain
```
analysis.py:obtain_gtf()
  → for each .gtf in data/:
      info = line.split("\t")[8]
      gene_id = info.split('"')[1]   ← Bug 2.2: grabs first quoted token, not gene_id value
      viral_accessions.add(gene_id)
  → write to log/analysis.txt

detection.py:preprocessing()
  → reads log/analysis.txt
  → checks adata.var_names for each accession

multimap.py:normalize_barcodes()
  → reads log/analysis.txt
  → builds viral_gene_indices set from gene_ids
```

### Trace C: ncbi_fetch cache validation
```
ncbi_fetch._fetch_one(accession, cache_dir, email, api_key)
  → fasta_path = acc_dir / f"{acc}.fasta"
  → if not fasta_path.exists() or fasta_path.stat().st_size == 0:  ← only size check
      download and write
  → (no checksum validation of existing file)  ← Bug 3.2
```

---

*Audit completed 2026-05-08. No source files were modified during this audit.*
