# Tripwire Audit — COVID ViralScan findings — 2026-07-15

**Mode**: Audit (no code runs; behavioral spec only)
**Source review**: `.living/outputs/reviews/2026-07-15-covid-viralscan-findings.md`
**Instrumentation detected**: 0 / 4 hooks (no checkpoint emission, no `--stop-after`, no
`analysis_labels.yml`, no drop ledger found in `covid_viralscan/`)
**Next step to run these tests**: scaffold mode — see "What's missing" section below.

---

## What tripwires apply here

This analysis is a **documentation/findings update** after a detection run, not a pipeline code
change. The applicable tripwire categories are:

- **Freshness checks**: does the findings document accurately reflect current primary artifacts?
- **Known-answer tests**: do the claimed negatives and positives hold under direct verification?
- **Fault-injection tests**: does the pipeline correctly reject/flag known-bad inputs?

Fault-injection and metamorphic tripwires (which perturb the pipeline's data flow) are deferred
until the pipeline re-run (Fix F3: dedup NC_002076.2, recommit, re-run `viralscan evidence`)
is complete, because the decisive artifact (`coverage.tsv`) is currently from an undocumented
remediation job.

---

## Tripwire T1 — Freshness: breadth figures in findings doc match coverage.tsv

**Links to**: F1 (3.41% vs 1.99% factual error), F3 (undocumented hf_align job)

**Perturbation**: read `coverage.tsv` produced by hf_align job 25181135 and extract the
`coverage` (breadth) column for KP343825.1 (x216), NC_001479.1 (both samples), and
MW455365.1 (x213).

**Expected outcome**: KP343825.1 x216 breadth = 3.41% (now stated in findings doc).
NC_001479.1 x213 breadth = 1.99107% (156/7835), x216 breadth = 1.99107% (same 156 bases).
MW455365.1 x213 breadth should be in the 2.15–2.5% range (value now implied but not
explicitly stated in the findings doc).

**Pass criterion**: findings doc figures match coverage.tsv to 2 significant figures.

**Failure signal**: if numbers diverge, one of (a) the findings doc was edited without
regenerating coverage.tsv, or (b) job 25181135 used different parameters than assumed.

**Implementation sketch**: a one-line shell command —
```bash
awk 'NR==1 || $1 ~ /KP343825|NC_001479|MW455365/' \
  covid_viralscan/results/hf_align/coverage.tsv
```
(adjust path to wherever job 25181135 wrote its output).

---

## Tripwire T2 — Freshness: pseudoalignment rates in findings doc match run_info.json

**Links to**: F11 (cited rates 6.4%/8.9% vs run_info.json 4.5%/7.1%)

**Perturbation**: read `run_info.json` from the corrected-whitelist quant run (job 25140008)
for both samples and extract `p_pseudoaligned`.

**Expected outcome**: the findings doc should cite whichever set of numbers reflects the
authoritative run. Currently the SURVEY_SUMMARY.md header states 6.4% (x213) / 8.9% (x216),
while the STARsolo run_info.json may show 4.5%/7.1%.

**Pass criterion**: a single authoritative source (run_info.json from job 25140008) is cited;
any discrepancy with earlier drafts is annotated.

**Failure signal**: two different numbers for the same sample from the same run → one was
from a different (possibly broken-whitelist) run.

**Implementation sketch**:
```bash
cat covid_viralscan/results/x213-g/run_info.json | python -c \
  "import json,sys; d=json.load(sys.stdin); print(d['p_pseudoaligned'])"
```

---

## Tripwire T3 — Known-answer: SARS-CoV-2 UMI = 0 in both samples

**Links to**: F7 (SARS-CoV-2 = 0 lacks clinical context)

**Perturbation**: re-extract the SARS-CoV-2 row from the corrected-whitelist count matrix
(job 25140008 output) for both x213-g and x216-g.

**Expected outcome**: 0 UMI in both samples at any threshold.

**Pass criterion**: sum of all barcodes' SARS-CoV-2 UMI = 0 for both samples.

**Failure signal**: any non-zero value → re-examine whether the reference panel includes
NC_045512.2 correctly and whether a new run introduced a different matrix.

**Implementation sketch**: already confirmed in SURVEY_SUMMARY.md section 3 (SARS-CoV-2
specificity control: 0/0). Re-confirmation via direct matrix query:
```bash
python - <<'EOF'
import scipy.io, numpy as np
# adjust paths to corrected-whitelist matrix dirs
for sample in ['x213-g', 'x216-g']:
    m = scipy.io.mmread(f'covid_viralscan/results/{sample}/counts_unfiltered/cells_x_genes.mtx')
    # need t2g to identify SARS-CoV-2 rows
    print(sample, 'done')
EOF
```
(requires t2g lookup for NC_045512.2 row index).

---

## Tripwire T4 — Known-answer: VARV (smallpox) UMI ≤ 1 as noise floor

**Links to**: per-virus-plausibility verdict (eradicated 1980 → any signal is noise)

**Perturbation**: extract the VARV row from the corrected-whitelist matrix post-STAR-filter
(the host-filtered UMI table).

**Expected outcome**: VARV UMI = 0–1 in both samples post-STAR-filter. Pre-filter VARV was
already very low; post-filter should be 0.

**Pass criterion**: post-STAR-filter VARV UMI ≤ 1 in both samples.

**Failure signal**: post-filter VARV > 1 → the STAR filter is not suppressing this signal,
which would undermine the filter's reliability claim for all other non-anellovirus signals.

**Implementation sketch**: extract from `per_cell_viral.tsv` (SURVEY_SUMMARY.md section 2
source) after re-running with `--host-filter starsolo`. Current evidence (1 UMI pre-filter,
likely 0 post) should be verified against the host-filtered matrix directly.

---

## Tripwire T5 — Fault injection: NC_002076.2 duplicate causes `viralscan evidence` crash

**Links to**: F3 (undocumented hf_align job; pipeline not reproducible)

**Perturbation**: run `viralscan evidence` on the CURRENT reference (before dedup fix) and
verify the crash reproduces as documented.

**Expected outcome**: `samtools sort` exits non-zero with a duplicate-header error for
NC_002076.2; `viralscan evidence` exits with a non-zero code; no `coverage.tsv` is written.

**Pass criterion**: the crash is observable and the error message matches the documented
`[E::sam_hrecs_update_hashes] Duplicate entry "NC_002076.2"` pattern.

**Post-fix tripwire**: after the NC_002076.2 dedup fix is applied, run `viralscan evidence`
again and verify (a) it exits 0 and (b) `coverage.tsv` is written with all expected contigs.

**Note**: this is the ONLY tripwire that requires a pipeline re-run, not just a read of
existing artifacts. It is the critical path to closing finding F3.

---

## Tripwire T6 — EBV per-cell resolution (unresolved signal)

**Links to**: per-virus-plausibility verdict (EBV = unresolved, 3–4 UMI post-filter)

**Perturbation**: extract the CB IDs of cells carrying ≥1 EBV (NC_007605.1) UMI from the
host-filtered matrix. Cross-reference with the host transcriptomic UMAP cell-type annotations
from CellRanger.

**Expected outcome (artifact hypothesis)**: EBV-UMI cells are distributed across all cell
types with no B cell enrichment. The 58–69 EBV bases covered cluster in the terminal repeat
region (highly repetitive, known homology trap).

**Expected outcome (genuine hypothesis)**: EBV-UMI cells are enriched in B cells (Fisher's
exact p < 0.05 after correcting for overall B cell fraction). Covered bases fall in the EBNA-2
coding region or other unique-sequence loci, not the terminal repeat.

**Pass criterion for artifact**: Fisher's OR < 2 for B cell enrichment, OR covered-base loci
= terminal repeat region (need locus-level inspection of the minimap2 BAM).

**Implementation sketch**:
```bash
# Step 1: extract EBV-UMI barcodes
python - <<'EOF'
import pandas as pd
df = pd.read_csv('covid_viralscan/results/x213-g/per_cell_viral.tsv', sep='\t')
ebv_cells = df[df['virus'].str.contains('NC_007605') & (df['UMI'] >= 1)]['barcode']
print(ebv_cells.to_list())
EOF

# Step 2: intersect with CellRanger cell-type annotations
# (CellRanger output at paths in covid_viralscan/README.md)

# Step 3: inspect which EBV bases are covered
samtools view covid_viralscan/results/hf_align/x213-g.bam NC_007605.1 | \
  awk '{print $4}' | sort -n | uniq -c  # start positions of covering reads
```

---

## What's missing for these tripwires to run

| Requirement | Status | Fix |
|---|---|---|
| Committed hf_align script + known output path for coverage.tsv | Missing | Fix F3: commit script, document path |
| NC_002076.2 dedup in viral reference | Pending | `viralscan build-ref` with dedup guard |
| `viralscan evidence` exits 0 on deduped reference | Pending | Requires fix above |
| Per-cell viral TSV from HOST-FILTERED matrix (not pre-filter) | Partial | Need to confirm job 25181135 produced per-cell data, not just coverage.tsv |
| CellRanger cell-type annotations for EBV enrichment test | Missing | Path in README.md; need to load filtered barcodes + cluster labels |
| VARV row index in t2g | Available | Look up `NC_001481.1` (VARV) row in panel t2g |

**Recommended first action**: fix T5 (the NC_002076.2 fault injection) — it unblocks T1 (freshness),
re-establishes reproducibility for F3, and produces the authoritative `coverage.tsv` that makes T1
a clean comparison against the current findings doc numbers.

---

## What was deliberately not tripled

- **D-list effectiveness (15% removal vs 95% STAR)**: the conclusion is already replicated by the
  two-STAR-method agreement; a third independent test would be the BLAST-annotation-of-panel
  approach, which is a build-time quality gate rather than a runtime tripwire.
- **STARsolo = 0 as cross-validation**: F14 explains this is a GTF artifact (no exon records for
  anellovirus), not independent evidence. No tripwire can rehabilitate it as a positive control
  without a GTF fix.
- **CeHV2 cross-mapping from HSV-1/VZV**: the mechanistic explanation (UL24/UL25 conservation)
  is well-established; testing it would require aligning CeHV2-mapped reads to HSV-1, which is
  informative but not a pipeline correctness test.
