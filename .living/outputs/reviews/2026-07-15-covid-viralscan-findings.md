# Review — COVID ViralScan findings — 2026-07-15

**Scope**: `.living/findings/covid-viralscan-no-sars2-anellovirus-dominant.md` + `covid_viralscan/results/SURVEY_SUMMARY.md`
**Focus questions**: (1) How to minimise the host-homology artefact? (2) Are the non-anellovirus detections (HHV-6, HHV-1, CeHV2, MPXV, EBV, molluscum, etc.) genuine or artefact?
**Sub-agents run**: 6 — bioinformatics, stats-causal, data-pipeline-leakage, narrative-consistency, remediation-strategies, per-virus-plausibility

---

## Key decisions in this analysis

- **Root-cause identification** — cDNA-only kallisto host reference fails to suppress intronic/intergenic reads that share non-coding GRCh38 homology with the anellovirus panel. This is the structural cause of all artifact signal.
- **Three-method convergence** — read-origin (STAR, 0/4.5M viral-primary), STAR host-filter (~95% removal), coverage breadth (≤3.4% across all contigs). Two of the three methods share STAR/GRCh38 as the underlying tool.
- **STAR host-filter as the validated fix** — `--host-filter starsolo` removes ~95% of artifact with no measured sensitivity loss on SARS-CoV-2 = 0.
- **Breadth-of-coverage as discriminator** — `viralscan evidence` (minimap2 + samtools coverage) distinguishes fixed-locus artifact (deep+narrow) from genuine infection (breadth spread). Threshold ≥5% breadth proposed but uncalibrated against a positive control.
- **Sweeping "no genuine viral infection" conclusion** — established for anellovirus; currently extended by implication to non-anellovirus signals that have not been independently subjected to the same three-method pipeline.
- **D-list rejected** — `--genome-dlist` (exact k-mer masking) shown to remove only ~15% of anellovirus artifact vs ~95% for STAR; rejected for B5 bulk re-run.

---

## Questions for the analyst

1. **Clinical metadata** — What is the COVID status, sampling timepoint, and tissue type for LUM-SJ-x213-g and x216-g? SARS-CoV-2 = 0 in PBMC is biologically expected even in active COVID-19; knowing the expectation changes how to frame the negative result.
2. **Positive control** — Was there a known-positive sample (e.g., CMV-viremic, HHV-6B-viremic, or synthetic spike-in) run through the same pipeline? Without it, the ≤5% breadth threshold and the "no real infection" conclusion are uncalibrated.
3. **EBV (EBNA-2) unresolved** — 4 UMI survive the STAR host-filter (EBV is not in GRCh38). Do these 4 cells fall in the B cell cluster from the host transcriptomic data? This is a 15-minute check.
4. **ciHHV-6 exposure risk** — Were the TINO patients known to be at risk for chromosomally integrated HHV-6B (ciHHV-6B)? The STAR filter suppresses ciHHV-6B reads (correctly, as host) — a ciHHV-6B carrier would be invisible to this pipeline.
5. **B5 bulk re-run** — Has `--genome-dlist` been replaced with `--host-filter starsolo` for the 99-sample bulk run (GSE128078)? The review confirmed d-list is inadequate; any B5 result citing anellovirus signal is currently artifact-contaminated.

---

## Findings

### Bioinformatics

#### Major

##### F1. Maximum coverage breadth is 3.41%, not 1.99% — factual error in findings document

`findings/covid-viralscan-no-sars2-anellovirus-dominant.md` (breadth verdict section)

```
"Max coverage breadth across EVERY viral contig: 1.99%"
"Breadth 1.99107% is identical to 5 decimals in both samples"
```

**Why it matters here**: The 1.99107% figure is specific to NC_001479.1 (EMCV, a picornavirus — not anellovirus). Actual maximum across all contigs is KP343825.1 (Torque teno Gammatorquevirus): 3.33% in x213, 3.41% in x216. The conclusion remains correct (3.4% is still far below the ≥10–20% breadth expected for real infection) but the stated statistic is wrong. A reviewer checking the coverage.tsv would find the discrepancy.

**Fix**: Correct to "maximum breadth 3.41% (KP343825.1, x216); most contigs ≤2.5%." Separately note that NC_001479.1 (EMCV) shows the most extreme depth-to-breadth ratio (841×/1621× depth at 1.99% breadth — a fixed homologous locus) as the clearest artifact signature.

---

##### F2. NC_001479.1 is EMCV (picornavirus), not an anellovirus — taxonomic misattribution

`findings/covid-viralscan-no-sars2-anellovirus-dominant.md` (breadth verdict section)

**Why it matters here**: The findings document uses NC_001479.1 as the anchor for the anellovirus breadth argument. EMCV (Encephalomyocarditis virus) is a cardiovirus, not related to anellovirus. Its extreme depth-at-narrow-breadth signature is strong evidence of artifact but should be attributed to a GRCh38 locus with IRES-homologous sequence (the EMCV IRES is used in expression vectors; this may represent a GRCh38 region with homology to the IRES), and discussed separately from the anellovirus finding.

**Fix**: Add "NC_001479.1 is EMCV (picornavirus), not anellovirus; its deep/narrow pile-up is an additional host-homology locus, consistent with GRCh38 IRES-region homology." Anellovirus contigs (KP343825.1, MW455365.1, etc.) should be cited separately with their own breadth values (2.15–3.41%).

---

##### F3. Coverage breadth evidence derives from an undocumented hf_align job (25181135) after pipeline crash — not reproducible from documented workflow

`covid_viralscan/logs/hf_evidence_25180994_0.err` / `covid_viralscan/results_hostfilter/*/evidence/coverage.tsv`

**Why it matters here**: `viralscan evidence` (job 25180994) crashed at `samtools sort` for both samples with `[E::sam_hrecs_update_hashes] Duplicate entry "NC_002076.2"` (exit 1). The authoritative `coverage.tsv` files were produced by a subsequent manually crafted `hf_align` job (25181135) whose script is not in `covid_viralscan/scripts/`. Anyone following the documented pipeline from RUNBOOK.md would find that `viralscan evidence` crashed and produced no BAM or coverage output — the decisive breadth evidence cannot be reproduced.

**Fix**: (1) Fix NC_002076.2 duplicate in the viral reference FASTA/GTF so `viralscan evidence` no longer crashes. (2) Commit the hf_align script to `covid_viralscan/scripts/` with a comment explaining its remediation role. (3) Add RUNBOOK.md note that job 25181135 produced the authoritative `coverage.tsv`. (4) Re-run `viralscan evidence` against the deduplicated reference so the coverage.tsv has fully traceable provenance.

---

##### F4. "No genuine viral infection" sweeps over non-anellovirus signals that lack independent coverage-breadth validation

`findings/covid-viralscan-no-sars2-anellovirus-dominant.md` (final verdict)

```
"FINAL VERDICT: No genuine viral infection is supported in either sample"
```

**Why it matters here**: The three-method pipeline was applied to anellovirus. The non-anellovirus signals (HHV-6 68→0 UMI, HHV-1 63→1–3 UMI, CeHV2 92→1–1.5 UMI, EBV 12→4 UMI, Molluscum 9→3 UMI, MPXV 31→0 UMI after STAR filter) have not been formally subjected to the same breadth analysis. The conclusion is likely correct for all of them (see F7–F10 below), but stating it as a global verdict without verifying each virus exceeds what the data formally support.

**Fix**: Either (a) explicitly run `viralscan evidence` on the host-filter survivors for the top 5 non-anellovirus entries and document their breadths (most are already in `coverage.tsv` — this is largely a documentation gap), or (b) qualify the verdict: "no genuine viral infection supported for anellovirus (three convergent methods); remaining signals collapse to <5 UMI post-STAR-filter with breadth signatures consistent with artifact, but have not been independently verified via the full three-method pipeline."

---

#### Minor

##### F5. ciHHV-6B is a genuine blind spot in the STAR host-filter strategy

**Why it matters here**: Chromosomally integrated HHV-6B (ciHHV-6B) is present in ~1% of the general population — every cell carries a copy integrated at telomeric ends of GRCh38. The STAR host-filter correctly classifies these reads as host-origin (they're in GRCh38) and removes them. This means the pipeline cannot detect ciHHV-6B, and an HHV-6B-positive result from STAR filter would require alignment to a combined GRCh38+HHV-6B reference where the viral chromosome is marked separately. The current "HHV-6B = 0 post-filter" result is therefore ambiguous between "no HHV-6B infection" and "ciHHV-6B reads correctly suppressed as host."

**Fix**: Document this limitation explicitly. If ciHHV-6 is clinically relevant for this cohort, run STAR against a combined GRCh38+HHV-6B reference and count reads that primary-align to NC_003663.2 specifically.

---

##### F6. CeHV2 (Herpes B virus) signal is alphaherpesvirus conserved-gene cross-mapping, not macaque virus detection

**Why it matters here**: CeHV2 UL24 and UL25 are the most conserved herpesvirus core genes across all alphaherpesviruses (HSV-1/2, VZV, pseudorabies). The near-identical UMI counts in both genes (23.5/22.5 UMI each) in the same cells, dropping to 1–1.5 post-STAR, and 0.12–0.15% breadth across a 150 kb genome confirm this is a cross-map from a human herpesvirus at conserved loci. Herpes B is a BSL-3/4 pathogen — any plausible-sounding detection would trigger biosafety concern. The panel should flag CeHV2 as cross-reactive with human alphaherpesviruses.

**Fix**: Add a note to the findings document. Consider adding CeHV2 to the sibling cross-mapping logic in `constants.py` (it currently covers HSV-1/HSV-2 and HHV-6A/B pairs).

---

### Statistics & causal inference

#### Major

##### F7. SARS-CoV-2 = 0 presented as "true negative" without clinical context or detection limit

`findings/covid-viralscan-no-sars2-anellovirus-dominant.md` (SARS-CoV-2 verdict)

```
"SARS-CoV-2 = 0 remains the only trustworthy result (true negative, no host homology)"
```

**Why it matters here**: "True negative" is interpretable only relative to an expected positive. PBMC is a low-viral-load compartment even in confirmed COVID-19; 5' GEX captures the 5' end of transcripts (sub-optimal for detecting SARS-CoV-2 subgenomic RNA structure). If these are convalescent samples, SARS-CoV-2 = 0 is expected regardless of past infection. Without clinical metadata (COVID status, sampling timepoint), this cannot be framed as a specificity validation — it is a methodological statement ("no reads matched NC_045512.2") not a biological one.

**Fix**: Add one sentence specifying sample type and clinical expectation. Reframe as: "SARS-CoV-2 yields 0 UMI across all pipeline modes — consistent with absence of active viremia in this sample type. Clinical interpretation depends on patient COVID status and sampling timepoint [unknown]."

---

##### F8. No positive control — the ≤5% breadth threshold is uncalibrated

**Why it matters here**: The 5% breadth threshold for "real detection" is proposed based on the observed artifact ceiling (~3.4%), not on a known-positive sample showing what a real viral detection produces in this assay. If the sensitivity of the pipeline to, say, latent HSV-1 in a ganglionic-origin blood cell is actually <5% breadth, any genuine latent herpesvirus would be excluded by this threshold. The "no real infection" conclusion requires knowing the threshold separates artifact from real — which requires a positive control.

**Fix**: Document that the threshold is provisional (calibrated only against the artifact floor, not against positive controls). Soften the conclusion to "below the threshold consistent with productive infection at detectable levels; threshold not validated against confirmed positive samples."

---

#### Minor

##### F9. "Three independent methods" — two share STAR/GRCh38, not truly independent

**Why it matters here**: The read-origin test (STAR + combined GRCh38+anellovirus genome, NH flag) and the STAR host-filter (STAR + GRCh38-only genome, unmapped reads forwarded) both depend on STAR alignment to GRCh38. A systematic STAR failure mode (e.g., a misassembled GRCh38 region) would affect both. Only the coverage breadth analysis (minimap2, viral-only reference) is genuinely independent.

**Fix**: Describe as "two converging STAR-based lines of evidence, independently confirmed by coverage breadth analysis via minimap2." 

---

##### F10. "Identical to 5 decimal places" is arithmetic artifact — 156/7835 = 1.99107% always

**Why it matters here**: The framing invites a coincidence-probability interpretation ("how likely to agree to 5 decimals by chance?") that is a category error. The breadth is exactly 156/7835 = 0.019910657... in both samples because the same 156 fixed positions are covered, not because of probabilistic coincidence. The correct interpretation is saturation: 841× depth (x213) → 1621× depth (x216) with no increase in covered bases — breadth is flat despite 2× more reads — which is the genuine evidence of a fixed homologous locus.

**Fix**: Replace "identical to 5 decimal places" with "the same 156 reference bases covered in both samples despite 2× depth difference — saturation at a fixed host-homology locus."

---

##### F11. Pseudoalignment rate discrepancy (6.4%/8.9% cited vs 4.5%/7.1% in run_info.json)

`covid_viralscan/RUNBOOK.md` vs `results_genomic/LUM-SJ-x213-g/run_info.json`

**Why it matters here**: If the 6.4%/8.9% figures came from an earlier run (wrong whitelist), citing them alongside results from the corrected run conflates two different runs' statistics. The 4.5%/7.1% from the correct-whitelist `run_info.json` should be used consistently.

**Fix**: Confirm which run produced the cited percentages; replace with figures from the authoritative job 25140008 run.

---

##### F12. SARS-CoV-1 = 0 provides no independent specificity evidence

**Why it matters here**: SARS-CoV-1 shares >80% sequence identity with SARS-CoV-2. If SARS-CoV-2 = 0 for biological reasons, SARS-CoV-1 = 0 follows automatically — this is not an independent negative control. The SARS-CoV-1 result is useful as a "cross-mapping guard" (if SARS-CoV-1 > 0 UMI, suspect cross-mapping or index contamination) but should not be cited alongside SARS-CoV-2 = 0 as independent specificity evidence.

**Fix**: Reframe: "SARS-CoV-1 = 0 serves as cross-mapping guard (not an independent control)."

---

### Data pipeline & reproducibility

#### Major

*(See F3 above — hf_align job 25181135 undocumented. This is the most critical pipeline integrity gap.)*

#### Minor

##### F13. Read-origin test covers x213-g only (5M read subsample), not x216-g

`covid_viralscan/scripts/diag_viral_read_origin.sh` — hardcodes x213-g R2, `head -n 20000000`

**Why it matters here**: x216-g has different cell count (17,209 vs 28,921) and higher post-filter anellovirus UMI (78k vs 57k). If the artifact mechanism were sample-specific, x213-g would be the lower-signal case to test. The read-origin test result (0/4.5M) is strongly directional but technically covers only one of two samples.

**Fix**: Run the NH test on x216-g R2 using the same script (~15 min SLURM job). Low-priority given the overwhelming consistency of the other evidence.

---

##### F14. STARsolo = 0 viral UMI is GTF artifact, not a biological cross-validation

**Why it matters here**: STARsolo's viral UMI = 0 is cited in some contexts as corroboration of the "no viral infection" conclusion. The actual cause is that the anellovirus GTF lacks `exon` records (only 2,292 of 4,650 gene_ids have exons), so STARsolo GeneFull counts no targets. This is a counting-schema artifact, not biological evidence.

**Fix**: Ensure STARsolo = 0 is never cited as biological confirmation of the SARS-CoV-2 = 0 result. Only cite it for cell-calling validation (cell count agreement with CellRanger).

---

### Per-virus verdicts: genuine or artefact?

#### Rank-ordered assessment (most → least artifact probability)

| Virus | Post-STAR UMI | Breadth on viral genome | Verdict |
|-------|--------------|------------------------|---------|
| VARVgp184 (smallpox) | 0 | 0 | **Certain artifact** — eradicated 1980; eliminated by d-list |
| HHV-2 (HHV2p01/p15) | 0 | 0 | **Certain artifact** — single sample, eliminated by d-list |
| Ydvgp129 (Yatapoxvirus) | 0 | 0 | **Certain artifact** — eliminated by d-list |
| MPXV_gp132/gp028 | 0 | 0.027–0.035% | **Certain artifact** — completely eliminated by d-list; mapQ=0; COVID-era predates 2022 mpox outbreak |
| HHV-6B (p23 + 6b genes) | 0 | 0.037–0.038% | **High probability artifact** — zero post-STAR; large d-list drop indicates GRCh38 homology. Note: genuine ciHHV-6B reads are correctly suppressed as host — that is not detectably wrong |
| CeHV2 (UL24/UL25) | 1–1.5 | 0.12–0.15% | **High probability artifact** — conserved core herpesvirus gene cross-mapping; macaque herpesvirus biologically impossible in routine PBMC |
| HHV-1 (gp00p39/p61) | 1–3 | 0.071–0.084% | **High probability artifact** — large d-list drop; residual likely VZV/HHV-2 cross-map at conserved locus |
| MOCVgp001 (Molluscum) | 3 (both samples) | 0.099% / 0.099% | **High probability artifact** — survives STAR but near-identical breadth (189 vs 188 bases at fixed locus); epithelial virus with no blood tropism |
| HHV-4/EBV (EBNA-2) | 4 (x213) / 3 (x216) | 0.034–0.040% | **Unresolved** — EBV not in GRCh38; 4 UMI fully survive STAR filter. Latency III EBNA-2 expression is biologically plausible in B cells. 0.034% breadth is below threshold but small-sample. Requires per-cell check. |

**Internal calibration anchor**: VARVgp184 = 1 UMI (pre-filter) in one sample. Smallpox is definitively eradicated; this establishes the noise floor at 1 UMI. Every other signal is within 1–2 orders of magnitude of a known-impossible detection. After STAR host-filter, all non-anellovirus signals collapse to 0–4 UMI with breadths 0.03–0.15%.

---

### What would distinguish real from spurious for remaining candidates

**EBV (priority 1 — do this first)**:
1. Extract the 4 cells with HHV4_EBNA-2 UMI from `per_cell_viral.tsv` (STAR host-filter run). Do these cells cluster in the B cell population in the host transcriptomic UMAP? Latency III EBV infects B cells — random distribution = artifact, B-cell enrichment = real candidate.
2. Inspect which 58–69 bases of NC_007605.1 are covered in the BAM. Terminal repeat regions in EBV are multi-copy within the genome — reads mapping there are ambiguous. EBNA-2 coding body (unique region) coverage = more specific.
3. For real Latency III: expect >10 UMI in the infected cell, co-detection of LMP1/LMP2A/EBNA-3 transcripts, B cell cluster enrichment. Current signal fails all three.

**HHV-6B (priority 3 — if clinically relevant)**:
- Build STAR genome containing GRCh38 + NC_003663.2 as a separate "viral chromosome." Primary alignments to NC_003663.2 = ciHHV-6B or active infection. Expected for ~1% of population if these patients carry ciHHV-6B. This is the only way to recover ciHHV-6B signal that is correctly suppressed by the current host-only filter.

---

## Q1 Summary — Artifact minimisation: top 3 recommendations

### Rank 1 — Promote `--host-filter starsolo` as the recommended default

Already implemented (`host_filter.py`). Empirically removes ~95% of artifact with no measured sensitivity cost on SARS-CoV-2 = 0. The only user-facing addition needed:
- `viralscan host-ref` convenience command to build the GRCh38 STAR index (~24 GB, ~75 min)
- Documentation update: recommend for all 5' GEX, 3' GEX, and ribo-depleted bulk RNA-seq
- For B5 bulk re-run (GSE128078): replace `--genome-dlist` with `--host-filter starsolo` immediately — the d-list fix is empirically shown to be ineffective

### Rank 2 — Integrate `viralscan evidence` breadth output as a standard per-run tier

Currently an optional post-hoc subcommand. Integrating coverage breadth (breadth %, mean depth, covered bases) into the main viral summary table (`viral_summary.tsv`) would make the artifact signature visible without a second pipeline invocation. Proposed tiers:
- `breadth ≥ 5%, UMI ≥ 500` → **candidate** (consistent with possible infection)
- `breadth 2–5%` → **inconclusive** (low signal or narrow transcriptional program)
- `breadth < 2%, depth > 100×` → **suspect artifact** (fixed-locus pile-up signature)

Note: thresholds should scale by genome size (5% of 2.8 kb anellovirus = 140 bp covered; 5% of 160 kb HHV-6 = 8 kb covered — very different evidential weights). Report covered-bases alongside breadth %.

### Rank 3 — BLAST-based panel annotation: flag high-homology sequences at `build-ref` time

Run BLAST of the viral panel against GRCh38 non-coding sequence (introns + intergenic) at reference build time. Sequences with ≥80% identity over ≥100 bp get a `high_homology_flag` in the panel manifest. Users see this flag in the panel documentation and know to apply `--host-filter starsolo` for flagged families. Do NOT auto-exclude flagged sequences — whole-genome exclusion removes all detection capability for that virus. Flagging preserves user agency.

---

## What was checked but is fine

- **SARS-CoV-2 = 0 as a result** — zero UMI across all pipeline modes is methodologically robust; the framing concern (F7) is about clinical interpretation, not about whether the 0 is correct.
- **STAR host-filter implementation** — the 95% removal rate is consistent with the expected mechanism; no evidence of over-removal affecting legitimate viral reads.
- **Coverage breadth calculation** — the minimap2 + samtools coverage approach is methodologically sound; the provenance issue (F3) is about documentation, not the values.
- **D-list rejection** — empirically validated (15% vs 95% removal); correctly deprecated for this use case.
- **Sibling cross-mapping logic** — the existing HSV-1/HSV-2 and HHV-6A/HHV-6B cross-mapping pairs address the most common herpesvirus ambiguity; CeHV2 is an extension candidate (F6).
- **Multiplicity** — 19 viruses tested without explicit multiple testing correction; this is appropriate for a viral survey where each detection is reported with UMI/breadth evidence rather than a p-value.

---

## Notes

**Cross-cutting**: The most important single action before citing any of these results in a manuscript is fixing the NC_002076.2 duplicate (so `viralscan evidence` doesn't crash) and committing the hf_align script (F3). The breadth evidence is the strongest pillar of the artifact verdict, and its provenance chain currently has a gap.

**Published precedent**: Viral-Track (Bost et al., *Cell Host & Microbe* 2020) uses genome-mode STAR alignment to a combined host+viral reference as its primary detection strategy — implicitly avoiding the cDNA-reference artifact. ViralScan's `--host-filter starsolo` approach achieves the same mismatch-tolerant host removal while keeping kallisto for UMI quantification. This is worth one sentence of methods comparison if the paper claims novelty over Viral-Track.

**ciHHV-6 blind spot** (F5) is a known limitation of any host-filter-based approach and should appear in the methods limitations section.

---

*Report generated by 6-agent parallel review (bioinformatics, stats-causal, data-pipeline-leakage, narrative-consistency, remediation-strategies, per-virus-plausibility) — 2026-07-15*
