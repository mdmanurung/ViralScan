# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-07-01] 64 GB references/ was untracked but NOT gitignored

**Category**: gotcha

**What happened**: During ingest of the reference set, `git check-ignore references` returned nothing — the 64 GB `references/` tree (STARsolo genome dirs) was untracked but not ignored. A stray `git add -A` / `git add references` would have tried to stage 64 GB.

**Why it matters**: Accidentally staging/committing multi-GB genome indices bloats the repo irreversibly (git history keeps them forever) and can hang or OOM the commit.

**Resolution**: Added `references/`, `starsolo_p22_6/`, `starsolo_p22_6b/` to `.gitignore`; registered the data via mycelium ingest with in-place pointer doc + provenance instead of committing bytes.

**Tags**: git, large-files, gitignore, ingest, bioinformatics, references

**mitigation_type**: convention

**structural_mitigation_candidate**: A pre-commit hook rejecting staged files > ~50 MB would structurally catch this class of error; not yet shipped.

### [2026-07-02] Controlling for a composite covariate that CONTAINS the tested feature is circular

**Category**: gotcha

**What happened**: A verification review flagged that the %mito control for MT-ND4L was circular — %mito was computed from a mito-gene set that INCLUDED MT-ND4L, so regressing MT-ND4L on %mito partly regresses it on itself (self-suppression), which would drop its significance even absent any real confound.

**Why it matters**: Composite covariates (%mito, %ribo, total counts, module scores, a cell-type signature) often contain the very feature being tested. Controlling for them then absorbs the feature's own variance and manufactures a null. The verdict is uninterpretable until the tested feature is excluded from the covariate.

**Resolution**: Recomputed %mito with MT-ND4L excluded (leave-one-out); MT-ND4L still dropped (FDR 4e-7→0.067), so the mito-QC confound is real and now non-circular. Fixed `go_enrichment.py`. General fix: exclude the tested feature from any composite covariate (leave-one-out) before adjustment.

**Tags**: confounding, circularity, covariate, scrna-seq, mito, review, statistics

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: When a covariate is a sum/score over features, build it leave-one-out per tested feature (or assert the tested feature ∉ the covariate's inputs).

### [2026-07-02] Mitochondrial genes need a %mito control before claiming them as biology

**Category**: gotcha

**What happened**: MT-ND4L was one of the 5 "depth-robust" EBV host-response genes (depth-adjusted FDR 4e-7). Adding a percent-mitochondrial covariate dropped it to FDR 0.07 — it was a mitochondrial-content/QC effect, not an EBV response. Across the full set, 217 of 518 depth-robust genes were lost once %mito was controlled.

**Why it matters**: In scRNA-seq, %mito tracks cell stress/quality and correlates with many conditions; a mitochondrial gene surviving depth adjustment can still be a QC artifact. Depth control is not enough — %mito is a separate confounder for mito-encoded genes (and for stress signatures generally).

**Resolution**: Added a %mito covariate alongside log-depth; reported only the depth-AND-mito-robust set (301 genes) and dropped MT-ND4L. See `go_enrichment.py`, finding F-001.

**Tags**: confounding, mitochondrial, scrna-seq, qc, gene-selection, causal-inference

**mitigation_type**: convention

**structural_mitigation_candidate**: Convention: "for scRNA-seq association/DE, control for %mito (and cell-cycle where relevant) in addition to depth; flag any mito-encoded hit for a %mito control." Fits alongside the depth-label convention.

### [2026-07-01] Raw-count positivity thresholds silently confound with sequencing depth

**Category**: gotcha

**What happened**: The EBV+ label (≥10 raw viral UMI) turned out to be strongly depth-dependent — EBV+ rate rose 27.5%→96.3% across host-depth quintiles, and depth alone predicted the label at AUC 0.80–0.97, beating the host-gene classifier. Class balancing + top-depth filtering did NOT fix it (they match class counts, not between-class depth distributions; the filter even widened the gap).

**Why it matters**: Any "positive/detected/expressed at ≥N raw counts" label makes the label a proxy for sequencing depth, so any depth-correlated feature looks predictive. This can turn a technical artifact into a headline biological result. It nearly did here.

**Resolution**: Ran a depth-confounder check (depth-alone AUC baseline, quintile rates, E-values); reported the headline as substantially confounded; recommended a depth-normalized label. See findings F-001 (contradicted) and F-003.

**Tags**: confounding, sequencing-depth, count-data, thresholding, label-definition, scrna-seq, causal-inference

**mitigation_type**: convention

**structural_mitigation_candidate**: A convention "define positive/detected labels depth-independently (CPM/fraction or depth-matched); always report depth-alone predictive baseline for any count-threshold label." Candidate for `.living/conventions.md`.

### [2026-07-01] Feature selection (HVG) before the CV split leaks into held-out metrics

**Category**: gotcha

**What happened**: Review of the EBV host-response classifier found `sc.pp.highly_variable_genes` fit on all 1906 cells before the per-seed train/test split. Because EBV+ cells drive variance in EBV-responsive genes, the HVG feature set is chosen using cells later used for evaluation — the held-out AUC/MCC are upward-biased. Documentation ("evaluated on held-out cells") was true for rows but not feature columns.

**Why it matters**: This is the single most common silent leak in scRNA-seq classifiers, and it inflated a headline metric that was about to go into a report. Any per-cell feature selection, normalization-parameter fit, or gene selection done before the split is suspect.

**Resolution**: Flagged as Major (F1) in `.living/outputs/reviews/2026-07-01-hostresponse-ebv-matched.md`; fix is to move HVG/gene selection inside the CV loop and re-run.

**Tags**: leakage, feature-selection, hvg, scrna-seq, cross-validation, classifier, review

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: A review-checklist item (and a possible convention) "feature selection / normalization params must be fit inside the CV split"; a test that asserts HVG is recomputed per fold. Candidate for promotion to `.living/conventions.md` if it recurs.

**Update (fixed 2026-07-01)**: After moving HVG inside the CV split, the corrected metrics differed from the pre-fix values by **less than one seed-SD on every metric** — i.e. unchanged within noise. HVG is *unsupervised* (variance-based, never sees labels), so this leak class is expected to be negligible, and it was. Two lessons: (1) don't assume "remove leak ⇒ lower number" — an unsupervised-selection leak can move the number either way or not at all; (2) more importantly, don't narrate a mechanism for a sub-SD shift (an earlier draft claimed per-fold HVG "tracks the EBV contrast better" — that was story-fitting noise; corrected per advisor). Fix the leak because it's methodologically correct, not because it changes the answer here.

**Category**: insight

**What happened**: Re-ran `_run_l2_regression` (fixed seeds) in a different conda env (bioenv, sklearn 1.7.1) than the original run to add MCC. Sensitivity/specificity/balanced-accuracy and the 15-gene stable set reproduced **exactly**; AUC differed in the 4th decimal (0.84456 → 0.84487).

**Why it matters**: AUC is rank-based on continuous `predict_proba` outputs, so it's sensitive to tiny BLAS/sklearn-version numerical differences; thresholded metrics (predictions at 0.5) are not. When reporting AUC across environments, quote ≤3 decimals or pin the env, or the "same" analysis looks irreproducible at high precision.

**Resolution**: Reported AUC as 0.845 (3 dp); documented the reproduction in the analysis provenance. Adopted the re-run outputs wholesale so all metrics come from one consistent run.

**Tags**: reproducibility, auc, sklearn, numerical, metrics, bioinformatics

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Pin sklearn/BLAS versions in the analysis env and record them in ENVIRONMENTS_INSTALLATIONS.md; assert AUC reproduces to 3 dp in a regression test.

### [2026-07-01] ACTIVE_CONVENTIONS.yaml is malformed after install_convention.py

**Category**: gotcha

**What happened**: After installing convention packs, `.living/conventions/ACTIVE_CONVENTIONS.yaml` contains `active_conventions: []` followed by dangling list entries. As YAML, `active_conventions` parses as an empty list and the entries below are orphaned — so a tool that reads `active_conventions` to see what's installed would see nothing.

**Why it matters**: The ingest protocol step 3 says "read ACTIVE_CONVENTIONS.yaml to see what's installed"; if parsed literally it reports no active conventions, so domain validation (bioinformatics) could be silently skipped.

**Resolution**: Worked around by reading the file's list entries directly (bioinformatics + skill-bridge are installed). Did not patch the mycelium-generated file. Consider filing a mycelium convention-gap issue.

**Tags**: mycelium, tooling, yaml, conventions, bug

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: install_convention.py should append entries under the `active_conventions:` key (or replace the `[]`), and a yaml.safe_load round-trip assertion in its test suite would catch the malformed output.

### [2026-07-02] kb ref / kallisto index have two silent FASTA-vs-GTF contracts

**Category**: gotcha

**What happened**: The covid_viralscan Stage 2 build hit two distinct failures on the
same combined host+viral reference. (1) `kb ref` **hung for 3.5 h at 0 bytes** because
`combined.gtf` carried Ensembl *chromosomal* seqnames (`1`, `2`, `X`) while the cDNA FASTA
headers are ENST transcript IDs — ngs_tools' genome-split step scanned 1.4 GB looking for
chromosome sequences that don't exist. (2) After switching to a cDNA-level GTF (seqname =
ENST) and letting the 2 h cDNA extraction finish, the final `kallisto index` **aborted**
with `Error: repeated name in FASTA file` because the source panel `viral_genome.fa`
contained accession `NC_002076.2` (Torque teno virus 1) twice (byte-identical, present in
both the anellovirus and Serratus sets), so ngs_tools extracted its 4 gene models once per
copy → 4 duplicate names.

**Why it matters**: `kb ref` exits 0 and logs "complete" even when the kallisto index step
failed with a 0-byte `index.idx` — you MUST verify by artifact (`kallisto inspect`, non-zero
size), never by exit code. Both failure modes are silent contracts: GTF seqnames must match
FASTA headers, and every target name in the FASTA must be unique.

**Resolution**: (1) generate a cDNA-level GTF via `gen_combined_cdna_gtf.py`; (2) keep-first
dedup of `cdna.fa`/`t2g.txt`/`combined.fa` (removed records byte-identical → equals a clean
build), resume `kallisto index` directly on deduped inputs, and add a dedup guard (Step 2.6)
to `slurm_build_ref.sh` so any future duplicate accession is handled automatically. Related:
the bulk `panel.idx` (P23) was separately found host-less because its June 24 build died at
the Ensembl `current_gtf` 404 (fixed 2026-07-01) — same "verify the artifact, not the log"
lesson.

**Tags**: kb-python, kallisto, kb-ref, ngs_tools, bioinformatics, reference-build, gotcha, verify-by-artifact
