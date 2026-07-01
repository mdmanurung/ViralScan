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

### [2026-07-01] Thresholded metrics reproduce bit-for-bit; AUC only to 3 decimals across envs

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
