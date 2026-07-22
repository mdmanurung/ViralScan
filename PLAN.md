# ViralScan 3.0 execution tracker

Status: **active**

Branch: `codex/viralscan-v3`

Last reconciled: 2026-07-22

Release target: `3.0.0rc1`, then `3.0.0`

Scientific target: methods-grade manuscript generated only from frozen v3 artifacts

This is the single operational tracker. The accepted scientific and product
contracts are in
[`docs/plans/2026-07-22-viralscan-v3-correctness-and-publication.md`](docs/plans/2026-07-22-viralscan-v3-correctness-and-publication.md).
The former pre-v3 tracker is retained only in the local governance archive and
excluded from Git and source distributions because it contains institutional
paths and private-analysis notes. It must not be used to establish v3
completion.

## Next action

**Do `SCI-03`: freeze whole-sample training/holdout partitions, the continuous
cell score and tier-calibration rule, numeric factors/seeds, LOD estimation, and
biological-sample uncertainty.** Estimated effort: 6-10 focused hours.
In parallel, `GOV-03` and `REL-01` are ready and do not depend on scientific
results.

### Do now

1. `SCI-03` — freeze partitions, calibration, metrics, LOD, and uncertainty
   without viewing holdout outcomes.
2. `GOV-03` — inventory claim-bearing artifacts with versions and SHA-256.
3. `REL-01` — make the pip tier installable without pretending it provides the
   full external-tool workflow.

## How to use this tracker

1. Take the first unchecked item whose dependencies are all `[x]`.
2. Add or update its test before changing production behavior.
3. Run the stated acceptance command and save the named artifact.
4. Mark `[x]` only with dated evidence; use `[~]` for partial and `[!]` for an
   external blocker.
5. Update **Next action**, the work-package row, and the claim registry in the
   same commit.

Status legend: `[x]` verified; `[~]` partial; `[ ]` ready/not started; `[!]`
blocked on a person, credential, private datum, external service, or unavailable
compute.

## Critical path

| Order | Work package | Status | Depends on | Exit gate |
|---:|---|:---:|---|---|
| 0 | Governance and legacy freeze | `[~]` | none | `G0` |
| 1 | Core software contracts | `[~]` | `G0` | `G1` |
| 2 | Install, lock, and artifact parity | `[~]` | `G1` | `G2` |
| 3 | Preregister scientific validation | `[~]` | `G0` | `G3` |
| 4 | Freeze production references | `[~]` | `G3` | `G4` |
| 5 | Build and validate the truth panel | `[ ]` | `G3`, `G4` | `G5a` |
| 6 | Run public positives and comparators | `[ ]` | `G2`, `G4`, `G5a` | `G5b` |
| 7 | Score, audit, and freeze results | `[ ]` | `G5b` | `G5` |
| 8 | Regenerate docs and claims | `[~]` | `G1`; numbers require `G5` | `G6a` |
| 9 | Regenerate and review manuscript | `[ ]` | `G5`, `G6a` | `G6` |
| 10 | Publish and test `3.0.0rc1` | `[ ]` | `G2`, `G5`, `G6a` | `G7` |
| 11 | Publish stable `3.0.0` and submit | `[ ]` | `G6`, `G7` | `G8` |
| 12 | Observe and maintain | `[ ]` | `G8` | ongoing |

Work packages 1-4 may proceed in parallel. Do not start truth-panel outcome
analysis before `G3`, regenerate quantitative documentation before `G5`, or tag
a release before its gate.

## Verified baseline

- [x] `BASE-01` — v3 branch and locked correctness roadmap exist.
- [x] `BASE-02` — molecule-safe streaming allocation, distinct-gene EC
  projection, five explicit methods, contracted H5AD layers, and exact mass
  audits pass synthetic/property tests.
- [x] `BASE-03` — retained EBV input completed with 103,145,071 BUS records,
  57,957,364 input molecules, and exact partition/matrix conservation. Evidence:
  [`analysis/v3_ebv_baseline/README.md`](analysis/v3_ebv_baseline/README.md).
- [x] `BASE-04` — exact-target evidence completes extraction, competitive
  alignment/BLAST, deduplication, QC, coverage, BAM indexing, and IGV on the tiny
  real-tool fixture.
- [x] `BASE-05` — the 2026-07-22 local gate passed 720 tests with 21 deselected;
  Ruff check/format, data governance, and protocol validation also passed. Rerun
  the full gate after each remaining software slice rather than treating this
  development-branch result as release evidence.

## WP0 — Governance and legacy freeze

Objective: prevent unsafe pre-v3 results, private data, or unsupported claims
from entering the release. Estimated remaining effort: 1-2 days.

- [x] `GOV-01` — mark pre-v3 corrected/combined counts scientifically
  incompatible and ineligible for v3 claims.
- [x] `GOV-02` — add a ship-scope data-governance check and exclude private
  clinical outputs from package/container/source-distribution scope.
- [ ] `GOV-03` — inventory every claim-bearing input, reference, intermediate,
  result, and scheduler record with software version and SHA-256. Write
  `analysis/v3_artifact_inventory.tsv`.
- [~] `GOV-04` — define the public ship-scope allowlist, then remove or quarantine
  institutional defaults from every included script and manifest.
- [~] `GOV-05` — expand `claims/registry.json` into a validated claim graph with
  source location, artifact digest, Git SHA, input/reference hashes, schema,
  layer, denominator, generation command, scope, and status.

`G0` passes when the governance scan is green, every public quantitative claim
is registered, all pre-v3 quantitative claims are rejected or historical, and
the inventory identifies enough retained BUS/reference/FASTQ material to rebuild
each eligible result.

Evidence to record:

```text
python3 scripts/check_data_governance.py
python3 scripts/validate_claim_registry.py --coverage
sha256sum analysis/v3_artifact_inventory.tsv claims/registry.json
```

## WP1 — Core software contracts

Objective: close remaining correctness and workflow-consistency gaps before
scientific-scale execution. Estimated remaining effort: 4-7 engineering days.

### WP1A — Schemas, assignments, and reruns

- [x] `SW-01` — move/package all v3 schemas inside the installed `viralscan`
  distribution, load them with `importlib.resources`, and fail closed when a
  required schema is missing. Add wheel and sdist tests.
- [ ] `SW-02` — enforce every public JSON/TSV/H5AD v3 schema at write and
  `validate-run` boundaries; remove generic silent skips.
- [ ] `SW-03` — add optional compressed molecule-assignment evidence containing
  CB, UMI, ECs, distinct genes, ambiguity class, method, weights, and exclusion
  reason without changing default matrix mass.
- [ ] `SW-04` — make `rerun-multimap` regenerate every method-dependent artifact
  in a new result tree: matrix/layers, count audit, summaries, evidence tiers,
  UMAPs, and host-response inputs.
- [ ] `SW-05` — add an integration test proving no stale artifact survives a
  method change and the source result remains untouched.

### WP1B — Workflow safety and architecture

- [~] `SW-06` — `doctor`, `validate-run`, fingerprints, resume/overwrite safety,
  and atomic manifests exist; finish whole-workflow staging and an atomic
  completion marker.
- [~] `SW-07` — exact-fragment STAR filtering and mate synchronization exist;
  emit a reason for every retained/removed fragment plus lost-truth and
  host-virus-ambiguous boundary counts.
- [x] `SW-08` — remove unsafe kallisto CB-UMI-wide host filtering from the stable
  CLI because exact fragment identifiers are unavailable.
- [ ] `SW-09` — split the oversized CLI into thin parsers plus importable service
  functions; convert Snakemake scripts to minimal wrappers without changing
  outputs.
- [ ] `SW-10` — run one tiny paired-end fixture through documented CLI commands:
  preflight, reference, quantification, molecule allocation, cell calling,
  summaries, evidence, BAM/BLAST/plots/IGV, and `validate-run`.
- [ ] `SW-11` — make production cell calling fail closed: caller exceptions,
  zero-match external lists, invalid barcode geometry, and canonical collisions
  must never silently turn every barcode into a cell; `none` remains explicit.

`G1` passes when all count invariants, schema checks, safety scenarios, rerun
consistency, and the full tiny workflow are green. No known correctness or data-
loss defect may remain.

Required gate:

```text
NUMBA_CACHE_DIR=/tmp/viralscan-numba-cache PYTHONPATH=src python3 -m pytest tests/ -q
python3 -m ruff check .
python3 -m ruff format --check .
PYTHONPATH=src python3 -m pytest -m "integration and not network" -q
```

## WP2 — Install, lock, and artifact parity

Objective: make the pip tier honest and make conda, OCI, and Apptainer execute
the same tested build. Estimated remaining effort: 5-8 engineering days plus
runner time.

### WP2A — Pip and locked environment

- [ ] `REL-01` — move full-workflow-only dependencies such as Snakemake out of
  mandatory pip runtime dependencies; make `doctor --profile pip` check only
  Python/API/reporting/validation capabilities.
- [ ] `REL-02` — define Python 3.11 on `linux-64` as the canonical full-workflow
  toolchain; advertise other platforms only after the same workflow passes.
- [ ] `REL-03` — generate and commit reproducible runtime/development lockfiles
  and an exact external-tool version manifest.
- [~] `REL-04` — build wheel and sdist once, run `twine check`, inspect packaged
  assets, install each in a clean environment, and write `SHA256SUMS`.
- [ ] `REL-05` — run `doctor --profile pip` and `validate-run` against a packaged
  v3 fixture from both clean installations.

### WP2B — Containers and parity

- [~] `REL-06` — Miniforge is digest-pinned; rebuild Docker from the committed
  lock and the exact tested wheel rather than from a mutable source install.
- [ ] `REL-07` — publish a versioned OCI candidate, record its digest and tool
  manifest, and never assign `latest` to an RC.
- [ ] `REL-08` — build Apptainer/SIF from that OCI digest rather than performing
  a separate dependency solve; record the SIF SHA-256.
- [ ] `REL-09` — run the identical tiny workflow in locked conda, Docker, and
  Apptainer and compare validated molecule counts and deterministic output hashes.
- [ ] `REL-10` — save `analysis/release_parity/parity_report.json` with commands,
  versions, digests, hashes, and explained nondeterministic files.

### WP2C — Supply chain and release workflow

- [~] `REL-11` — security CI exists; audit the locked product environment rather
  than the scanner job, add OCI scanning, a full SBOM, dependency/data licence
  report, and reviewed vulnerability exceptions.
- [ ] `REL-12` — pin GitHub Actions by commit SHA and generate provenance/
  attestations for wheel, sdist, OCI, SIF, locks, and checksums.
- [ ] `REL-13` — make release publication depend on green CI for the exact tagged
  SHA; build distributions once and make all downstream jobs consume them.
- [ ] `REL-14` — add protected TestPyPI/PyPI environments, stable approval, GitHub
  release assets, and separate RC/stable container tag behavior.
- [!] `REL-15` — maintainers must provide final author, licence-holder,
  maintainer/contact, ORCID, affiliation, support-window, and trusted-publisher
  metadata.

`G2` passes when a clean lock-created environment plus Docker and Apptainer run
the same tiny workflow with identical validated counts; installable schemas and
assets are present; security/licence reports contain no unresolved actionable
finding; and the exact-SHA release workflow is dry-run verified without
publishing.

## WP3 — Preregister scientific validation

Objective: freeze outcome-independent decisions before producing new v3
scientific results. Estimated effort: 1-2 days plus reviewer sign-off.

- [x] `SCI-01` — create schema-valid `analysis/v3_validation/protocol.yaml` and
  `README.md` with hypotheses, supported scope, datasets, factors, seeds,
  chemistry, input/reference hashes, and exclusions.
- [x] `SCI-02` — freeze the barcode universe/cell-calling rule, feature
  intersection, primary unique-only cross-tool layer, secondary ambiguity-aware
  layers, denominators, and failure-reporting policy.
  Frozen subsection SHA-256:
  `699854b71169222e74d26c2119f1c9d7742d5962ac4ed8c02dfa9f288e3339dc`;
  dependent hypotheses/endpoints/exclusions digest:
  `fd818661a8d06f28c5c23e78bfb1c5b48b151023f5f749453260d023adbff333`.
- [ ] `SCI-03` — freeze training/holdout partitions, evidence-tier calibration
  rules, metrics, limit-of-detection method, and biological-sample bootstrap unit.
- [ ] `SCI-04` — enumerate every ViralScan, STARsolo, traditional alignment,
  Venus, and Viral-Track row with exact environment/reference requirements.
- [ ] `SCI-05` — obtain an independent protocol review, resolve findings, then
  record the protocol SHA-256 and Git SHA before outcome-generating runs.

`G3` passes when the protocol validates against its schema, has no unresolved
review finding, is hashed, and outcome-generating jobs have not preceded its
freeze commit.

## WP4 — Freeze production references

Objective: make reference contents reproducible and calibrate host-homology
safeguards without holdout leakage. Estimated effort: 3-5 engineering days plus
about 8 cluster hours per full GRCh38 build.

### WP4A — Profiles and provenance

- [~] `REF-01` — profile names and opt-in expanded anellovirus behavior exist;
  freeze accession lists for `curated`, `broad-discovery`,
  `anellovirus-representative`, and `anellovirus-expanded`.
- [~] `REF-02` — fail-closed fetches and manifests exist; complete accession
  version, taxonomy, snapshot, retrieval date, SHA-256, length, licence, cluster,
  representative status, rationale, and missing-accession fields.
- [ ] `REF-03` — apply identical masking, duplicate-ID/sequence validation, and
  manifest generation to dedicated and combined build paths.
- [ ] `REF-04` — make frozen inputs rebuild byte-identical panel FASTA/GTF/t2g
  contents and save a reproducibility audit.
- [ ] `REF-05` — replace vague source-data licence text with reviewed terms for
  every redistributed or fetched reference source.

### WP4B — Full-genome competition and anellovirus

- [~] `REF-06` — the public `--genome-dlist` path and raw annotations exist;
  build the production curated human-plus-virus index with a checksum-pinned
  GRCh38 D-list.
- [ ] `REF-07` — save the reference manifest, host-homology/low-complexity table,
  index/t2g/GTF, commands, versions, checksums, and build resource accounting.
- [ ] `REF-08` — calibrate homology/complexity exclusion thresholds using only
  preregistered training controls and freeze `thresholds.json`.
- [ ] `REF-09` — prove planted human-homology reads cannot reach probable/strong
  evidence on holdout; retain raw measurements and all excluded calls.
- [!] `REF-10` — an orthogonally confirmed anellovirus-positive sample is needed
  for real sensitivity claims; without it, ship screening support only and label
  every result accordingly.

`G4` passes when panel contents reproduce byte-for-byte, manifests validate,
duplicates are absent, GRCh38 competition is operational, and holdout host-
homology negatives cannot become probable/strong calls.

## WP5 — Build and validate the truth panel

Objective: create deterministic read/molecule truth across supported chemistry
and ambiguity regimes. Estimated effort: 1-2 engineering weeks plus compute.

### WP5A — Generator

- [ ] `VAL-01` — implement a seeded generator spanning viral abundance, infected-
  cell fraction, host homology, sibling viruses, low complexity, PCR duplication,
  CB/UMI collisions, ambient/index hopping, and 10x v2/v3/Drop-seq geometry.
- [ ] `VAL-02` — emit paired FASTQs, `truth_manifest.tsv`, read/molecule truth
  tables, barcode/chemistry metadata, input hashes, and a run manifest.
- [ ] `VAL-03` — add synthetic host-only and adversarial GRCh38-homology
  conditions plus reagent/empty-droplet controls when available.
- [ ] `VAL-04` — add a checksum-pinned 10x healthy-donor PBMC v3 presumed-negative
  observational control; never label it absolute ground truth.
- [ ] `VAL-05` — add planted target and sibling-virus positives for EBV, HHV-6,
  HSV-1/2, KSHV, and a separately scored anellovirus panel.

### WP5B — Scorer and tiny gate

- [ ] `VAL-06` — implement molecule/cell precision, recall, F1, AUPRC, sibling
  confusion, host-homology false positives, burden concordance, calibration, and
  limit-of-detection scoring.
- [ ] `VAL-07` — validate the scorer against hand-computed fixtures and reject
  denominator, feature, barcode, or count-layer mismatch.
- [ ] `VAL-08` — run the golden tiny panel end to end and prove exact planted-
  molecule recovery, count conservation, determinism, row-order invariance, and
  chunk-size invariance.
- [ ] `VAL-09` — freeze generator/scorer version, seeds, manifests, and expected
  tiny outputs before cluster-scale execution.
- [ ] `VAL-10` — run the full training and untouched holdout panels, retaining
  failures rather than silently dropping conditions.

`G5a` passes when the golden tiny panel is exact and deterministic, the full
panel is manifest-complete, and every condition has a result or reproducible
failure record.

## WP6 — Run public positives and comparators

Objective: regenerate all eligible evidence with v3 and make cross-tool
comparisons denominator-, annotation-, barcode-, and layer-matched. Estimated
elapsed time: 2-4 weeks, dominated by downloads, queues, and comparator setup.

### WP6A — Public v3 reruns

- [ ] `RUN-01` — HHV-6B `SRR20710641` under measured 10x v3 geometry.
- [ ] `RUN-02` — EBV `SRR12682296` under measured 10x v2 geometry.
- [ ] `RUN-03` — HSV-1 `SRR8315713` under measured Drop-seq geometry.
- [ ] `RUN-04` — KSHV latent/lytic series `GSE190558`.
- [ ] `RUN-05` — KSHV plus EBV replicates `GSE154900`.

Each sample must produce a validated output tree plus input/reference hashes,
command manifest, tool versions, scheduler accounting, and failure log. Expected
virus recovery is contextual evidence, not ground truth.

### WP6B — Harmonized workflow matrix

- [ ] `CMP-00` — replace historical output-derived barcode unions/intersections
  with the frozen pre-outcome anchor/feature manifests and audited count-layer
  adapters; use `counts_unique`, retain structural zeros only from declared
  complete barcode domains, and record every missing/extra mapping.
- [ ] `CMP-01` — ViralScan combined and exact-fragment STAR two-step; keep
  kallisto two-step excluded unless exact fragment lineage becomes available.
- [ ] `CMP-02` — STARsolo combined and STAR host-filter/two-step.
- [ ] `CMP-03` — traditional host-genome subtraction followed by viral alignment,
  with and without CB/UMI retention.
- [ ] `CMP-04` — Venus in an isolated version/digest-pinned environment with no
  silent algorithm patch.
- [ ] `CMP-05` — Viral-Track in an isolated version/digest-pinned environment with
  no silent algorithm patch.

Primary comparison rules: identical FASTQ hashes and viral sequences, same host
release, outcome-independent cell anchor, audited feature intersection, and
unique-only molecule parity. Recommended ambiguity-aware outputs are a separate
secondary analysis. Run at least EBV, HHV-6B, and HSV-1 for both dedicated
comparators.

`G5b` passes when every preregistered row is complete or transparently failed,
raw outputs and commands are retained, and an independent audit finds no
outcome-selected barcodes or mismatched denominator, annotation, feature, or
count layer.

## WP7 — Score, audit, and freeze results

Objective: turn completed runs into a single immutable scientific result bundle.
Estimated effort: 3-5 days after all runs finish.

- [ ] `RES-01` — compute all preregistered molecule/cell metrics, runtime, peak
  RAM, temporary storage, output size, calibration, and failure rates.
- [ ] `RES-02` — estimate uncertainty by bootstrapping biological samples, never
  by treating cells as independent experimental replicates.
- [ ] `RES-03` — calibrate thresholds on training data, freeze them, then evaluate
  the untouched holdout once; record every deviation.
- [ ] `RES-04` — independently audit inputs, references, barcodes, features,
  denominators, layers, manifests, scheduler records, and failed rows.
- [ ] `RES-05` — freeze `analysis/v3_results_bundle/` with a manifest containing
  every file digest, generation command, environment/tool version, Git SHA, and
  validation report.

`G5` passes only with 100% count-invariant compliance, complete ambiguous/
unresolved reporting, no probable/strong call in synthetic host-only or planted
host-homology negatives, positive recovery with sample-level uncertainty, and no
general superiority claim based on one virus or sample.

## WP8 — Regenerate documentation and claim evidence

Objective: make public instructions and claims match the actual v3 CLI, schemas,
and frozen results. Estimated effort: 4-7 days; quantitative pages wait for `G5`.

### WP8A — User documentation

- [~] `DOC-01` — reconcile README, installation, quickstart, CLI, outputs, API,
  FAQ, reference-panel, support, security, and migration docs with v3 contracts.
- [~] `DOC-02` — generate CLI/default tables from the parser and test exact
  defaults; remove plain `em`, removed primary-call modes, silent knee fallback,
  v2.5 container commands, and raw-UMI language for fractional estimates.
- [~] `DOC-03` — clearly label every output as observation, model estimate,
  evidence tier, diagnostic flag, or biological interpretation.
- [~] `DOC-04` — document combined versus two-step information loss, anellovirus
  screening limits, pip/full-workflow tiers, and legacy rebuild-only migration.
- [ ] `DOC-05` — execute clean-install quickstart commands and validate the
  resulting run using only documented steps.

### WP8B — Vignettes and claims

- [~] `DOC-06` — eight notebooks exist but contain pre-v3 calls/values; rebuild
  them with negative and ambiguous examples and only v3 APIs/artifacts.
- [ ] `DOC-07` — execute six lightweight notebooks in CI and the reference/full-
  workflow notebooks in the locked scheduled workflow; save logs and hashes.
- [ ] `DOC-08` — add a balanced five-workflow pros/cons table generated from the
  harmonized benchmark rather than rhetorical claims.
- [ ] `DOC-09` — add a claim-registry schema, validator, stale-hash detection, and
  coverage check for README, docs, manuscript, tables, figures, and captions.
- [ ] `DOC-10` — require each quantitative, comparative, validated, performance,
  specificity, and installation claim to resolve to a `validated_v3` artifact.

`G6a` passes when Sphinx warnings are errors, links and shell examples pass,
all notebooks execute in their declared tier, parser/default parity is exact,
and claim coverage contains no missing or stale artifact.

Documentation gate:

```text
python3 -m sphinx -W -b html docs docs/_build/html
python3 -m sphinx -W -b linkcheck docs docs/_build/linkcheck
PYTHONPATH=src python3 -m pytest tests/test_docs_consistency.py -q
```

## WP9 — Regenerate and review the manuscript

Objective: create a submission package exclusively from the frozen v3 bundle.
Estimated effort: 1-2 writing weeks after `G5`, excluding author review.

- [ ] `MS-01` — replace the historical draft and hard-coded legacy figure script
  with generators for registered values, tables, figures, captions, and
  supplement sourced only from `analysis/v3_results_bundle/`.
- [~] `MS-02` — frame molecule-aware host-virus ambiguity, combined/two-step
  evaluation, read-level specificity/QC, reference provenance, and honest
  evidence tiers; remove every legacy count and unsupported superiority claim.
- [ ] `MS-03` — rerun EBV host response using v3 labels with depth and
  mitochondrial controls; retain only if it replicates. Keep TTV solely as a
  host-homology case study absent orthogonal validation.
- [!] `MS-04` — authors must supply author order, affiliations, ORCIDs, CRediT,
  lead contact, funding, conflicts, acknowledgments, and ethics/consent/data-
  access text. Exclude private COVID libraries unless all requirements are met.
- [ ] `MS-05` — verify every citation, run independent methods/statistics and
  claim-to-artifact reviews, resolve all findings, and generate the final
  submission/availability package.

`G6` passes when every number/table/figure resolves to a frozen artifact digest,
there are no placeholders or ineligible private results, citations support their
sentences, and independent review has no unresolved correctness, denominator,
annotation, uncertainty, ethics, or unsupported-claim finding.

## WP10 — Publish and test `3.0.0rc1`

Objective: expose the exact candidate artifacts to real users without promoting
them as stable. Estimated release work: 1-2 days; testing window: 2-4 weeks.

- [ ] `RC-01` — freeze `3.0.0rc1` version, changelog, migration notes, CFF, locks,
  reference manifest, checksums, SBOM, licence report, and release notes on one
  green protected-main SHA.
- [!] `RC-02` — maintainers configure/approve TestPyPI or PyPI trusted publishing,
  GHCR/GitHub release permissions, and public artifact hosting.
- [ ] `RC-03` — tag `v3.0.0rc1` only after exact-SHA gates; publish prerelease
  wheel/sdist, versioned OCI, SIF, locks, checksums, attestations, reference
  archive, and GitHub prerelease. Do not update `latest`.
- [!] `RC-04` — recruit at least three external laboratories covering supported
  chemistries; give each a tester packet for clean install, tiny workflow,
  `validate-run`, and one supported real workflow.
- [ ] `RC-05` — track every correctness, data-loss, installation, documentation,
  and reproducibility failure to closure; rerun all gates after fixes.

`G7` passes when all three testers submit complete manifests, each supported
path works from public artifacts, and no release-blocking issue remains.

## WP11 — Publish stable `3.0.0` and submit

Objective: publish one exact reviewed build everywhere, archive it, then submit
the methods manuscript. Estimated release work: 1-3 days after approvals.

- [ ] `STB-01` — freeze final truth-panel/comparator/manuscript bundles and prove
  the release candidate's fixes did not change scientific outputs unexpectedly.
- [!] `STB-02` — reserve the Zenodo software DOI and benchmark/archive DOI, keep
  them distinct from the reference-data DOI, and provide them for metadata.
- [ ] `STB-03` — update CFF/metadata, pass exact-SHA gates, tag `v3.0.0`, and
  publish wheel, sdist, versioned plus `latest` OCI, SIF, locks, reference archive,
  release notes, migration guide, checksums, SBOM, licences, and attestations.
- [ ] `STB-04` — after the real PyPI sdist exists, insert its SHA-256 into the
  Bioconda recipe, lint/build/test it, submit the PR, and verify the installed
  package on the tiny workflow.
- [ ] `STB-05` — archive software and benchmark artifacts, verify every public
  DOI/link/download, update `CITATION.cff`, then submit the manuscript.

`G8` passes when the tagged SHA is the reviewed green SHA, public distribution
formats reproduce validated counts, all links and DOIs resolve, the software and
benchmark archives are durable, and the manuscript has been submitted. Journal
acceptance timing is external and is not a software completion condition.

## WP12 — Observe and maintain

- [~] `OPS-01` — publish supported-version, patch-release, deprecation, and
  security-response policies with named contact routes.
- [ ] `OPS-02` — establish reference-update cadence and manifest compatibility
  rules; never silently change a frozen reference under an existing version.
- [ ] `OPS-03` — triage false positives, false negatives, data loss, and reference
  drift as scientific incidents with reproducible packets.
- [ ] `OPS-04` — monitor install/usage failures during the first 90 days and ship
  patch releases from the same gate process.
- [ ] `OPS-05` — schedule independent truth-panel refreshes without changing v3
  thresholds retrospectively.

## Gate dashboard

| Gate | State | Required proof |
|---|:---:|---|
| `G0` governance | `[~]` | inventory, ship-scope scan, complete claim graph |
| `G1` software | `[~]` | full unit/property/safety/tiny-workflow suite |
| `G2` distribution | `[~]` | clean installs, locks, Docker/Apptainer parity, supply-chain reports |
| `G3` preregistration | `[ ]` | reviewed, schema-valid, hashed protocol frozen before outcomes |
| `G4` references | `[~]` | byte-rebuild, GRCh38 D-list, calibrated holdout safeguard |
| `G5` science | `[ ]` | truth holdout, public positives, all comparators, uncertainty, frozen bundle |
| `G6a` docs | `[~]` | parser parity, executable docs, complete validated claim graph |
| `G6` manuscript | `[ ]` | artifact-generated manuscript plus independent audits |
| `G7` release candidate | `[ ]` | public RC artifacts and three-laboratory closure |
| `G8` stable/publication | `[ ]` | stable artifacts, archives/DOIs, Bioconda, manuscript submission |

## External inputs and authority

These items cannot be invented or completed by an implementation agent:

| Input | Needed by | Owner action |
|---|---|---|
| Author/ethics/funding/conflict metadata | `REL-15`, `MS-04` | authors approve final factual text |
| Publishing credentials and protected environments | `RC-02` | repository owner configures services |
| Three external laboratories | `RC-04` | maintainers recruit and coordinate testers |
| Zenodo DOI reservations | `STB-02` | archive owner reserves distinct DOIs |
| Orthogonal anellovirus-positive sample | `REF-10` | collaborator supplies lawful validated data, or claim stays screening-only |

## Stop rules

- Never tag or publish from a dirty, unreviewed, or non-green SHA.
- Never migrate or reinterpret a pre-v3 H5AD value as a v3 molecule count.
- Never select shared barcodes using ViralScan-positive outcomes.
- Never turn ambiguity/QC flags into biological conclusions automatically.
- Never claim real-anellovirus sensitivity, general superiority, or formal
  infection from a nonzero molecule without the required validation evidence.

## Evidence log

Append one line after each completed item:

```text
YYYY-MM-DD ITEM — command/result; artifact path(s); Git SHA; reviewer if required
```

- 2026-07-22 `BASE-03` — EBV v3 baseline conserved exact molecule mass;
  `analysis/v3_ebv_baseline/host_conservative.json`; SLURM 25316336.
- 2026-07-22 `BASE-04` — real-tool evidence integration passed;
  `tests/integration/test_exact_lineage.py` and
  `tests/integration/test_evidence_chain.py`; Git SHA `a146050`.
- 2026-07-22 `BASE-05` — full non-network suite: 720 passed, 21 deselected;
  Ruff check/format, data-governance check, and draft protocol validation passed.
- 2026-07-22 `SW-01` — wheel and sdist contain all six byte-matched v3 schemas;
  fresh wheel install loaded 6/6, missing-schema paths fail closed; 8 focused
  tests passed; Git SHA `6e1fe66`.
- 2026-07-22 `REL-04` partial — isolated wheel/sdist build and fresh-wheel schema
  smoke test passed; `twine check`, sdist clean-install test, and `SHA256SUMS`
  remain.
- 2026-07-22 `DOC-01`–`DOC-04` partial — public docs and tracked manuscript
  placeholder reject legacy quantitative/private claims and use v3 molecule and
  candidate-evidence terminology; 94 docs/CLI tests passed; Git SHA `3bb7d1b`.
- 2026-07-22 `SCI-01` — schema-valid, non-executable protocol draft with separate
  training/holdout gates; 16 focused validation tests and 4 docs-consistency
  tests passed; data-governance check passed; independent implementation review
  resolved four blockers and passed; `analysis/v3_validation/protocol.yaml`,
  `schemas/v3/validation_protocol.schema.json`, and
  `scripts/validate_v3_protocol.py`; Git SHA `d941261`.
- 2026-07-22 `SCI-02` — outcome-independent harmonization frozen with 24
  endpoint-mapped denominator contracts, canonical and dependent-field digests,
  49 focused adversarial tests, schema/CLI fail-closed checks, data-governance
  pass, and three independent reviewer passes; protocol remains non-executable;
  Git SHA `d941261`.
