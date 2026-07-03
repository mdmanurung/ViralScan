# ViralScan 99 Percent Publication Readiness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring ViralScan to a release-clean, install-tested, reproducible, scientifically honest state suitable for a software release and near-ready methods-manuscript submission.

**Architecture:** Split readiness into five independent tracks: package/release, CI/installability, documentation/API consistency, benchmark provenance, and manuscript/scientific evidence. Each track has its own tests and commit point so software-release work can land before longer manuscript validation. The manuscript track must not publish claims unless the supporting artifacts are tracked, reproducible, and internally consistent.

**Tech Stack:** Python packaging with setuptools/pyproject, GitHub Actions, pytest, Sphinx/MyST docs, Snakemake/kb-python/STARsolo workflows, SLURM benchmark runs, Markdown manuscript artifacts.

---

## Scope And Readiness Definition

This plan makes two readiness tracks explicit:

1. **Track A: software release readiness** means a clean tagged source tree, installed wheel smoke tests, packaged runtime resources, consistent docs, and release metadata prepared for PyPI/container/bioconda follow-up.
2. **Track B: methods-manuscript readiness** means all manuscript claims match tracked artifacts, the host-response depth confound is corrected in text and figures, incomplete benchmark evidence is either finished or removed, and a dedicated viral scRNA-seq tool comparison is added or the target venue is narrowed.

The repository is not 99 percent ready until all Track A tasks and either Task 8A or Task 8B are complete. For a strong methods paper, Task 8A is required. For a software-only release plus honest preprint draft, Task 8B is acceptable.

## File Map

- `pyproject.toml`: package metadata and package-data rules.
- `.github/workflows/ci.yml`: installed-package test gate and documented conda environment gate.
- `.github/workflows/release.yml`: wheel/sdist install smoke test before PyPI publish.
- `tests/test_cellcalling.py`: runtime resource regression for `emptydrops.R`.
- `tests/test_docs_consistency.py`: documentation/default/version consistency checks.
- `docs/installation.md`, `docs/quickstart.md`, `docs/cli_reference.md`, `docs/output_reference.md`, `docs/api.md`: user-facing documentation fixes.
- `src/viralscan/menu.py`: root help summary if evidence/check-whitelist are omitted.
- `CHANGELOG.md`, `CITATION.cff`, `conda-recipe/meta.yaml`: release metadata.
- `scripts/audit_reference_strategy.py`, `scripts/fetch_reference_strategy_fastqs.py`, `scripts/package_starsolo_viral_references.py`, `scripts/prepare_reference_strategy_benchmark.py`, `scripts/slurm_build_starsolo_reference_strategy_refs.sh`, `scripts/summarize_reference_strategy.py`: benchmark tooling that must be tracked if cited.
- `analysis/reference_strategy_benchmark/`: tracked lightweight benchmark provenance and final status.
- `docs/manuscript_draft.md`, `docs/figures/figure2_benchmark.png`, `docs/figures/figure2_benchmark.pdf`, `scripts/make_manuscript_figures.py`: manuscript claim and figure updates.
- `results/hostresponse_ebv_matched/*`: depth-confound and deconfounded host-response artifacts.
- `docs/PUBLICATION_READINESS.md`: final readiness ledger.

---

### Task 0: Freeze Baseline And Work Safely

**Files:**
- Read: whole repository
- Create: optional execution branch/worktree outside this plan
- Modify: none

- [ ] **Step 1: Record current branch and dirty state**

Run:

```bash
git branch --show-current
git status --short --untracked-files=all
git log --oneline -5
```

Expected: branch is `claude/multimap-memory-and-showcase`; dirty state includes manuscript and untracked readiness scripts.

- [ ] **Step 2: Decide execution isolation**

Run one of these before editing:

```bash
# Preferred for implementation work
git worktree add ../ViralScan-publication-ready HEAD
cd ../ViralScan-publication-ready
git switch -c chore/publication-readiness-99
```

or, if staying in the current checkout:

```bash
git switch -c chore/publication-readiness-99
```

Expected: a dedicated branch exists and user changes are not reverted.

- [ ] **Step 3: Create a readiness checklist issue file**

Create `docs/PUBLICATION_READINESS_EXECUTION.md` with:

```markdown
# Publication Readiness Execution Log

Started: 2026-07-03
Branch: chore/publication-readiness-99

## Gates

- [ ] Wheel includes every runtime resource.
- [ ] CI tests installed package, not only PYTHONPATH.
- [ ] Release workflow install-tests wheel before publish.
- [ ] Public docs agree with runtime defaults.
- [ ] Manuscript host-response AUROC is deconfounded or explicitly labelled confounded.
- [ ] Incomplete benchmark claims are removed or benchmark is complete.
- [ ] Clean clone can regenerate every cited lightweight result.
```

- [ ] **Step 4: Commit the log scaffold**

Run:

```bash
git add docs/PUBLICATION_READINESS_EXECUTION.md
git commit -m "docs: add publication readiness execution log"
```

Expected: one small docs commit.

---

### Task 1: Package The Missing `emptydrops.R` Runtime Resource

**Files:**
- Modify: `pyproject.toml`
- Modify: `tests/test_cellcalling.py`
- Verify: `src/viralscan/scripts/emptydrops.R`

- [ ] **Step 1: Write a failing resource test**

Append to `tests/test_cellcalling.py`:

```python
def test_emptydrops_r_script_is_next_to_cellcalling_module():
    from pathlib import Path

    import viralscan.scripts.cellcalling as cellcalling

    script = Path(cellcalling.__file__).with_name("emptydrops.R")
    assert script.exists()
    assert "emptyDrops" in script.read_text()
```

- [ ] **Step 2: Run the focused test**

Run:

```bash
python -m pytest tests/test_cellcalling.py::test_emptydrops_r_script_is_next_to_cellcalling_module -q
```

Expected before packaging fix in source checkout: PASS if the source file exists. This source-level test protects accidental deletion, not wheel packaging.

- [ ] **Step 3: Add package-data rule**

Edit `[tool.setuptools.package-data]` in `pyproject.toml` to include:

```toml
"viralscan.scripts" = ["emptydrops.R"]
```

Keep the existing lines:

```toml
"viralscan" = ["Snakefile", "templates/*.j2"]
"viralscan.data" = ["anellovirus_accessions.tsv"]
```

- [ ] **Step 4: Rebuild distribution artifacts**

Run in a Python 3.11 or 3.12 environment with build tooling:

```bash
rm -rf build dist src/ViralScan.egg-info
python -m pip install build twine
python -m build
python -m twine check dist/*
```

Expected: `twine check` reports both artifacts `PASSED`.

- [ ] **Step 5: Verify wheel and sdist contain the R script**

Run:

```bash
python -m zipfile -l dist/viralscan-2.5.0-py3-none-any.whl | grep 'viralscan/scripts/emptydrops.R'
tar -tzf dist/viralscan-2.5.0.tar.gz | grep 'src/viralscan/scripts/emptydrops.R'
```

Expected: both commands print exactly one matching path.

- [ ] **Step 6: Commit package resource fix**

Run:

```bash
git add pyproject.toml tests/test_cellcalling.py
git commit -m "fix(package): include emptyDrops runtime script"
```

Expected: commit contains only package-data and focused test changes.

---

### Task 2: Add Installed-Package CI And Release Gates

**Files:**
- Modify: `.github/workflows/ci.yml`
- Modify: `.github/workflows/release.yml`

- [ ] **Step 1: Update unit-test job to install the package**

In `.github/workflows/ci.yml`, replace the unit-test install/run pattern with:

```yaml
      - name: Install test dependencies
        run: |
          python -m pip install --upgrade pip
          python -m pip install pytest pytest-cov responses numpy pandas scipy scikit-learn anndata scanpy matplotlib seaborn plotly jinja2 pyyaml "requests>=2.28" pyfiglet
          python -m pip install --no-deps -e .

      - name: Run unit tests
        run: pytest tests/ -v --cov=viralscan --cov-report=term-missing --cov-report=xml --cov-fail-under=60

      - name: CLI smoke test (installed console script)
        run: |
          viralscan --version
          viralscan --help
          viralscan data fetch --help
          viralscan build-ref --help
          viralscan evidence --help
          viralscan check-whitelist --help
```

Remove `PYTHONPATH: src` from this job.

- [ ] **Step 2: Update integration job to install the package**

In `.github/workflows/ci.yml`, after micromamba dependency setup and before pytest, add:

```yaml
      - name: Install ViralScan package into integration environment
        run: python -m pip install --no-deps -e .
```

Then remove `PYTHONPATH: src` from the integration test step.

- [ ] **Step 3: Add documented environment validation job**

Add a CI job:

```yaml
  environment-file:
    name: Validate documented conda environment
    runs-on: ubuntu-latest
    defaults:
      run:
        shell: bash -el {0}
    steps:
      - uses: actions/checkout@v4
      - uses: mamba-org/setup-micromamba@v1
        with:
          environment-file: environment.yml
          environment-name: viralscan
          condarc: |
            channels:
              - conda-forge
              - bioconda
      - name: Install package and smoke-test documented commands
        run: |
          python -m pip install --no-deps .
          viralscan --version
          viralscan --help
          kb --version
          snakemake --version
```

Expected: CI validates the same `environment.yml` flow shown in docs.

- [ ] **Step 4: Add release artifact install-test before publish**

In `.github/workflows/release.yml`, after `python -m build`, add:

```yaml
      - name: Install-test built wheel
        run: |
          python -m venv /tmp/viralscan-wheel-test
          /tmp/viralscan-wheel-test/bin/python -m pip install --upgrade pip
          /tmp/viralscan-wheel-test/bin/python -m pip install pyyaml pyfiglet
          /tmp/viralscan-wheel-test/bin/python -m pip install --no-deps dist/*.whl
          /tmp/viralscan-wheel-test/bin/viralscan --version
          /tmp/viralscan-wheel-test/bin/viralscan --help
          /tmp/viralscan-wheel-test/bin/viralscan data fetch --help
          test -f /tmp/viralscan-wheel-test/lib/python*/site-packages/viralscan/scripts/emptydrops.R
```

Expected: publish cannot run if the wheel omits runtime resources or console scripts fail.

- [ ] **Step 5: Run local syntax validation for workflows**

Run:

```bash
python - <<'PY'
from pathlib import Path
import yaml
for path in [Path(".github/workflows/ci.yml"), Path(".github/workflows/release.yml")]:
    yaml.safe_load(path.read_text())
    print(f"valid yaml: {path}")
PY
```

Expected: both workflow files parse.

- [ ] **Step 6: Commit CI gates**

Run:

```bash
git add .github/workflows/ci.yml .github/workflows/release.yml
git commit -m "ci: test installed package and release artifacts"
```

Expected: CI/release changes are isolated.

---

### Task 3: Fix Documentation Defaults, API Drift, And Version Drift

**Files:**
- Create: `tests/test_docs_consistency.py`
- Modify: `docs/installation.md`
- Modify: `docs/quickstart.md`
- Modify: `docs/cli_reference.md`
- Modify: `docs/output_reference.md`
- Modify: `docs/api.md`
- Modify: `src/viralscan/menu.py`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Add documentation consistency tests**

Create `tests/test_docs_consistency.py`:

```python
from pathlib import Path

from viralscan import __version__
from viralscan.defaults import DEFAULT_MULTIMAP_METHOD


DOCS = [
    Path("README.md"),
    Path("docs/installation.md"),
    Path("docs/quickstart.md"),
    Path("docs/cli_reference.md"),
    Path("docs/output_reference.md"),
    Path("docs/api.md"),
]


def test_docs_do_not_reference_stale_container_version():
    for path in DOCS:
        text = path.read_text()
        assert "2.3.0" not in text, f"{path} still references 2.3.0"
        assert "viralscan_2.3.0" not in text, f"{path} still references old SIF name"


def test_runtime_multimap_default_is_documented_consistently():
    assert DEFAULT_MULTIMAP_METHOD == "equal"
    stale_phrases = [
        "`--multimap-method METHOD` | | `host-conservative`",
        "default multimapping method is `host-conservative`",
        "host-conservative is the default",
    ]
    for path in DOCS:
        text = path.read_text()
        for phrase in stale_phrases:
            assert phrase not in text, f"{path} contains stale default phrase: {phrase}"


def test_cli_reference_mentions_current_version():
    text = Path("docs/cli_reference.md").read_text()
    assert f"**{__version__}**" in text


def test_output_reference_documents_called_cell_column():
    text = Path("docs/output_reference.md").read_text()
    assert "is_called_cell" in text
```

- [ ] **Step 2: Run tests and confirm failures**

Run:

```bash
python -m pytest tests/test_docs_consistency.py -q
```

Expected: FAIL before docs are fixed, at least on stale defaults/version/schema.

- [ ] **Step 3: Fix container version examples**

In `docs/installation.md`, replace:

```bash
docker build -t viralscan:2.3.0 .
docker run --rm -it -v "$PWD:/data" viralscan:2.3.0 --help
singularity build viralscan_2.3.0.sif Singularity.def
singularity exec viralscan_2.3.0.sif viralscan --help
```

with:

```bash
docker build -t viralscan:2.5.0 .
docker run --rm -it -v "$PWD:/data" viralscan:2.5.0 --help
singularity build viralscan_2.5.0.sif Singularity.def
singularity exec viralscan_2.5.0.sif viralscan --help
```

- [ ] **Step 4: Fix multimap default docs**

In `docs/quickstart.md`, `docs/cli_reference.md`, and `docs/output_reference.md`, make all default descriptions match:

```markdown
The default `--multimap-method` is `equal` for fast first-pass reporting. For combined host+virus references where host-virus cross-homology matters, rerun or start with `--multimap-method host-conservative`.
```

In `docs/cli_reference.md`, set the table row to:

```markdown
| `--multimap-method METHOD` | | `equal` | Multimapper allocation: `equal`, `host-conservative`, `unique-weighted`, or `em` |
```

- [ ] **Step 5: Document live `build-ref` flags**

In `docs/cli_reference.md` build-ref table, add:

```markdown
| `--anellovirus` / `--no-anellovirus` | | on | Include the packaged Anelloviridae accession table in host+virus references |
| `--reference-panel anellovirus` | | *(none)* | Build a predefined Anelloviridae panel |
| `--no-mask` | | off | Skip low-complexity masking for anellovirus references |
| `--cluster` | | off | Cluster anellovirus sequences at 95% identity after masking |
```

Add this warning below the table:

```markdown
By default, `build-ref` includes the packaged Anelloviridae accession table. Pass `--no-anellovirus` when you want a reference containing only the host and the explicit `--virus-accessions`.
```

- [ ] **Step 6: Fix API examples**

In `docs/api.md`, replace dict-based `cell_type_enrichment()` examples with:

```python
from types import SimpleNamespace
from viralscan.enrichment import cell_type_enrichment

cfg = SimpleNamespace(cell_types="cell_types.csv")
cell_type_enrichment(
    config=cfg,
    stats=stats,
    per_cell=per_cell,
    output_dir="output/sample/results",
)
```

In the `build_combined_reference()` API section, document:

```markdown
`include_anellovirus=True` is the default. Set `include_anellovirus=False` to build only the host plus explicitly requested viral accessions.
```

- [ ] **Step 7: Fix output schema docs**

In `docs/output_reference.md`, add `is_called_cell` to the `per_cell_viral.tsv` schema:

```markdown
| `is_called_cell` | Boolean flag indicating whether the barcode is in the primary called-cell denominator |
```

Update AnnData layer docs to state:

```markdown
`adata.X` contains the active corrected count matrix. The original unique counts and additional multimapping layers are stored in `.layers`, including `counts_original`, `counts_corrected`, `counts_combined`, and method-specific layers such as `counts_multimap_equal`, `counts_multimap_host_conservative`, `counts_multimap_unique_weighted`, and `counts_viral_ambiguous_upper` when available.
```

- [ ] **Step 8: Fix root help subcommand summary**

In `src/viralscan/menu.py`, update the subcommand summary block to include:

```text
  evidence        Extract/read-trace evidence behind viral calls.
  check-whitelist Check barcode whitelist/chemistry compatibility.
```

- [ ] **Step 9: Fix changelog compare links**

At the bottom of `CHANGELOG.md`, replace stale refs with:

```markdown
[Unreleased]: https://github.com/mdmanurung/ViralScan/compare/v2.5.0...HEAD
[2.5.0]: https://github.com/mdmanurung/ViralScan/compare/v2.4.0...v2.5.0
[2.4.0]: https://github.com/mdmanurung/ViralScan/compare/v2.3.0...v2.4.0
[2.3.0]: https://github.com/mdmanurung/ViralScan/compare/v2.2.0...v2.3.0
[2.2.0]: https://github.com/mdmanurung/ViralScan/compare/v2.1.0...v2.2.0
[2.1.0]: https://github.com/mdmanurung/ViralScan/releases/tag/v2.1.0
```

- [ ] **Step 10: Run docs consistency tests**

Run:

```bash
python -m pytest tests/test_docs_consistency.py tests/test_cli.py::TestDefaults::test_multimap_method_default -q
```

Expected: PASS.

- [ ] **Step 11: Commit docs consistency**

Run:

```bash
git add tests/test_docs_consistency.py docs/installation.md docs/quickstart.md docs/cli_reference.md docs/output_reference.md docs/api.md src/viralscan/menu.py CHANGELOG.md
git commit -m "docs: align public docs with runtime defaults"
```

Expected: one documentation/API consistency commit.

---

### Task 4: Make Release Metadata And Distribution Source Complete

**Files:**
- Modify: `pyproject.toml`
- Modify: `CITATION.cff`
- Modify: `conda-recipe/meta.yaml`
- Modify: `MANIFEST.in` if using a manifest file
- Verify: `dist/*`

- [ ] **Step 1: Add source distribution manifest**

Create `MANIFEST.in`:

```text
include LICENSE
include README.md
include CHANGELOG.md
include CITATION.cff
include environment.yml
include Dockerfile
include Singularity.def
recursive-include conda-recipe *.yaml *.md
recursive-include docs *.md *.py *.txt
recursive-include src/viralscan/templates *.j2
recursive-include src/viralscan/data *.tsv
recursive-include src/viralscan/scripts emptydrops.R
```

- [ ] **Step 2: Rebuild and inspect sdist contents**

Run:

```bash
rm -rf build dist src/ViralScan.egg-info
python -m build
tar -tzf dist/viralscan-2.5.0.tar.gz | grep -E 'CITATION.cff|CHANGELOG.md|environment.yml|Dockerfile|Singularity.def|conda-recipe/meta.yaml|emptydrops.R'
```

Expected: all listed files appear in the source distribution.

- [ ] **Step 3: Keep bioconda sha256 intentionally pending until PyPI publish**

Leave `conda-recipe/meta.yaml` sha256 as zeros before PyPI publish, but add this comment above it:

```yaml
  # Release procedure: replace this after PyPI upload with the sha256 of
  # https://pypi.org/project/ViralScan/2.5.0/#files sdist.
```

Expected: recipe remains clearly pre-publication, not accidentally submit-ready.

- [ ] **Step 4: Add software DOI placeholder language without fake DOI**

In `CITATION.cff`, keep GitHub URL identifier and add no DOI until Zenodo archives the release. After Zenodo archive, add:

```yaml
  - type: doi
    value: 10.5281/zenodo.REPLACE_WITH_RELEASE_DOI
    description: Archived software release
```

Do not commit a fake DOI.

- [ ] **Step 5: Commit release metadata**

Run:

```bash
git add MANIFEST.in pyproject.toml CITATION.cff conda-recipe/meta.yaml
git commit -m "chore(release): complete source distribution metadata"
```

Expected: source distribution now carries citation, changelog, environment, and container recipes.

---

### Task 5: Promote Reproducibility Provenance Out Of Ignored Directories

**Files:**
- Create: `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/README.md`
- Create: `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/commands.jsonl`
- Create: `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/fastq_audit.tsv`
- Create: `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/reference_audit.tsv`
- Create: `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/failure_summary.tsv`
- Modify: `data/DATA_MANIFEST.md`
- Modify: `data/raw/reference_strategy_fastqs/REFERENCE_STRATEGY_FASTQS.md`
- Modify: `data/metadata/reference_strategy_fastqs/provenance.md`
- Modify: `data/raw/reference_strategy_refs/REFERENCE_STRATEGY_REFS.md`
- Modify: `reference_manifest.json`

- [ ] **Step 1: Track benchmark helper scripts**

Run:

```bash
git add scripts/audit_reference_strategy.py scripts/fetch_reference_strategy_fastqs.py scripts/package_starsolo_viral_references.py scripts/prepare_reference_strategy_benchmark.py scripts/slurm_build_starsolo_reference_strategy_refs.sh scripts/summarize_reference_strategy.py
```

Expected: helper scripts are staged for a future commit.

- [ ] **Step 2: Create tracked run packet directory**

Run:

```bash
mkdir -p analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12
cp benchmark_runs/reference_strategy_2026-06-28_fresh12/commands.jsonl analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/commands.jsonl
cp benchmark_runs/reference_strategy_2026-06-28_fresh12/fastq_audit.tsv analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/fastq_audit.tsv
cp benchmark_runs/reference_strategy_2026-06-28_fresh12/reference_audit.tsv analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/reference_audit.tsv
```

Expected: lightweight provenance is tracked; bulky logs/results remain ignored.

- [ ] **Step 3: Create failure summary**

Create `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/failure_summary.tsv`:

```text
row_id	status	evidence_path	reason
hhv6b__viralscan__two_step	incomplete	benchmark_runs/reference_strategy_2026-06-28_fresh12/logs/vs_ref_strategy_25102849_3.err	rerun began host-filtering but no downstream viral_summary.tsv or per_cell_viral.tsv was produced
ebv__viralscan__two_step	incomplete	benchmark_runs/reference_strategy_2026-06-28_fresh12/logs/vs_ref_strategy_25102849_7.err	rerun began host-filtering but no downstream viral_summary.tsv or per_cell_viral.tsv was produced
hsv1__viralscan__two_step	incomplete	benchmark_runs/reference_strategy_2026-06-28_fresh12/logs/vs_ref_strategy_25102849_11.err	rerun began host-filtering but no downstream viral_summary.tsv or per_cell_viral.tsv was produced
```

- [ ] **Step 4: Create run packet README**

Create `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/README.md`:

```markdown
# Reference Strategy Run Packet: 2026-06-28 fresh12

This packet preserves lightweight provenance for the 12-row reference-strategy benchmark attempt.

Tracked files:

- `commands.jsonl`: command manifest used by the SLURM array.
- `fastq_audit.tsv`: raw-byte FASTQ audit with sizes, hashes, record counts, and source URLs.
- `reference_audit.tsv`: reference artifact checksums and feature counts.
- `failure_summary.tsv`: concise status of rows that did not produce final ViralScan outputs.

Bulky files intentionally remain ignored under `benchmark_runs/`: raw logs, intermediate matrices, STAR indices, kallisto outputs, and FASTQs.

Publication status: this run packet is provenance only. It is not a complete benchmark result and must not be cited as final performance evidence.
```

- [ ] **Step 5: Correct data manifest statements**

In `data/DATA_MANIFEST.md`, replace any claim that every sample has a tracked `source_urls.tsv` with:

```markdown
Only `SRR12682296` currently has a tracked `benchmark_inputs/.../source_urls.tsv`. The full 2026-06-28 run packet records source URLs and raw-byte hashes for SRR20710641, SRR12682296, and SRR8315713 under `analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12/fastq_audit.tsv`.
```

- [ ] **Step 6: Correct reference provenance language**

In `data/raw/reference_strategy_refs/REFERENCE_STRATEGY_REFS.md`, add:

```markdown
The current benchmark references were built on institutional storage paths and are documented for provenance, not redistributed here. A clean public regeneration path must use `viralscan build-ref` or publish the exact reference artifacts separately before claims based on these references are considered reproducible.
```

- [ ] **Step 7: Update root reference manifest**

Replace stale “preflight” language in `reference_manifest.json` with a `"status"` field:

```json
"status": "audited_internal_paths_not_publicly_redistributed"
```

Keep checksums and paths factual; do not remove private paths if they are needed for provenance.

- [ ] **Step 8: Commit provenance promotion**

Run:

```bash
git add scripts/audit_reference_strategy.py scripts/fetch_reference_strategy_fastqs.py scripts/package_starsolo_viral_references.py scripts/prepare_reference_strategy_benchmark.py scripts/slurm_build_starsolo_reference_strategy_refs.sh scripts/summarize_reference_strategy.py analysis/reference_strategy_benchmark/run_packets/2026-06-28_fresh12 data/DATA_MANIFEST.md data/raw/reference_strategy_fastqs/REFERENCE_STRATEGY_FASTQS.md data/metadata/reference_strategy_fastqs/provenance.md data/raw/reference_strategy_refs/REFERENCE_STRATEGY_REFS.md reference_manifest.json
git commit -m "docs: preserve benchmark provenance for clean clones"
```

Expected: clean-clone users can inspect commands, audits, and failure status without ignored run directories.

---

### Task 6: Remove Incomplete Reference-Strategy Benchmark From Manuscript Evidence

**Files:**
- Modify: `docs/manuscript_draft.md`
- Modify: `analysis/reference_strategy_benchmark/REFERENCE_STRATEGY_BENCHMARK.md`
- Modify: `docs/PUBLICATION_READINESS.md`

- [ ] **Step 1: Mark reference-strategy benchmark as excluded from submission claims**

In `analysis/reference_strategy_benchmark/REFERENCE_STRATEGY_BENCHMARK.md`, add immediately after the status table:

```markdown
## Publication Use

This benchmark is excluded from current manuscript claims because the 12-row design did not complete. It can be cited only as an incomplete provenance/audit record until all planned rows produce final outputs and `scripts/summarize_reference_strategy.py` exits successfully.
```

- [ ] **Step 2: Ensure manuscript does not cite incomplete benchmark**

Run:

```bash
grep -nE "reference-strategy|Selectivity Index|fresh12|reference_strategy_benchmark" docs/manuscript_draft.md
```

Expected: no matches. If matches exist, replace those sentences with:

```markdown
A broader reference-strategy benchmark is in progress and is not used for the present claims.
```

- [ ] **Step 3: Update readiness doc**

In `docs/PUBLICATION_READINESS.md`, move reference-strategy benchmark from “paper evidence” to “future work unless completed”:

```markdown
The incomplete reference-strategy benchmark is excluded from the current manuscript evidence package. Reintroduce it only after all 12 rows complete and the final tracked summary contains no failed, blocked, incomplete, or running rows.
```

- [ ] **Step 4: Commit benchmark exclusion**

Run:

```bash
git add docs/manuscript_draft.md analysis/reference_strategy_benchmark/REFERENCE_STRATEGY_BENCHMARK.md docs/PUBLICATION_READINESS.md
git commit -m "docs: exclude incomplete reference-strategy benchmark from claims"
```

Expected: manuscript no longer relies on incomplete benchmark rows.

---

### Task 7: Correct Host-Response Manuscript Text And Regenerate Figure 2

**Files:**
- Modify: `docs/manuscript_draft.md`
- Modify: `scripts/make_manuscript_figures.py`
- Modify: `docs/figures/figure2_benchmark.png`
- Modify: `docs/figures/figure2_benchmark.pdf`
- Verify: `results/hostresponse_ebv_matched/depth_confounder.txt`
- Verify: `results/hostresponse_ebv_matched/cpm_label_crosscheck.txt`
- Verify: `results/hostresponse_ebv_matched/depth_matched_reanalysis.txt`

- [ ] **Step 1: Replace the host-response subsection heading**

In `docs/manuscript_draft.md`, replace:

```markdown
### Viral burden supports leakage-controlled host-response modelling
```

with:

```markdown
### Host-response modelling requires explicit depth controls
```

- [ ] **Step 2: Replace the raw AUROC claim**

Replace the paragraph reporting `AUROC 0.845 +/- 0.032` with:

```markdown
Using the raw `>=10 corrected UMI` label, host-expression models predicted EBV-high status with AUROC 0.845 +/- 0.032. However, this raw label was strongly depth-confounded: sequencing depth alone reached AUROC 0.967 in the same matched evaluation design. We therefore treat the raw AUROC as a diagnostic rather than biological evidence. Depth-controlled analyses gave more conservative estimates: a CPM-normalized prevalence-matched label yielded AUROC 0.636, and a depth-matched cohort yielded AUROC 0.718. These depth-controlled values define the interpretable host-response signal reported here.
```

- [ ] **Step 3: Replace the Ensembl-only stable-feature paragraph**

Replace the paragraph beginning `Randomised Lasso stability selection` with:

```markdown
Randomised Lasso stability selection identified a small set of stable host features, but biological interpretation remains limited until gene-symbol annotation and pathway enrichment are finalized for the depth-controlled labels. The current result supports ViralScan's ability to expose depth confounding and run matched host-response analyses; it does not support interpreting the raw AUROC as a depth-independent host-response signature.
```

- [ ] **Step 4: Add discussion limitation**

In the limitations paragraph, add:

```markdown
Host-response modelling is sensitive to label definitions because raw viral UMI thresholds can track library size. ViralScan now reports a depth-alone baseline and supports CPM labels and depth-matched designs; manuscripts should report these baselines alongside model AUROC.
```

- [ ] **Step 5: Update Figure 2 footer generation**

In `scripts/make_manuscript_figures.py`, change the Figure 2 host-response footer text to:

```python
host_response_footer = (
    "EBV host-response: raw AUROC 0.845 is depth-confounded "
    "(depth-alone AUROC 0.967); depth-controlled AUROC 0.636 CPM-label "
    "and 0.718 depth-matched."
)
```

Use the script’s existing variable/style names if they already exist.

- [ ] **Step 6: Regenerate figures**

Run:

```bash
python scripts/make_manuscript_figures.py
```

Expected: `docs/figures/figure2_benchmark.png` and `.pdf` are updated.

- [ ] **Step 7: Verify stale AUROC framing is gone**

Run:

```bash
grep -n "AUROC 0.845" docs/manuscript_draft.md
grep -n "depth-alone" docs/manuscript_draft.md
```

Expected: first command appears only in a sentence labelling the raw value as depth-confounded; second command finds the new depth-baseline language.

- [ ] **Step 8: Commit host-response correction**

Run:

```bash
git add docs/manuscript_draft.md scripts/make_manuscript_figures.py docs/figures/figure2_benchmark.png docs/figures/figure2_benchmark.pdf
git commit -m "docs: correct depth-confounded host-response claim"
```

Expected: manuscript and figure now match host-response artifacts.

---

### Task 8A: Add Dedicated Viral scRNA-seq Tool Comparison For Methods-Manuscript Readiness

**Files:**
- Create: `analysis/dedicated_tool_comparison/README.md`
- Create: `analysis/dedicated_tool_comparison/run_venus_or_viraltrack.sh`
- Create: `analysis/dedicated_tool_comparison/summarize_dedicated_tool_comparison.py`
- Create: `results/dedicated_tool_comparison.tsv`
- Modify: `docs/manuscript_draft.md`
- Modify: `docs/PUBLICATION_READINESS.md`

- [ ] **Step 1: Define comparison scope**

Create `analysis/dedicated_tool_comparison/README.md`:

```markdown
# Dedicated Viral scRNA-seq Tool Comparison

Scope: compare ViralScan against one dedicated viral single-cell RNA-seq detector on SRR12682296 EBV LCL matched barcodes.

Primary matched cell set: 1,906 GSM4796271 LCL_777_B958 barcodes used in `results/matched_barcode_comparison.tsv`.

Primary metrics:

- EBV-positive cells at >=1 viral UMI or tool-equivalent positive call.
- EBV-high cells at >=10 viral UMI or tool-equivalent high-confidence call.
- Lytic-marker-positive cells using BZLF1, BRLF1, or BHRF1 where the tool exposes gene-level calls.
- Per-cell viral-burden rank concordance where both tools output continuous burden.
- Runtime and peak memory from scheduler accounting.

Acceptance rule: the manuscript may claim dedicated-tool comparison only after `results/dedicated_tool_comparison.tsv` has non-empty rows for both ViralScan and the comparator on the same cell set.
```

- [ ] **Step 2: Add runner wrapper**

Create `analysis/dedicated_tool_comparison/run_venus_or_viraltrack.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

echo "Comparator execution is intentionally explicit."
echo "Install the selected comparator in a separate conda environment."
echo "Write raw comparator outputs under analysis/dedicated_tool_comparison/raw/"
echo "Do not overwrite ViralScan benchmark artifacts."
exit 2
```

This wrapper intentionally fails until a specific comparator command is inserted after tool installation is confirmed.

- [ ] **Step 3: Add summarizer contract**

Create `analysis/dedicated_tool_comparison/summarize_dedicated_tool_comparison.py`:

```python
from __future__ import annotations

from pathlib import Path

import pandas as pd


def main() -> None:
    raw_dir = Path("analysis/dedicated_tool_comparison/raw")
    matched = Path("results/matched_barcode_comparison.tsv")
    if not raw_dir.exists():
        raise SystemExit("Missing comparator raw output directory")
    if not matched.exists():
        raise SystemExit("Missing ViralScan/STARsolo matched comparison table")

    # This file is deliberately small until the comparator format is pinned.
    out = Path("results/dedicated_tool_comparison.tsv")
    rows = [
        {
            "tool": "ViralScan",
            "dataset": "SRR12682296",
            "cell_set": "GSM4796271 matched 1906 cells",
            "status": "available_in_results_matched_barcode_comparison",
        },
        {
            "tool": "dedicated_comparator",
            "dataset": "SRR12682296",
            "cell_set": "GSM4796271 matched 1906 cells",
            "status": "raw_output_required",
        },
    ]
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    raise SystemExit("Comparator raw output required before manuscript use")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run and confirm it blocks manuscript use**

Run:

```bash
python analysis/dedicated_tool_comparison/summarize_dedicated_tool_comparison.py
```

Expected: nonzero exit with `Comparator raw output required before manuscript use`.

- [ ] **Step 5: Execute comparator on SLURM**

After selecting and installing Venus or ViralTrack, replace the intentional `exit 2` in `run_venus_or_viraltrack.sh` with the exact comparator command and submit:

```bash
sbatch --parsable analysis/dedicated_tool_comparison/run_venus_or_viraltrack.sh
```

Expected: one SLURM job ID and raw comparator outputs under `analysis/dedicated_tool_comparison/raw/`.

- [ ] **Step 6: Update summarizer for real comparator output**

Modify `summarize_dedicated_tool_comparison.py` so it reads the comparator’s actual output columns and writes:

```text
tool	dataset	cell_set	positive_cells	high_confidence_cells	lytic_marker_positive_cells	spearman_vs_viralscan	runtime_minutes	max_rss_gb	status
```

Expected: no `raw_output_required` status remains.

- [ ] **Step 7: Add manuscript paragraph**

Add to `docs/manuscript_draft.md` after the STARsolo comparison:

```markdown
We additionally compared ViralScan with a dedicated viral single-cell RNA-seq detector on the same 1,906 matched EBV LCL cells. This comparison uses identical barcode anchors and reports positive-cell count, high-confidence viral burden, lytic-marker calls, rank concordance, runtime, and peak memory. The dedicated-tool comparison is summarized in `results/dedicated_tool_comparison.tsv`.
```

Only add numerical claims after the TSV contains final values.

- [ ] **Step 8: Commit dedicated comparison scaffold or final comparison**

Run:

```bash
git add analysis/dedicated_tool_comparison results/dedicated_tool_comparison.tsv docs/manuscript_draft.md docs/PUBLICATION_READINESS.md
git commit -m "analysis: add dedicated viral tool comparison"
```

Expected: if final values are absent, docs must state scaffold only and manuscript must not claim results.

---

### Task 8B: Narrow Manuscript Scope If Dedicated Comparison Is Deferred

**Files:**
- Modify: `docs/manuscript_draft.md`
- Modify: `docs/PUBLICATION_READINESS.md`

- [ ] **Step 1: Narrow target wording**

In `docs/manuscript_draft.md`, replace any methods-paper competitive framing with:

```markdown
This manuscript presents ViralScan as an open-source workflow and validation case study rather than a comprehensive benchmark against all dedicated viral single-cell detectors.
```

- [ ] **Step 2: Add explicit comparison limitation**

In the Discussion, add:

```markdown
We did not complete a head-to-head comparison against dedicated viral single-cell tools such as Venus or ViralTrack in the present version. The STARsolo comparison tests a splice-aware general aligner and should not be interpreted as a comprehensive dedicated-tool benchmark.
```

- [ ] **Step 3: Update readiness doc**

In `docs/PUBLICATION_READINESS.md`, state:

```markdown
Dedicated-tool comparison deferred. With this deferral, the current manuscript should target a software/resource or preprint venue, not a methods venue that requires comprehensive head-to-head benchmarking.
```

- [ ] **Step 4: Commit scope narrowing**

Run:

```bash
git add docs/manuscript_draft.md docs/PUBLICATION_READINESS.md
git commit -m "docs: narrow manuscript claims pending dedicated-tool benchmark"
```

Expected: manuscript no longer overclaims methods-paper competitiveness.

---

### Task 9: Add Minimal Truth-Panel Specificity/Sensitivity Framing

**Files:**
- Create: `analysis/truth_panel/README.md`
- Modify: `README.md`
- Modify: `docs/manuscript_draft.md`

- [ ] **Step 1: Add truth-panel design note**

Create `analysis/truth_panel/README.md`:

```markdown
# Truth Panel Plan For ViralScan Specificity/Sensitivity

The current release validates against published public scRNA-seq datasets but does not estimate formal false-positive or false-negative rates against known single-cell infection truth.

Minimum acceptable truth panel:

1. A negative-control host-only scRNA-seq dataset processed with the same host+virus reference.
2. A synthetic spike-in or read-mixing dataset with known viral read/barcode assignments.
3. A positive-control infected dataset with an orthogonal published infection label or accepted threshold.

Until these are complete, README and manuscript claims must say FPR/FNR are not formally characterized.
```

- [ ] **Step 2: Keep README limitation explicit**

Ensure `README.md` Limitations retains:

```markdown
False-positive / false-negative rates are not yet characterized on an independent gold-standard single-cell truth panel.
```

- [ ] **Step 3: Add manuscript limitation**

In `docs/manuscript_draft.md` Discussion, add:

```markdown
The present validation uses published infection-rate ranges and matched tool comparisons rather than a gold-standard single-cell truth panel. Formal false-positive and false-negative rates require negative controls and planted viral-read simulations and are left as a release-gated validation extension.
```

- [ ] **Step 4: Commit truth-panel framing**

Run:

```bash
git add analysis/truth_panel/README.md README.md docs/manuscript_draft.md
git commit -m "docs: define viral detection truth-panel limits"
```

Expected: limitation is explicit and reviewers are not misled.

---

### Task 10: Complete Submission Metadata And Citation Hygiene

**Files:**
- Modify: `docs/manuscript_draft.md`
- Modify: `CITATION.cff`
- Modify: `docs/PUBLICATION_READINESS.md`

- [ ] **Step 1: Replace manuscript placeholders**

Replace:

```markdown
**Authors:** [Author list to be supplied]
**Affiliations:** [Affiliations to be supplied]
**Lead contact:** [Lead contact to be supplied]
```

with author-approved values. This must be done by the project owner or corresponding author.

- [ ] **Step 2: Replace acknowledgments and declarations**

Replace:

```markdown
[Acknowledgments to be supplied by authors.]
[Author contributions to be supplied by authors.]
[Declaration of interests to be supplied by authors.]
```

with author-approved text.

- [ ] **Step 3: Resolve VIRTUS2 citation**

Run:

```bash
grep -n "VIRTUS2" docs/manuscript_draft.md
```

Then choose one of these concrete outcomes:

```markdown
8. [Remove this reference and refer only to VIRTUS if no citable VIRTUS2 record is available.]
```

or replace with a verified DOI-backed citation. Do not leave bracketed citation placeholders.

- [ ] **Step 4: Commit submission metadata**

Run:

```bash
git add docs/manuscript_draft.md CITATION.cff docs/PUBLICATION_READINESS.md
git commit -m "docs: complete manuscript submission metadata"
```

Expected: no square-bracket placeholders remain except explanatory editorial comments intentionally excluded from public docs.

---

### Task 11: Final Local Verification Matrix

**Files:**
- Read: whole repository
- Modify: `docs/PUBLICATION_READINESS_EXECUTION.md`

- [ ] **Step 1: Run source syntax check**

Run:

```bash
python -m compileall -q src tests scripts
```

Expected: exit code 0.

- [ ] **Step 2: Run unit tests in the correct environment**

Run in Python 3.11 or 3.12 environment with dependencies:

```bash
python -m pytest -q
```

Expected: default test suite passes with `not network and not integration and not research` marker behavior.

- [ ] **Step 3: Run integration tests**

Run in an environment with `kb`, `snakemake`, `kallisto`, `bustools`, and required scientific Python deps:

```bash
python -m pytest -m "integration and not network" -q
```

Expected: integration subset passes.

- [ ] **Step 4: Build and install-test artifacts**

Run:

```bash
rm -rf build dist src/ViralScan.egg-info
python -m build
python -m twine check dist/*
python -m venv /tmp/viralscan-wheel-test
/tmp/viralscan-wheel-test/bin/python -m pip install --upgrade pip
/tmp/viralscan-wheel-test/bin/python -m pip install pyyaml pyfiglet
/tmp/viralscan-wheel-test/bin/python -m pip install --no-deps dist/*.whl
/tmp/viralscan-wheel-test/bin/viralscan --version
/tmp/viralscan-wheel-test/bin/viralscan --help
test -f /tmp/viralscan-wheel-test/lib/python*/site-packages/viralscan/scripts/emptydrops.R
```

Expected: wheel install smoke passes and `emptydrops.R` exists.

- [ ] **Step 5: Build docs**

Run:

```bash
python -m pip install -r docs/requirements.txt
sphinx-build -W -b html docs docs/_build/html
```

Expected: docs build with warnings treated as errors.

- [ ] **Step 6: Verify no stale or unsupported claims remain**

Run:

```bash
grep -RIn "2.3.0\\|host-conservative.*default\\|AUROC 0.845 +/- 0.032\\|to be supplied\\|TBD\\|citation to be verified" README.md docs analysis results scripts src tests CHANGELOG.md CITATION.cff
```

Expected: no matches that represent public-facing stale claims. Matches inside this plan file are acceptable.

- [ ] **Step 7: Update execution log**

Append to `docs/PUBLICATION_READINESS_EXECUTION.md`:

```markdown
## Final Verification

- `python -m compileall -q src tests scripts`: PASS
- `python -m pytest -q`: PASS
- `python -m pytest -m "integration and not network" -q`: PASS
- `python -m build`: PASS
- `python -m twine check dist/*`: PASS
- wheel install smoke: PASS
- `sphinx-build -W -b html docs docs/_build/html`: PASS
- stale-claim grep: PASS
```

- [ ] **Step 8: Commit verification log**

Run:

```bash
git add docs/PUBLICATION_READINESS_EXECUTION.md
git commit -m "docs: record publication readiness verification"
```

Expected: final verification evidence is tracked.

---

### Task 12: Release And Post-Release Steps

**Files:**
- Modify after PyPI publish: `conda-recipe/meta.yaml`
- Modify after Zenodo archive: `CITATION.cff`
- Read: GitHub Actions release logs

- [ ] **Step 1: Confirm clean tree**

Run:

```bash
git status --short --untracked-files=all
```

Expected: no uncommitted files except intentionally ignored local caches/build artifacts.

- [ ] **Step 2: Open PR and wait for CI**

Run:

```bash
git push -u origin chore/publication-readiness-99
```

Open a PR to `main`. Expected: CI includes installed-package, release-artifact smoke, environment-file validation, docs consistency, and integration jobs.

- [ ] **Step 3: Tag after merge**

After merge to `main`:

```bash
git checkout main
git pull --ff-only
git tag v2.5.0
git push origin v2.5.0
```

Expected: GitHub release workflow builds wheel/sdist, install-tests wheel, publishes to PyPI, and pushes GHCR image.

- [ ] **Step 4: Fill bioconda sha256 after PyPI publish**

Run:

```bash
python - <<'PY'
import hashlib
from pathlib import Path
path = Path("dist/viralscan-2.5.0.tar.gz")
print(hashlib.sha256(path.read_bytes()).hexdigest())
PY
```

Replace zeros in `conda-recipe/meta.yaml` with the PyPI sdist hash. Commit:

```bash
git add conda-recipe/meta.yaml
git commit -m "chore(bioconda): fill v2.5.0 sdist sha256"
```

- [ ] **Step 5: Archive GitHub release on Zenodo**

After Zenodo creates the software DOI, add the DOI identifier to `CITATION.cff` and commit:

```bash
git add CITATION.cff
git commit -m "docs: add Zenodo software DOI"
```

Expected: software DOI and data DOI are distinct and both citation paths are clear.

---

## Final Self-Review Checklist

- [ ] Every runtime file referenced by installed code is included in wheel and sdist.
- [ ] CI tests the installed package and console script.
- [ ] Release workflow install-tests wheel before publish.
- [ ] Docs agree with runtime defaults: `--multimap-method equal` is default; `host-conservative` is recommended for cross-homology.
- [ ] Docs do not show stale `2.3.0` examples.
- [ ] Output/API docs match live code.
- [ ] Manuscript does not report confounded AUROC as biological evidence.
- [ ] Figure 2 footer matches deconfounded host-response artifacts.
- [ ] Incomplete reference-strategy benchmark is either finished or explicitly excluded from manuscript claims.
- [ ] Dedicated viral-tool comparison is complete, or manuscript target/scope is narrowed.
- [ ] All tracked scripts needed for cited benchmarks are committed.
- [ ] Clean clone contains lightweight provenance for all cited benchmark artifacts.
- [ ] No submission placeholders remain in manuscript.
- [ ] `docs/PUBLICATION_READINESS.md` accurately states remaining post-tag actions only.
