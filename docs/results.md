# Results

After `viralscan validate-run --run-dir output/` reports `"ok": true`, open
`output/<sample>/report.html` for the overview and then read the tables below.
The first half of this page says what the main results mean; the second half
gives five checks that move a `detected` row from "the run is technically sound"
to "what this candidate supports". Exact column definitions are in the
[Output reference](output_reference.md). All `results/` paths are under
`output/<sample>/`.

## Which file answers which question

| Question | File | First fields to read |
|---|---|---|
| What did this run find? | `results/viral_summary.tsv` | `virus_name`, `viral_molecules_total_est`, `infected_called`, `n_called_cells`, `pct_infected_called`, `detected` |
| Which cells carry a candidate? | `results/per_cell_viral.tsv` | `barcode`, `virus_name`, `viral_molecules_total_est`, `viral_fraction`, `is_called_cell` |
| What is the per-cell denominator? | `results/called_cells.tsv` | `barcode` (one row per called cell) |
| Could a negative be real? | `results/sensitivity.tsv` | `virus_name`, `observed_molecules`, `detected`, `capture`, `lod95_molecules`, `informative_negative`, `negative_blockers` |
| What does the report show? | `report.html` | Run metadata, QC summary, per-virus detection table, multimapping evidence, histograms |
| What did the run actually do? | `run_manifest.json`, `run_complete.json` | Options, input hashes, software identity, completion marker and artifact hashes |
| Which reference was used? | `results/reference_provenance.json` | Index, `t2g`, GTF, technology, multimapping settings, viral accessions |

`report.html` is self-contained and works offline. Rows with `detected = false`
are excluded from its virus overlays. `summary.txt` gives the same per-virus
figures as plain text and marks a candidate that misses the called-cell gate as
`CANDIDATE ONLY`.

## Units: molecules, not reads

ViralScan counts corrected **(cell barcode, UMI) molecules**, not reads. PCR and
sequencing duplicate reads for the same molecule are collapsed before counting;
BUS `count` is read multiplicity and is audit-only.

`viral_molecules_total_est` has two components:

- **unique support**: the molecule's compatible gene set had exactly one gene;
- **allocated ambiguous support**: the molecule was compatible with several
  genes and the selected multimapping method allocated it. The allocation is a
  model estimate and moves if the method changes.

In the count matrix, `counts_unique` and `counts_ambiguous_allocated` are
disjoint and sum to `adata.X`. `viral_molecules_per_10k_est` is
`viral_molecules_total_est / molecules_total_est × 10 000`, both corrected
molecule estimates; in a two-step (`--host-filter starsolo`) run it is NA, with
the reason in `denominator_note`.

## `detected` is a configured threshold result

`detected = true` means the virus passed the gate configured for this run:

- at least `--detection-min-cells` called cells are infected with the virus
  (default 1), **and**
- the called cells hold at least `--detection-threshold` of its molecules in
  total (default 1).

A cell is infected for a virus when its summed allocated molecules reach
`--infected-min-molecules` (default 1.0). Weaker cells are counted in
`infected_subthreshold_cells`, never dropped.

`detected = false` means *candidate only*: the virus had signal in the
all-barcode screen but missed the called-cell gate. It is not proof of absence.
Passing the gate is candidate evidence, not a confirmed infection.

## Denominators: all barcodes versus called cells

| Field | Numerator / denominator | Use |
|---|---|---|
| `pct_infected` | infected cells over **all barcodes** | Screening view; diluted by empty droplets |
| `pct_infected_called` | infected called cells over `n_called_cells` | Within-run per-cell candidate-support rate |
| `pct_infected_comparable` | Infected cells within the comparable set over called barcodes clearing an absolute host-UMI floor (200 host molecules) | Comparing runs that used different host-filtering strategies |

`pct_infected_called` divides by the called cells of *this run*. In a two-step
host-subtraction run, cell calling happens after filtering, so the denominator
moves with the strategy: host filtering can change called cells and viral
molecules by different factors. Do not compare `pct_infected_called` across
host-filtering strategies.

`pct_infected_comparable` uses a denominator that host filtering cannot move. Its
denominator, numerator and percentage are **blank (unavailable), not zero**, when
no cell clears the host-depth floor, when host depth is unavailable, and for
virus-only matrices or external overrides without a host matrix. Cross-arm
comparisons should use the same frozen reference cell set for every arm.

None of these fields is proof of infection. The field names are retained for
schema compatibility.

## Reading one row

The values below are **illustrative**, not from a real sample.

| virus_name | viral_molecules_total_est | infected_called | n_called_cells | pct_infected_called | background_umis | detected |
|---|---:|---:|---:|---:|---:|---|
| Epstein-Barr virus | 412.0 | 37 | 9,850 | 0.38 | 19.5 | true |
| Human herpesvirus 6B | 6.0 | 0 | 9,850 | 0.00 | 6.0 | false |

Read the first row as: 412 estimated viral molecules in total, of which 37
called cells carry support, over 9,850 called cells (0.38 %), with 19.5 molecules
outside the called set. It passed the gate, so `detected` is true.

Read the second row as: 6 estimated molecules, all of them in barcodes outside
the called set (`background_umis = 6.0`), and no infected called cell. It is a
candidate that failed the called-cell gate. That pattern is consistent with
ambient RNA, barcode swapping or empty droplets; it does not establish that the
virus is absent.

For each row you intend to use:

1. Filter on `detected` to separate the detected set from candidates.
2. Read `n_called_cells` and `background_umis` before any percentage.
3. Check `role`, `reference_risk_flags`, `artifact_risk` and
   `sibling_crossmap_note`. These are labels, never filters.
4. Check `host_homology_status`; `not_measured` is not zero homology.
5. Before treating a call as biology, run the five checks below.

## The AnnData file

For downstream analysis, load the corrected matrix:

```python
import scanpy as sc

adata = sc.read_h5ad(
    "output/sample/kb-python/counts_unfiltered/adata_multimap.h5ad"
)
print(adata)
```

`adata.X` is the complete selected-method matrix. `obs` carries per-cell
metadata and `layers` keeps the components, including `counts_unique`,
`counts_ambiguous_allocated`, the alternative allocation layers, and the
host–virus ambiguity diagnostics. The [Output reference](output_reference.md)
lists every layer and column.

An empty field anywhere means the measurement was unavailable, not that the
value was zero or that the result was clean. Examples: `pct_infected_comparable`
when host depth cannot be established, `capture` when no control certifies the
row, and the anellovirus alignment columns when that branch did not run. Report
these as unavailable.

## Assess a candidate: five checks

A detected row is a candidate, not a conclusion. Run the checks in order; each
names the artifact, what it can and cannot show, and the next action. There is
no universal breadth, identity or molecule cutoff that separates a real
infection from an artefact. The evidence tiers are rule-based categories, not
calibrated probabilities.

### 1. Run integrity

**Artifacts.** `viralscan validate-run --run-dir output/` and its JSON report;
`output/run_complete.json`; `output/run_manifest.json`;
`results/reference_provenance.json`.

Validation confirms the manifest schema, the run fingerprint and input hashes,
the completion marker, the artifact hashes recorded in that marker, and the
count invariants in the H5AD. A missing `run_complete.json` means the run was
interrupted or is still in progress. Read these manifest fields before trusting
a number:

- `software_identity`: package content digest and, from a checkout, the Git
  commit and dirty-source status.
- `input_fingerprints` / `input_paths`: the exact FASTQ, index, `t2g`, GTF and
  external-list bytes.
- `options`: the run's configured thresholds and methods.
- `defaults_status`: `provisional` on every `3.0.0.devN` build, because the
  detection defaults are not yet validated against the truth panel. Numbers from
  such a run are provisional; the field becomes `final` from the release
  candidate.
- `chemistry_sanity`: post-count library checks, including failed or
  unavailable checks.
- `reference_catalogue` / `catalogue_drift`: whether the packaged virus
  catalogue still matches the index build. Drift is a warning, not a failure.

**Next action.** If validation reports an error, fix that before reading any
result: restore the input, rebuild the reference, or start a fresh output
directory. Validation cannot show that the reference covers your virus or that
a call is infection.

### 2. Cell support

**Artifacts.** `results/called_cells.tsv`, `results/viral_summary.tsv`
(`n_called_cells`, `infected_called`, `pct_infected_called`, `background_umis`),
`results/chemistry_sanity.json`, and `results/cell_calling_input.json` for
two-step runs.

The called-cell list is the denominator. Molecules in barcodes outside it are
`background_umis` and never count toward the detection gate. That signal is
consistent with ambient RNA, barcode swapping or empty droplets, and the run
does not say which. Check that the denominator matches the experiment:

- `--cell-calling auto` uses an external list when `--called-cells-file` is
  supplied and otherwise runs EmptyDrops on the host transcriptome only, never on
  viral reads. `external` requires the list.
- `knee` is sensitivity-only. On real libraries it lands at the 10-molecule floor
  and calls empty droplets as cells, so its numbers are not reportable.
- `none` reports over all barcodes and is for deliberate toy data only. It is
  not a workaround for failed cell calling.
- For a two-step run, EmptyDrops uses the complete STARsolo `GeneFull/raw` host
  matrix. A missing or incomplete matrix is an error; the virus-only matrix is
  not a substitute.
- An external list sets the full denominator, including cells absent from the
  viral matrix. `results/host_called_cells.tsv` records that universe;
  `results/called_cells.tsv` records its intersection with the viral matrix.

Check chemistry before trusting an all-empty or all-background pattern.
`viralscan check-whitelist` reports the barcode–whitelist match rate; the main
run also detects chemistry from the first 100k R1 reads and stops when the reads
fit no single chemistry, when samples in one run disagree, or when an explicit
`--technology` contradicts the reads. `--require-chemistry-sanity` turns an
error or an unavailable post-count check into a stop.

**Next action.** If the denominator is wrong, fix the external list or the
chemistry and re-run. If the signal is only in `background_umis`, treat the row
as a candidate and continue to checks 3 to 5. The
[cell-calling and denominators](vignettes/cell_calling_denominators.ipynb)
notebook shows an external-list example and a wrong-list failure on synthetic
data.

### 3. Mapping support

**Artifacts.** `results/viral_summary.tsv` (`host_viral_ambig_fraction`,
`accession_breadth`, `sibling_crossmap_note`, `host_homology_status`,
`host_homology_max_identity`, `artifact_risk`, `role`, `reference_risk_flags`),
`results/multimap_evidence.tsv`, and the H5AD count layers.

Separate the evidence by where the molecules came from:

- **Unique support** is unambiguous molecule evidence.
- **Allocated ambiguous support** is a model estimate. For a combined host +
  virus reference, the default `host-conservative` method sends mixed
  host–virus molecule mass only to compatible host genes, keeping it out of
  primary viral estimates. `host_viral_ambig_fraction` reports how much of the
  virus's signal was host–virus ambiguous.
- `accession_breadth` is the fraction of the virus's index genes with at least
  one molecule. An endogenous-element artefact concentrates on one or two loci;
  a genuine infection spreads across the genes. This is gene-level breadth, not
  genome coverage.
- `sibling_crossmap_note` warns when near-identical viruses can exchange
  molecules (for example HHV-6A/6B or EBV types). Allocation cannot prove which
  sibling is present.
- `host_homology_status` says whether the reference genomes were measured
  against the host genome. `not_measured` is **not** zero homology. The maximum
  identity and coverage fields are observations, with no threshold applied.
- `artifact_risk` labels families prone to low-complexity read artefacts, today
  Anelloviridae only. It is a label, never a filter: low-complexity evidence
  bounds a real infection, it does not exclude one. `role` and
  `reference_risk_flags` likewise qualify a call without suppressing it.

**Next action.** For a call that matters, re-run allocation with
`viralscan rerun-multimap --run-dir output/ -o output_equal/ --multimap-method equal`
as an explicit comparison, and inspect reads for the exact target with check 4.
If reference coverage is the problem, rebuild with the host genome or the
expanded anellovirus panel; see [Reference panel](reference_panel.md). The
[multimapping correction](vignettes/multimapping_correction.ipynb) notebook
works through an ambiguous host–virus case on synthetic data.

### 4. Read inspection

**Artifacts.** A `viralscan evidence` output directory. The command needs a
completed run, an exact target, and, for competitive alignment, both the viral
and host FASTA. Competitive alignment needs `minimap2` and `samtools` on `PATH`,
and `--blast` also needs BLAST+ (`blastn`, `makeblastdb`):

```bash
viralscan evidence \
  --run-dir output/sample/ \
  --viral-fasta EBV.fa \
  --host-fasta GRCh38.fa \
  --virus EBV \
  --blast --read-start-profile --cell-tags \
  -o output/sample/evidence/
```

`--viral-fasta` requires `--host-fasta`; alignment and BLAST are competitive, so
the target is judged against host sequence rather than in isolation. `--virus`
accepts an exact accession, a registered alias, or a canonical detected call;
substring matching is forbidden. Reusing a non-empty evidence directory requires
a matching `evidence_manifest.json`; use a fresh directory for older outputs.

The outputs trace the corrected molecules behind the call:

- `read_lineage.tsv.gz`: exact read ID, CB, UB, equivalence class, compatible
  genes, ambiguity class, assigned weight, method, tier and exclusion reason.
- `molecule_verdicts.tsv.gz` / `molecule_verdict_summary.tsv`: one row per
  corrected `(CB, UMI)` molecule, with the best viral and host alignment score
  and a verdict (`virus_best`, `host_best`, `tie`, `unaligned`). The default
  `--dedup umi` keeps one deterministic representative per corrected
  `(CB, UMI, reference)` molecule; `markdup` is positional deduplication, not
  molecule evidence.
- `alignment_qc.tsv`, `per_cell_alignment_qc.tsv`, `coverage.*`, and the
  competitive BAMs: coverage, identity, breadth, duplication, complexity, host
  competition, read-start hotspots and cell support.
- `interpretation_flags.tsv`: host-homology, low-complexity, ambiguity and
  hotspot diagnostics. All are diagnostic only.

A verdict compares the scores of emitted alignments; it is not a probability
and not a host-contamination call. Candidate search is limited to the
alignments the aligner emitted, so a single primary alignment cannot prove
specificity.

**Next action.** If the molecules are `host_best`, dominated by low-complexity or
reagent reads, or supported by one hotspot, treat the call as weaker. If they
are `virus_best` with broad coverage and complex bodies, the evidence for the
candidate is stronger. Neither outcome becomes a diagnosis. The
[QC and read evidence](vignettes/qc_and_read_evidence.ipynb) notebook runs this
evidence path on synthetic data.

### 5. Conclusion limits

`detected`, `detection_source`, `claim_scope` and the alignment labels are
evidence tiers assigned by fixed rules. They are not calibrated probabilities.
`claim_scope = screening_only` (Anelloviridae by default) means the call is a
screen, never a confirmed infection or a sensitivity claim.

Successful validation and threshold passage do not establish infection. A
candidate can be explained by ambient RNA, barcode swapping, doublets, uptake of
infected material, or host homology. ViralScan has no doublet caller, and a
single run cannot measure barcode swapping; the cross-run `scripts/swap_check.py`
(source checkout only) is a screening check, and without a pool manifest
swapping stays `not_assessed`.

**Negatives.** A negative depends on reference coverage and sensitivity, not only
on depth. Start with `results/sensitivity.tsv`:

- `informative_negative` is true only when the virus was not detected, depth is
  sufficient **and** a capture term was measured for that row. It is false on
  almost every run, deliberately; `negative_blockers` says which condition
  failed.
- `capture` is measured only by a valid positive control whose scope covers the
  row. Without one, `capture` stays empty and `capture_status` says why; the
  depth-only floor in `lod95_per_10k` is not a capture measurement.
- A different reference can change the answer. The bundled 20-genome
  Anelloviridae panel shares little 31-mer space with the full human anellovirus
  set, so a negative from that panel is not a statement about the sample. Build
  with `--anellovirus` when anellovirus sensitivity matters.

For a certifiable negative, plant a control at a known abundance and require it:
`--positive-control-gene`, `--positive-control-molecules`,
`--positive-control-scope exact_sequence`, `--positive-control-virus-key`,
`--positive-control-receipt` and `--require-positive-control`. The recovered
fraction is the capture term, and it certifies only the named target. The FAQ
entry [ViralScan found nothing](faq.md#viralscan-found-nothing-is-there-really-nothing-there)
gives the depth, capture and allocation terms behind this.

**Virus- and chemistry-specific limits.** Keep these with the candidate:

- Anelloviruses are a commensal virome; a detection is expected background. Use
  `viralscan evidence` on the specific accession before reporting one.
- HPV: presence and the E6/E7-versus-L1/L2 contrast are supported; per-genotype
  attribution from L1 or other late-region reads is not, because L1 cross-maps
  between genotypes.
- Latency: the assay is transcriptomic. A silent provirus, a
  transcriptionally silent integrated HPV genome and the HBV cccDNA pool produce
  no reads, so they are invisible, not latent.
- 10x 5′ libraries need `reverse` or `unstranded` strand handling; a 5′ library
  run as forward loses most reads.

**What you can and cannot say**

| Statement | Supported when |
|---|---|
| "This run produced candidate viral molecule support in N called cells." | Validation passed and the denominator is the one you intend |
| "The candidate passed the configured detection gate." | `detected = true` in this run |
| "Reads for this target compete with host sequence and carry complex bodies." | `viralscan evidence` ran with host and target FASTA and shows it |
| "The virus was absent." | Almost never from one run; requires an informative negative with measured capture |
| "These cells are infected." | Not from ViralScan alone; requires orthogonal confirmation |

Record which checks you ran, the reference digest, the run fingerprint, and any
unavailable field as unavailable. Then state the candidate, not the conclusion.
