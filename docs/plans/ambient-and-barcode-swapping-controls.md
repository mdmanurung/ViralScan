# Ambient RNA and barcode-swapping controls: design note (PLAN AUD-12, finding M-03)

Status: design only, written 2026-10-12. No code, no tests, nothing run: `src/` and `tests/` are
untouched. Claims about the current code cite the file and function they come from (read, not run).
Numbers that I could not verify are marked "unverified". Literature is deliberately not cited: no
reference in this note has been checked, so none is named beyond what the repository already cites.

Source of the finding: Biomni review, `analysis/audit_2026-10/findings.tsv` row `M-03` ("ambient
RNA and barcode/index swapping need separate controls; keep empty-droplet evidence separate from
called-cell denominator; test doublets/uptake"; proposed fix "run-specific background estimate,
negative-sample controls, report background FP rate"). Raised by Biomni, not independently
audited. Related: `AUD-07` (the called-cell detection gate, whose PLAN row says it "needs a
background-aware design (`M-03`), not a bare N-cells cutoff").

## 1. Problem

The detection gate (`AUD-07`) is deliberately permissive: a virus is `detected` when at least
`--detection-min-cells` (default 1) called cells each hold at least `--infected-min-molecules`
(default 1.0) of its molecules and the called cells hold at least `--detection-threshold` (default
1) in total (`detection.compute_stats`, `defaults.DEFAULTS`). One viral molecule in one called cell
passes. The package already says in three places that this is candidate evidence that "does not
rule out ambient RNA or barcode swapping" (`detection._gate_note`, `docs/faq.md`, the
`background_umis` row of `docs/output_reference.md`), but it offers no way to *assess* either
explanation. Three distinct mechanisms produce low-level viral signal that is not infection of the
cell carrying it:

1. **Ambient (cell-free) RNA.** Viral RNA in the suspension (virions, lysed infected cells) is
   captured in every droplet, so it appears in called cells and in empty droplets alike. Its
   expectation scales with the droplet's total molecules.
2. **Barcode or index swapping between libraries.** A molecule from library A is read out under
   library B's sample index (index hopping on a shared sequencing pool), or a chimeric read joins a
   barcode and an insert from different molecules. Observable signature: the same `(CB, UMI)` key
   appears in two libraries that share a whitelist and a sequencing run. I assert only the
   signature, not the mechanism, because the outputs cannot tell hopping from chimeras or from
   cross-contamination at library preparation.
3. **Doublets and uptake.** A droplet holding an infected and an uninfected cell, or a cell that
   engulfed infected material, carries real viral RNA that is not replicating infection in that
   cell. RNA counts alone cannot separate this from infection (section 9).

Where a *swapped* molecule lands matters, and it is easy to get wrong. A hopped molecule keeps its
cell barcode and UMI; only the sample index changes. In the recipient library it therefore sits at
a barcode chosen by the *donor* library, which is a called cell in the recipient only if the two
barcode sets happen to coincide. Expected consequence (to be measured, not assumed): swapped
molecules mostly inflate non-called barcodes of the recipient (so they show up in `background_umis`
and in rows of `per_cell_viral.tsv` with `is_called_cell` false) and only occasionally a called
cell. Ambient RNA is the mechanism that reaches called cells in bulk. The two therefore need
different tests, which is why the finding asks for separate controls.

## 2. What ViralScan already writes (verified by reading)

| Output | Holds | Use here | Limit |
|---|---|---|---|
| `results/viral_summary.tsv` `background_umis` | Viral molecules in non-called barcodes: `total_umi_raw - called_umi` (`detection.compute_stats`) | Size of the empty-droplet viral signal per virus | One scalar. Not per barcode, not a fraction, and it mixes soup, swapping and low-depth cells. `--cell-calling none` makes it 0 by construction (`output_reference.md`) |
| `viral_summary.tsv` `infected_called`, `n_called_cells`, `infected_subthreshold_*`, `accession_breadth`, `host_viral_ambig_fraction` | Called-cell counts and breadth | Numerators of any enrichment | `infected_called` is `>= infected_min_molecules` only |
| `results/per_cell_viral.tsv` | One row per infected barcode x virus, **including non-called barcodes** (`is_called_cell`), with `viral_molecules_total_est`, `molecules_total_est`, `viral_fraction` | Cross-run barcode test (3.2); per-barcode distribution | Only barcodes at or above `infected_min_molecules`; zero-viral cells are absent, so no full distribution |
| `results/called_cells.tsv` (`cellcalling.write_called_cells`) | Called barcodes (header `barcode`) | The donor set in the cross-run test | Two-step runs: the called set comes from the STARsolo host matrix |
| `kb-python/counts_unfiltered/adata.h5ad` (`kb_outputs.KbCountOutputs.adata`), `cells_x_genes.barcodes.txt` | Every observed barcode with its viral and host counts | The ambient estimate (3.1); the null for 3.2 | Unfiltered matrix, so empty droplets are present. A two-step virus-only matrix has no host totals |
| `kb-python/output.resolved.sorted.bus(.txt)` (`KbCountOutputs.resolved_bus`, `resolved_bus_txt`) | Corrected, CB-UMI-sorted BUS | Key-level swap test (3.3) | Large; needs a streaming reader |
| `viral_summary.tsv` is computed from `adata.X` selected layer (`multimapping.select_detection_matrix`) | Fractional weights possible | Same layer as `detected`, so flags are consistent with the gate | Fractional mass in empties is noisy |
| `evidence.interpretation_flags` | Static `contamination` row: status `not_assessed`, basis "requires negative controls and ambient model" | The natural place for the new flag to land | Today it is a constant, nothing feeds it |

Also verified: there is no doublet code in `src/` (grep for scrublet or doublet returns nothing);
`report.html.j2` suggests SoupX/DecontX for values below 0.1 % but ViralScan runs neither; the
`umap.py` nearest-neighbour permutation test exists, and `report.html.j2` states that significant
clustering "suggests true infection rather than ambient contamination". That sentence goes beyond
what a permutation test on neighbours can show (uptake and doublets cluster too) and should be
reviewed by the docs owner.

The synthetic generator already models both mechanisms, which is the validation route, with one
gap (section 7): `scripts/generate_truth_panel.py` `viral_reads` relocates a fraction of each
virus's molecules to `lay["uninfected"]` (uninfected cells plus empty droplets, line 601) with
`origin = ambient`, and `hop` copies molecules into another same-chemistry, same-partition sample
with `origin = hopped`. The protocol factor is `ambient_index_hopping` (levels 0, 0.01, 0.05,
`status: pending`, frozen at G3).

## 3. What can be measured without new data

### 3.1 Single run: the soup bound (tier 0)

All inputs are in the unfiltered matrix and the called mask. Let `A` be the ambient-range barcodes
(total molecules at or below `emptydrops_lower`, default 100). `emptydrops.R` passes `lower` to
`DropletUtils::emptyDrops`; my understanding is that `emptyDrops` treats barcodes at or below
`lower` as empty and estimates the ambient profile from them (unverified against the installed
version). `emptydrops.R` also labels every barcode at or above the knee a cell regardless of FDR
(its header comment), so a very viral-rich droplet is never excluded by the ambient test. Define for
virus `v`:

- `V_A`, `T_A`: viral and total molecules in `A`; `f_hat = V_A / T_A`.
- `f_U`: one-sided 95 % upper confidence bound on `f`, treating `V_A` as Poisson with exposure
  `T_A`. For `V_A = 0` this is `2.996 / T_A`, the same constant `sensitivity.py` uses for its LOD.
  The bound, not the point estimate, goes into the test: a clean run with zero empty-droplet viral
  molecules would otherwise flag every called-cell molecule as "above soup".

Under a soup-only hypothesis with contamination fraction `rho` in `[0, 1]`, a called cell `i` with
`N_i` molecules receives `Poisson(rho * f * N_i)` viral molecules. Setting `rho = 1` (the whole
droplet is soup) gives the largest soup expectation any `rho` can produce:

- `E1 = sum_i P(Poisson(f_U * N_i) >= t)` over called cells, `t = infected_min_molecules`, the
  expected number of called cells that soup alone could push over the infected rule;
- `E2`, the same with `>= 2` molecules (soup rarely puts two molecules in one cell, infection
  often does);
- `mu = f_U * sum_i N_i`, the expected soup molecules in called cells.

Compare `infected_called` with `E1` (and the called-cell molecules with `mu`) by an upper Poisson
tail. If the observation is above the tail at the declared level, soup alone cannot explain it at
any `rho <= 1`. This needs no calibration, which is its only virtue: it is a **weak, necessary
screen**. It cannot show that a signal is soup; passing it does not show infection; and for a
signal that is enriched over soup by a small factor it will not discriminate. A tighter test needs
an estimate of `rho`, which is the calibration problem of section 6.

Comparing each cell's `viral_fraction` with `f_hat` is not a test and is not proposed: with one
molecule in a cell of a few thousand, almost every infected cell "exceeds" the ambient fraction by
chance. The comparison is aggregate, against a Poisson expectation, with `f_U`.

Per-barcode descriptors from the same vectors (computed where full vectors exist, in
`compute_stats`, because `per_cell_viral.tsv` omits zeros):

- Spearman correlation of viral molecules with total molecules over called cells. Soup scales with
  droplet size; a confound is that infected cells may be larger.
- Share of called-cell viral molecules in the top cell, and the count of called cells with `>= 2`
  molecules against `E2`.

Two-step runs: the matrix handed to detection holds only barcodes with non-host reads, so `T_A`
cannot come from it. The host totals exist per barcode in STARsolo `GeneFull/raw`
(`cellcalling.read_host_matrix_totals`), but `compute_stats` receives only `host_cells`, which
carries the `called` and `comparable` barcode sets, no totals. Either a new field or a re-read is
needed (decision 5). Virus-only indices (`index_kind == virus_only`, no host depth) are
`not_assessed`, as `_comparable_called_cells` already treats them.

### 3.2 Two runs: barcode overlap (tier 1, no BUS needed)

Applies to 10x runs on the same sequencing pool and the same whitelist; not to Drop-seq (random
bead barcodes, no shared whitelist). For an ordered pair of runs `A`, `B`:

- `S_B`: barcodes in `B`'s `per_cell_viral.tsv` with `is_called_cell` false, for virus `v`.
- `D_A`: barcodes that carry `v` in `A` (called or not, as listed in `A`'s `per_cell_viral.tsv`),
  and the broader `A` called set from `A`'s `called_cells.tsv` as a secondary reference.
- Null: among *all observed non-called barcodes of `B`* (rows of `B`'s
  `cells_x_genes.barcodes.txt` minus the called set), the fraction `q` that lie in `D_A`.
- Test: one-sided hypergeometric (Fisher) of `|S_B ∩ D_A|` against `q`.

A hopped molecule carries the donor's barcode, so overlap far above `q` is a signature that needs
no UMI information and no BUS reading. The expected overlap by chance is of the order
`|D_A| / |whitelist|`; I have not verified whitelist sizes here, so compute `q` from the run's own
barcode lists, never from recalled constants.

Caveats: rows exist only at or above `infected_min_molecules` (fractional weights below it are
invisible); barcode suffixes (`-1`) must be stripped the same way on both sides
(`cellcalling._strip_suffix`, and audit finding `B-10` shows the failure when they are not);
`A` and `B` must be different libraries. Two lanes of one library share barcodes legitimately, and
a declared technical replicate would read as 100 % swapping.

### 3.3 Two runs: key-level swap and library hop rate (tier 2, BUS needed)

`output.resolved.sorted.bus.txt` is CB-UMI sorted, so two runs can be merge-joined in one
streaming pass with constant memory.

- **Viral key sharing.** Viral molecules `(CB, UMI)` in `B` that also occur as a key in `A`
  (any equivalence class), against a chance rate derived from the number of keys in `A` and the
  size of the key space. Decisive when positive: independent molecules do not share a barcode and
  a UMI.
- **Library hop rate.** The same join over *host* molecules gives `h(A->B)`, the shared-key
  fraction after subtracting chance. The prediction for virus `v` is then `h * (viral molecules in
  A)` shared keys in `B`; the observation is compared with it. This separates "this library leaks
  at rate `h`" from "this virus is special".

The existing record: PLAN `EXPL-HPV16` notes that the normal-raft HPV16 molecules look like hopping,
"4 of its 9 clean HPV16 CB-UMI pairs also occur in the HPV16 raft", done by hand on one pair and
not re-verified here. It is the one real-data observation in the repository of this signature, and
it is also an example of tier 2 detecting something the default gate would call positive (section 4
of the M-04 note covers what that means for scoring).

### 3.4 What none of this measures

Cell-free viral RNA that is *not* soup-like (for example enriched at specific barcodes because of
droplet-level carry-over) and swapping that preserves barcode and UMI inside one library (barcode
collisions, `cb_umi_collisions` in the generator) are not covered by 3.1 to 3.3.

## 4. What needs extra controls (tier 3)

These are wet or design inputs, so they are user-owned. The protocol already lists two as optional
datasets without a source: `reagent_control` and `empty_droplet_control`
(`analysis/v3_validation/protocol.yaml`, `required: false`, "report the row as unsupported if
unavailable at freeze").

| Control | What it measures | What it cannot |
|---|---|---|
| Same-pool virus-negative library (different donor or cell line with the virus known absent, on the same lane or flowcell) | The empirical cross-library leak rate into a library that cannot be infected: any viral key it holds came from elsewhere. Gives `h` with a known-zero truth | Needs a pool designed with it, and it must share the whitelist |
| Reagent control (no cells) | Kit or index contamination by the virus | Not ambient soup of the real sample |
| Empty-droplet-only comparison (the existing `A` set) | The sample's own soup composition | It is not independent of the cells' own lysis |
| Cell-free spike (a sequence absent from host and panel, added to the suspension at a known amount) | Soup capture and `rho`, directly | Recipe unspecified; not available today. The `--positive-control-gene` machinery plants *molecules at the read level*, and its `over-recovered` state (`detection.measure_positive_control`: "reads leaked from another cell") already treats leakage as a confound, which shows the control is not spike-in-specific |
| Species-mixed (barnyard) library | Soup and doublet rates in one experiment | No such dataset is registered in the protocol |
| External doublet score | A per-barcode doublet call from the host matrix | None exists in ViralScan; see 9 |

## 5. Proposed columns and flags

Names are proposals, to be fixed at implementation. All are descriptive. None changes `detected`,
`pct_infected_called` or any denominator, so empty-droplet evidence stays separate from the
called-cell denominator, as the finding asks. Outputs documented in `docs/output_reference.md`
(pinned) need a separate docs commit at implementation.

`viral_summary.tsv`, per virus:

| Column | Definition | Applicable when |
|---|---|---|
| `ambient_viral_fraction` | `V_A / T_A` | An ambient-range set exists and host depth is known; else `not_applicable` |
| `ambient_viral_fraction_ub95` | `f_U` | Same |
| `n_ambient_barcodes_with_virus` | Barcodes in `A` with at least one viral molecule | Same |
| `soup_expected_cells_max` | `E1` | Same |
| `soup_expected_cells_ge2_max` | `E2` | Same |
| `soup_tail_p_cells` | Upper Poisson tail of `infected_called` given `E1` | Same |
| `viral_depth_spearman` | Viral vs total molecules over called cells | At least 3 infected called cells, else blank |
| `ambient_assessment` | `not_assessed`, `soup_sufficient` (observation at or below the tail), `soup_rejected_rho1` | `not_assessed` for `--cell-calling none`, `knee`, no empties, virus-only without host depth |
| `ambient_basis` | One sentence naming the reason, in the style of the evidence flag table | Always |

Blank or `not_applicable`, never 0, where the quantity is undefined (the convention of
`pct_infected_comparable` and of `n_cells_reactivating`).

Cross-run, from a new sample-sheet-driven step (not part of a single run), one row per
`(recipient run, donor run, virus)`: `n_nonCalled_virus_barcodes`, `overlap_with_donor`,
`overlap_expected`, `overlap_p`; with BUS: `shared_viral_keys`, `shared_keys_expected`,
`library_hop_rate`. The sample sheet must declare which runs shared a pool and a whitelist; ViralScan
cannot infer that.

`evidence.interpretation_flags`: replace the static `contamination` row by one derived from
`ambient_assessment` and, when present, the cross-run table. Status stays `not_assessed` where the
inputs do not exist; it never reads "clean".

## 6. Thresholds that must be calibrated, and on what

The first release proposed here has **no tuned threshold**: the ambient flag uses only the
calibration-free `rho = 1` bound and a conventional tail level. That is on purpose. The protocol's
threshold rule is restrictive (`calibration.threshold_search`: cutoffs come from a pre-declared
grid, training partition only, "no value outside this grid may be adopted"), so any numeric cutoff
that gates a call would need a protocol amendment and a deviation record. The following would need
calibration before they could gate anything:

| Quantity | Why it is not free | Calibrate on |
|---|---|---|
| Contamination fraction `rho` and its spread across cells | Decides how weak the `rho = 1` bound is; varies by cell type and size | Real presumed-negative libraries (`pbmc_10x_healthy_v3`, `VAL-04`; the tonsil-atlas negative subset, `DSR-07`) for the host side; planted truth for the viral side. No value is available in the repository |
| Upper-bound level (0.95) | Convention only | Planted `origin` labels (see 7) |
| Ambient-range cutoff (`emptydrops_lower`, default 100) | A depth-dependent proxy; chemistry and depth shift it | Sweep over the same real libraries |
| Minimum ambient molecules for `soup_sufficient` to be meaningful | With 0 observed the bound is `2.996 / T_A`, which may exceed the signal | Same |
| Hop-rate cutoff for a cross-run flag | No distribution exists | Same-pool negative libraries (section 4); planted CB-preserving hops (section 7) |
| Any enrichment cutoff on `viral_depth_spearman` or `E2` | Cell type confounds depth | Planted data with `origin` in `{cell, ambient, hopped}` |

Calibration must follow the protocol's separation: training partition only, one holdout look,
outcome-independent cell anchors, the library as the resampling unit.

## 7. Validation, and a generator gap

The synthetic panel can score the flags because it records the origin of each planted molecule
(`origin` in `{cell, ambient, hopped, artefact:<class>}`, a required `truth_manifest.tsv` column).
Two limits, both verified in `scripts/generate_truth_panel.py`:

1. **`hop()` does not preserve the barcode.** The copy gets `rng.choice(self.layout(target)
   ["uninfected"])` as its cell barcode and keeps only the UMI (and sequence). That is a model of
   "a molecule appears in another sample", not of index hopping, which keeps both barcode and UMI.
   Consequently the generator cannot validate tier 1 or tier 2, whose signature is the shared
   barcode or key. A CB-preserving mode is needed: hopped molecules keep the donor `(CB, UMI)` and
   land in the recipient's barcode space (called cell if the barcode is in the recipient's cell
   pool, else a non-called barcode that the recipient's matrix then holds). The default output
   must stay byte-identical (the AUD-18 hash test sets the precedent). The generator's barcodes
   are random and not on a 10x on-list (its docstring), so the barcode-coincidence rate would also
   have to be modelled, or the generator given `--onlist` (it already has one, AUD-18).
2. **The ambient factor relocates a fraction of each virus's molecules.** It is not a soup
   concentration. Mapping its levels (0, 0.01, 0.05) to a contamination fraction needs an
   explicit statement, and empty droplets carry "ambient-only content, log-normal UMI total,
   median 50" (generator parameters, not measurements; `docs/plans/2026-10-05-val01-generator-design.md`
   section 4).

Acceptance for the validation (to be fixed before running): on planted data, the flag must read
`soup_rejected_rho1` for planted `cell`-origin molecules above a stated abundance and must not
read it for `ambient`-origin molecules at the ambient levels; and cross-run tests must detect
planted CB-preserving hops without flagging independent pairs. Outcome runs stay blocked until the
protocol freeze (`freeze.training_execution_allowed: false`); this note authorises none.

## 8. Failure modes

1. **Soup is not homogeneous.** Lysis of infected cells makes droplets near them richer in virus;
   `rho` varies by cell type and size; `f` from the empties need not equal `f` in the cell's own
   soup share.
2. **In a genuinely infected sample the soup is real signal.** A virus with high `f` from lysed
   infected cells can still be cell-associated in other cells. The screen then returns
   `soup_sufficient`: that is "not separable", not "not infected".
3. **Zero or tiny empty-droplet counts.** `f_hat = 0` gives no information; hence `f_U`. If `T_A` is
   small the bound is large and nothing is rejected.
4. **Cell calling interacts with the estimate.** `--cell-calling none` removes the empties; `knee`
   is sensitivity-only and calls empties as cells (`cellcalling` header); two-step runs call cells
   on the host matrix. A different caller moves `A` and the called mask together.
5. **Layer choice.** Fractional EM or multimap weights in empties add noise; use the same layer as
   `detected` and say so.
6. **Swapping versus other shared-key causes.** Key sharing is a signature, not a mechanism: pooled
   amplification chimeras, library-prep cross-contamination and technical replicates all produce
   it. Same-library pairs must be excluded by declaration.
7. **A hop into a called cell is rare by construction but not impossible.** Tier 1 has power only
   when `|S_B|` is large enough; a run with five non-called viral barcodes has little.
8. **Barcode-suffix and whitelist bookkeeping errors** turn overlap into a silent zero (audit
   `B-10`).
9. **Thin Poisson null.** The null treats soup counts as independent across cells; PCR-duplicate
   UMIs and overdispersed soup make the tail optimistic (too many soup rejections), the unsafe
   direction. The tail level is therefore reported, not hidden.
10. **Report wording.** A flag column must not reach `summary.txt` as "clean".

## 9. What a result may and may not support

| Result | May support | May not support |
|---|---|---|
| `soup_rejected_rho1` | Viral molecules in called cells exceed what soup at any `rho <= 1` would give at this confidence | Infection; replication; that the cells are the source (doublets, uptake) |
| `soup_sufficient` | The signal is not distinguishable from soup in this run | That the virus is absent or the signal is soup |
| Significant tier-1/2 overlap | These molecules (or barcodes) have a shared origin with another library on the pool; quantifies leak into this library | Which mechanism; that the leak is the only source in this run |
| No overlap with the one donor tested | Nothing about other donors; no leak from that library | Absence of swapping from untested libraries |
| `ambient_assessment = not_assessed` | That the run cannot be assessed (reason in `ambient_basis`) | That the run is clean |

Doublets and uptake: RNA alone cannot tell uptake from infection, so the limit is stated, not
engineered away. Evidence of replication (spliced or lytic transcripts, strand, programme state:
`gene_programs.py`, `evidence.py`) is a different axis and is not claimed by this note. A doublet
score would be an optional external per-barcode table, reported as the fraction of infected called
cells flagged, never as a filter; ViralScan should not grow its own doublet caller (decision 7).

## 10. Recommendation

Deliver in this order, stopping when the data stop justifying the next step:

1. Tier 0: the per-virus ambient columns and the `ambient_assessment` flag, descriptive, feeding the
   `contamination` flag. Pure functions over vectors already inside `compute_stats`; testable
   offline on tiny arrays. No new input.
2. Generator: CB-preserving hop mode (section 7), because tiers 1 and 2 cannot be validated without
   it.
3. Tier 1 as a cross-run step over existing outputs, then tier 2 only if tier 1 on real pooled runs
   shows a signal worth the BUS reader.
4. Calibrate any numeric gate only on frozen training data, after the protocol freeze.

The first step is small enough to implement and test offline; it is not implemented here because
the task restricts this change to design notes.

## Documentation rows implied

For the docs owner (pinned files, not edited here): `docs/output_reference.md` rows for the new
`viral_summary.tsv` columns when implemented; the `background_umis` row should point to them; the
`report.html.j2` sentence about clustering and ambient contamination (section 2) needs review; the
`evidence.interpretation_flags` `contamination` basis text changes with step 1.

## Decisions needed from the user

1. Descriptive only for now (columns and a flag, `detected` unchanged), or may a future flag demote
   `detected` or enter a tier? Recommended: descriptive; anything that gates needs a protocol
   amendment and calibration on frozen training data.
2. Approve the calibration-free `rho = 1` bound with the Poisson upper confidence bound as the only
   flag logic in the first release, deferring an estimate of `rho`.
3. A cross-run swap step: a new subcommand or a `scripts/` tool, and who declares which libraries
   shared a pool and a whitelist (ViralScan cannot infer it)?
4. Approve a CB-preserving hop option in `generate_truth_panel.py` (default output byte-identical;
   the generator is a pinned file, so a re-pin follows).
5. Two-step runs: carry per-barcode host totals into `compute_stats` (extend `host_cells`) or
   re-read `GeneFull/raw`?
6. Which wet controls to request or select for the optional protocol rows: same-pool negative
   library, reagent control, cell-free spike, species mix. Without one, `contamination` stays
   `not_assessed` for pooled-run swapping.
7. Doublets and uptake: out of scope (claim limit only), or an optional external per-barcode table?
