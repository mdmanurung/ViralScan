# ViralScan archived v2 versus v3 identical-BUS diagnostic

## Decision summary

The archived ViralScan 2.2.0 results are reproducible from their retained
matrices and summaries. The fork also completed molecule-aware re-counting for
all 44 frozen BUS rows. The two outputs are not numerically interchangeable,
and the reason is mechanical rather than a matter of convention.

The archived value is the sum of two quantities in different units. Its
unique-mapping term, `counts_original`, comes from `kb count` and is
UMI-deduplicated. Its multimapping term, `counts_corrected`, is built in
`multimap.py` by distributing each BUS record's read `count` fractionally across
the genes of its equivalence class; the UMI column is parsed and then never
used, so no deduplication occurs on that path. The multimap term is also derived
from the raw `output.bus` rather than the corrected BUS that produced
`counts_original`, so its barcodes are not whitelist-corrected. The summary then
reports the multimap-only total under the label `Total viral UMIs (corrected)`,
a quantity that contains no UMI collapse.

v3 reports corrected-cell-barcode/UMI molecules throughout, with explicit
ambiguous-allocation layers. Because the archive's mixing ratio between its two
terms varies per sample with duplication rate and multimapper fraction, no
rescaling puts the two on a common scale. This is a mechanistic incompatibility,
not evidence that either is more accurate.

The strongest supported conclusions are:

1. the retained v2 result tables are internally reproducible;
2. v3 is deterministic on the repeated pilot and conserves molecule mass on
   the completed cohort;
3. all three EBV positive controls remain qualitatively positive in v3;
4. expected HIV-1 recovery is not demonstrated because the frozen 106-virus
   v3 feature universe has no HIV-labelled feature;
5. low-mass unexpected skin signals persist under some v3 allocation
   endpoints and must be treated as specificity concerns, not confirmed
   infections; and
6. a matched-FASTQ, fresh-v2-versus-fresh-v3 comparison is still required to
   separate historical pipeline effects from counting changes.

This is a diagnostic-development audit of outcomes that were already viewed.
It is not confirmatory evidence and does not establish sensitivity,
specificity, biological truth, viral absence, disease association, package
superiority, or release/publication readiness.

## Cohort and execution

- Frozen denominator: 44 technical rows representing 42 logical inputs.
- Classes: 3 EBV controls, 3 HIV-control technical rows representing 2
  logical inputs, and 38 skin rows.
- Chemistry is confounded with class: all six control rows use 10x v2
  whitelists and all 38 skin rows use 10x v3 whitelists.
- All 44 rows reached wrapper status `success`.
- The aggregate contains 4,664 run-virus rows, exactly 44 rows by the frozen
  106-virus universe.
- The aggregate failure table has zero data rows.

The first full-array attempt was canceled and retained after an output-layout
defect was detected. The corrected full-wrapper arrays used the frozen
environment, source, wheel, references, inputs, and scientific parameters.

## Archived v2 reconstruction

The audit contains 468 checks: 356 gene values, 68 virus values, and 44 total
viral-load values. All 459 reported numeric values match their reconstruction
within the frozen absolute tolerance of `1e-6`; the maximum absolute delta is
`1.397e-9`. The other nine rows reconstruct to zero and correspond to archived
summaries that omitted a no-call total.

This establishes internal reproduction from retained v2 artifacts. It is not
yet a fresh execution of Emma's original package from the original FASTQs.

## V3 molecule accounting

Across all 44 technical rows:

| Quantity | Total |
|---|---:|
| Unique molecule mass | 1,295,627,351 |
| Allocated ambiguous mass | 27,282,599 |
| Selected matrix mass | 1,322,909,950 |

For every row, selected matrix mass equals unique molecule mass plus allocated
ambiguous mass exactly. Ambiguous allocation contributes 2.062% of the pooled
selected mass; its per-row share ranges from 1.162% to 2.652%.

These whole-matrix totals are dominated by host features. They must not be
compared as though they were viral loads from v2.

## Expected-control findings

### EBV

Epstein-Barr virus is nonzero in all three EBV controls at unique, equal, and
host-conservative endpoints:

| Run | Legacy EBV value | V3 unique | V3 equal | V3 host-conservative |
|---|---:|---:|---:|---:|
| SRR12682296 | 1,792,860 | 256,487 | 874,039 | 874,039 |
| SRR12682297 | 262,394 | 36,287 | 121,749 | 121,749 |
| SRR12682298 | 31,011.3333 | 4,501 | 18,345 | 18,345 |

This supports qualitative positive-control recovery. Magnitudes are not
directly commensurate because the endpoints use different units and counting
rules.

### HIV-1

The frozen v3 aggregate universe contains no HIV-labelled feature. Expected
HIV-1 recovery is therefore not demonstrated. Under the frozen interpretation
rule this is a sensitivity-failure description, not evidence that HIV is
biologically absent. The three HIV-control rows contain only low-mass
non-target viral endpoints: totals of 17 unique, 20.667 equal, and 20
host-conservative molecules.

## Skin findings

The 38 skin rows have a summed archived viral load of 118.8014, with median
2.0211 and range 0 to 9. V3 produces:

| Endpoint | Viral mass | Nonzero run-virus entries | Runs with any nonzero entry |
|---|---:|---:|---:|
| Unique | 36 | 33 | 22 |
| Equal | 89.0559 | 94 | 38 |
| Host-conservative | 76 | 59 | 32 |

The largest host-conservative skin totals are Cercopithecine herpesvirus
(37), human herpesvirus 8 (10), human herpesvirus 6b (8), Epstein-Barr virus
(7), molluscum contagiosum virus (7), and human herpesvirus 1 (3); every other
virus total is at most 1.

Of 50 legacy-reported skin run-virus entries, unique retains 33 and loses 17;
equal and host-conservative each retain 43 and lose 7. Host-conservative adds
16 nonzero entries outside the legacy-reported set, while equal adds 51.
Cercopithecine herpesvirus is the only host-conservative label repeated across
multiple archived lanes of the same skin sample: it recurs in both observed
lanes for seven sample identifiers.

These are nonzero endpoint entries, not validated positive calls. The
unexpected labels and low mass are specificity concerns that require
read-level and external-truth follow-up. No calibrated positivity threshold
was available, so sensitivity, specificity, and call-set accuracy were not
estimated.

## Allocation effects

Across all 4,664 run-virus entries, unique is nonzero in 43 entries, equal in
109, and host-conservative in 74. Equal allocation exceeds unique in 71
entries; host-conservative exceeds unique in 36. Equal exceeds
host-conservative in 61 mostly low-mass entries totaling 13.7226, while no
host-conservative entry exceeds equal.

These differences describe how ambiguous molecules are assigned. They do not
measure method accuracy without truth labels and aligned decision thresholds.

## Cell-level concordance

The legacy and v3 barcode universes match exactly in all 44 rows (barcode
Jaccard 1.0). Using only a nonzero viral-mass indicator, not a diagnostic
positivity threshold, median legacy-to-v3 cell-set Jaccard is 0.052 for unique,
0.909 for equal, and 0.128 for host-conservative. The class-specific median
host-conservative Jaccard is 0.975 for EBV controls, 0.455 for HIV controls,
and 0.118 for skin.

This shows that barcode selection did not drive the observed differences.
Ambiguous allocation strongly changes which cells have nonzero viral mass,
especially in skin. It does not show which endpoint is biologically correct.

## Technical repeats

Both declared repeat pairs have identical aggregate run metrics and identical
106-virus vectors at all three v3 endpoints (L1 distance zero). This supports
repeat-level reproducibility for the archived inputs. The pairs are technical
only and are not independent biological replication.

## Remaining work

The completed identical-BUS comparison isolates counting behavior only
partially. The following steps remain blocked on a separately approved,
immutable control-FASTQ acquisition:

1. freeze identical FASTQs for the three EBV and two HIV logical controls;
2. execute fresh ViralScan 2.2.0 and frozen v3 on those same reads and
   references;
3. audit processed reads, pseudoalignment, barcodes, matrices, target recovery,
   and exact-read evidence; and
4. retain all failures without changing thresholds or scientific parameters.

Until those steps and external truth-panel/comparator/calibration work are
complete, this report remains diagnostic only.
