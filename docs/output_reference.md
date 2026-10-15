# Output Reference

Runs record `software_identity` in the root `run_manifest.json` and in
`results/reference_provenance.json`: the executing package's content digest, plus Git
commit and dirty-source status when run from its own checkout.
This descriptive field does not alter the resume fingerprint. Post-count library checks
write `results/chemistry_sanity.json` and the manifest's per-sample `chemistry_sanity` block,
including failed or unavailable checks. `--require-chemistry-sanity` stops on errors or
unavailable checks after retaining this receipt; warnings alone remain advisory.
When the index has a build manifest, the manifest's per-sample `catalogue_drift` block records
whether the packaged virus catalogue is the one the index was built against (`status` is
`match`, `drift`, `not_pinned` or `unavailable`, with both SHA-256 values); a drift also logs a
warning, and virus names still come from the packaged catalogue.

EmptyDrops calls cells on the host transcriptome only, never on viral reads, so the called set
is the same for every analysis of a library. Two-step EmptyDrops uses the complete STARsolo
`GeneFull/raw` host matrix. Combined mode uses the host gene columns of the kb matrix (viral
columns are dropped, `count_layer` `kb.host_genes`); the `total` column of
`kb-python/counts_unfiltered/emptydrops_cells.tsv` is then host UMIs. A combined run whose viral
gene partition cannot be determined stops with an error rather than calling cells on every
column; a matrix column the Virus Identity table does not list stops it too. Both modes record the
input matrix and called-cell-list digests in `results/cell_calling_input.json`, together with the
seed and the `effective_seed` R used (a seed above 2147483647 is reduced modulo 2147483647;
smaller seeds are used as given), the R and DropletUtils versions and the code revision
(`code_sha`, null outside a git checkout, and `build_sha256`); combined mode also records the
partition file and its SHA-256 (`partition_source`, `partition_sha256`). A re-run deletes the
previous `emptydrops_cells.tsv` and receipt first, so a failed run leaves neither, and a run
that writes no receipt (`--cell-calling none` or `knee`, or `external` in combined mode) removes
the old one. The combined-mode host columns depend on the reference's
D-list (swapping `cat42b` for `cat42d` moved EmptyDrops calls by -22 % on the HPV16 rafts,
`REF-06`); STARsolo host `GeneFull` counts (two-step) and an external list are independent of
the reference.
From `CELLS-02` the `knee`, `inflection` and `total` columns of a combined run's
`emptydrops_cells.tsv` are host-only, so the dossier's `median_umi_called`, `knee` and
`inflection` (`dsr05_dossier.py`) are host UMIs there, not host plus viral. A
`called_cells.tsv` or `emptydrops_cells.tsv` written before `CELLS-02` carries host plus viral
semantics and must not be reused as a cell anchor.
An explicit external cell list supplies the full called-cell denominator, including cells
with no viral molecules. Comparable metrics use host counts at the existing host-depth
floor; they remain unavailable when the host matrix or gene partition is unavailable.

Grouped host-response evaluation writes `hostresponse_manifest.json` alongside its tables.
It records input and output digests, software identity, raw-depth source, fold memberships
and training features, seeds, cell-type strata, and the declared permutation settings.
Failed strata and permutation replicates remain visible; these software diagnostics do
not establish scientific validity. The existing cell-split panel metric retains its
documented descriptive interpretation; group-split panel selection uses training folds.
The default `--cv auto` resolves to `group` or `cell` before any fit; the manifest's `cv_mode`
is the concrete mode, `cv_mode_requested` is what was asked for, and `cv_auto_fallback` is
the reason when an auto-selected group column could not be used (otherwise null).

ViralScan writes one sample directory under the path passed to `--output / -o`.
The sample directory is inferred from the R1 FASTQ filename before the first
underscore. For `sample_R1.fastq.gz`, the run directory is `output/sample/`.

---

## Directory layout

```
output/
└── sample/
    ├── config.yaml
    ├── summary.txt
    ├── report.html
    ├── log/
    │   ├── analysis.txt
    │   ├── detection.done
    │   ├── kb.done
    │   ├── kb_count.log
    │   ├── multimap.done
    │   └── umap.done
    ├── kb-python/
    │   ├── counts_unfiltered/
    │   │   ├── adata.h5ad
    │   │   └── adata_multimap.h5ad
    │   ├── output.bus
    │   ├── output.resolved.sorted.bus
    │   ├── output.resolved.sorted.bus.txt
    │   ├── run_info.json
    │   └── ...
    ├── plots/
    │   ├── <virus>_histogram.png
    │   ├── SuperExpressor_<virus>.png
    │   ├── umap_binary.html
    │   └── umap_continuous.html
    └── results/
        ├── virus_identity.tsv
        ├── viral_summary.tsv
        ├── sensitivity.tsv
        ├── positive_control.json
        ├── gene_program_summary.tsv
        ├── gene_program_cells.tsv
        ├── per_cell_viral.tsv
        ├── called_cells.tsv
        ├── multimap_evidence.tsv
        ├── molecule_assignments.tsv.gz
        ├── anello_alignment_by_accession.tsv
        ├── anello_alignment_by_virus.tsv
        ├── cell_type_enrichment.tsv
        └── reference_provenance.json
```

`reference_provenance.json` records the viral reference used (index/t2g/GTF,
technology, multimap settings, and the viral accessions in the reference and
detected) so results are traceable to their annotation.
`gene_program_*.tsv` are present only when `--gene-programs` is supplied; see
`results/gene_program_summary.tsv` below.
`called_cells.tsv` lists the barcodes cell calling kept (header `barcode`), the
set every `*_called` rate is over and the set layer 2 scores (since 2026-10-03,
PLAN `PROG-17`).
`cell_type_enrichment.tsv` is present only when `--cell-types` is supplied.
`multimap_evidence.tsv` is present only when multimapping is enabled.
`molecule_assignments.tsv.gz` is present only with `--multimap-molecule-assignments`.
`log/analysis.txt` lists the viral gene IDs of the index, one per line, sorted, with version
suffixes as stored. They come from the index build manifest (`<index>.build_manifest.json`,
`viral_gene_ids`) when it exists, plus any `-gtf` gene IDs and the anellovirus candidate IDs;
only an index with no build manifest falls back to the bundled GTF panel
(`viralscan data fetch`). The manifest's optional `catalogue_sha256` is the sha256 of the
packaged `virus_catalog.tsv` at build time, and a run warns when the packaged catalogue no
longer matches it; manifests without the field remain valid.
The run manifest records this comparison in `reference_catalogue`: the
`expected_sha256` from the index build, `observed_sha256` from the running package,
and `matches` (`true`/`false`, or `null` when either hash is unavailable).
This provenance is outside the resume fingerprint; drift remains a warning.
`results/evidence/<genus>/` and `log/auto_evidence.done` are present only with
`--auto-evidence`; see `results/evidence/<genus>/` below.
`anello_alignment_by_accession.tsv` and `anello_alignment_by_virus.tsv` are
present only when the anellovirus alignment branch ran (see below).
UMAP files are present only when `--umap` is supplied. `host_filtered/` is
present only when `--host-filter` is supplied.
`log/kb_count.log` retains the count command's output. Counting writes into a
unique temporary directory and checks for a non-empty `counts_unfiltered/adata.h5ad` before replacing
`kb-python/`. A failed count preserves the previous count directory, sample
configuration, and user notes.

`virus_identity.tsv` is written by every run; see
`results/virus_identity.tsv` below.

The raw `output.bus` is not a v3 counting input. ViralScan barcode-corrects it
when a whitelist is supplied, sorts it with bustools, and retains the resolved
sorted BUS plus text representation used for molecule counting and read lineage.

The v3 STAR filter also writes `host_filter_audit.tsv` with input, retained,
and removed fragment totals, plus `fragment_lineage.tsv.gz` with one row per
input fragment, in input order. Its columns are `read_id`, `filter_decision`
and `reason`:

- `retained` / `host_unmapped`: STAR left the pair unmapped to the host, so it
  reached viral quantification.
- `removed` / `host_mapped`: the pair is absent from STAR's unmapped output,
  meaning it was host-aligned or ambiguous. STAR runs with `--outSAMtype None`,
  so no finer subclass is available. Removed rows are derived as input IDs
  minus retained IDs.

ViralScan validates mate synchronization before and after filtering. For
validation, the Python-only `viralscan.scripts.host_filter.lost_truth_counts`
counts truth-viral fragments and molecules removed at this boundary (protocol
endpoints D15 and D16). Given the `read_id`s of `mixed_host_virus_manifest.tsv`
(`mixed_manifest_tsv=`) it also returns `mixed_fragments` and `removed_host_virus_mixed`,
the host-virus-ambiguous boundary count: the manifest's fragments the lineage saw, and
those the boundary removed. Both are `None` without the manifest, and mixed reads the
lineage never saw are not counted.

**Kind column.** Every column table below labels each field with one of five
kinds, so an observed count is never read as a model output or a conclusion:

| Kind | Meaning |
|------|---------|
| observation | A count, identifier or measurement taken directly from the inputs or the reference, with no model between it and the data |
| model estimate | A value that depends on a model choice (multimap allocation, capture model, statistical test, classifier) and moves if that choice changes |
| evidence tier | A categorical level assigned from molecule-level evidence by a fixed rule; never a calibrated probability |
| diagnostic flag | A QC, provenance or caveat field that qualifies other fields; it is not itself a result |
| biological interpretation | A statement about infection, latency or biology; v3 emits none as a column and requires orthogonal confirmation |

---

## `run_complete.json`

Written at the run root (beside `run_manifest.json`) only after every sample has finished, and again after `rerun-multimap`. It is the completion marker: a run without it was interrupted or is still in progress. Fields: `schema_version`, `run_fingerprint` (must equal the manifest's), `viralscan_version`, `samples`, `completed_at`, and `artifacts`, a map of `<sample>/<path>` to sha256 for the narrow set `results/viral_summary.tsv`, `results/virus_identity.tsv`, `results/multimap_evidence.tsv` and `kb-python/counts_unfiltered/adata_multimap.h5ad` (each only if present).

The marker is deleted whenever a run starts, resumes or overwrites, and re-computed by in-place commands (`rerun-programs`, `hostresponse`) when one existed. New runs set `completion_marker: true` in `run_manifest.json`; `viralscan validate-run` then treats a missing marker, a fingerprint mismatch or an artifact hash mismatch as an error. Runs whose manifest lacks the field only get a warning for a missing marker. Schema: `schemas/v3/run_complete.schema.json`.

`run_manifest.json` also carries `defaults_status`. It is `provisional` on every `3.0.0.devN` build: the detection defaults are not yet validated against the truth panel (gate G5), so any number from such a run is provisional. It becomes `final` from the release candidate. Like `completion_marker`, the field is outside `run_fingerprint`, so `--resume` still matches older manifests.

Two more `run_manifest.json` fields are also outside `run_fingerprint`. `input_paths` maps each hashed input (`<option>:<index>`, as in `input_fingerprints`) to its resolved absolute path, so `viralscan validate-run` finds the inputs from any working directory; manifests without it fall back to the option value. `fetched_reference` is written by an `-acc` run once the NCBI FASTA and GTF have been fetched or read from the cache: the accessions and the sha256 of both files. On `--resume` a different hash is refused (the accession changed upstream); use `--overwrite`.

---

## `virus_identity.tsv`

Tab-separated, one row per gene of the index: whether it is viral and which
virus it belongs to. The `analysis` step builds it from the index t2g, whose
fifth column is each transcript's genome accession, joined to the packaged
virus catalogue. On every cell-calling path, detection and multimapping stop when a kb matrix
column is not listed in the table (another run's table; re-run the analysis step to rewrite it)
and warn once when viral table genes are absent from the matrix (they cannot carry counts).

| Column | Description | Kind |
|--------|-------------|------|
| `gene_id` | Gene ID as in the index | observation |
| `genome_accession` | t2g column 5: the genome (GTF seqname) of the gene; empty for a pre-v3 t2g | observation |
| `status` | `catalogued` (accession in the catalogue), `uncatalogued` (in the viral GTF set, not catalogued; the run warns), `host`, or `legacy_prefix` (pre-v3 t2g, named by prefix maps) | observation |
| `viral` | `true` / `false` | observation |
| `virus_key` | One virus: `taxid:<n>`; `taxid:<n>\|<strain>` for one isolate of a segmented virus whose taxid several isolates share; `genus:<genus>` for anelloviruses; `accession:<acc>` for uncatalogued genomes; `name:<name>` for the legacy fallback | observation |
| `virus_name` | Display name of the virus: the curated name, else the NCBI organism | observation |
| `organism`, `taxid`, `species`, `genus`, `family` | NCBI organism and taxonomy of the genome | observation |
| `segment`, `strain` | Segment and strain/isolate of the genome, when recorded | observation |
| `sibling_group` | Near-identical viruses that allocation can move molecules between (e.g. `HHV-4` for EBV types 1 and 2) | observation |
| `risk_class` | `eve` for families with endogenous or ubiquitous-commensal risk | diagnostic flag |
| `role` | `decoy` for lab-contaminant decoy genomes | diagnostic flag |

A run with no viral gene in its index stops at this step.

## `viral_summary.tsv`

Tab-separated, one row per candidate virus: every virus with signal in the
all-barcode screen, whether or not it passes the called-cell detection gate. The
`detected` column carries the gate result; filter on it for the detected set.
Since ANDET-09 an anellovirus that only the alignment branch saw also gets a row
(`detection_source = alignment_only`, kallisto counts 0); see the anellovirus
alignment branch below.

| Column | Description | Kind |
|--------|-------------|------|
| `virus_name` | Human-readable virus name | observation |
| `role` | Curated reference role from `virus_identity.tsv`: `target`, `contaminant`, `endogenous` or `decoy`. Empty means uncurated/unknown, including legacy runs without role metadata. Shown in the HTML report as Reference role. This label does not filter counts or establish infection; contaminant, endogenous and decoy signals require interpretation in that context. Conflicting nonempty roles within one virus are rejected | diagnostic flag |
| `reference_risk_flags` | Pipe-separated reference-context cautions for endogenous sequences, vector/reagent contaminants and possible chromosomally integrated HHV-6. Flags do not filter counts or diagnose integration | diagnostic flag |
| `viral_molecules_total_est` | Unique viral molecules plus allocated ambiguous molecule mass, over all barcodes | model estimate |
| `infected_cells` | Legacy-named schema field: cells (all barcodes) whose summed allocated molecules across the virus's genes reach `--infected-min-molecules` (default 1.0); not confirmed infection. Cells with less than that are in `infected_subthreshold_cells` | model estimate |
| `total_cells` | Total cells in the count matrix (**all** barcodes) | observation |
| `pct_infected` | `infected_cells / total_cells × 100` (all-barcode denominator; secondary to `pct_infected_called`). **Blank (NA) in a two-step (`--host-filter starsolo`) run**: the viral matrix holds only barcodes with a non-host read, so there is no usable all-barcode denominator; see `denominator_note` | model estimate |
| `infected_subthreshold_cells` | Cells (all barcodes) with `0 <` summed molecules `<` `--infected-min-molecules`: fractional multimapper or EM weight that does not make a cell infected. Counted, never dropped | model estimate |
| `infected_subthreshold_called` | The same, restricted to called cells | model estimate |
| `background_umis` | Viral molecules in barcodes that are **not** called cells (empty droplets), i.e. `viral_molecules_total_est` minus the molecules in called cells. They never count toward the detection gate. Molecules outside called cells can come from ambient RNA, barcode swapping or empty droplets; this column does not tell which. With `--cell-calling none` every barcode is called and it is 0. The `ambient_*` and `soup_*` columns below size and test this signal | model estimate |
| `detected` | `true` when the virus passes the called-cell gate: at least `--detection-min-cells` (default 1) called cells are infected, **and** the called cells hold at least `--detection-threshold` of its molecules in total (compared with a relative tolerance of 1e-6, the same as the infected-cell rule, so a molecule split over several genes whose float shares sum to just under 1 still counts as 1). `false` marks a candidate that has signal in the all-barcode screen but misses the called-cell gate (for example, signal only in empty droplets). Passing is candidate evidence, not a confirmed infection, and does not exclude ambient RNA | evidence tier |
| `viral_molecules_per_10k_est` | Viral molecule estimate divided by full-matrix molecule estimate × 10,000. **Blank (NA) in a two-step run** (the full matrix is the virus-only matrix, so the ratio would read 10,000 whatever the load); see `denominator_note` | model estimate |
| `denominator_note` | Why `pct_infected` and `viral_molecules_per_10k_est` are blank, in a two-step run; empty otherwise (combined-index runs, values unchanged) | diagnostic flag |
| `n_called_cells` | Number of **called** cells (real, non-empty droplets) — see cell-calling below | observation |
| `infected_called` | Infected cells (same `--infected-min-molecules` rule) restricted to the called-cell set | model estimate |
| `pct_infected_called` | `infected_called / n_called_cells × 100`; a called-cell candidate-support rate, not a biological infection rate. **Within-run only** — see below | model estimate |
| `infected_comparable` | Candidate-support cells over the strategy-independent denominator | model estimate |
| `n_comparable_cells` | Barcodes clearing an absolute host-UMI floor (200 molecules), intersected with the called set | observation |
| `pct_infected_comparable` | `infected_comparable / n_comparable_cells × 100`. Use this to compare runs that used different host-filtering strategies. **Blank (unavailable), not 0, when `n_comparable_cells` is 0 or blank**: no cell cleared the floor, or host depth is unavailable, so the rate is undefined | model estimate |
| `sibling_crossmap_note` | Text note when the virus has near-identical siblings that allocation can move molecules between; empty otherwise | diagnostic flag |
| `accession_breadth` | Fraction of the virus's **index** genes with ≥1 molecule in any barcode (endogenous-element artefact check: an EVE concentrates on 1–2 loci, a genuine infection spreads across its genes). Before 2026-10-04 it was computed over the detected genes only and was always 1.0 (`ANDET-01`). Gene-level breadth, not genome-coverage breadth | diagnostic flag |
| `n_features_supported`, `n_index_genes` | Tier measure F and its ceiling (`calibration.freeze_scope`, Q135): the number of the virus's index genes whose summed `adata.X` over the **called cells** is ≥1 molecule (fractions kept; ≥1, not >0, so 0.5 + 0.5 over two cells counts and a 0.5/0.5 split over two genes does not), and the number of index genes. A virus with fewer index genes than the cutoff `f` cannot reach `probable` | observation |
| `host_viral_ambig_fraction` | Fraction of the virus's signal that is host-virus ambiguous | diagnostic flag |
| `eve_risk` | Flag copied from the identity table's `risk_class` (`eve` families: germline endogenous viral elements). No shipped family has it: Anelloviridae were cleared on 2026-10-03, because no germline human anellovirus EVE is known (F-019 revised the original basis) | diagnostic flag |
| `artifact_risk` | `low_complexity` when the virus's family is prone to low-complexity read artefacts (poly-G no-signal reads, poly-A sinks; F-019, F-021), otherwise empty. Today this is Anelloviridae only. It is a label, **not a filter**: anelloviruses are commensal and a real call is expected. Check reads with `viralscan evidence` | diagnostic flag |
| `claim_scope` | What the call may claim. `screening_only` for Anelloviridae (by catalogue family; genus-name fallback on legacy indexes): no orthogonally confirmed anellovirus-positive sample exists (`REF-10`), so a call is a screen, never a confirmed infection or a sensitivity claim. Empty means no scope restriction | diagnostic flag |
| `host_homology_status` | Genome-level host homology of the virus's reference genomes, read from `reference_manifest.json` beside the index (`viralscan build-ref --genome-dlist`): `measured` (every genome), `partial` (some) or `not_measured` (none, or no manifest). `not_measured` is not zero homology; indexes built without a host genome, or by another builder, always report it (the build now logs a warning when it records this status, AUD3-05) | observation |
| `host_homology_max_identity` | Highest minimap2 `asm10` identity of any alignment of the virus's genomes to the host genome (maximum over genomes). Empty when `not_measured`. No threshold is applied; where reads land is the evidence route's question (`ANDET-02`, `REF-08`) | observation |
| `host_homology_max_query_coverage` | Fraction of the genome covered by its longest host alignment (maximum over genomes). Empty when `not_measured` | observation |
| `host_homology_max_aligned_bases` | Length in bases of the longest host alignment (maximum over genomes). Empty when `not_measured` | observation |
| `ambient_viral_fraction` | Viral molecules over total molecules in the **ambient barcodes**: barcodes that are not called cells and hold at most `--emptydrops-lower` (default 100) molecules. Total molecules are host plus viral (a two-step run takes the host part from STARsolo `GeneFull/raw`). Blank, never 0, when the screen is `not_assessed` | model estimate |
| `ambient_viral_fraction_ub95` | One-sided 95 % exact Poisson upper bound on that fraction (`f_U`; with no viral molecule in the ambient barcodes it is `2.996 / total ambient molecules`). The bound, not the point estimate, enters the test, so a clean run does not look enriched by default | model estimate |
| `n_ambient_barcodes_with_virus` | Ambient barcodes with any viral molecule mass (`> 0` on the detection layer) | observation |
| `soup_expected_cells_max` | `E1`: the expected number of called cells that ambient RNA alone could push over the infected rule (`--infected-min-molecules`), taking the whole droplet to be soup (contamination fraction `rho = 1`, the largest soup any `rho <= 1` gives): the sum over called cells of `P(Poisson(f_U × N) ≥ ceil(rule))`, `N` the cell's total molecules. Called cells absent from the viral matrix (two-step runs) are included | model estimate |
| `soup_expected_cells_ge2_max` | `E2`: the same with at least two molecules (soup rarely puts two in one cell, infection often does) | model estimate |
| `soup_tail_p_cells` | Upper Poisson tail `P(X ≥ infected_called)` with mean `E1`. A small value means ambient RNA at any `rho <= 1` would rarely give this many infected called cells. It assumes independent soup counts, so over-dispersed soup or PCR duplicates make it optimistic | model estimate |
| `viral_depth_spearman` | Spearman correlation of viral molecules with total molecules over the called cells (zeros included). Soup scales with droplet size, but infected cells may also be larger. Blank below 3 infected called cells | diagnostic flag |
| `ambient_assessment` | `soup_rejected_rho1` (`soup_tail_p_cells` below 0.05), `soup_sufficient` (at or above 0.05, including when no called cell is infected) or `not_assessed`. **Descriptive: `detected` and every tier ignore it** (decision 2026-10-12, M-03). It is a weak necessary screen: `soup_rejected_rho1` does not show infection, replication, or that the cells are the source (doublets, uptake), and `soup_sufficient` means "not separable from soup", not "absent" or "soup". `not_assessed` for `--cell-calling none`/`knee`, no ambient barcode, no called cell, a virus-only index with no host depth, a two-step run without STARsolo `GeneFull/raw` totals, and alignment-only rows | diagnostic flag |
| `ambient_basis` | One sentence: the level and tail used, or the reason it is `not_assessed`. It always says that swapping between libraries is not assessed in a single run (see `scripts/swap_check.py` below) | diagnostic flag |
| `detection_source` | `kallisto`, `kallisto+alignment` (anellovirus with ≥1 genus-unique alignment molecule) or `alignment_only` (the alignment branch saw it, kallisto did not) | evidence tier |
| `alignment_status` | Anellovirus rows only: `ok`, `disabled` (`--no-anello-align`), `skipped_no_host_filter` or `skipped_no_anello_index`. Empty for other viruses | diagnostic flag |
| `alignment_reads` | Host-unmapped reads STARsolo aligned to the virus's genomes | observation |
| `alignment_unique_reads` | Of those, reads with a single placement (NH = 1) | observation |
| `alignment_molecules_unique` | Distinct (CB, UB) among reads whose every placement is this virus; STARsolo-corrected barcodes, homopolymer UMIs dropped | observation |
| `alignment_cells_unique` | Barcodes holding those molecules (all barcodes, not only called cells) | observation |
| `alignment_median_identity` | Median over accessions of the per-accession median `1 − NM / aligned bases` (edit-distance identity, indels included). **Read it beside `alignment_median_query_coverage`**: normalised by *aligned* length, so a soft-clipped 34-nt perfect match scores 1.0 exactly like a full-length one | observation |
| `alignment_median_query_coverage` | Median fraction of the read that took part in the alignment (aligned bases / read length). The soft-clipped remainder is where falsifying evidence hides — in the covid case the whole TSO sat in it | diagnostic flag |
| `alignment_accessions` | Genomes of this virus with ≥1 aligned read | observation |
| `alignment_start_sites` | Distinct read start positions, summed over accessions; one hotspot gives a low number | diagnostic flag |
| `alignment_homopolymer_fraction` | Read-weighted fraction of aligned reads carrying a ≥15-nt homopolymer (the F-019/F-021 artefact class). **Not diagnostic on its own**: a genuine 3′-end read carries a poly-A tail by definition | diagnostic flag |
| `alignment_complex_body_fraction` | Read-weighted fraction of aligned reads whose *body* — the sequence after any leading TSO/TruSeq reagent and 5′ of the first ≥15-nt homopolymer run — is ≥20 nt and either has dinucleotide entropy ≥2.0 bits or has ≥25 nt (the whole body, if shorter) inside the alignment. A genuine 3′-end read is `[complex viral sequence][untemplated poly-A]` and scores ≈1; a TSO/poly-A chimera has no body and scores ≈0. **Measures complexity, not viral origin**: G/C/A mosaic chimeric bodies (2.1–2.9 bits) and Alu or other non-viral bodies also score 1. Read it with `alignment_reagent_fraction` and `alignment_median_query_coverage` | diagnostic flag |
| `alignment_reagent_fraction` | Read-weighted fraction of aligned reads carrying reagent in a shape no genuine molecule has: forward TruSeq R1 (an adapter chimera), or the 10x TSO (either orientation, ≤2 mismatches, or a read starting inside it) **not** followed by a complex body — the `TSO|poly-T` zero-length insert. A TSO followed by a complex body is not counted: 10x documents TSO at the R2 start on full-length short molecules. Reverse-complement TruSeq after the poly-A is read-through on a short genuine insert and is not counted either | diagnostic flag |
| `alignment_r1_tso_fraction` | Read-weighted fraction of aligned reads whose R1 — raw barcode + UMI (`CR`/`UR`, else `CB`/`UB`) — contains TSO sequence (any 15-mer, ≤1 mismatch). The bead oligo defines those positions, so this is a chimera by construction; empty when no read carried an R1 tag | diagnostic flag |
| `alignment_splice_reads` | Aligned reads with an `N` CIGAR operation (spliced) | observation |
| `alignment_breadth` | Read-weighted mean over the virus's genomes of the fraction of genome bases covered by any aligned record (`breadth` in `anello_alignment_by_accession.tsv`, secondary placements included). **A label, never a gate**: 5′ capture concentrates true reads near the start site, so a low value does not argue against a call, and the F-005 ≤3.41 % line is not applied. Anellovirus alignment branch only (off by default); empty otherwise and for rows without breadth | diagnostic flag |
| `alignment_breadth_unique` | The same, counting only single-placement (NH = 1) records, so secondary placements cannot make related genomes look covered. Same caveats; not a gate | diagnostic flag |
| `alignment_breadth_max` | Breadth of the single best-covered genome of this virus (maximum of the per-accession `breadth`, circular-safe: positions are taken modulo the genome length). The read-weighted mean in `alignment_breadth` dilutes a well-covered genome among many sparsely covered relatives; this does not. Reporting only: not a gate, no call changes. Anellovirus alignment branch only; empty otherwise | diagnostic flag |
| `alignment_breadth_unique_max` | The same, counting only single-placement (NH = 1) records (maximum of `breadth_unique`). Reporting only | diagnostic flag |

**Reading the artefact columns (F-019, update 2026-10-04).** The pileup
position, the homopolymer fraction and the identity are each equally consistent
with a genuine 3′-end read, so none of them settles anything alone: 10x 3′
chemistry sees only transcript ends, and an artefactual poly-A read lands on the
genome's longest templated A-tract, which in an anellovirus sits at the polyA
site. Judge a call on `alignment_complex_body_fraction`,
`alignment_reagent_fraction` and `alignment_r1_tso_fraction`, with
`alignment_median_query_coverage` beside the identity. A call with
complex-body ≈ 1, both reagent columns ≈ 0 and coverage ≈ 1 is real evidence;
complex-body ≈ 0 with reagent > 0 is the chimera class.

**These columns are filled only when `--anello-align` is on**, which is not the
default (PLAN `ANDET-09e`). On a default run an anellovirus row carries
`alignment_status = disabled` and these columns are empty, and the only
artefact signal is `artifact_risk`. Empty is "not measured", never "clean".
Note also that the branch's own `--outFilterMatchNminOverLread 0.80` already
rejects most chimeras before they are counted, so a branch run reports few of
them — the columns describe what survived that filter, not the raw library.

All of these are **labels, never filters**. Anelloviruses are a commensal
virome, so a call is biologically expected, and every homopolymer-based filter
removes precisely the genuine 3′-end reads it would take to prove one — such
evidence *bounds* a real infection, it never excludes one (user, 2026-10-04).

**Two denominators.** The all-barcode `pct_infected` field is diluted by empty
droplets; `pct_infected_called` uses only the declared called-cell set and is the
appropriate denominator for a per-cell candidate-support rate. The field names
are retained for schema compatibility and do not establish infection. The v3
`auto` default uses an external called-cell list when one is supplied and
otherwise runs DropletUtils EmptyDrops. The approximate knee caller is explicit
only and is never a fallback. It is sensitivity-only and must not supply
reported numbers: on real libraries it calls empty droplets as cells. With `cell_calling=none`, the `*_called` columns
equal the all-barcode values.

**Two-step runs.** With `--host-filter starsolo` the kallisto matrix holds only barcodes
that carry a non-host read (host-filtered barcodes), so dividing by it gives 100 %
and 10,000 molecules per 10k for any virus present. `pct_infected` and
`viral_molecules_per_10k_est` are therefore written as NA (blank in the TSV, `NA` in
`summary.txt` and the HTML report), with the reason in `denominator_note`. Use
`pct_infected_called` (denominator: the host called cells) or `pct_infected_comparable`.
The same data gave `0.04` / `0.0156` with a combined index and `100.0` / `10000.0` in two-step
mode before this change, while both called and comparable rates agreed (`0.2`).

**Report headline.** `report.html` leads with the called-cell rate
(`pct_infected_called`: the run summary and the main table, with the infected-called and
called-cell counts), and its "high" badge (above 1 %) uses that rate; the comparable rate
stands in when there are no called cells. The all-barcode rate is shown beside it,
labelled secondary.

**Three denominators, and why the middle one is not comparable across runs.**
`pct_infected_called` divides by the called cells *of this run*, and cell calling
runs after host subtraction, so the denominator moves with the strategy. Host filtering changes the number of called cells and the number of viral molecules by different factors, so the reported prevalence can rise while viral molecules fall; that is why `pct_infected_called` inverts across strategies. Note that the published covid *Alphatorquevirus* signal was an artefact of poly-G no-signal reads (≈90 %) and host-homologous reads (≈10 %) (F-005, F-019); low-level divergent anellovirus is not excluded. `pct_infected_comparable` uses an absolute host-UMI floor that host filtering cannot move, so it is the field to use across strategies; `pct_infected_called` remains the within-run primary.

**`summary.txt` wording.** The per-virus lines lead with the called-cell figures
(`<n>/<n_called_cells> called cells with viral reads (<pct>%)`) and give the all-barcode
figure second, labelled as such; they also report the sub-threshold and background counts
when nonzero. The gene table is headed "Candidate viral Gene IDs" because it is the
all-barcode screen. A virus that fails the called-cell gate is marked `CANDIDATE ONLY`
(with "no infected called cell" when `infected_called` is 0), and a closing "Detection gate:"
paragraph states the rule in force (`--detection-min-cells`, `--infected-min-molecules`,
`--detection-threshold`) and lists the candidate-only viruses. "Official name of viral load
detected" lists only viruses with `detected` true. A parser of the old "cells infected (x%)"
phrasing must be updated.

**Count layer.** V3 summaries use `adata.X`, the complete selected-method
molecule matrix. `counts_unique` and `counts_ambiguous_allocated` are disjoint
and sum to `X`. Nonzero molecule support is candidate evidence; biological
interpretation requires calibrated evidence and may still require orthogonal
confirmation. The separate read-level workflow supplies diagnostics rather than
an automatic infection call.

---

## `sensitivity.tsv`

Tab-separated, one row per virus **in the index**, written on **every** run;
an undetected virus gets a row (`MECH-B`,
2026-10-04; before then only detected viruses had rows). Answers the question
a zero otherwise cannot: *is there nothing there, or did we not look hard enough?*

| Column | Description | Kind |
|--------|-------------|------|
| `virus_name` | Virus the row describes | observation |
| `observed_molecules` | Molecules attributed to the virus. For an undetected virus: the count in the detection matrix below the gate (`0` when there is none), so a nonzero value here with `detected = false` is a count below the gate, not an inferred absence | observation |
| `detection_threshold` | The sample-level gate that decided the call, applied to the virus's molecules summed over all its genes and, for the final call, over the **called cells** only (`--detection-min-cells` also applies; empty-droplet molecules are `background_umis` in `viral_summary.tsv`). Before 2026-10-04 it was applied per gene, which could miss a virus spread thinly across genes | observation |
| `detected` | The actual detection decision: the virus cleared the gate. Not inferred from the count, so a nonzero `observed_molecules` with `detected = false` is a count below the gate (a negative, not an absence) | observation |
| `capture` | Measured recovery of the positive control, only on rows that control certifies; **empty** otherwise (the depth-only floor is not printed as a measurement) | observation |
| `capture_measured` | `true` only if the capture term came from a valid positive control whose scope covers this virus: never a default, another virus's control, a failed control or an over-recovered one | diagnostic flag |
| `capture_status` | Why: `measured`, `not-measured` (no control configured), `control-out-of-scope` (a control was measured but its scope does not cover this row, including every `panel_mechanics` control), `control-failed`, `control-over-recovered`, `control-gene-missing` | diagnostic flag |
| `lod95_per_10k` | True viral molecules per 10k host molecules at which the virus would be reported with at least 95 % probability, **adjusted by the effective capture** (1.0 when capture is not measured, so the depth-only floor; a measured capture above 1 is clamped to 1) | model estimate |
| `lod95_molecules` | **Required true** molecules at that limit: `2.996 x detection_threshold / capture`. For `depth > 0` it equals `depth x lod95_per_10k / 1e4`. At capture 0.5 and threshold 1 it is 5.991, not 2.996 | model estimate |
| `lod95_expected_observed_molecules` | **Expected observed** molecules at that limit: `2.996 x detection_threshold`, independent of capture. `lod95_molecules x capture` equals this | model estimate |
| `lod_interpretation` | `informative` / `adequate` / `shallow` / `insufficient-depth`, from the capture-adjusted `lod95_per_10k` | diagnostic flag |
| `depth_only_sufficient` | The LOD95 at capture 1.0 (depth alone) is in the `informative` band (at most 0.01 per 10k) | diagnostic flag |
| `depth_sufficient` | The capture-adjusted `lod95_per_10k` is in the `informative` band. Legacy column, defined this way: it equals `depth_only_sufficient` for every row without a measured capture, and is stricter where capture costs depth | diagnostic flag |
| `sensitivity_eligible` | `depth_sufficient` **and** `capture_measured`: an adequate, measured sensitivity estimate exists for the row. Says nothing about the result | diagnostic flag |
| `informative_negative` | `sensitivity_eligible` **and not** `detected`. Always `false` for a called virus | diagnostic flag |
| `negative_blockers` | `;`-joined reasons `informative_negative` is false: `detected`, `capture-not-measured`, `depth-insufficient`. Empty when it is true | diagnostic flag |
| `expected_molecules_at_1_per_10k` | Expected observed molecules for a virus at 1 estimated molecule per 10k host molecules (cf. `viral_molecules_per_10k_est`) | model estimate |
| `p_detect_at_1_per_10k` | Probability of clearing the threshold at that abundance | model estimate |
| `p_zero_at_lod95` | ≈ 0.05 by construction; reported so the arithmetic is checkable | model estimate |
| `notes` | Why the LOD is a bound rather than an estimate, when a row is zero, etc. | diagnostic flag |

Two molecule counts, deliberately separate: **required true molecules**
(`lod95_molecules`) is how many viral molecules must be present, **expected observed
molecules** (`lod95_expected_observed_molecules`) is how many of them survive k-mer
capture and are counted. They coincide only at capture 1.

All four outputs (`sensitivity.tsv`, `positive_control.json`, `summary.txt` and
`report.html`) are rendered from the same per-target records, so for every target they
agree on whether capture was measured, whether it applies, whether the negative is
informative and why. The "LOD95" statements are probabilities under stated assumptions
(Poisson sampling of quantified molecules in the selected count layer, the stated
capture, a fixed gate), not a guarantee that a virus above the limit would have been
reported.

**`informative_negative` is the column that matters, and it is false almost
always.** Depth is not the limiting term in practice: the three real covid
configurations produced LOD95 values of 0.0003–0.0056 estimated viral molecules per 10k host molecules, all in
the `informative` band, and the covid samples called SARS-CoV-2 = 0 at 21.6 M
quantified molecules. What cannot be measured from inside a run is the k-mer
**capture** term. As a **substitution-only heuristic** (exact for i.i.d. substitutions
against one target; 90 bp, k=31; it ignores indels, sequencing error and competition
from the host and other panel members, so it is not a measurement) capture falls to
0.708 at 5 % sequence divergence, 0.254 at 10 %, 0.063 at 15 % and 0.013 at 20 %. Without
a measured capture term a negative cannot be certified at any depth — which is why
`informative_negative` requires an adequate depth, a measured capture **and** an
undetected row.

Depth is the sum of the count matrix, **not** raw reads: only quantified
molecules can be detected, and the two differ substantially — in this repo's own
benchmarks, pseudoaligned-read counts exceed quantified-molecule counts by
roughly 2x, so substituting reads for molecules would understate the LOD95 by
about the same factor. The Poisson floor was verified by molecule-level
downsampling of the bundled EBV LCL run: P(detect) stayed 1.0000 down to 1,270
downsampled reads and first reached 0 at 127.

## `positive_control.json`

Written on every run. Present so a negative can be audited.

| Field | Description | Kind |
|-------|-------------|------|
| `status` | `not-configured` / `measured` / `failed` / `over-recovered` / `gene-not-in-reference` | diagnostic flag |
| `certifies_negatives` | `true` only when `informative_negative_targets` is non-empty, i.e. at least one undetected row has a measured, in-scope capture **and** adequate adjusted depth. A run with a valid control on a shallow library, or on a virus that was called, is `false` | diagnostic flag |
| `scope`, `target` | What the control certifies: `exact_sequence`, `virus_key` (rejected until a transfer calibration exists) or `panel_mechanics` (certifies no virus); `target` names the virus. A control with no scope is treated as `panel_mechanics` and warns | diagnostic flag |
| `certified_targets` | The sensitivity rows whose capture was measured in the control's scope (`capture_measured` in `sensitivity.tsv`); never panel-wide. A called or shallow row can be listed here without being an informative negative | diagnostic flag |
| `informative_negative_targets` | The subset of `certified_targets` that is also an informative negative (`informative_negative` in `sensitivity.tsv`) | diagnostic flag |
| `capture_used_for_sensitivity` | The control's raw recovery ratio when at least one row used it, else `null` | observation |
| `scope_note` | Plain-language statement of what was and was not certified | diagnostic flag |
| `gene`, `expected_molecules`, `observed_molecules`, `capture` | The recovery ratio | observation |
| `implied_divergence` | Per-base divergence whose substitution-only heuristic capture matches, by bisection; `null` when unidentifiable | model estimate |
| `implied_divergence_note` | Labels `implied_divergence` as the substitution-only heuristic, never a measurement | diagnostic flag |

### Receipt-verified `exact_sequence` controls

An `exact_sequence` control certifies only with a checksum-bound receipt
(`--positive-control-receipt JSON`). The receipt has these required fields: `schema_version`
(`3.0.0`), `gene`, `accession_version`, `sequence_sha256`, `index_sha256`,
`reference_fasta_sha256`, `expected_molecules` (> 0), `count_layer`, `count_strategy`
(`kb_count` or the multimap method), `reference_manifest{path,sha256}` and
`evidence{path,sha256}`; relative paths resolve against the receipt's own directory.
`positive_control.json` then adds the fields below. Failed or unavailable verification forces
an empty `certified_targets` and `certifies_negatives = false`.

| Field | Description | Kind |
|-------|-------------|------|
| `verification_status` | `verified`, `failed` or `unavailable` | diagnostic flag |
| `verification_detail` | Which check failed or why verification was unavailable | diagnostic flag |
| `certification_status` | `single_accession_group`, or `group_contains_other_accessions` when the target's group has other accessions (no certification) | diagnostic flag |
| `resolved_target`, `accession_version`, `sequence_sha256` | The target the receipt resolved to through the virus identity table | observation |
| `receipt_sha256`, `receipt` | Digest and verbatim copy of the receipt | observation |
| `measurement_archive{path,sha256}` | The immutable copy `results/control_receipts/<sha256>.json`, written by exclusive create | observation |
| `archive_error` | Set when the archive could not be written; the run then fails | diagnostic flag |

Supply a control with `--positive-control-gene` and `--positive-control-molecules`
(both required together), and state what it covers with `--positive-control-scope`
(and `--positive-control-virus-key`). `failed` means the planted control was invisible, which
makes every negative in that run uninterpretable. `over-recovered` means more was
recovered than planted — the control is not spike-in-specific — and is treated as
no measurable loss, which is the optimistic direction.

Add `--require-positive-control` to fail the run outright when nothing is
detected and no capture term could be measured.

---

## `gene_program_summary.tsv` and `gene_program_cells.tsv`

Second layer, over the viruses the first layer detected. Present only when
`--gene-programs` is supplied.

### What it is for, and why a naive version is wrong

The obvious implementation is to sum the latent genes, sum the lytic genes, and
compare. On the bundled EBV LCL run (`SRR12682296`), which is latently
infected by construction, that gives:

| | UMI |
|---|---:|
| LATENT (10 genes) | 236,342 |
| LYTIC (12 genes) | 247,633 |

A latent:lytic aggregate ratio of **1.15** in a cell line defined by latency is
not biology. `EBNA-1` — expressed from the latent episome in latency I-III, so a
latency I-III cell line should show it well above a lytic gene — is 920 UMI, ~155× below `BHLF1`. The cause is pervasive
overlapping-ORF cross-mapping: EBV's latent transcripts are transcribed from a
region densely packed with nested and antisense lytic ORFs, so reads cross-map
in both directions. **Per-gene aggregate totals are uninformative**, and no
amount of care applied to them recovers an answer.

> **Correction (2026-09-27).** `BARF1` is latent only in epithelial cancers (NPC, EBV-gastric; PMID 32708965) and `BaRF1.1` is the lytic ribonucleotide reductase, so neither is a latency marker in a B-cell line. Both have been removed from the catalogue. The per-marker figures and the 2,240 / 1,277 cell counts in this section were measured with them included. Re-measured without them on the same run: **895** cells latent on the uniquely-placing layer versus **856** on the allocated layer, so most of the apparent latent-sensitivity gain came from `BARF1.2` (PLAN `PROG-07`).

> **Correction (2026-10-02).** `BNLF2a`, `BNLF2b` and `BHRF1` are early lytic genes (CAGE kinetic class *early*, PMID 29864140), not latency markers, and have been removed from the catalogue: `BNLF2a`/`BNLF2b` lie inside `LMP-1`'s overlap group, and `BHRF1`-locus reads can come from latent `EBNA-LP` transcripts (PMID 28950226). The `BNLF2a`/`BNLF2b` rows below were labelled `latent` when measured. HHV-6A `U90`/`U86`, HHV-6B `U95` and HHV-7 `U90` are immediate-early, so they are `productive` (PMID 12706083, 33627386, 10573164). HHV-6A and HHV-7 are now `partial` with latency not observable: no remaining latent marker separates latency from productive infection.

> **Re-measured (2026-10-03, PLAN `PROG-07`, `PROG-17`).** The figures above came from an older index and a run without barcode correction (`SW-13`). On the current panel (`CAT-42`), after `PROG-11`, the same LCL run gives these counts over the **932** called cells with EBV marker evidence (of the run's 2,763 emptyDrops-called cells). A first re-measure the same day scored 1,679 barcodes, empty droplets included, because layer 2 did not yet restrict itself to called cells (`PROG-17`); those numbers are superseded.
> - **Uniquely-placing layer:** 526 latent, 67 productive, 184 mixed, 155 indeterminate.
> - **Allocated layer:** 339 latent, 12 productive.
>
> Still **0** cells go from latent on the unique layer to productive on the allocated layer. Of the 526 latent cells, 319 stay latent on the allocated layer and 207 become `mixed` there. So in this run cross-mapping mostly adds lytic evidence; it does not drain latent evidence to `indeterminate`.

The per-marker breakdown on that same run shows the mechanism directly:

| marker | programme | uniquely-placing | multimap-allocated |
|---|---|---:|---:|
| `BARF1.2` | latent | **13,668** | 0 |
| `BNLF2a` | latent (since removed) | **0** | 49,662 |
| `BNLF2b` | latent (since removed) | **0** | 45,108 |
| `BZLF1` | productive | **0** | 9,308 |
| `BMRF1` | productive | **0** | 45,005 |
| `BcLF1` | productive | **6,609** | 59 |

The uniquely-placing layer puts **zero** molecules on `BZLF1` — the canonical
lytic marker — which is the correct answer for a latent cell line, and puts
13,668 on `BARF1.2` where the allocated layer put none.

Two defences, and it is worth being precise about what each buys:

1. **Breadth in overlap groups, never a gene count.** Overlap groups are
   computed in `extras/build_gene_programs.py` from shared 31-mers, not from
   exon-interval intersection: two genes are linked when at least 50 % of the
   smaller gene's canonical 31-mers (from the normalised exon sequences) occur in
   the other, and a group is a greedy maximal clique of that link graph, never a
   transitive chain. The groups are maximal, not a minimum partition, so breadth
   can still over-count linked genes that land in different cliques (the sensitivity
   table `analysis/audit_2026-10/aud01_sensitivity.tsv` counts such pairs). In EBV, `EBNA-1` is its own latent group,
   `EBNA-2` and `EBNA-LP` share one, and `BTRF1` and `BcLF1` are separate groups; a
   productive marker such as `BMRF1` is no longer pulled into a latent group by
   a chain of overlapping genes. Genes in one group cannot give independent
   evidence for each other. This is what removes the need for an aggregate
   comparison at all.
2. **Uniquely-placing evidence.** Per-cell breadth calling is directionally
   consistent on *both* layers: on the real run the two never disagree in the
   dangerous direction (**0** cells go latent-on-unique to
   productive-on-allocated). What the unique layer buys is **sensitivity** —
   **2,240** cells called latent versus **1,277** on the allocated layer,
   because 1,263 cells fall to `indeterminate` there once cross-mapping has
   drained their latent signal onto lytic ORFs. The failure mode is lost
   sensitivity, not an inverted call.

A third defence is about honesty rather than arithmetic:

3. **No absence claim from a thin anchor set.** A virus whose latency side is a
   single transcript (HSV-1's `LAT`) cannot support "not detected, therefore
   latent", so `latency_observable_in_rna` is false and `latent` is unreachable.

4. **No absence claim from a marker that may not be captured.** The shipped catalogue
   (`gene_programs.tsv`) carries `polya_capture` (`yes`, `no` or `unknown`) and `capture_basis` per
   marker: whether oligo-dT priming (10x 3' and Drop-seq chemistries) is expected to capture the
   transcript. `unknown` is the default and is never a guess; `yes` means only that the RefSeq
   annotation describes an mRNA or a protein-coding gene with a CDS, so it is annotation-implied and
   not measured per transcript. HSV-1 `LAT` is `unknown`: the stable LAT introns that dominate 3'
   scRNA-seq signal are not polyadenylated and are seen only through off-target priming
   (PMID 38445887), and the indexed record is the spliced transcript. A virus with such a marker
   gets a note in the `caveat` column of `gene_program_summary.tsv` ("no reads on them is not evidence
   that the transcript is absent"). If the marker is *productive*, the absence-based states
   (`reactivating`, `latent`) are not called and read `indeterminate`, the same rule as for an
   unresolved productive marker. A catalogue without the columns behaves as all `yes`.

### States

| State | Meaning |
|---|---|
| `productive` | ≥ `min_breadth` distinct non-overlapping productive overlap groups carry uniquely-placing molecules |
| `latent` | ≥ `latent_min_breadth` (default 1) latent overlap group(s), productive below `min_breadth`, **and** `latency_observable_in_rna` |
| `reactivating` | productive evidence in which every productive marker carrying molecules is immediate-early (`kinetic_class`), and the cell is not `productive`, `mixed` or `latent`. It only replaces what would otherwise be `indeterminate`. It is reachable only for a virus whose catalogue rows set `latency_observable_in_rna` true **and** that has an immediate-early productive marker: today that is EBV alone. VZV, HHV-6A/6B, HHV-7, HSV-1/2 and HCMV never emit it (`not_applicable`: there is no latency anchor, and for HCMV latency mirrors late-lytic expression at low levels while the IE genes UL122/UL123 are detectable in latent cells); KSHV has no immediate-early productive marker. It is also not called when a productive catalogue marker is unresolved against the index (see `indeterminate`) |
| `mixed` | both of the above in the same cell |
| `indeterminate` | the virus was detected but no programme met its threshold, and the cell is not `reactivating` — not a negative. It is also what a cell reads when a productive catalogue marker is unresolved against this index: both absence-based states (`latent`, which asserts productive evidence is below threshold, and `reactivating`) degrade to `indeterminate`, because the markers that would show the absence are missing. Positive calls (`productive`, `mixed`) are unaffected. `n_markers_unresolved` and `caveat` record it |
| `not_applicable` | no programme model exists for this virus (summary rows only) |

### `gene_program_cells.tsv`

One row per called cell and detected virus that carries any marker evidence on either
layer: `barcode`, `virus_name`, `state`, `productive_breadth`, `latent_breadth`
(uniquely-placing layer), `selected_state`, `selected_productive_breadth`,
`selected_latent_breadth` (multimap-allocated layer), `evidence_layer` and
`latency_not_observable`. A cell is omitted only when all four breadths are 0. A cell
with allocated-only evidence (unique breadth 0, selected breadth above 0) **is** listed,
so `n_cells_*_selected_layer` in the summary counts it. Before 2026-10-11 such cells were
dropped and the selected-layer counts were understated. The never-filled `n_cells`
column was removed.

### Coverage and honesty fields

| Column | Description | Kind |
|---|---|------|
| `panel_completeness` | `complete` (EBV and KSHV — a real latency *and* reactivation split; KSHV since 2026-10-03, `PROG-08`), `partial` (CMV, HSV-1/2, HHV-6A/6B, HHV-7, VZV), `not_applicable`. CMV is partial because single-cell latency mirrors a low-level late-lytic programme (PMID 29535194), so no marker's presence separates the states | diagnostic flag |
| `latency_observable_in_rna` | Whether a `latent` call is reachable. `false` ⇒ `n_cells_latent` and `n_cells_mixed` are 0 **by construction** | diagnostic flag |
| `evidence_layer` | Always `counts_unique_viral` | observation |
| `min_breadth` | The `--programme-min-breadth` used | observation |
| `latent_min_breadth` | The `--programme-latent-min-breadth` used (default 1) | observation |
| `marker_min_umi` | The `--programme-min-umi` used: a marker counts only at ≥ this many molecules (default 0 = any nonzero value, so the allocated layer's fractional values still count) | observation |
| `n_cells_reactivating` | Cells called `reactivating` on the unique layer; the string `not_applicable` where the state is unreachable (every virus but EBV, and any virus with an unresolved productive marker), so a number here is never a structural zero | model estimate |
| `productive_breadth_median` / `latent_breadth_median` | Breadth on the unique layer | model estimate |
| `n_cells_latent_selected_layer` / `n_cells_productive_selected_layer` | What the same rule would have called on the multimap-allocated layer — the honest comparison | model estimate |
| `selected_*_breadth_median` | Breadth on the allocated layer, for reference | model estimate |
| `layer1_molecules` | Layer 1's molecule total, so the two layers need not be joined by hand | observation |
| `n_called_cells` | Cells scored: the called set from `called_cells.tsv`. Before 2026-10-03 (`PROG-17`) layer 2 scored every barcode in the multimap H5AD, so older `n_cells_*` counts include empty droplets | observation |
| `n_markers_resolved` | Catalogue markers found in this index's gene set and used for the calls (`PROG-14`); empty in summaries written before it | observation |
| `n_markers_unresolved` | Catalogue markers absent from this index; non-zero means the calls use a reduced marker set and `panel_completeness` describes the catalogue, not this index | observation |
| `caveat` | Why a row is weaker than it looks | diagnostic flag |

### Scope

Only viruses **detected by layer 1** appear: those with `detected` true in
`viral_summary.tsv`. A virus layer 1 did not call gets no programme row, and a
candidate that only passed the all-barcode screen (`detected` false, see
`viral_summary.tsv`) is skipped too (a summary written by an older run without the
`detected` column is read as all detected). Not calling a virus is a
detection-limit problem (`sensitivity.tsv`), not something layer 2 can repair. This is a transcriptomic assay throughout:
DNA-level latency (a silent HIV provirus, a transcriptionally silent integrated
HPV genome, the HBV cccDNA pool) produces no reads and is **invisible, not
latent**.

---

## `per_cell_viral.tsv`

Tab-separated, one row per infected cell × virus combination, for every candidate
virus of `viral_summary.tsv` (`detected` true or false). A cell is listed for a virus only
when its summed allocated molecules reach `--infected-min-molecules` (the rule behind
`infected_cells`); sub-threshold cells are counted in `viral_summary.tsv`
(`infected_subthreshold_cells`) and are not listed here. Rows with `is_called_cell` false
are barcodes outside the called-cell set.

| Column | Description | Kind |
|--------|-------------|------|
| `barcode` | Cell barcode | observation |
| `virus_name` | Virus name | observation |
| `viral_molecules_total_est` | Viral molecule estimate for this cell from the selected-method matrix | model estimate |
| `molecules_total_est` | Full selected-method molecule estimate for this cell | model estimate |
| `viral_fraction` | Viral molecule estimate divided by full molecule estimate | model estimate |
| `is_called_cell` | Boolean flag: whether the barcode is in the primary called-cell denominator (see cell-calling above) | diagnostic flag |

---

## Ambient RNA and barcode swapping (`M-03`)

Two mechanisms give low-level viral signal that is not infection of the cell carrying it, and they
need different controls. **Ambient (cell-free) RNA** appears in every droplet, so the single-run
screen compares called cells with the empty droplets of the same run (the `ambient_*` and `soup_*`
columns of `viral_summary.tsv`). **Barcode or index swapping** keeps the donor library's barcode,
so a swapped molecule mostly sits in a barcode the recipient did not call; it can only be seen by
comparing libraries that shared a sequencing pool and a barcode whitelist, which ViralScan cannot
infer. Both are descriptive: `detected`, `pct_infected_*` and every denominator are unchanged, and
empty-droplet evidence stays separate from the called-cell denominator. The `contamination` row of
`interpretation_flags.tsv` (`viralscan evidence`) shows the screen's `soup_tail_p_cells` in `metric`
and the assessment in `basis`, and its `status` stays `not_assessed`: one run never measures swapping.

### `scripts/swap_check.py` (cross-run, tier 1)

```bash
python scripts/swap_check.py pool_manifest.tsv --out swap_check.tsv
```

`pool_manifest.tsv` is written by the user, one row per **library**, tab-separated, header required.
`library_id` and `run_dir` identify a library and its finished ViralScan run directory (relative paths
resolve against the manifest's folder). `pool_id` and `whitelist` are free labels: two libraries are
compared only when both match and `library_id` differs. List the lanes or technical replicates of one
library once, because their shared barcodes are legitimate and would read as 100 % swapping. Drop-seq
and other chemistries with random per-library barcodes have no shared whitelist and do not apply.

```text
library_id	run_dir	pool_id	whitelist
donor_lib	runs/donor_lib/out	flowcell_A_lane1	10x-3M-february-2018
recipient_lib	runs/recipient_lib/out	flowcell_A_lane1	10x-3M-february-2018
neg_lib	runs/neg_lib/out	flowcell_A_lane1	10x-3M-february-2018
```

One row per ordered (recipient, donor, virus) with at least one non-called viral barcode in the
recipient. With `S` the recipient's non-called barcodes carrying the virus, `D` the donor's barcodes
carrying it (called or not) and `U` every barcode the recipient observed that is not a called cell
(`kb-python/counts_unfiltered/cells_x_genes.barcodes.txt` minus `results/called_cells.tsv`; `-1`
suffixes stripped on both sides):

| Column | Description | Kind |
|--------|-------------|------|
| `pool_id`, `whitelist` | The manifest labels the pair shares | observation |
| `recipient_library`, `donor_library` | The ordered pair | observation |
| `virus_name` | Virus | observation |
| `n_noncalled_virus_barcodes` | `\|S\|` | observation |
| `n_donor_virus_barcodes` | `\|D\|` | observation |
| `n_noncalled_observed` | `\|U\|` | observation |
| `n_noncalled_in_donor` | `\|U ∩ D\|`, the barcodes that could overlap | observation |
| `overlap_with_donor` | `\|S ∩ D\|` | observation |
| `overlap_expected` | `\|S\| × \|U ∩ D\| / \|U\|`, the overlap expected by chance | model estimate |
| `overlap_p` | Upper hypergeometric tail `P(X ≥ overlap_with_donor)` | model estimate |
| `swap_assessment` | `overlap_above_chance` (`overlap_p` below 0.05), `not_above_chance`, or `not_assessed` (no donor barcode is observed unclaimed in the recipient, so there is no power). Descriptive; several donors and viruses are several tests | diagnostic flag |

Overlap above chance is a signature of a shared origin (index hopping, a chimeric read, carry-over, a
technical replicate declared as two libraries), not a mechanism. No overlap with the donors listed
says nothing about libraries not listed. `per_cell_viral.tsv` lists only barcodes at or above
`--infected-min-molecules`, so weaker swapping is invisible. A shared `(CB, UMI)` key test over the BUS
files (tier 2) is not implemented. Without a pool manifest, swapping stays `not_assessed`.

### What these controls do not support

- **Doublets and uptake.** A droplet holding an infected and an uninfected cell, or a cell that took up
  infected material, carries real viral RNA that RNA counts cannot separate from infection. ViralScan
  has no doublet caller and none of the columns above addresses this; a called-cell hit is
  candidate evidence, not infection or replication in that cell. Replication evidence (spliced or lytic
  transcripts, strand, programme state) is a different axis.
- **Calibration.** The screen uses only the calibration-free `rho = 1` bound and conventional levels
  (0.95 bound, 0.05 tail). `rho`, the ambient-range ceiling and any hop-rate cutoff are not estimated
  or tuned, so nothing here can gate a call before a protocol amendment and calibration on frozen
  training data.
- **Soup that is not homogeneous.** Lysis of infected cells enriches nearby droplets, so the empty
  droplets need not match the soup inside a cell, and in a genuinely infected sample the soup is itself
  signal (`soup_sufficient` is then "not separable").

---

## `report.html`

A self-contained HTML file with:

- Run metadata (date, parameters, sample names)
- QC summary table
- Per-virus detection table with normalised metrics
- Multimapping evidence table with unique and ambiguous viral support
- Embedded histogram plots (base64 PNG)
- Interpretation guidance

UMAP clustering annotations are written to the interactive UMAP HTML files
when `--umap` is supplied. The UMAP is viral-inclusive (viral genes are in the embedding); the
annotated p-value comes from a nearest-neighbour permutation test on a separate viral-free host
graph and is descriptive only: it does not distinguish infection from ambient RNA, uptake or
doublets. Virus overlays use the downstream detected-row filter.
Rows with `detected = false` are excluded. A legacy summary without `detected`
warns and retains all rows. Malformed flags in a newer summary stop processing.
Blank flags are retained only for `detection_source = alignment_only` rows.

Open in any modern browser — no internet connection required.

---

## `cell_type_enrichment.tsv`

Tab-separated, one row per detected virus and labeled cell type. Written only
when `--cell-types cell_types.csv` is supplied.

| Column | Description | Kind |
|--------|-------------|------|
| `virus` | Virus name | observation |
| `cell_type` | Cell-type label from the CSV | observation |
| `n_infected` | Legacy-named field: labeled cells with candidate molecule support in this cell type, using the selected-method matrix | model estimate |
| `n_total` | Total labeled cells of this type | observation |
| `pct` | `n_infected / n_total × 100` | model estimate |
| `OR` | One-sided Fisher exact odds ratio | model estimate |
| `pvalue` | Raw Fisher exact p-value | model estimate |
| `padj` | Benjamini-Hochberg adjusted p-value | model estimate |

Input CSV requirements: identical duplicate barcode labels are accepted.
Conflicting labels for one barcode stop processing.

```text
barcode,cell_type
AAACCCAAGAGT-1,T cell
AAACCCAGTGCA-1,Monocyte
```

Barcodes must match `adata.obs_names`. If no barcodes overlap, ViralScan skips
the enrichment table and logs a warning.

---

## `count_audit.tsv` and `molecule_audit`

One row per sample (also stored as `adata.uns["molecule_audit"]`), validated against `schemas/v3/count_audit.schema.json`. The molecule fields (`input_molecules`, `resolved_molecules`, `unique_molecules`, `ambiguous_molecules`, `unresolved_molecules`, `ignored_read_multiplicity`, `allocated_ambiguous_mass`) count only what survives barcode correction. The off-list drop fields (PLAN `MECH-F`) record what correction removed:

| Column | Description | Kind |
|--------|-------------|------|
| `barcode_correction` | Who corrected: `viralscan` (user on-list, `bustools correct` re-run), `kb` (kb's own corrected BUS), or `none` | observation |
| `bus_records_raw`, `bus_reads_raw` | Records / reads in kb's raw `output.bus` | observation |
| `bus_records_after_correction`, `bus_reads_after_correction` | Records / reads in the corrected BUS that molecule resolution used | observation |
| `offlist_dropped_records`, `offlist_dropped_reads` | Raw minus corrected: barcodes off the on-list with no unique one-mismatch neighbour. Reads is the figure comparable to kb's "reads on the whitelist" in `inspect.json` | observation |
| `bus_totals_source` | `bustools_inspect`, or `not_corrected` (then drops are 0 and the raw/corrected totals are absent) | observation |

The new fields are optional in the schema, so audits written before this change still validate.

## `multimap_evidence.tsv`

Tab-separated, one row per viral gene in the reference. This table is additive:
it does not replace `viral_summary.tsv` or change its default schema.

| Column | Description | Kind |
|--------|-------------|------|
| `virus_name` | Human-readable virus name or accession fallback | observation |
| `gene_id` | Viral gene/accession ID | observation |
| `viral_molecules_unique` | Integer molecules resolving only to this viral gene | observation |
| `viral_molecules_ambiguous_allocated` | Fractional ambiguous molecule mass assigned here | model estimate |
| `host_virus_ambiguous_molecules` | Molecules compatible with host and viral genes | observation |
| `viral_molecules_total_est` | Unique molecules plus allocated ambiguous mass | model estimate |
| `viral_molecules_upper_bound` | Unique signal plus all viral-compatible ambiguous mass | model estimate |
| `n_unique_viral_cells` | Cells with unique viral signal | observation |
| `n_ambiguous_viral_cells` | Cells with viral-compatible ambiguous signal | observation |
| `multimap_method` | `equal`, `host-conservative`, `unique-weighted`, `em-global`, or `em-cell` | observation |
| `evidence_tier` | `candidate_unique`, `candidate_virus_ambiguous`, `candidate_host_virus_ambiguous`, or `not_detected` | evidence tier |
| `locus_id` | Independent locus of this gene (LOCK-01 D1): genes whose annotated intervals share at least one base on either strand form one locus; empty when `results/viral_loci.tsv` is absent | observation |
| `n_unique_positive_cells_called` | Per virus (repeated on each of its rows): called barcodes with at least one unique-layer molecule at any gene of the virus, the union over genes and not the maximum (LOCK-01 D2) | observation |

The default `multimap_method` is `host-conservative` (keeps host-virus ambiguous
mass off viral genes); use `equal` for an equal-allocation comparison. These are
molecule-evidence tiers only. `probable` and `strong` require calibrated
read-level evidence and are never assigned from molecule counts alone.

(evidence-tiers)=
### Probable and strong tiers (`viralscan.evidence_tiers`, PLAN `TIER-01`)

The `evidence_tier` column above is never `probable` or `strong`. Those two tiers come from
`viralscan.evidence_tiers`, a pure module (no file I/O, not yet wired into `detection` or the Snakefile)
built to the frozen `calibration.threshold_search` contract in `analysis/v3_validation/protocol.yaml`:

- the cutoff triple is (distinct supporting molecules per cell, distinct supporting features per sample,
  exact-fragment support fraction), searched over the protocol's 45-point grid, read by
  `grid_from_protocol`; no cutoff off the grid is accepted;
- a grid point is eligible only with zero probable-or-strong calls in synthetic host-only and in planted
  host-homology training negatives; among eligible points within one standard error of the best training
  F1, the most conservative wins (largest molecule cutoff, then feature cutoff, then support fraction);
- no eligible point means calibration failed: the tier is `no_tier`, never a relaxed point;
- without a supplied calibration result the tier is `uncalibrated`. Comparisons are inclusive (`>=`);
  an unmeasured support fraction never reaches a tier. The molecule cut-off carries the same 1e-6 relative
  slack as `infected_masks` (`INFECTED_RELATIVE_TOLERANCE`), so 1.999999 molecules meets a cut-off of 2;
  `assign_tier` has no production caller yet, so no shipped output is affected;
- the grid defines only *probable-or-stronger*. `strong` exists only as a caller-supplied grid point at least as
  strict as `probable` on every cutoff; without one, `strong` is never emitted.

Package tiers and the study claim levels (`scripts/adjudicate_candidates.py`, `level`) are separate scales. The
adjudicator's last column `package_tier_analogue` is the one explicit mapping between them. It names the analogue
by position, not an equivalence and not the tier the package assigned to the row: a package tier is a calibrated
molecule-count rule, and a study level adds read-level evidence on top of its package analogue.

| Study level | Package tier analogue |
|-------------|-----------------------|
| `L0` | `candidate` |
| `L1` | `probable` |
| `L2` | `strong` |
| `L3` | `strong` (the gene-programme gate has no package counterpart) |

(molecule-assignments)=
## `molecule_assignments.tsv.gz`

Optional (`--multimap-molecule-assignments`, off by default; PLAN `SW-03`). One gzip TSV
row per corrected (CB, UMI) molecule of an on-list barcode, written in a second pass over
the same BUS text that the count layers come from. It changes no count. The weights are
those of the *selected* `--multimap-method`: summed per cell and gene they equal
`counts_ambiguous_allocated`, and the `unique` rows equal `counts_unique`. Records of
barcodes absent from the barcode list are not molecules of the run and are not written.
`rerun-multimap` to another method runs a full multimap pass instead of a layer swap
when this file exists, so it never describes a different method than its tree.

| Column | Description | Kind |
|--------|-------------|------|
| `barcode` | Corrected cell barcode (trailing `-1` removed) | observation |
| `umi` | UMI of the molecule | observation |
| `ec_ids` | Equivalence-class IDs of its BUS records, `;`-joined | observation |
| `status` | `unique` (one compatible gene), `ambiguous` (several), or `unresolved` (the records' gene sets are disjoint) | observation |
| `genes` | Compatible gene IDs, `;`-joined; empty when `unresolved` | observation |
| `weights` | Fraction of the molecule given to each gene, aligned with `genes` and summing to 1; `1` for `unique`; empty when `unresolved` | model estimate |
| `host_virus_boundary` | `true` when an ambiguous molecule is compatible with both host and viral genes; blank when `unresolved` | observation |

---

## `host_homology_loci.tsv` (PLAN `ANDET-02`)

One row per host-homology locus that `viralscan build-ref --genome-dlist` retained for a viral
genome of the index (read from `reference_manifest.json` beside the index). The file is always
written, with the header only when the index carries no measurement. It records **where** the
homology lies; no call, `detected` value or `eve_risk` changes because of it (demotion waits for
the `REF-08` calibration).

| Column | Description | Kind |
|--------|-------------|------|
| `virus_name` | Display name of the virus the genome belongs to | observation |
| `accession_version` | Viral genome accession with version | observation |
| `viral_start`, `viral_end` | Aligned interval on the viral genome (0-based, half-open) | observation |
| `host_contig`, `host_start`, `host_end` | Aligned host interval (0-based, half-open) | observation |
| `strand` | PAF strand: `+` when the host contig aligns to the viral genome in the same orientation, `-` for the reverse complement | observation |
| `identity` | Matching bases divided by alignment block length (0 to 1) | observation |
| `aligned_bases` | Alignment block length of the locus, in bases | observation |
| `mapq` | minimap2 mapping quality | observation |
| `eve_risk` | The virus's `eve_risk` flag from `viral_summary.tsv` | diagnostic flag |

---

## Anellovirus alignment branch (`anello_align/`, PLAN `ANDET-09`)

kallisto needs an exact 31-mer, and the panel's anellovirus set was clustered
at 95 % identity, so a divergent strain can share almost no k-mer with its
nearest genome (F-013). `viralscan evidence` only re-aligns reads that kallisto
already placed. This branch aligns **every** host-unmapped read to the panel's
own anellovirus genomes with STARsolo and reports what it finds. That covers
reads kallisto never placed.

It runs when all three hold: `--anello-align` (the default), `--host-filter
starsolo`, and an `anello_star/` index next to the kb index (built by
`scripts/build_bundled_panel_ref.py`, or `--anello-star-only` for an existing
build). Otherwise `alignment_status` in `viral_summary.tsv` says which was
missing.

STAR settings (`viralscan.anello_align.ALIGN_ARGS`):
- ≥80 % of the read aligned, ≤8 % mismatches;
- up to 100 placements, all written to the BAM (`NH` is the true count);
- two-pass splice discovery;
- unstranded counting (F-020);
- barcodes corrected against the same on-list kb uses.

These are starting values, calibrated by the `ANDET-09e` plant.

- `anello_align/Aligned.sortedByCoord.out.bam` (+ `.bai`), `Solo.out/`,
  `SJ.out.tab`: for read review.
- `results/anello_alignment_by_accession.tsv`: one row per genome with ≥1
  aligned read. A `# {json}` first line records the STAR version, the index
  and the barcode list. Columns:
  - `reads`, `unique_reads`;
  - `weighted_reads` (Σ 1/NH; describes ambiguity, it is not EM);
  - `median_nh`, `median_identity`;
  - `breadth` (all placements) and `breadth_unique` (NH = 1 only), so 100
    secondary placements cannot make every related genome look covered;
  - `start_sites`, `homopolymer_fraction`, `splice_reads`;
  - `sense_fraction` (reads on the record's forward strand).
- `results/anello_alignment_by_virus.tsv`: the virus-level columns merged into
  `viral_summary.tsv`.

All of it is a **label, never a filter**. Anelloviruses are commensal, so a real
call is expected. No row is a confirmed infection (`REF-10`).

---

<!-- viralscan-claim:v3-evidence-workflow status=provenance_incomplete -->
## `viralscan evidence` output

Evidence is generated in the explicit `--output` directory for one exact
accession, registered alias, or canonical call. With `--viral-fasta` it also
requires a full `--host-fasta`; alignment and BLAST are competitive rather
than virus-only.

Reusing a non-empty evidence directory requires a matching canonical target,
gene set, allocation method, and non-empty run fingerprint in its receipt.
Use a fresh directory for older outputs without that provenance.
Batch extraction disambiguates colliding target names with deterministic hash suffixes.
Reads compatible with several genes of the same virus remain `unique_target`
when no other virus or host gene competes.

The implementation and test sources cover this workflow, but no immutable
successful execution receipt is registered; its execution claim therefore
remains provenance-incomplete.

| Output | Description | Kind |
|--------|-------------|------|
| `read_lineage.tsv.gz` | Exact read ID, CB, UB, EC, compatible genes, ambiguity class, assigned weight, method, tier, and exclusion reason | observation |
| `evidence_manifest.json` | Schema version, exact target, allocation method, run fingerprint, reference hashes, the `replay` block (`strand` from `config.yaml`, null when `kb count` was given no flag; `kallisto_strand_flag`; `technology`; `kallisto_technology`, the `-x` string kallisto received, where `10xv4` becomes `10XV3`), and output hashes. Manifests written before this block lack `replay` | observation |
| `competitive_reads.raw.bam` | Sorted/indexed raw competitive host-plus-target alignments | observation |
| `competitive_reads.<mode>_dedup.bam` | Separate UMI, coordinate-marked, or non-deduplicated evidence BAM | observation |
| `competitive_reads.deduplicated.tagged.bam` | Optional indexed BAM with CB/UB tags for per-cell IGV grouping | observation |
| `coverage.raw.tsv`, `coverage.deduplicated.tsv` | Raw and deduplicated coverage summaries | observation |
| `coverage.raw_vs_deduplicated.png` | Direct depth-track comparison | observation |
| `alignment_qc.tsv` | Per-reference breadth at 1x/3x/10x, depth, intervals, hotspots (`read_start_entropy`, `read_start_gini`, `max_50bp_window_fraction`, the maximum fraction of starts in any 50 consecutive reference bases; `read_start_gini` is the Gini coefficient of read-start counts over **all** reference positions, positions without a start counting as zero, so a single-position pile-up scores near 1 rather than 0), strands, mapping quality, identity, host competition, cells, molecules, duplicate fraction, and the read-side artefact labels `complex_body_fraction` (primary reads with a complex non-poly-A body) and `reagent_fraction` (primary reads carrying TSO or TruSeq adapter sequence); both are measured in sequencing orientation and are labels, never filters | observation |
| `per_cell_alignment_qc.tsv` | Per-cell host/virus competitive reads, molecules, strands, mapping quality, and identity | observation |
| `blast_identity.tsv` | Best host and viral hit with identity, query coverage, E-value, bit score, score difference, and sequence-complexity flag. A blank `viral_minus_host_bitscore` means the read was not assessed (no hit), and `interpretation_flags.tsv` does not count it as host-preferred | observation |
| `blast_sampling.json` | Deterministic sampling strategy, seed, counts, and fraction | observation |
| `interpretation_flags.tsv` | Transparent host-homology, low-complexity, ambiguity, and hotspot diagnostics; all are diagnostic only | diagnostic flag |
| `viralscan_evidence.igv.xml` | IGV session containing raw, deduplicated, and optional CB/UB-tagged BAMs | observation |
| `molecule_verdicts.tsv.gz` | One row per corrected (CB, UMI) molecule: `reads`, `counted` (at least one read with `assigned_weight > 0`), best viral and host alignment score (`virus_as`, `host_as`; blank when that side scored below the minimap2 secondary cutoff), and `verdict`. Written only without `--competitor-fasta` | observation |
| `molecule_support.tsv`, `read_support_summary.tsv` | Q137 provisional read-support **flags** per molecule (one row per corrected (CB, UMI); the representative read is the one with the best viral alignment) and one call-level summary row. Gates: viral aligned length >= 60 (M/=/X) and >= 60% of the trimmed query (`aligned_bases`, `query_coverage`); >= 20 DUST-unmasked aligned bases (`dust_unmasked_bases`, symmetric DUST w64 l20), with `aligned_AT_warning` at max(A,T) >= 0.70; competitive human alignment (`human_equal_or_better` when human AS >= viral AS); `hotspot_only` and the DNAJB6 `host_decoy` from `src/viralscan/data/read_support_loci.tsv`. `known_artefact` (`human_DNAJB6`, `low_complexity_hotspot`) is set per molecule and, in the summary, only when every molecule fails for that reason (`known_artefact_mixed` when every molecule fails for one of several). `n_molecules`, `n_pass_all_gates`. Flags remove nothing and move no threshold. Written only without `--competitor-fasta` | observation |
| `viral_reads.umi_dedup.fasta` | One representative extracted read per corrected (CB, UMI) molecule, chosen after the competitive alignment; the BLAST sample and the interpretation-flag fractions are drawn from it | observation |
| `molecule_verdict_summary.tsv` | Molecules per `scope` (`all`, `counted`) and `verdict`; every verdict is listed, zeros included. `scripts/dsr02_verdicts.py` decides on the `counted` scope | observation |

Read-start hotspot flags are uncalibrated diagnostics. With two or fewer reads,
or unavailable read counts, hotspot assessment is `not_assessed`.

`verdict` is `virus_best` when the best viral `AS` exceeds the best host `AS`,
`host_best` for the reverse, `tie` when equal, and `unaligned` when neither side
aligned. When only one side has an alignment the other scored below the
secondary cutoff, so the present side wins. A verdict compares scores of the
emitted alignments; it is not a probability and not a host-contamination call.

### `results/evidence/<genus>/` and `log/auto_evidence.done`

Written only with `--auto-evidence` (off by default; PLAN `ANDET-04`), which needs
`--host-fasta` and a viral-only FASTA (`--viral-fasta`, default `viral.fa` next to the index).
After detection, `viralscan evidence` runs with host confirmation once for each
detected Anelloviridae genus in `results/viral_summary.tsv`, skipping `detection_source =
alignment_only` rows, whose reads the kallisto replay cannot see. Each `<genus>/` holds the
normal `viralscan evidence` outputs, including `evidence_manifest.json`; the evidence is
diagnostic and does not by itself confirm infection. `log/auto_evidence.done` lists the
genera processed, or `none`.

### Competitor mode (`--competitor-fasta`)

With `--competitor-fasta` and `--competitor-manifest`, the competitive FASTA
holds target, host, and alternate/related/decoy references, so a two-sided
host-versus-virus verdict is undefined. `molecule_verdicts.tsv.gz` and
`molecule_verdict_summary.tsv` are **not written**, and the manifest records
`competition.molecule_verdict = "not_produced_in_competitor_mode"`. The run
needs a fresh, empty output directory. Everything is diagnostic only.

| Output | Description | Kind |
|--------|-------------|------|
| `molecule_competition.tsv` | One selected (CB, UMI, reference) representative per molecule: `selected_reference_id`, `selected_class`, `n_candidate_alignments`, `candidate_classes` (semicolon-joined), `mapping_quality`, `identity`, `alignment_span`, `ambiguity_class` (`emitted_multiclass`, `candidate_search_not_exhaustive`, or `unresolved`), `status` (`selected`, `unresolved`, or `failed`), `candidate_search`. Molecules with no alignment are kept as `unresolved` with blank selection fields | observation |
| `competitor_summary.tsv` | One row: molecule and read counts (`n_molecules`, `n_reads_total`, `n_reads_sampled_blast`, `n_unsampled_*`, `n_no_hit_reads`, `n_failed_reads`), `fraction_<class>` per BLAST diagnostic class, `median_target_minus_related_bitscore`, `median_target_minus_host_bitscore`, `blast_status`, `competitor_manifest_sha256`. Fractions use the **sampled UMI representatives** (one read per molecule, `viral_reads.umi_dedup.fasta`) as the denominator (`blast_fraction_denominator`), so they are molecule fractions and are blank when BLAST failed; molecule counts come from the exact lineage | diagnostic flag |
| `blast_identity.tsv` | Per sampled read: `diagnostic_class`, `top_class`, tied hits (`tie_classes`, `top_tied_hits`, `n_top_tied_hits`), `search_status`, and best hit, bitscore, identity and coverage per class (`target`, `same_group`, `related`, `host`, `decoy`). Reads with no hit are kept as `no_hits` rows | diagnostic flag |
| `competition_failure.txt` | Present only when competitive alignment or BLAST failed; lineage and status artifacts are kept and `competition.qc_status` is `failed` | diagnostic flag |

Candidate search is limited to the alignments minimap2 emitted
(`-ax sr`, default secondary policy: score ratio 0.8, limit 5), so
`candidate_classes` is evidence of competition, never an exhaustive search.
`unresolved` molecules are not counted for or against any reference.

The `competition` block of `evidence_manifest.json` records:

| Field | Description | Kind |
|-------|-------------|------|
| `mode`, `manifest_sha256` | `explicit_manifest` and the SHA-256 of the competitor manifest file | observation |
| `sequence_hash_normalization` | How FASTA records were hashed for the manifest check (`uppercase_remove_whitespace_DNA_IUPAC`) | observation |
| `tie_delta`, `sampling_seed`, `max_sampled_reads` | BLAST tie tolerance, sampling seed, and the sampled-read cap (500) | observation |
| `molecule_representative`, `candidate_search` | Representative rule (`corrected_CB_UMI_reference_rank_v1`) and search scope (`emitted_alignments_only`) | observation |
| `minimap2_options`, `minimap2_candidate_policy` | Aligner options and the statement that the secondary policy is not exhaustive | observation |
| `blast_status`, `qc_status` | `not_assessed`, `assessed`, or `failed` | observation |
| `extraction_status`, `n_extracted_reads`, `n_candidate_molecules` | `extracted` or `empty`, with the extracted read count; `n_candidate_molecules` is written only on the `empty` path (0) | observation |
| `interpretation` | Always `diagnostic_only` | diagnostic flag |
| `molecule_verdict` | `not_produced_in_competitor_mode` (see above) | observation |
| `summary`, `failure` | Copy of the `competitor_summary.tsv` row; failure message when `qc_status` is `failed` | diagnostic flag |

`r1_tso_fraction` is deliberately not reported here: exact-lineage reads carry
corrected on-list barcodes only, so a TSO-in-R1 fraction would be about zero by
construction.

Contamination, expected 3-prime bias, and subgenomic-RNA-like patterns remain
`not_assessed` unless a suitable negative-control or target-specific model is
available. No diagnostic flag is an automatic biological conclusion.

---

## `hostresponse/` — host-response outputs

Written by `viralscan hostresponse` (or a main run with `--host-h5ad`) under
`<output>/hostresponse/`. See the CLI reference for the flags; the key point is
that the raw candidate-support label can be **depth-confounded**, so interpret
the model AUC together with the depth baseline, cohort balance, and sensitivity
metrics.

Each virus is one feature: the viral gene columns of the run are pooled (summed) per
virus before anything else, using the run's `virus_identity.tsv` (or, for a legacy
`log/analysis.txt`, the name rule), so another gene of the same virus never counts as a
negative. `virus` below is therefore the pooled virus **display name** (for example
`Human herpesvirus 1`), not a gene or accession ID, and the per-virus file names use it
with characters other than letters, digits, `_` and `-` replaced by `_`
(`Human_herpesvirus_1_gene_weights.csv`). Two viruses whose names sanitise to the same
string are refused rather than overwritten.

**`hostresponse_metrics.csv`** — one row per virus:

| Column | Description | Kind |
|--------|-------------|------|
| `virus`, `n_positive` | Virus display name (pooled over its genes); number of candidate-positive cells under the chosen label. With the default `raw` label a cell is positive when its pooled viral molecule estimate is at least `--hostresponse-min-umi` (default 10; separate from `--detection-threshold`, the sample gate) | observation |
| `cv_mode` | The concrete evaluation scheme the row was computed with, `group` or `cell` (`--cv auto`, the default, is resolved before the fit; see the manifest's `cv_mode_requested`). Cell-mode and group-mode numbers are different estimands | observation |
| `label`, `depth_matched`, `mito_controlled` | Which de-confounding design was used (`raw`/`cpm`/`fraction`; depth-matched cohort; %mito covariate) | observation |
| `auc_mean` / `sensitivity_*` / `specificity_*` / `balanced_acc_*` / `mcc_*` | Held-out model metrics (mean/SD over seeds) | model estimate |
| `depth_alone_auc_mean` | **AUC from sequencing depth ALONE** under the identical split. If this ≈ `auc_mean`, the model AUC is a depth artifact | diagnostic flag |
| `n_stable_genes`, `n_genes_evalue_ge2` | Genes at or above `--stab-min-prob` (uncalibrated selection frequency; the column name is kept), and how many survive depth adjustment with an E-value ≥ 2 | model estimate |
| `n_differential_fdr05` | (with `--differential`) genes significant at FDR < 0.05 in the selected host-gene family. Counts testable genes only: a gene with no residual variance has a NaN `fdr` and is not counted | model estimate |
| `observed_cv_auc` | (group CV with `--hostresponse-permutations N`) the observed statistic the permutation null is compared with: the mean over seeds of the fold-mean held-out AUC (`auc_statistic` in `hostresponse_manifest.json`). It equals `auc_mean` of the same row. Called `observed_group_auc` before 2026-10-11 although it is not group-specific | model estimate |
| `empirical_p_auc` | (group CV with `--hostresponse-permutations N`) the permutation p-value of `observed_cv_auc`: `(1 + k) / (1 + n)` over `n` replicates, `k` of them at least as large. Absent unless every replicate was valid; see `permutation_status` | model estimate |
| `null_auc_mean`, `null_auc_sd`, `observed_minus_null_auc` | (group CV with `--hostresponse-permutations N`) mean and SD of the null AUCs, and `observed_cv_auc` minus the null mean | model estimate |

**`<virus>_gene_weights.csv`** — per-fold-HVG L2 coefficients (`weight_mean/sd`,
`n_folds_selected`). **`<virus>_stability.csv`** — per-gene randomized-Lasso
**uncalibrated selection frequency** (`stab_prob`; `stable` = `stab_prob >= --stab-min-prob`, a
ranking flag, not error control: C is fixed, so noise genes pass it with many positives, and
a warning is logged above 300 positives). Both gain a `symbol` column with
`--gene-symbols`.
A failed requested selection-frequency fit stops the analysis; failed fits are not
silently removed from the selection-frequency denominator.

**`<virus>_depth_diagnostics.csv`** — per above-threshold gene, the depth- (and, by
default, %mito-) adjusted odds ratio (`adj_OR`) and its `E_value` (the confounder
strength needed to explain the association away under the stated sensitivity
model). It does not prove absence of residual or unmeasured confounding.

**`<virus>_differential.csv`** (with `--differential`) — a
depth/%mito-adjusted differential table over selected host HVGs, or all host
genes when HVG selection is disabled: `partial_r`,
`p_value`, `fdr` (Benjamini-Hochberg), `direction` (`up`/`down`/`none`).
`fdr_family` names that selection (`selected_host_hvg` or `all_host_genes`);
`n_genes_tested` and `n_finite_pvalues` report the tested and finite-p-value
family sizes. `inference_unit = cell` and `donor_adjusted = false` make the
inferential scope explicit. These p-values do not provide donor-level inference.
The t-test has `n - 2 - k` degrees of freedom (`k` residualised covariates: log depth,
and %mito when used, plus the intercept and the label), as in an OLS fit. A gene with no
residual variance after residualisation (or a label with none) is untestable: `partial_r`,
`p_value` and `fdr` are NaN and `direction` is `none`. Benjamini-Hochberg runs over the
testable genes only, so untestable genes do not change the other genes' `fdr`.

`hostresponse_manifest.json` records `barcode_alignment` with policy
`host_cell_universe_zero_fill`. The host AnnData defines the cell set. Host cells
absent from the viral matrix receive zero viral counts; viral barcodes absent
from the host matrix are dropped. The record includes `n_host`, `n_viral`,
`n_shared`, `n_viral_zero_filled`, `frac_viral_zero_filled`,
`n_viral_dropped_not_host_called`, and `frac_viral_dropped_not_host_called`.
An empty host cell set or a nonempty viral matrix with no matching host barcode
raises an error. Barcode strings must match exactly.

---

## AnnData files (`.h5ad`)

The AnnData objects can be loaded with [scanpy](https://scanpy.readthedocs.io/):

```python
import scanpy as sc

adata = sc.read_h5ad("output/sample/kb-python/counts_unfiltered/adata_multimap.h5ad")
print(adata)
# Layers: counts_unique, counts_ambiguous_allocated, plus diagnostic layers
```

Key layers:

| Layer | Description | Kind |
|-------|-------------|------|
| `counts_unique` | Bustools-resolved one-gene molecule counts | observation |
| `counts_ambiguous_allocated` | Selected-method ambiguous molecule allocation | model estimate |
| `counts_multimap_equal` | Equal-split ambiguous molecule allocation | model estimate |
| `counts_multimap_host_conservative` | Allocation of mixed host-virus molecule mass only among compatible host genes | model estimate |
| `counts_multimap_unique_weighted` | Heuristic allocation weighted by unique-gene evidence plus pseudocount | model estimate |
| `counts_unique_viral` | Unambiguous viral molecule evidence retained for auditing | observation |
| `counts_host_viral_ambiguous` | Viral-compatible host-virus ambiguous signal | observation |
| `counts_host_viral_selected` | Portion of the selected allocation assigned to viral genes from host-virus ambiguous ECs | model estimate |

`adata.X` = `counts_unique + counts_ambiguous_allocated` (complete selected-method matrix).

---

## Multiple samples

When `--sample1` and `--sample2` contain comma-separated FASTQ lists,
ViralScan processes each pair separately:

```bash
viralscan \
  -t t2g.txt -i index.idx -o output/ \
  -s1 A_R1.fastq.gz,B_R1.fastq.gz \
  -s2 A_R2.fastq.gz,B_R2.fastq.gz
```

Expected directories:

```text
output/A/
output/B/
```
