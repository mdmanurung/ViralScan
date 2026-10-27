# Kraken2 taxonomic discovery (`viralscan discover`)

`viralscan discover` classifies the reads a completed run could not place on the host
genome with [Kraken2](https://github.com/DerrickWood/kraken2), then groups them by cell
barcode and UMI. It reports which taxa those reads resemble in which cells. It is a
separate command with its own output directory and its own fingerprint. It never writes
into the source run, the H5AD layers, `detected` or the evidence tiers.

## What it is, and is not

**A Kraken2 classification is not an infection call.** A read is assigned to the taxon
whose database k-mers it shares. That says the read resembles the database sequence. It
does not say the organism is present, replicating, or in that cell.

- Use it to see what else is in the host-unmapped reads and to find taxa the targeted
  reference did not contain.
- Do not use it to replace the primary viral load. The `viral_summary.tsv` molecule
  counts and tiers are unchanged by it.
- Every row of `candidates.tsv` has `state = observation`: a Kraken2 assignment alone never
  goes higher. `candidate_states.tsv` (see "Candidate states") can promote a candidate only
  with orthogonal evidence. No calibrated threshold or tier exists yet for discovery; the
  flags and cut-offs are documented screens.

## Prerequisites

- A completed ViralScan run with host filtering (`host_filtered/` present) and
  `results/called_cells.tsv`. A run without `host_filtered/` needs
  `--host-subtract starsolo --host-index <STAR genomeDir>`; those reads go under `-o`.
- `--called-cells FILE` replaces the run's own `called_cells.tsv`. A two-step run calls cells on barcodes
  that carry viral reads, so an uninfected library can have only a handful; pass the EmptyDrops
  `called_cells.tsv` of a full-matrix run of the same library instead. Its barcode space comes from the
  `cell_calling_input.json` beside that file (`unknown` without one); its path and sha256 are recorded in
  `provenance.json` and the sha256 enters the discovery fingerprint.
- Kraken2 2.17.2 and `bustools` on `PATH`. Install the optional environment and check it:

  ```bash
  conda env create -f environment.discovery.yml   # kraken2 only; layer it on the main env
  viralscan doctor --profile discovery
  ```

- A Kraken2 database directory (`hash.k2d`, `opts.k2d`, `taxo.k2d`). ViralScan never
  downloads or modifies it. Download one yourself, for example the Standard plus
  protozoa/fungi/plant database capped at 16 GB (PlusPF-16) from the
  [Kraken2 index listing](https://benlangmead.github.io/aws-indexes/k2), and check it
  against the provider's md5 before use:

  ```bash
  mkdir -p k2_pluspf_16 && cd k2_pluspf_16
  tar -xzf ../k2_pluspf_16gb_<date>.tar.gz
  # compare with the md5 the provider lists for that archive
  md5sum ../k2_pluspf_16gb_<date>.tar.gz
  ls hash.k2d opts.k2d taxo.k2d
  ```

  If the directory holds an `*.md5` or `md5sum.txt` file, `discover` re-checks the listed
  files and records the result in `provenance.json`.

## Command

```bash
viralscan discover --run-dir run/ --sample SAMPLE \
    --kraken2-db k2_pluspf_16/ --threads 8
```

### From `viralscan run`

```bash
viralscan run ... --host-filter starsolo --host-index star_hg38/ \
    --discover --kraken2-db k2_pluspf_16/
```

runs the same discovery on every sample after the run completes, with the defaults below,
into `<output>/posthoc/discovery/<sample>/`. The primary results are written and
`run_complete.json` exists before discovery starts; a discovery failure exits non-zero but
leaves the run valid. For other settings run `viralscan discover` afterwards.

`--sample` is required only when the run holds more than one. Outputs go to
`run/posthoc/discovery/SAMPLE` (layout version 2; `validate-run` ignores `posthoc/`), or to
`-o DIR` outside the run. Either must be empty, absent or a previous publication. The defaults are
`--confidence 0.05` and `--minimum-hit-groups 4` (the Cellsnake documented defaults, not
calibrated here; Kraken2's own are 0.0 and 2) and `--polya-trim` on. See the
[CLI reference](cli_reference.md) for every flag.

Re-running with the same inputs and settings exits 0 and does nothing. Changing the
database, an input or a setting changes the fingerprint; the command then refuses to
overwrite unless you pass `--force`. A `.partial/` directory left by an interrupted run is
cleared and the run restarts.

What the fingerprint does and does not cover:

- `--threads`, `--memory-mapping` and file paths are not in it, so changing them is still
  a no-op.
- The ViralScan version is in it, so after an upgrade a rerun needs `--force` or a new
  `-o`.
- With `--host-subtract`, the STAR index is identified by its path, not its contents. If
  you rebuild an index in place, use a new `-o`.
- `--force` deletes the previous outputs before the new run starts. If that run fails, the
  old results are gone.
- A crash during the final move of files into `-o` leaves outputs without
  `discovery_complete.json`. The next run refuses that directory, so empty it by hand.

### Slurm

Kraken2 loads the database into memory unless `--memory-mapping` is set, so request
about the database size plus headroom. Use a partition without a per-user CPU/memory cap
that is smaller than this; on the LUMC cluster that is `all`, not `medium` (its QOS caps a
user at 4 CPUs and 8 GB).

```bash
#!/bin/bash
#SBATCH --job-name=viralscan-discover
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00

viralscan discover --run-dir run/ --sample SAMPLE -o discovery/ \
    --kraken2-db k2_pluspf_16/ --threads "$SLURM_CPUS_PER_TASK"
```

For scale, one EBV library (17.76 M host-unmapped reads, PlusPF-16, 8 CPUs) took 35 to 37
minutes and 15.2 GB peak memory per run. A rerun that finds the same fingerprint still
hashes the database (about 8 minutes, 0.2 GB). These are single observations to size a
job, not performance claims. Use `--memory-mapping` only on fast local storage.

## Real-data check: one EBV library

EBV-positive 10x 3' v2 library SRR12682296 (PlusPF-16 2026-06-26, 2,524 called cells), run
three ways with the source run's files byte-identical before and after. Read level and one
library only; this describes the settings, it does not measure accuracy.

| | default 0.05 / 4, trim | 0.0 / 2, trim | 0.05 / 4, no trim |
|---|---:|---:|---:|
| Reads unclassified | 93.4% | 71.3% | 92.9% |
| Reads human (clade) | 4.1% | 20.9% | 4.5% |
| Reads EBV (clade) | 2.3% | 6.8% | 2.3% |
| EBV UMI groups in called cells | 246,233 | 623,096 | 246,882 |
| Called cells with EBV (ViralScan: 2,524) | 2,224 | 2,471 | 2,227 |
| Per-cell Spearman with ViralScan EBV molecules | 0.890 | 0.966 | 0.891 |
| Median discovery / ViralScan EBV per cell | 0.26 | 0.82 | 0.27 |
| Candidate taxa (flagged low minimizer diversity) | 47 (26) | 1,906 (1,680) | 56 (34) |

The default settings cut the human-labelled reads by a factor of five and the candidate
list by a factor of forty, and they also drop about 60% of the EBV signal. Use
`--confidence 0.0 --minimum-hit-groups 2` when sensitivity matters more than a short
candidate list. PolyA/T trimming changed 2.03 M reads (11%) but moved the human fraction
only from 4.5% to 4.1%.

## Outputs

All files are staged in `-o/.partial/` and moved into place; `discovery_complete.json` is
moved last, so its presence means the publication is whole.

| File | Content |
|------|---------|
| `reads.tsv.gz` | One row per fragment, keyed by `ordinal`: identity, `r2_trim` (`none`, `polyA`, `polyT`, `both`) and `r2_len_classified` |
| `kraken_input.fa.gz` | The R2 sequences sent to Kraken2, named by ordinal |
| `kraken.out.gz`, `kraken_report.txt` | Raw Kraken2 output and its report (`--report-minimizer-data`) |
| `classifications.tsv.gz` | The per-read calls in ordinal order |
| `barcodes.tsv.gz` | Corrected barcode, UMI and barcode status per fragment |
| `umi_groups.tsv.gz` | One row per (`corrected_cb`, `umi`): read counts, `taxids`, `resolved_taxid`, `group_class`, `lineage_consistent`, `deepest_consistent_taxid`, `conflict`, `called` |
| `taxonomy.tsv` | The database's own taxonomy: `taxid`, `parent_taxid`, `rank`, `name`, `depth`, `host_relation` |
| `cell_taxon.h5ad` | Called cells x resolved taxa; `X` counts resolved UMI groups; `var` has name, rank, parent, `is_host_lineage` |
| `offcell_taxon.tsv.gz` | Resolved groups on barcodes that are not called cells |
| `candidates.tsv` | One row per taxon with a resolved group (below) |
| `qc.json`, `provenance.json`, `discovery_complete.json` | Counts, provenance, and output sha256s with the discovery fingerprint |

`group_class` is `resolved` (one taxon), `ambiguous` (the result is above species level,
a Kraken2 `R*` rank) or `unclassified`. When every classified read of a group lies on one
root-to-leaf path, the group takes the deepest taxon; when they conflict it takes their
lowest common ancestor. The audit columns stay in the table. Only called-cell groups
enter `cell_taxon.h5ad`; ambiguous and unclassified groups never do.

`candidates.tsv` columns: `taxid`, `name`, `rank`, `groups_called`, `groups_offcell`,
`n_cells`, `reads`, `distinct_minimizers`, `minimizer_per_read`, `host_relation`,
`db_has_9606`, `state`, and the 0/1 flags `low_minimizer_diversity` (fewer than 20
distinct minimizers), `offcell_dominant` (more than half the groups off called cells),
`host_lineage` and `viral_only_db` (taxid 9606 absent from the database). The cut-offs are
documented constants, not calibrated.

### Candidate states

```bash
viralscan discover-states run/posthoc/discovery/S1 --run-dir run/
```

writes `candidate_states.tsv` beside `candidates.tsv` (`taxid`, `name`, `rank`, `state`,
`corroboration`, `read_support`, `reason`) and `candidate_states.provenance.json` (sources
used, `k`, the constants, `calibrated: false`). Route (i) takes one source by availability:
the run's pruned `star_hostvirus/viral.bam` (`--star-confirm`), else
`results/molecule_assignments.tsv.gz` (`--multimap-molecule-assignments`), else
`not_assessed`. Route (ii) runs only with `--minimap2 --host-fasta F --candidate-fasta
TAXID=FASTA` (at most 3). Without `--run-dir` no evidence is produced and every taxon is
`observation` or `prioritized`. States:

- `observation`: Kraken2 only, or evidence that did not pass.
- `prioritized`: Kraken2 only, at least `--k` resolved called-cell groups (default 3, the
  DISC-12 blind rule; reported, not calibrated), no host lineage, not a viral-only database.
  Evidence that came back negative does not demote a `prioritized` taxon; its explanation is
  appended to `reason` (for example `prioritized(k=3): Kraken2 only; corroboration: 0/10 ...`).
- `corroborated` (route i): the same taxon is supported by the kb (panel) arm or the STAR
  host+virus arm. The molecules were competitively adjudicated (not better explained by
  host or another virus) and fall in cells that are also called-cell groups of the taxon.
- `read_supported_observation` (route ii, kb-independent, for off-panel taxa): reads align
  competitively to the hit genome versus host + panel, with shared barcodes and UMIs.
  It is still an observation, never an infection claim.

`corroboration` and `read_support` are `not_assessed` (no evidence row), `negative`
(evidence supplied, rule not met) or `positive`; `not_assessed` is never `negative`.
Host-lineage taxa and a viral-only database stay at `observation`.

The package does not yet produce a per-taxon alignment summary, so the adapter reads a
documented tab-separated table: `taxid`, `route` (`corroboration` or `read_support`),
`source` (`kb` or `star` for route i; a free label for route ii), `n_adjudicated` (units
tested), `n_taxon_best` (units that scored best for the taxon alone; ties and wins by
another reference are not counted), `n_shared_cells`, `n_shared_umis` (cells and distinct
(CB, UMI) molecules also resolved called-cell groups of the taxon). The caller builds it
from existing kb or STAR outputs, or from its own competitive alignment. A route is positive
with at least 3 taxon-best units, at least 80% of the adjudicated units, and at least one
shared cell (route ii also one shared molecule). These cut-offs are documented constants
(`MIN_TAXON_BEST`, `MIN_TAXON_BEST_FRACTION`), not calibrated. Malformed, duplicate or
unmatched evidence rows abort.

`build_kb_evidence(discovery_dir, sample_dir)` (what `discover-states` runs for `source=kb`) writes the route (i)
`source=kb` table `evidence_kb.tsv` from the primary run's
`results/molecule_assignments.tsv.gz` (needs `--multimap-molecule-assignments`) and
`results/virus_identity.tsv`. A viral panel gene belongs to a candidate when its `taxid` equals
the candidate or sits on one root-to-leaf lineage with it in `taxonomy.tsv` (strain gene for a
species candidate, or the reverse); host genes never match. Per distinct (CB, UMI) touching the
taxon, `n_adjudicated` counts it whatever else it is compatible with, and `n_taxon_best` only when
every gene is the taxon, a descendant of it, or (taxon below species) a species-or-below ancestor.
A gene at an ancestor taxid above species rank (genus, family, root bins; ranks from `taxonomy.tsv`)
still counts in `n_adjudicated` (it competes) but never makes a molecule taxon-best. `n_adjudicated`
and `n_taxon_best` cover all barcodes; the shared counts cover called cells only. Shared cells and
(CB, UMI) are those of the taxon-best molecules that are resolved called-cell groups of exactly that taxid. A candidate with no panel
gene gets no row (`not_assessed`); one with panel genes but no molecule gets a zero row
(`negative`). An empty molecule file, or one with no molecule in a called cell, aborts. It refuses unless the called cells are in kb barcode space (a STARsolo or unknown
space, or kb molecules sharing no barcode with the called cells, aborts). 
`build_star_evidence(discovery_dir, sample_dir)` (run by `discover-states` when the STAR BAM exists) is the `source=star`
counterpart (`evidence_star.tsv`). It reads the run's pruned `star_hostvirus/viral.bam`
(`VIRUS|` records only), maps `VIRUS|<accession>` to taxids through `results/virus_identity.tsv`
(same lineage rule as above) and counts molecules as distinct (CB, UB) with the `-N` suffix of the
CB stripped; records with CB or UB missing or `-` are not molecules, and since STARsolo Unique mode
gives multimappers no UB, counts are a lower bound. Host placements are invisible in the pruned BAM,
so a read is taxon-best for a taxon only when its visible records on matching contigs equal its `NH`
(above-species ancestors never win); a molecule is taxon-best when all its reads are. No called-cell
overlap, a missing BAM or no molecule aborts.

`build_read_support_evidence(discovery_dir, sample_dir, r1, r2, host_fasta, {taxid: candidate_fasta})`
(run by `discover-states --minimap2`; no genome download) writes the route (ii) `source=minimap2` table
`evidence_read_support.tsv` for at most `max_candidates` (default 3) taxa. `r1`/`r2` must be the
resolved pair discovery classified (read IDs are checked). The R2 reads of the taxon's resolved
called-cell groups are aligned competitively to host, the candidate FASTA and optionally a panel
FASTA. `n_adjudicated` counts every extracted read, unmapped ones included; `n_taxon_best` only reads
whose best alignment score on the candidate is strictly above any host or panel score (ties and host
or panel wins are not). A taxon with no resolved called-cell read gets no row (`not_assessed`).

To count at species, genus or family level, run
`viralscan discover-collapse OUT --rank {species,genus,family}`. It sums columns into the
rank ancestor without reclassifying; taxa above the rank stay as themselves with
`above_rank = True`.

### Reading `qc.json`

`qc.json` carries every count the run checks on itself, so a reader can see where reads
went.

- `stream`: fragments read from the pair, bad geometry, barcodes or UMIs with `N`, and
  reads and bases trimmed by `--polya-trim`.
- `kraken2`: `n_classified` and `n_unclassified`, which add to `stream.n_input`.
- `barcodes`: how many fragments had an on-list barcode, a corrected one, or none
  (`no_onlist`, for chemistries without an on-list).
- `groups`: resolved, ambiguous and unclassified UMI groups and the reads they hold.
- `cells`: `n_called_cells`, how many of them appear in the reads
  (`n_called_cells_seen`), resolved groups in and off called cells, `n_candidates`.
- `called_cells_overlap`: `n_called_cells_seen / n_called_cells`. A low value means the
  called-cell barcodes and the reads do not share a barcode space; check it first.
- `fractions`: classified and unclassified fractions of the input reads.

The counts must add up, or the run stops rather than publish. The called-cell set is the
primary `results/called_cells.tsv`, which can differ from the `n_called_cells` in
`viral_summary.tsv` when that column comes from the host-only calling.

`provenance.json` holds the primary run fingerprint, input and database hashes (the three
`.k2d` files and the taxonomy), Kraken2 and bustools versions, settings, and the on-list
and called-cells hashes. `discovery_fingerprint` is the hash of the part that decides
whether a rerun is the same.

## Limits

- **Residual host reads.** STAR leaves reads it could not align well, so Kraken2 still
  labels many reads human. In the EBV library above, 4.1% (default settings) to 20.9%
  (0.0 / 2) of host-unmapped reads were human, and a relaxed realignment of a 5,000-read sample put 97.8% on GRCh38, mostly as
  partial alignments (median 71% of the read). PolyA/T trimming removes only a minority
  (17.5% polyA and 5.3% polyT in that sample). Host-lineage taxa are labelled
  (`host_relation`) and kept out of candidates, but the human fraction is not zero.
- **The confidence denominator.** Kraken2 2.17.2 divides by all k-mers of the read,
  including those with an ambiguous base (`N`), although its manual says they are
  excluded. A read with 11 k-mers on a taxon and 45 ambiguous ones is unclassified at
  `--confidence 0.5`. A read with `N` bases can fail the threshold for reasons unrelated to
  the organism. `--confidence` is a k-mer fraction, not a probability.
- **Read level, not molecule level.** Each read is classified alone. A UMI group is the
  unit counted in `cell_taxon.h5ad`, so counts are not comparable one to one with the
  molecule counts of the primary run, which resolve multimapping against its own
  reference. See the EBV table above for how the two compare on one library.
- **Database size.** A 16 GB-capped database keeps only part of the full library, so a
  taxon can be missed or placed at a higher rank than a larger database would. Another
  database or another build date gives different results; the database hashes are
  recorded for that reason.
- **Low-complexity and repeated sequence.** A taxon supported by few distinct minimizers
  can be a single repeated locus. In the EBV library most bacterial candidates looked like
  this; `low_minimizer_diversity` flags them.
- **A negative result is not absence.** No resolved group in a cell, or no taxon in the
  table, does not show the organism is absent: it may be missing from the database,
  lost to host subtraction or the capture chemistry, or below the thresholds.
- **Not calibrated.** The confidence and hit-group defaults and the flag cut-offs are
  not tuned against a benchmark, and no discovery performance claim is made.
