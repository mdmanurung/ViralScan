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
- Every row of `candidates.tsv` has `state = observation`. No calibrated threshold or
  tier exists yet for discovery; the flags are documented screens.

## Prerequisites

- A completed ViralScan run with host filtering (`host_filtered/` present) and
  `results/called_cells.tsv`. A run without `host_filtered/` needs
  `--host-subtract starsolo --host-index <STAR genomeDir>`; those reads go under `-o`.
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
viralscan discover --run-dir run/ --sample SAMPLE -o discovery/ \
    --kraken2-db k2_pluspf_16/ --threads 8
```

`--sample` is required only when the run holds more than one. `-o` must be outside
`--run-dir`, and empty, absent or a previous publication. The defaults are
`--confidence 0.05` and `--minimum-hit-groups 4` (the Cellsnake documented defaults, not
calibrated here; Kraken2's own are 0.0 and 2) and `--polya-trim` on. See the
[CLI reference](cli_reference.md) for every flag.

Re-running with the same inputs and settings exits 0 and does nothing. Changing the
database, an input or a setting changes the fingerprint; the command then refuses to
overwrite unless you pass `--force`. A `.partial/` directory left by an interrupted run is
cleared and the run restarts.

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

For scale, the first EBV run (one library, 17.76 M host-unmapped reads, PlusPF-16, 8 CPUs,
confidence 0.0 and 2 hit groups) took 26 minutes and 15.9 GB peak memory. That is one
observation to size a job, not a performance claim. Use `--memory-mapping` only on fast
local storage.

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
  labels many reads human. In the first EBV run about 23% of host-unmapped reads were
  human, and a relaxed realignment of a 5,000-read sample put 97.8% on GRCh38, mostly as
  partial alignments (median 71% of the read). PolyA/T trimming removes only a minority
  (17.5% polyA and 5.3% polyT in that sample). Host-lineage taxa are labelled
  (`host_relation`) and kept out of candidates, but the human fraction is not zero.
- **The confidence denominator.** Kraken2's `--confidence` divides by the k-mers of the
  read that fall in the database taxonomy, not by all k-mers. A read with few usable
  k-mers can fail or pass the threshold for reasons unrelated to the organism.
  `--confidence` is a k-mer fraction, not a probability.
- **Read level, not molecule level.** Each read is classified alone. A UMI group is the
  unit counted in `cell_taxon.h5ad`, so counts are not comparable one to one with the
  molecule counts of the primary run, which resolve multimapping against its own
  reference. In the first EBV run (read level, one library, not a performance claim) the
  per-cell EBV counts of the two methods had Spearman 0.966 across called cells, and
  discovery found EBV in 2,474 of 2,524 called cells against 2,524 for ViralScan.
- **Database size.** A 16 GB-capped database keeps only part of the full library, so a
  taxon can be missed or placed at a higher rank than a larger database would. Another
  database or another build date gives different results; the database hashes are
  recorded for that reason.
- **Low-complexity and repeated sequence.** A taxon supported by few distinct minimizers
  can be a single repeated locus. In the first EBV run the bacterial hits looked like
  this; `low_minimizer_diversity` flags them.
- **A negative result is not absence.** No resolved group in a cell, or no taxon in the
  table, does not show the organism is absent: it may be missing from the database,
  lost to host subtraction or the capture chemistry, or below the thresholds.
- **Not calibrated.** The confidence and hit-group defaults and the flag cut-offs are
  not tuned against a benchmark, and no discovery performance claim is made.
