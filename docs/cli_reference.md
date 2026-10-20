# CLI Reference

Run `viralscan --help` or `viralscan <subcommand> --help` to see options for
the installed version. This page documents the public CLI as of ViralScan
**3.0.0.dev1**. Available subcommands: `build-ref`, `data fetch`, `evidence`, `evidence-plot`,
`rerun-multimap`, `doctor`, `validate-run`, `hostresponse`, and
`check-whitelist`.

---

## `viralscan run` — viral quantification

```
viralscan run [OPTIONS]
```

A bare `viralscan [OPTIONS]` (no subcommand) still runs `run`, with a deprecation warning.

Use this command with paired FASTQ files. ViralScan validates inputs, prepares
or reuses the reference, then dispatches the Snakemake workflow.

For host-aware viral detection, the recommended workflow is to build a combined
host+virus reference with `viralscan build-ref`, then pass its `index.idx` and
`t2g.txt` to this command. The default multimapping method is
`conservative` (canonical name `host-conservative`), which keeps host-virus ambiguous EC mass out of primary
viral estimates — the conservative choice for combined host+virus references where
host-virus cross-homology can create candidate evidence. Use
`--multimap-method equal` for an explicit equal-allocation comparison, or
`--multimap-method em` for iterated, cell-local allocation. `--multimap-method` shows three options
(`conservative`, `em`, `equal`). They are CLI aliases: the parser resolves them to the canonical engine
names `host-conservative`, `equal` and `em-cell`, and run manifests, `config.yaml`, h5ad `uns` and outputs
record only the canonical name. The expert values `unique-weighted`, `sibling-weighted`, `em-global` (and the
canonical names themselves) are still accepted but not listed. `em-global` can leave a small residue on a gene with zero unique support at extreme count ratios and its tolerance matters; see
[the counting contract](v3_counting_contract.md) before using it.

### Input / output

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--output PATH` | `-o` | *(required)* | Root output directory. Each FASTQ pair is written to a sample subdirectory. |
| `--r1 PATHS` | | *(required)* | R1 FASTQ path, or comma-separated R1 paths for multiple samples. |
| `--r2 PATHS` | | *(required)* | R2 FASTQ path, or comma-separated R2 paths matching `--r1`. |

### Reference (choose one of three modes)

| Flag | Short | Description |
|------|-------|-------------|
| `--index PATH` | `-i` | Pre-built kallisto index from `kb ref`. Required in pre-built mode. |
| `--transcripts PATH` | `-t` | `t2g.txt` transcript-to-gene map from `kb ref`. Required in pre-built mode. |
| `--f1 PATH` | | Optional cDNA FASTA passed through for kb workflows that need it. |
| `--reference` | | Build index from `--fasta` + `--gtf` |
| `--fasta PATH` | | FASTA file(s), comma-separated |
| `--gtf PATH` | | GTF file(s), comma-separated |
| `--ncbi-accession ACC` | | NCBI accession(s), comma-separated; fetch + build |
| `--ncbi-email EMAIL` | | Contact email for NCBI (or `$NCBI_EMAIL`); needed only when a sequence must be downloaded, not for a cache hit |
| `--cache-dir PATH` | | Viral annotation cache root for the bundled panel (`PATH/data/`); read only for an index without a `<index>.build_manifest.json` |

Reference modes are mutually exclusive:

- Pre-built: use `-i` and `-t`.
- FASTA/GTF: use `--reference --fasta ... --gtf ...`.
- NCBI: use `--ncbi-accession ...` and provide an NCBI email.

`--reference`, `--visual` and `--multimapping` take no value: use the `--no-` form to turn
one off. The original single-dash and value forms (`-ref True`, `-v False`, `-mm False`) were
removed in 3.0 and now fail with the replacement named; see [migration.md](migration.md).

### Analysis parameters

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--technology STRING` | `-x` | `10xv3` | Single-cell technology (`kb --list` for all) |
| `--whitelist PATH` | `-w` | *(bundled)* | Barcode whitelist file |
| `--strand {forward,reverse,unstranded,auto}` | | *(kb default)* | Read strandedness passed to `kb count --strand`. 10x 5′ libraries need `reverse` or `unstranded` (F-020). `auto` (opt-in) pilots all three on the first 1M read pairs per sample and picks `reverse` or `forward` if it keeps ≥0.8 of the unstranded pseudoalignment rate (larger ratio wins, tie → `unstranded`), else `unstranded`; recorded as `strand_inference` in `run_manifest.json` and reused on `--resume` |
| `--cores N` | `-c` | `6` | CPU cores |
| `--multimapping` / `--no-multimapping` | | on | Multimapping correction |
| `--multimap-method METHOD` | | `host-conservative` | Multimapper allocation: `conservative`, `em`, or `equal` |
| `--multimap-pseudocount FLOAT` | | `1.0` | Positive pseudocount for `unique-weighted` |
| `--multimap-primary-call MODE` | | `selected-method` | V3 fixed contract: summaries use the complete selected-method molecule matrix. |
| `--multimap-em-max-iter N` | | `100` | Maximum iterations for `em-global` or `em-cell` |
| `--multimap-em-tol FLOAT` | | `1e-6` | Convergence tolerance for `em-global` or `em-cell` |
| `--cell-calling METHOD` | | `auto` | `auto`, `emptydrops`, `external`, `knee`, or `none`; `auto` uses an external list when supplied and otherwise EmptyDrops on the host transcriptome only (needs R with DropletUtils; `Rscript` on `PATH`). `knee` is sensitivity-only and not for reported numbers: on real libraries it lands at the 10-molecule floor and calls empty droplets as cells |
| `--called-cells-file PATH` | | *(none)* | External called-cell barcodes used by `auto` or required by `external` |
| `--umap` | | off | Generate UMAP plot |
| `--visual` / `--no-visual` | | on | Generate visualisations |
| `--host-filter ALIGNER` | | *(none)* | Optional irreversible host subtraction before quantification. V3 supports `starsolo` only. |
| `--host-index PATH` | | *(none)* | Required with `--host-filter`; a full-host-genome STAR index directory. |
| `--host-filter-star-params SET` | | `pinned` | STAR filter parameters for `--host-filter starsolo`: `pinned` or `star-default` (STAR's own defaults). |
| `--read-filter MODE` | | `off` | Read-artefact filter before `kb count` (DEF-01): `off` or `artefact`. Drops TSO/reagent chimeras and reads with no complex body; also drops low-complexity viral reads. |

### Detection thresholds

| Flag | Default | Description |
|------|---------|-------------|
| `--detection-threshold N` | `1` | Min estimated viral molecules summed over the **called cells** for a virus to pass the detection gate; molecules in non-called barcodes are reported as `background_umis` |
| `--detection-min-cells N` | `1` | Distinct called cells that must be infected for a virus to pass the gate. Viruses seen only in the all-barcode screen stay in `viral_summary.tsv` as candidates (`detected` = False). With `--cell-calling none` every barcode counts as called |
| `--infected-min-molecules X` | `1.0` | A cell is infected for a virus when its summed allocated molecules reach X (one rule for detection, enrichment and the UMAP). Cells with 0 < weight < X are counted in `infected_subthreshold_cells`, not dropped |
| `--hostresponse-min-umi N` | `10` | Per-cell positive label for host-response (the sample gate stays `--detection-threshold`) |
| `--gene-programs` / `--no-gene-programs` | off | Second layer: for viruses already detected, infer the viral gene programme (latent vs productive) per cell and write `results/gene_program_summary.tsv` + `results/gene_program_cells.tsv` |
| `--programme-min-breadth N` | `2` | Distinct non-overlapping overlap groups required before a programme is called. Counted in overlap groups rather than genes because herpesvirus latent and lytic ORFs share exonic sequence; a per-gene comparison gives a latent:lytic ratio of 1.15 in a latently-infected cell line. Must be >= 1 |
| `--se-threshold N` | `10` | Legacy-named display threshold for high candidate molecule support; not a biological classification |
| `--cell-type-labels PATH` | *(none)* | CSV with `barcode,cell_type` columns for per-virus cell-type enrichment |

### Host-response analysis (optional)

Associate viral presence with host gene expression by providing a host h5ad.
When `--host-h5ad` is supplied, ViralScan trains per-virus L2 logistic
regression models (Luebbert et al. 2025) and runs a randomized Lasso selection-frequency
ranking (uncalibrated); results are written to `<output>/hostresponse/`. Pathway enrichment
requires `pip install "viralscan[enrichment]"`.

| Flag | Default | Description |
|------|---------|-------------|
| `--host-h5ad PATH` | *(none)* | Host gene-expression h5ad (cells × genes). Triggers host-response analysis when provided. |
| `--hostresponse-n-seeds N` | `6` | Number of random seeds for the multi-seed L2 regression |
| `--hostresponse-n-stab-iter N` | `100` | Iterations for randomized Lasso stability selection |
| `--hostresponse-use-hvg` / `--no-hostresponse-use-hvg` | on | Use highly variable genes as features (disable to use all genes) |
| `--hostresponse-stab-min-prob PROB` | `0.6` | Min selection probability to call a gene stably associated |
| `--hostresponse-top-n-genes N` | `50` | Top N stable genes to pass to pathway enrichment |
| `--enrichment` | off | Run pathway enrichment via gget (requires `viralscan[enrichment]`) |
| `--enrichment-db DB` | `GO_Biological_Process_2023` | gget.enrichr database |

### UMAP / QC parameters

These flags affect only the optional `--umap` workflow.

| Flag | Default | Description |
|------|---------|-------------|
| `--min-counts N` | `1000` | Min selected-method molecule estimate per cell (UMAP QC) |
| `--min-genes N` | `200` | Min detected genes per cell (UMAP QC) |
| `--hvg-min-mean X` | `0.0125` | Scanpy highly-variable-gene `min_mean` |
| `--hvg-max-mean X` | `3.0` | Scanpy highly-variable-gene `max_mean` |
| `--hvg-min-disp X` | `0.5` | Scanpy highly-variable-gene `min_disp` |
| `--umap-n-neighbors N` | `15` | Neighbor count for Scanpy graph construction |

### Verbosity

| Flag | Description |
|------|-------------|
| `--verbose` | Enable DEBUG-level logging |
| `--quiet` | Suppress INFO; show warnings and errors only |

### Output reuse and safety

The default mode refuses a non-empty output directory. Reuse must be explicit:

| Flag | Behavior |
|------|----------|
| `--resume` | Resume only when the existing `run_manifest.json` fingerprint exactly matches the current invocation and input bytes |
| `--overwrite` | Replace only the resolved output target after an explicit confirmation |
| `--yes`, `-y` | Answer the `--overwrite` confirmation, and the `--reference` prompt that asks whether to continue when `-i`/`-t`/`--f1` are also given (they are ignored; the reference is built into the output directory). It does not imply `--overwrite` or `--resume`. Without a terminal and without `--yes`, that prompt is an error, not a hang |

`--resume` and `--overwrite` are mutually exclusive. A missing manifest or
fingerprint mismatch makes resume fail without changing the output directory, and
the error lists the options and input files that differ (for example
`Options that differ: --strand.`).

The fingerprint covers the input bytes, not just their paths: the FASTQs, `-t`, `-i`,
`--f1`, `-w`, `--gtf`, `--fasta`, `--called-cells-file`, `--positive-control-receipt`, and also the
STAR host index directory (`--host-index`, hashed as a directory), `--host-h5ad`,
`--cell-type-labels`, `--host-fasta` and `--viral-fasta`. For `--ncbi-accession` runs the sha256 of the
fetched NCBI FASTA and GTF is recorded in the manifest after the fetch, and `--resume`
refuses a run whose accession now returns different bytes. A run that used any of these
newer inputs and was started on an earlier version has no matching fingerprint: finish it
with `--overwrite` rather than `--resume`. Plain runs that use none of them keep their
fingerprint. `--ncbi-email` and `--ncbi-api-key` are not recorded in the manifest.

---

## `viralscan build-ref` — reference builder

```
viralscan build-ref [OPTIONS]   # --cores N bounds kb ref / minimap2 / cd-hit-est (unset: legacy 8 / 3 / all CPUs)
```

Build a combined host + virus kallisto reference without running a full
analysis. The command downloads the host cDNA FASTA/GTF from Ensembl, viral
FASTA/GTF from NCBI, concatenates them, and runs `kb ref` unless
`--no-kb-ref` is supplied. This is the preferred host-aware setup for routine
runs.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--host SPECIES` | | *(none)* | Host species (see `--list-species`) |
| `--virus-accessions ACC [ACC ...]` | | *(none)* | One or more NCBI accessions separated by spaces |
| `--profile PROFILE` | | `curated` | Frozen profile label: `curated` or explicitly requested `broad-discovery` |
| `--output PATH` | `-o` | `viralscan_ref` | Output directory |
| `--ensembl-release N` | | *(current release)* | Pin the Ensembl release of the host cDNA FASTA/GTF (for example `116`). Without it the current release is read from `ftp.ensembl.org/pub/VERSION`; if that cannot be reached the build stops and asks for a pin instead of silently using an old release. The cache is keyed by release (`<cache>/ensembl/<species>/release-N/`), so two releases never share a file |
| `--ncbi-email EMAIL` | | *(none)* | Contact email for NCBI (needed only when a sequence has to be downloaded, not for a cache hit) |
| `--ncbi-api-key KEY` | | *(none)* | NCBI API key |
| `--cache-dir PATH` | | `~/.cache/viralscan/` | Download cache root |
| `--no-kb-ref` | | off | Stop after writing FASTA + GTF; skip `kb ref` |
| `--genome-dlist FASTA` | | *(none)* | Full host genome used as kallisto D-list and for raw viral host-homology measurements; requires minimap2. **See the warning below — measured to be the weakest of the three host-control options** |
| `--anellovirus` / `--no-anellovirus` | | off | Explicitly include the expanded packaged Anelloviridae table in a host+virus reference |
| `--allow-partial-panel` | | off | Permit an incomplete expanded panel and write the complete missing-accession report; default fails closed |
| `--reference-panel anellovirus` | | *(none)* | Build a predefined Anelloviridae panel (bundled FASTA if cached, else NCBI download) |
| `--no-mask` | | off | Skip dustmasker low-complexity masking of the viral sequence (every build path); the k-mer gate then allows a few low-complexity k-mers instead of none |
| `--cluster` | | off | Cluster anellovirus sequences at 95% identity (cd-hit-est) after masking |
| `--list-species` | | off | Print supported host species and exit |
| `--verbose` | | off | DEBUG-level logging |
| `--quiet` | | off | Warnings + errors only |

The default curated reference contains the host and explicit
`--virus-accessions`; expanded anellovirus inclusion is opt-in. Requested
masking, clustering, and index construction fail if their tools are missing.
Every successful build writes `reference_manifest.json` with a combined hash
and per-sequence identifier, source, retrieval time, digest, and length.
`kb ref` is always run with `--overwrite`, so rebuilding into a non-empty output
directory replaces the index and t2g instead of silently keeping the old ones. A host
cDNA cache written by an earlier version (without the `release-N/` level) is
downloaded again once.

### ⚠ `--genome-dlist` is the weakest host-control option, not the strongest

Measured on the same covid PBMC sample (`LUM-SJ-x213-g`, 1,203,332,091 reads),
comparing three ways of removing host signal before viral quantification:

| | no filter | `--genome-dlist` | `--host-filter starsolo` |
|---|---:|---:|---:|
| reads reaching viral quant | 1,203,332,091 | 1,203,332,091 | 120,405,696 |
| `p_pseudoaligned` | 6.4 % | 4.5 % | 5.8 % |
| **`p_unique`** | **2.1 %** | **0.6 %** | **4.4 %** |
| quantified molecules | 21,613,840 | 8,404,326 | 5,308,302 |
| called cells | 143,243 | 28,922 | 28,921 |

D-list masking left most of the anellovirus artifact in place (the published covid *Alphatorquevirus* signal was an artefact of poly-G no-signal reads (≈90 %) and host-homologous reads (≈10 %) (F-005, F-019); low-level divergent anellovirus is not excluded) while cutting
`p_unique` **3.5x** (2.1 % → 0.6 %), i.e. it destroyed uniquely-placed molecules
to achieve very little. The reason is mechanical: a kallisto D-list masks shared
k-mers by *exact match*, so it cannot see host sequence that has diverged even
slightly. Divergent host sequence is exactly the case that produces the artifact
in the first place.

In the table above, `--host-filter starsolo` gave the highest `p_unique`
(4.4 %) and `--genome-dlist` the lowest (0.6 %). The table does not include a
combined host+virus *genome* reference arm, so this page does not rank that
option against the other two; no ordering beyond the measured three is claimed.
Prefer `--host-filter starsolo` over `--genome-dlist` for host subtraction, and
use `--genome-dlist` only when the host genome FASTA is already on hand and its
host-homology measurements are the actual goal.
With `--genome-dlist`, `host_homology_annotations.tsv` retains maximum identity,
query coverage, aligned bases, and best host target for every viral sequence;
the genome path and SHA-256 are frozen in the manifest.

Example:

```bash
viralscan build-ref \
  --host human \
  --virus-accessions NC_045512.2 NC_002021.3 \
  --output ref_human_virus/ \
  --ncbi-email you@example.org
```

Expected outputs when `kb ref` succeeds:

| File | Use |
|------|-----|
| `combined.fa` | Concatenated host + virus FASTA |
| `combined.gtf` | Concatenated host + virus GTF |
| `index.idx` | Pass to `viralscan run -i` |
| `t2g.txt` | Pass to `viralscan run -t` |
| `reference_manifest.json` | Frozen profile and per-sequence provenance/digests |
| `cdna.fa` | cDNA FASTA produced by `kb ref -f1` |

### Supported host species

Run `viralscan build-ref --list-species` for the current list.
Examples: `human`, `mouse`, `rat`, `zebrafish`, `chicken`, `macaque`, `pig`.

---

## `viralscan data fetch` — bundled viral panel

```
viralscan data fetch [OPTIONS]
```

Download the bundled viral GTF panel from Zenodo, verify checksums, and unpack
the GTF files into the local ViralScan cache. A run needs the panel only for an index
without a `<index>.build_manifest.json` (such as one not built by `viralscan build-ref`);
with a manifest the viral gene list comes from the manifest.

| Flag | Default | Description |
|------|---------|-------------|
| `--cache-dir PATH` | `~/.cache/viralscan/` | Cache root. GTF files are written to `PATH/data/`. |
| `--url URL` | Zenodo archive URL | Override archive URL, mainly for tests or mirrors. |
| `--sha256 DIGEST` | *(none)* | Optional expected SHA-256 digest for the downloaded archive. |
| `--force` | off | Re-download and replace cached GTF files. |
| `--verbose` | off | DEBUG-level logging. |
| `--quiet` | off | Warnings + errors only. |

Typical use:

```bash
viralscan data fetch
```

---

## `viralscan evidence` — trace reads behind viral calls

```
viralscan evidence [OPTIONS]
```

Trace the reads whose corrected (barcode, UMI) molecule was assigned to viral
genes in a completed ViralScan run. With exact host and target FASTA inputs,
ViralScan aligns reads competitively, writes IGV assets, and reports coverage,
identity, duplication, complexity, and host-competition diagnostics. These
outputs strengthen or weaken candidate evidence; they do not by themselves
confirm infection.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--run-dir PATH` | | *(required)* | A completed ViralScan run output directory |
| `--output PATH` | `-o` | *(required)* | Directory for evidence outputs |
| `--viral-fasta PATH` | | *(none)* | Exact target-virus FASTA; enables competitive BAM/coverage/QC and requires `--host-fasta` |
| `--host-fasta PATH` | | *(none)* | Full host-genome FASTA required with `--viral-fasta` |
| `--virus STRING` | | *(required)* | Exact accession/gene, registered alias, or canonical detected call; substring matching is forbidden |
| `--blast` | | off | Competitively BLAST a deterministic read sample against host plus target (requires `blast+`) |
| `--dedup MODE` | | `umi` | `umi`, `markdup`, or `none`; raw and selected deduplicated BAMs remain separate. `umi` keeps one deterministic representative per corrected (CB, UMI, reference) molecule regardless of start, strand or CIGAR (rank: MAPQ, AS, NM, aligned query span, then lexical); `markdup` is positional alignment deduplication, not molecule evidence. Primary mapped alignments only (FLAG 0xF04 excluded in coverage, depth, QC and start tallies) |
| `--read-start-profile` | | off | Write a per-position 5-prime read-start profile |
| `--cell-tags` | | off | Write an indexed CB/UB-tagged BAM and add it to the IGV session |
| `--sampling-seed N` | | `42` | Seed for order-independent deterministic BLAST sampling |
| `--cores N` | `-c` | `4` | Threads for minimap2/samtools/blast |
| `--verbose` | | off | Enable DEBUG-level logging |
| `--quiet` | | off | Suppress INFO messages |

Example:

```bash
viralscan evidence \
  --run-dir output/sample/ \
  --viral-fasta EBV.fa \
  --host-fasta GRCh38.fa \
  --virus EBV \
  --blast --read-start-profile --cell-tags \
  -o output/sample/evidence/
```

---

## `viralscan evidence-plot` — molecule-level host versus virus alignment

Draws `competitive_alignment.png` and `.pdf` for each `viralscan evidence` directory (the evidence files are
only read). Counts default to one representative read per corrected cell barcode and UMI (the read with the best viral
alignment), selected **after kb alignment and competitive alignment**, retaining that read's primary, secondary and
supplementary hits.
Panel A shows molecule lanes (virus row, host row, clipped bases, poly-A/T runs, 10x TSO, verdict and
both alignment scores), panel B the virus-versus-host score scatter, panel C the human loci hit (1 kb bins).
Panel D is a histogram of viral aligned length (M/=/X, no clips) per molecule with the 60-base Q137 gate marked. The lane block (panel A) is drawn only up to `--max-lanes` molecules (default 60) and is left out above that. `--layer raw` shows
explicitly labelled raw-read diagnostics. `viralscan evidence`
writes the figure itself when the competitive BAM exists. With `--cross-call PREFIX`, the alignments of all given
directories are pooled into `PREFIX.tsv` and `PREFIX.png` (host loci in 1 kb bins, viral positions in 100 bp
bins, with distinct molecules, calls and libraries per locus, recurrent loci first). The cross-call TSV adds
`count_layer` and `n_molecules`; its legacy `n_reads` column counts representative reads in the UMI layer.

`viralscan evidence` also writes `viral_reads.umi_dedup.fasta` for BLAST sampling, `molecule_support.tsv` and
`read_support_summary.tsv` (Q137 provisional flags per molecule and per call; columns in `docs/output_reference.md`;
loci in `src/viralscan/data/read_support_loci.tsv`; `python -m viralscan.read_support <evidence dir> -o <dir>` recomputes
them read-only). Flags remove no reads and change no detection threshold. Raw BAMs and raw-read verdict
tables remain diagnostic artifacts. DSR postprocessing writes separate `molecule_read_verdicts.tsv.gz` and
`molecule_read_verdict_summary.tsv`; call-level complexity and genome structure use UMI-deduplicated QC.
Legacy evidence with only raw QC must be regenerated before it can support those checks.

```bash
viralscan evidence-plot --evidence-dir ev/callA ev/callB -o plots/ --cross-call plots/cross_call_loci
```

---

## `viralscan rerun-multimap` — swap multimapping method

```
viralscan rerun-multimap [OPTIONS]
```

Re-run the multimapping step with a different algorithm for an existing run,
skipping the expensive pseudoalignment (`kb count`). The source tree is copied
to the required new output directory and is never modified. For `equal`,
`host-conservative`, and `unique-weighted`, a schema-valid v3 H5AD can use its
pre-stored allocation layers without BUS reprocessing. `em-global`, `em-cell`,
and older/incomplete H5AD files reprocess the retained BUS input.

**Development limitation:** the current command refreshes multimapping,
detection, and UMAP workflow checkpoints, but complete invalidation and
regeneration of every downstream method-dependent artifact is not yet release-
gated. In particular, copied host-response outputs may be stale and must not be
used without an explicit rerun and audit. This remains tracked as `SW-04` and
`SW-05` in `PLAN.md`.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--run-dir PATH` | | *(required)* | Completed source run; it is never modified. |
| `--output PATH` | `-o` | *(required)* | New result directory; it must not already be non-empty. |
| `--multimap-method METHOD` | | *(required)* | Multimapping resolution method: `conservative`, `em`, or `equal` |
| `--cores N` | `-c` | `6` | Total CPU budget; native workers (kb, bustools, STAR, ...) are bounded to it |
| `--multimap-em-max-iter N` | | *(preserve)* | (EM methods) Maximum iterations; omit to keep the existing config value |
| `--multimap-em-tol TOL` | | *(preserve)* | (EM methods) Convergence tolerance; omit to keep the existing config value |
| `--verbose` | | off | Enable DEBUG-level logging |
| `--quiet` | | off | Suppress INFO messages |

Examples:

```bash
viralscan rerun-multimap --run-dir out/ -o out_hc/ --multimap-method conservative
viralscan rerun-multimap --run-dir out/ -o out_em/ --multimap-method em --cores 8
```

**Tip:** the default (`conservative`) is the safe choice; use
`equal` for an equal-allocation comparison or `em` for iterated cell-local
allocation. Always validate the new result tree before analysis.

---

## `viralscan hostresponse` — host-response analysis

```
viralscan hostresponse [OPTIONS]
```

Associate candidate viral molecule support with host gene expression via logistic regression
(Luebbert et al. 2025 approach) on a completed viralscan run. For each
detected virus, ViralScan trains an L2 logistic regression model predicting
candidate-positive vs. candidate-negative labels from host gene expression, then runs
randomized Lasso selection frequency to rank associated host genes (uncalibrated: no error
rate is controlled, and with many positives noise genes pass the threshold).
Optional gget pathway enrichment (`--enrichment`) identifies enriched biological
processes in the above-threshold gene set.

The viralscan output directory must already contain `config.yaml` and a
completed multimap step (`log/multimap.done`). Results are written to
`<output>/hostresponse/`.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--output PATH` | `-o` | *(required)* | Existing viralscan sample output directory (contains `config.yaml`) |
| `--host-h5ad PATH` | | *(required)* | Host gene-expression h5ad (cells × genes, matched barcodes) |
| `--n-seeds N` | | from config or `6` | Random seeds for multi-seed L2 regression |
| `--n-stab-iter N` | | from config or `100` | Stability-selection iterations |
| `--use-hvg` / `--no-use-hvg` | | on | Use highly variable genes as features; `--no-use-hvg` uses all genes |
| `--stab-min-prob P` | | from config or `0.6` | Min selection probability to call a gene stably associated |
| `--top-n-genes N` | | from config or `50` | Top N stable genes to pass to pathway enrichment |
| `--min-umi N` | | from config or `10` | Min selected-method viral molecule estimate for the per-cell raw candidate-support label |
| `--label {raw,cpm,fraction}` | | `raw` | Candidate-support label (see **Depth-confound controls** below) |
| `--depth-match` | | off | Restrict each virus to a depth-matched cohort to reduce measured depth imbalance; residual confounding can remain |
| `--mito-control` / `--no-mito-control` | | on | Add per-cell %mito as a covariate to the per-gene E-values |
| `--gene-symbols` | | off | Annotate output CSVs with HGNC symbols from Ensembl IDs via mygene.info (network) |
| `--differential` | | off | Write a genome-wide depth/%mito-adjusted differential table (`<virus>_differential.csv`) |
| `--enrichment` | | off | Run pathway enrichment via gget (requires `pip install "viralscan[enrichment]"`) |
| `--enrichment-db DB` | | `GO_Biological_Process_2023` | gget.enrichr database |
| `--verbose` | | off | Enable DEBUG-level logging |
| `--quiet` | | off | Suppress INFO messages |

### Cell alignment (host cells are the universe)

The host h5ad defines the analysed cells. A two-step run's viral matrix holds only
barcodes with non-host reads, so a host cell missing from it was observed with host
reads and has **zero viral reads**: it is zero-filled and counts as a negative. A
viral-matrix barcode that is not in the host h5ad has no host observation: it is
**dropped** (unavailable, never a zero). `hostresponse_manifest.json` records both as
`barcode_alignment` (`n_viral_zero_filled`/`frac_viral_zero_filled` over host cells,
`n_viral_dropped_not_host_called`/`frac_viral_dropped_not_host_called` over viral
barcodes). The run errors, rather than guessing a universe, when the host h5ad has no
cells or when no viral barcode matches any host barcode (typically mismatched barcode
strings, which would otherwise zero-fill every cell). A grouping or cell-type column
that exists only in the viral matrix is missing for zero-filled cells; put it in the host
h5ad obs.

### Depth-confound controls (important)

The default candidate-support label (`--label raw`, `counts >= min-umi`)
**tracks sequencing depth**: deeper cells carry more viral *and* more host counts,
so a naive host-gene AUC can partly reflect library size rather than biology.
Every run reports a depth-confound baseline. Two opt-in controls change the
label or cohort to reduce measured depth imbalance, but neither proves that
technical or biological confounding has been eliminated:

- **Always reported** (no flag needed): `hostresponse_metrics.csv` includes
  `depth_alone_auc_mean` (the AUC from sequencing depth alone, under the identical
  split) next to the model AUC, and `<virus>_depth_diagnostics.csv` gives a
  depth-adjusted odds ratio + [E-value](https://doi.org/10.7326/M16-2607) per stable
  gene. If `depth_alone_auc` ≈ the model AUC, treat the headline as depth-confounded.
- **`--label cpm`** (or `fraction`): a depth-normalized, prevalence-matched label —
  candidate-positive cells have the highest viral molecule estimate per total
  molecule estimate, keeping the same prevalence for comparison. The ratio can
  still correlate with depth or other technical factors.
- **`--depth-match`**: builds a coarsened-exact depth-matched case/control cohort so
  the observed depth distributions are closer. Inspect balance diagnostics and
  the depth-only baseline; residual and unmeasured confounding can remain.

Examples:

```bash
# Minimal — run with defaults from config (raw label; still reports the depth baseline)
viralscan hostresponse -o output/sample/ --host-h5ad host_genes.h5ad

# Depth-normalized label + %mito control + gene symbols
viralscan hostresponse \
  -o output/sample/ \
  --host-h5ad host_genes.h5ad \
  --label cpm --gene-symbols

# Depth-matched cohort + genome-wide differential table
viralscan hostresponse \
  -o output/sample/ \
  --host-h5ad host_genes.h5ad \
  --depth-match --differential

# With stability tuning and pathway enrichment
viralscan hostresponse \
  -o output/sample/ \
  --host-h5ad host_genes.h5ad \
  --n-stab-iter 200 \
  --enrichment
```

> **Pathway enrichment**: requires the optional `gget` dependency.
> Install it with `pip install "viralscan[enrichment]"` before using `--enrichment`.

---

## `viralscan check-whitelist` — chemistry/whitelist preflight

```
viralscan check-whitelist --r1 R1.fastq.gz -w WHITELIST -x TECHNOLOGY
```

Sample the first reads of R1, extract the cell barcode using the technology's
barcode geometry, and report the fraction that match the whitelist. A low match
rate means `--technology` / `--whitelist` do **not** match the library chemistry —
in which case `bustools` silently discards most reads and `kb count` produces an
all-empty-droplet matrix with no error. (In one real case a GEM-X 5′ library
mislabeled as `10xv3` lost 96.5% of its reads before this was caught.)

Run this before a long quantification job when you are unsure of the chemistry.
The main `viralscan` run also performs this check automatically as a **warning**
whenever an explicit `--whitelist` is supplied.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--r1 R1` | | *(required)* | R1 FASTQ (the barcode read; optionally `.gz`) |
| `--whitelist PATH` | `-w` | *(required)* | Barcode whitelist, one barcode per line (optionally `.gz`) |
| `--technology X` | `-x` | `10xv3` | Single-cell technology (sets the barcode length/offset) |
| `--min-match-rate F` | | `0.5` | Minimum acceptable barcode match rate |
| `--n-sample N` | | `100000` | Number of R1 reads to sample |

Exits `0` when the match rate is acceptable and `1` on a likely mismatch, so it
can gate a pipeline. Example:

```bash
viralscan check-whitelist \
  --r1 sample_R1.fastq.gz \
  -w 3M-february-2018.txt \
  -x 10xv3
```

---

## `viralscan check-chemistry` — library chemistry and strand diagnosis

```
viralscan check-chemistry --r1 R1.fastq.gz [--r2 R2.fastq.gz -i index.idx -t t2g.txt] [-x TECH] [-w LIST] [--json PATH]
```

Read-only check before a run. From the first R1 reads it reports bulk vs single-cell, the
on-list/chemistry call and UMI length; with `--r2`, `-i` and `-t` it also runs the 1M-pair strand
pilot, which tells 3' from 5' and picks `--strand`. It prints the `-x`/`-w`/`--strand` to use and
exits 1 when the chemistry is unresolved or the evidence conflicts. Flags: see the generated
table below.

---

## `viralscan doctor` — dependency/profile check

```
viralscan doctor [--profile {pip,full}] [--json PATH] [--verbose | --quiet]
```

Check whether the current installation provides the dependencies required by a
declared installation tier. This is a diagnostic command; a green pip profile
does not imply that native full-workflow tools are installed.

| Flag | Default | Description |
|------|---------|-------------|
| `--profile {pip,full}` | `full` | Check the Python/reporting/validation tier or the complete external-tool workflow tier |
| `--json PATH` | *(none)* | Also write the report as JSON (it is always printed) |

Example:

```bash
viralscan doctor --profile full --json doctor_report.json
```

---

## `viralscan validate-run` — schema and count-invariant check

```
viralscan validate-run --run-dir RUN_DIR [--no-verify-inputs] [--json PATH]
```

Validate a v3 output tree against its run manifest, schemas, fingerprints, and
available count invariants. Validation reports defects; it does not convert or
reinterpret a pre-v3 H5AD file.

| Argument/flag | Default | Description |
|---------------|---------|-------------|
| `--run-dir RUN_DIR` | *(required)* | Output directory containing `run_manifest.json` |
| `--no-verify-inputs` | off | Skip re-hashing manifested input files; use only when those inputs are intentionally unavailable |
| `--json PATH` | *(none)* | Write the validation report as JSON |

Example:

```bash
viralscan validate-run --run-dir output/ --json output/validation_report.json
```

Missing required schemas, incompatible schema versions, fingerprint failures,
or broken count invariants make validation fail rather than silently downgrade.

## `viralscan rerun-programs` — re-score gene programmes

Run layer 2 over a completed run, **in place**:

```bash
viralscan rerun-programs --run-dir output/sample/ --programme-min-breadth 2
```

Unlike `rerun-multimap` this does not copy the run and does not need a separate
`--output`. Layer 2 only reads the count matrices and layer 1's
`viral_summary.tsv` and adds two files, so there is no reason for the two
directories to be able to disagree. It requires `log/detection.done` in the run
directory; without it the command refuses rather than emitting empty tables.

---

## Complete flag index (generated)

These tables are generated from the argument parser by
`scripts/gen_cli_reference.py` (`PYTHONPATH=src python scripts/gen_cli_reference.py`;
`--check` fails if they are stale). They list every flag with its parser default and help text;
the sections above add context. Do not edit between the markers.

### viralscan run

<!-- BEGIN GENERATED: run -->

| Flag | Default | Help |
|------|---------|------|
| `--output, -o OUTPUT` | *required* | The path to the output directory. |
| `--r1 R1` | *required* | The path to the forward FASTQ sample (gunzipped is preferred). |
| `--r2 R2` | *required* | The path to the backward FASTQ sample (gunzipped is preferred). |
| `--transcripts, -t TRANSCRIPTS` | *(none)* | The path to the transcripts (t2g) file produced by kb ref. |
| `--index, -i INDEX` | *(none)* | The path to the reference index created by kb ref. |
| `--cores, -c CORES` | `6` | The amount of cores the workflow can use. Default: 6. |
| `--reference, --no-reference` | `False` | Build a kb ref index from --fasta and --gtf into the output directory. |
| `--gtf GTF` | *(none)* | Path to GTF files (comma-delimited, without space in-between). |
| `--fasta FASTA` | *(none)* | Path to FASTA files (comma-delimited, without space in-between). |
| `--f1 F1` | *(none)* | Path to the cDNA FASTA (lamanno, nucleus) or mismatch FASTA (kite) to be generated |
| `--visual, --no-visual` | `True` | Add visualizations to the output. Use --no-visual to disable. Default: True. |
| `--technology, -x TECHNOLOGY` | *(none)* | Single-cell technology (`kb --list` to view). Default: detected from the first 100k R1 reads of each sample (on-list match, TSO/poly-T position). An explicit -x that the reads contradict, or reads that fit no single chemistry, stop the run. GEM-X 5' needs its on-list via -w. |
| `--force-technology` | `False` | Run with the explicit -x even when the chemistry check disagrees or cannot decide. The detection is still logged. |
| `--require-chemistry-sanity` | `False` | Stop after kb count if library sanity reports an error or cannot be checked. Diagnostics are retained; warnings alone do not stop the run. |
| `--whitelist, -w WHITELIST` | *(none)* | Path to file of whitelisted barcodes. If absent, kb-python's bundled whitelist is used. |
| `--strand STRAND` | *(none)* | Read strandedness passed to `kb count --strand`. Default (no flag): with the pinned kallisto 0.52.0 the effective behaviour is `unstranded` for every chemistry — kallisto 0.52.0 marks the 10x technologies strand-specific yet leaves the strand unset, so the legacy per-technology forward default no longer engages (upstream regression; kallisto <=0.51.1 defaulted 10x to forward). Set the flag explicitly for stranded quantification. 10x 5' libraries need `reverse` or `unstranded`: forward pseudoaligns only 6.5-8.9% of reads (F-020). `auto` (opt-in) pilots forward/reverse/unstranded on the first 1M read pairs of each sample and picks one; the choice is recorded in run_manifest.json. |
| `--multimapping, --no-multimapping` | `True` | Take multimapping into account. Use --no-multimapping to disable. Default: True. |
| `--umap` | `False` | Generate a UMAP plot. Significantly increases runtime. Default: off. |
| `--ncbi-accession NCBI_ACCESSION` | *(none)* | One or more NCBI nucleotide accessions (e.g. 'NC_002021.3'), comma-separated. ViralScan will download FASTA + GTF for each and build the index. Mutually exclusive with --reference / --fasta / --gtf. |
| `--ncbi-email NCBI_EMAIL` | *(none)* | Contact email for NCBI E-utilities. Falls back to $NCBI_EMAIL. |
| `--ncbi-api-key NCBI_API_KEY` | *(none)* | NCBI API key for higher E-utilities request rates with --ncbi-accession. Falls back to $NCBI_API_KEY. |
| `--cache-dir PATH` | *(none)* | Root cache directory for ViralScan's fetched viral annotation panel. Default: $VIRALSCAN_CACHE or ~/.cache/viralscan/. |
| `--se-threshold SE_THRESHOLD` | `10` | Selected-method viral molecule estimate above which a cell is flagged as a 'super-expressor'. Default: 10. |
| `--detection-threshold DETECTION_THRESHOLD` | `1` | Minimum selected-method viral molecule estimate, summed over the called cells, required for a virus to pass the called-cell detection gate (see --detection-min-cells). Molecules in non-called barcodes (empty droplets) do not count and are reported as background_umis. Default: 1. |
| `--detection-min-cells N` | *(none)* | Minimum number of distinct called cells that must be infected (see --infected-min-molecules) for a virus to pass the detection gate. A virus seen only in the all-barcode screen is listed as a candidate, never dropped. With --cell-calling none every barcode counts as a called cell. Default: 1. |
| `--infected-min-molecules X` | *(none)* | A cell is infected for a virus when its summed allocated molecules across that virus's genes reach X. One rule for detection, cell-type enrichment and the UMAP. Cells with 0 < weight < X are counted as sub-threshold (infected_subthreshold_cells), not dropped. Fractional EM or multimapper weights make a bare 'greater than 0' rule too permissive. Default: 1.0. |
| `--positive-control-gene GENE_ID` | *(none)* | Gene ID of a spike-in planted at a known molecule count. It is the only way to measure the k-mer capture term, and therefore the only way to turn 'no virus detected' into a certifiable negative rather than a sampling statement. Must be given together with --positive-control-molecules. |
| `--positive-control-molecules N` | *(none)* | Molecules of --positive-control-gene planted in the library. Capture is measured as observed/N and reported in results/positive_control.json; capture=1.0 means no loss was measurable and implies nothing about sequence divergence. |
| `--positive-control-scope POSITIVE_CONTROL_SCOPE` | *(none)* | What the positive control may certify. exact_sequence: only the virus row named by --positive-control-virus-key, and only that exact sequence. panel_mechanics: the pipeline recovers a planted molecule; certifies no virus. virus_key: needs an approved transfer calibration and is currently rejected. Unset on a configured control behaves as panel_mechanics and warns once. |
| `--positive-control-virus-key VIRUS_KEY` | *(none)* | Virus row (the name in results/sensitivity.tsv) an exact_sequence control was measured on. Required with --positive-control-scope exact_sequence. |
| `--positive-control-receipt JSON` | *(none)* | Pinned control identity/measurement receipt required for exact-sequence certification. |
| `--require-positive-control, --no-require-positive-control` | `False` | Fail the run when nothing is detected and no positive control could measure a capture term. A mechanics/recovery gate: passing it does not certify any virus outside the control's declared scope. Default: False. |
| `--anellovirus-gene-ids, --no-anellovirus-gene-ids` | `True` | Treat the expanded anellovirus panel's {accession}_geneN IDs as viral. Required for any reference built with `viralscan build-ref --reference-panel anellovirus` (or the bundled-panel builder), because those GTFs are materialized into the index rather than the panel directory. Off means 2,021 of 2,041 anellovirus genomes are counted but never reported. Default: True. |
| `--gene-programs, --no-gene-programs` | `False` | Second layer: for viruses the detection rule already called, infer the viral gene programme (latent vs productive) per cell from uniquely-placing molecules, and write results/gene_program_summary.tsv and results/gene_program_cells.tsv. Only nine viruses have a programme model; for the rest a 'not_applicable' row is emitted so silence is not read as 'programme not detected'. Off by default. Default: False. |
| `--anello-align, --no-anello-align` | `False` | Anellovirus alignment branch (PLAN ANDET-09): align every host-unmapped read to the panel's anellovirus genomes with STARsolo, independent of kallisto, and add alignment_* evidence columns plus detection_source to viral_summary.tsv. Runs only with --host-filter starsolo and an anello_star/ index next to the kb index (built by scripts/build_bundled_panel_ref.py); otherwise alignment_status records why it was skipped. The columns are labels, never filters. Default: False. |
| `--auto-evidence, --no-auto-evidence` | *(none)* | After detection, run `viralscan evidence` with host confirmation (minimap2 against host + viral) for every detected Anelloviridae genus, writing results/evidence/<genus>/. Needs --host-fasta, a viral-only FASTA (--viral-fasta, or viral.fa next to the index, as written by scripts/build_bundled_panel_ref.py), minimap2 and samtools, and multimapping. Re-reads the FASTQs once per genus. Opt-in, and stays so until ANELLO-PRIOR.3 closes. Default: off. |
| `--host-capture, --no-host-capture` | *(none)* | Split each virus's reads by host compatibility: replay the quantified FASTQs through `kallisto bus -n`, capture the reads that touch a viral transcript, then `bustools capture --complement` the host transcripts. Writes results/host_capture.tsv with reads_virus_only (no host transcript in the read's equivalence class), reads_also_host and reads_total per virus. Needs an index that holds host transcripts (a combined host + virus index); on a virus-only index every count is 'not_applicable', never 0. Counts reads, not molecules. Needs kallisto and bustools. Opt-in. Default: off. |
| `--host-fasta FASTA` | *(none)* | Full host-genome FASTA for --auto-evidence's competitive alignment. |
| `--viral-fasta FASTA` | *(none)* | Viral-only FASTA for --auto-evidence. Default: viral.fa next to the kb index. |
| `--programme-min-breadth N` | `2` | Distinct non-overlapping overlap groups required before a programme is called. Counted in overlap groups rather than genes because EBV's latent and lytic ORFs share exonic sequence: on the EBV LCL run a naive per-gene comparison gives a latent:lytic ratio of 1.15 in a cell line defined by latency. Must be >= 1. Default: 2. |
| `--programme-latent-min-breadth N` | `1` | Distinct latent overlap groups required before 'latent' is called. Must be >= 1. Default: 1. |
| `--programme-min-umi X` | `0.0` | A marker counts as evidence only at >= X molecules; 0 means any nonzero value. Must be >= 0. Default: 0.0. |
| `--min-counts MIN_COUNTS` | `1000` | Minimum total molecule count per cell (for UMAP QC). Default: 1000. |
| `--min-genes MIN_GENES` | `200` | Minimum detected genes per cell (for UMAP QC filter). Default: 200. |
| `--hvg-min-mean HVG_MIN_MEAN` | `0.0125` | Scanpy highly-variable-gene min_mean parameter. Default: 0.0125. |
| `--hvg-max-mean HVG_MAX_MEAN` | `3.0` | Scanpy highly-variable-gene max_mean parameter. Default: 3.0. |
| `--hvg-min-disp HVG_MIN_DISP` | `0.5` | Scanpy highly-variable-gene min_disp parameter. Default: 0.5. |
| `--umap-n-neighbors UMAP_N_NEIGHBORS` | `15` | Number of neighbors for Scanpy graph construction before UMAP. Default: 15. |
| `--cell-calling CELL_CALLING` | `auto` | How to identify real (non-empty-droplet) cells so viral rates are reported over called cells (primary) as well as all barcodes (secondary). 'auto' uses --called-cells-file when supplied, otherwise EmptyDrops; 'external' requires that called-cell list; 'emptydrops' runs DropletUtils::emptyDrops (needs R) on the host transcriptome only, never on viral reads; 'knee' is a sensitivity-only barcode-rank approximation, not for reported numbers (it calls empty droplets as cells on real libraries); 'none' = all barcodes. Default: auto. |
| `--called-cells-file CALLED_CELLS_FILE` | *(none)* | Path to an external called-cell barcode list (one per line, optional -1 suffix) used when --cell-calling external. Typically a CellRanger/STARsolo filtered barcodes.tsv(.gz). |
| `--emptydrops-seed EMPTYDROPS_SEED` | `100` | Random seed for DropletUtils::emptyDrops. It is a Monte-Carlo test, so this decides which borderline barcodes are called. Set it from the preregistered seed when running under a frozen protocol. Default: 100. |
| `--emptydrops-niters EMPTYDROPS_NITERS` | `10000` | Monte-Carlo iterations for DropletUtils::emptyDrops. Default: 10000. |
| `--multimap-method {conservative,em,equal}` | `conservative` | How to allocate multi-gene EC counts. 'conservative' excludes host-virus ambiguous EC mass from viral genes (recommended for combined host+virus references); 'em' iterates the allocation using cell-local evidence; 'equal' splits reads equally (fast, unbiased first pass). Use 'viralscan rerun-multimap' to switch methods after the run without redoing the pseudoalignment. Default: conservative. |
| `--multimap-pseudocount MULTIMAP_PSEUDOCOUNT` | `1.0` | Positive pseudocount used by --multimap-method unique-weighted. Default: 1.0. |
| `--multimap-primary-call MULTIMAP_PRIMARY_CALL` | `selected-method` | V3 detection-count contract. The only supported value is 'selected-method': summaries use the complete molecule matrix for --multimap-method, while unique and ambiguous evidence remain separately reported. Default: selected-method. |
| `--multimap-em-max-iter MULTIMAP_EM_MAX_ITER` | `100` | Maximum EM iterations for --multimap-method em. Default: 100. |
| `--multimap-em-tol MULTIMAP_EM_TOL` | `1e-06` | EM convergence tolerance for --multimap-method em. Default: 1e-06. |
| `--multimap-molecule-assignments, --no-multimap-molecule-assignments` | `False` | Write results/molecule_assignments.tsv.gz: one row per corrected CB-UMI molecule with its EC ids, status (unique/ambiguous/unresolved), compatible genes and the selected --multimap-method's per-gene weights. Optional, off by default, and it changes no count. Needs one extra pass over the BUS text; the file is large on deep samples. Default: False. |
| `--cell-type-labels CELL_TYPES` | *(none)* | Path to a CSV (barcode,cell_type) providing cell-type labels for per-type viral enrichment in the HTML report. Optional. |
| `--host-filter ALIGNER` | *(none)* | Optional advanced host-subtraction pre-step before viral quantification. Removes reads that align to the host genome, reducing false positives. V3 supports 'starsolo' (full-genome STAR alignment) because it preserves exact fragment identity. Requires --host-index. Usually not needed when using a combined host+virus reference. |
| `--read-filter READ_FILTER` | *(none)* | Read-artefact filter before kb count (after --host-filter, if any). 'off' (the default) counts every read. 'artefact' drops pairs whose R1 carries the 10x TSO in the barcode/UMI span, whose R2 is reagent (TSO\|poly-T, TruSeq chimera), or whose R2 has no complex body before its first >=15 nt homopolymer in either orientation. It also drops genuine low-complexity viral reads. 'tso-trim' keeps every pair but cuts the 10x TSO (30 nt) from the start of R2, because TSO-led host reads that STAR's host filter misses pseudoalign to viral CAG tracts (F-028). Writes read_filtered/read_filter_audit.tsv and fragment_lineage.tsv.gz. |
| `--host-filter-star-params HOST_FILTER_STAR_PARAMS` | *(none)* | STAR filter parameter set for --host-filter starsolo. 'pinned' (the default) removes only near-identical, near-full-length host alignments (at most 4 mismatches, 90% of the read aligned, up to 20 loci). 'star-default' uses STAR's own defaults (10 mismatches, 66% aligned, up to 10 loci). Both are passed explicitly and recorded in host_filter_audit.tsv. |
| `--host-index PATH` | *(none)* | Path to the STAR host-genome directory required by --host-filter, built with STAR --runMode genomeGenerate. |
| `--host-h5ad PATH` | *(none)* | Path to a host gene-expression h5ad file (cells × host genes, log-normalised or raw). When provided, ViralScan trains per-virus logistic regression models predicting virus presence from host gene expression (Luebbert et al. 2025 approach) and writes results to <output>/hostresponse/. |
| `--hostresponse-cv HOSTRESPONSE_CV` | *(none)* | CV scheme (default: auto). auto = group-disjoint CV on --groups, or on the one donor/sample-like obs column if there is exactly one, else cell-level CV with a leakage warning (also if group CV cannot run on that auto-selected column); cell = legacy cell-level CV (a donor can span train and test); group = group-disjoint CV, requires --groups, fails closed. |
| `--hostresponse-groups HOSTRESPONSE_GROUPS` | *(none)* | obs column declaring donor/sample groups; used by --cv group and --cv auto (ignored, with a warning, under --cv cell). |
| `--hostresponse-cv-folds HOSTRESPONSE_CV_FOLDS` | *(none)* | Number of group-disjoint folds (default: 5). |
| `--hostresponse-cell-type-column HOSTRESPONSE_CELL_TYPE_COLUMN` | *(none)* | obs column used for cell-type baselines and strata. |
| `--hostresponse-cell-types HOSTRESPONSE_CELL_TYPES` | *(none)* | Declared cell types to analyze separately (space- or comma-separated). |
| `--hostresponse-panel-in-fold, --no-hostresponse-panel-in-fold` | *(none)* | Fit the selection-frequency panel within each training fold. |
| `--hostresponse-permutations HOSTRESPONSE_PERMUTATIONS` | *(none)* | Structured-null replicates (default: 0). |
| `--hostresponse-permutation-unit HOSTRESPONSE_PERMUTATION_UNIT` | *(none)* | Explicit exchangeability unit for the null. |
| `--hostresponse-permutation-block HOSTRESPONSE_PERMUTATION_BLOCK` | *(none)* | obs column restricting label exchangeability. |
| `--hostresponse-min-negative-cells HOSTRESPONSE_MIN_NEGATIVE_CELLS` | *(none)* | Minimum negative cells per stratum (default: 10). |
| `--hostresponse-min-groups HOSTRESPONSE_MIN_GROUPS` | *(none)* | Minimum groups per stratum (default: 2). |
| `--hostresponse-n-seeds N` | *(none)* | Number of random seeds for the multi-seed L2 logistic regression, 1-6 (default: 6). |
| `--hostresponse-n-stab-iter N` | *(none)* | Iterations for randomized Lasso selection frequency (uncalibrated; default: 100). |
| `--hostresponse-min-umi N` | *(none)* | Per-cell positive label for host-response: a cell is virus-positive when its viral molecule estimate is >= N. Separate from --detection-threshold, which is the sample-level gate (default: 10). |
| `--hostresponse-use-hvg, --no-hostresponse-use-hvg` | `True` | Use highly variable genes as features (default: on). --no-hostresponse-use-hvg uses all genes. |
| `--hostresponse-stab-min-prob PROB` | *(none)* | Minimum uncalibrated selection frequency to flag a gene; a ranking threshold, not an error rate (default: 0.6). |
| `--hostresponse-top-n-genes N` | *(none)* | Top N above-threshold genes to pass to pathway enrichment (default: 50). |
| `--hostresponse-label HOSTRESPONSE_LABEL` | *(none)* | Virus-support labeling strategy for hostresponse (default: raw = selected-method molecule estimate >= detection_threshold). 'cpm'/'fraction' normalize by measured host molecule depth but do not remove all depth confounding. |
| `--hostresponse-depth-match` | `False` | Restrict hostresponse to a depth-matched cohort (removes depth as a design-level confounder). Recommended when depth_alone_auc is close to model_auc. |
| `--hostresponse-mito-control, --no-hostresponse-mito-control` | `True` | Include %mito as a covariate in hostresponse depth-adjusted E-values (default: on). Disable with --no-hostresponse-mito-control if the host h5ad has no mitochondrial genes. |
| `--hostresponse-differential` | `False` | Run a genome-wide depth-and-mito-adjusted differential expression test alongside the selection-frequency model (written to <output>/hostresponse/<virus>_differential.csv). |
| `--enrichment` | `False` | Run pathway enrichment on above-threshold host genes via gget.enrichr (requires gget; install with pip install 'ViralScan[enrichment]'). |
| `--enrichment-db DB` | *(none)* | Enrichment database for gget.enrichr (default: GO_Biological_Process_2023). |
| `--resume` | `False` | Resume only when run_manifest.json exactly matches this invocation. |
| `--overwrite` | `False` | Explicitly replace a non-empty output directory after confirmation. |
| `--yes, -y` | `False` | Answer the --overwrite confirmation; does not imply overwrite or resume. |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages; only show warnings and errors. |

<!-- END GENERATED -->

### viralscan data fetch

<!-- BEGIN GENERATED: data fetch -->

| Flag | Default | Help |
|------|---------|------|
| `--cache-dir CACHE_DIR` | *(none)* | Root cache directory. Default: ~/.cache/viralscan/ |
| `--url URL` | *(none)* | Override archive URL. Intended for tests or mirrors; defaults to Zenodo. |
| `--sha256 SHA256` | *(none)* | Optional expected SHA-256 digest for the downloaded archive. |
| `--force` | `False` | Re-download and replace cached GTF files even if data already exists. |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->

### viralscan build-ref

<!-- BEGIN GENERATED: build-ref -->

| Flag | Default | Help |
|------|---------|------|
| `--host HOST` | *(none)* | Host species, e.g. 'human', 'mouse'. Run --list-species for all options. |
| `--virus-accessions ACCESSION` | *(none)* | One or more NCBI nucleotide accessions, e.g. NC_045512.2. |
| `--profile PROFILE` | `curated` | Frozen combined-reference profile. Default: curated. |
| `--output, -o OUTPUT` | `viralscan_ref` | Output directory for reference files. Default: viralscan_ref/ |
| `--ncbi-email NCBI_EMAIL` | *(none)* | Contact e-mail for NCBI E-utilities (avoids throttling). |
| `--ncbi-api-key NCBI_API_KEY` | *(none)* | NCBI API key for higher request rates. |
| `--cache-dir CACHE_DIR` | *(none)* | Root directory for download cache. Default: ~/.cache/viralscan/ |
| `--ensembl-release N` | *(none)* | Pin the Ensembl release used for the host cDNA/GTF (e.g. 116). Default: the current release from ftp.ensembl.org/pub/VERSION (an error if unreachable). Cached per release and recorded in the reference manifest. |
| `--no-kb-ref` | `False` | Skip running 'kb ref'; only produce concatenated FASTA and GTF. |
| `--genome-dlist FASTA` | *(none)* | Full host-genome FASTA passed to kallisto as a D-list and used for raw viral host-homology annotation. Requires minimap2 and the full tool profile. |
| `--strict-reference, --no-strict-reference` | `False` | Fail the build when two reference records are near-copies of each other on either strand (under 5% of one record's 31-mers lie outside the other) and are not declared catalogue siblings. Default: log a warning and continue. |
| `--anellovirus, --no-anellovirus` | `False` | Include the full packaged Anelloviridae accession table (2,041 accessions) in the combined host+viral reference (explicit opt-in; default: off). When --reference-panel anellovirus is used instead, builds an Anelloviridae-only reference without a host transcriptome; combine with --no-mask / --cluster for masking/clustering options. |
| `--allow-partial-panel` | `False` | Allow an incomplete expanded panel and write missing_accessions.tsv; default fails closed. |
| `--no-mask` | `False` | Skip dustmasker hard-masking of low-complexity viral sequence; the k-mer gate then allows a few. |
| `--cluster` | `False` | (--anellovirus) Run cd-hit-est clustering at 95% identity after masking. |
| `--reference-panel PANEL` | *(none)* | Build a pre-defined reference panel. Currently supported: 'anellovirus'. Uses the bundled FASTA from `viralscan data fetch` when available, otherwise falls back to NCBI accession download (same as --anellovirus). Combine with --no-mask / --cluster for masking/clustering options. |
| `--list-species` | `False` | Print all supported host species and exit. |
| `--cores, -c N` | *(none)* | Total CPU budget for kb ref, minimap2 homology and cd-hit-est. Native workers are bounded to N (kb ref min(8, N), minimap2 min(3, N), cd-hit-est -T N). Unset keeps the legacy requests (8 / 3 / all CPUs). |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->

### viralscan evidence

<!-- BEGIN GENERATED: evidence -->

| Flag | Default | Help |
|------|---------|------|
| `--review-candidates` | `False` | Review all Anelloviridae candidates, including below-gate and STAR-only rows; write results/evidence_review.tsv. |
| `--run-dir RUN_DIR` | *required* | A completed ViralScan run output directory. |
| `--output, -o OUTPUT` | *required* | Directory for evidence outputs. |
| `--viral-fasta VIRAL_FASTA` | *(none)* | Viral genome FASTA to align extracted reads against (enables BAM/coverage/BLAST). Omit for read-extraction only. |
| `--host-fasta HOST_FASTA` | *(none)* | Full host-genome FASTA required with --viral-fasta for competitive v3 alignment and BLAST. |
| `--competitor-fasta COMPETITOR_FASTA` | *(none)* | Explicit alternate/related/decoy FASTA; requires --competitor-manifest. |
| `--competitor-manifest COMPETITOR_MANIFEST` | *(none)* | Checksum-pinned TSV covering target, host and competitor FASTA records. |
| `--blast-tie-delta BLAST_TIE_DELTA` | `0.0` | Nonnegative bitscore distance retaining top BLAST ties (default: 0). |
| `--virus VIRUS` | *required* | Required exact target: accession/gene ID, canonical detected virus label, or registered alias (for example EBV or HHV6B). Substring matching is not used. |
| `--blast` | `False` | BLAST a sample of extracted reads against the viral reference (requires blast+). |
| `--read-start-profile` | `False` | Write a per-position 5' read-start distribution along the viral genome (read_start_profile.tsv). Requires --viral-fasta. |
| `--cell-tags` | `False` | Write viral_reads.tagged.bam with CB/UB cell-barcode tags (from read names) for per-cell IGV inspection (group by tag CB). Requires --viral-fasta. |
| `--dedup DEDUP` | `umi` | PCR-duplicate handling for --read-start-profile: umi (collapse per CB+UMI; default), markdup (samtools markdup), or none. Default: umi. |
| `--bin-size BIN_SIZE` | `1` | Bin width (bp) for the read-start profile. Default: 1. |
| `--sampling-seed SAMPLING_SEED` | `42` | Seed for deterministic BLAST read sampling. Default: 42. |
| `--cores, -c CORES` | `4` | Threads for minimap2/samtools/blast. |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->

### viralscan evidence-plot

<!-- BEGIN GENERATED: evidence-plot -->

| Flag | Default | Help |
|------|---------|------|
| `--evidence-dir EVIDENCE_DIR` | *required* | One or more 'viralscan evidence' output directories. |
| `--output, -o OUTPUT` | *(none)* | Write <output>/<evidence dir name>/competitive_alignment.{png,pdf}. Default: inside each evidence directory. |
| `--cross-call PREFIX` | *(none)* | Also write PREFIX.tsv and PREFIX.png pooling all --evidence-dir directories. |
| `--layer LAYER` | `umi` | Count post-alignment UMI representatives (default) or diagnostic raw reads. |
| `--max-lanes MAX_LANES` | `60` | Draw the lane block (panel A) only up to this many molecules; the aligned-length histogram (panel D) is always drawn. Default: 60. |
| `--coverage-map` | `False` | Draw coverage and optional gene track from retained BAMs. |
| `--target TARGET` | *(none)* | Exact BAM reference; optional only for one explicit viral reference. |
| `--annotation ANNOTATION` | *(none)* | GTF annotation; seqid = target, or the bare accession for viral CLASS\|ACC targets. |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->

### viralscan rerun-multimap

<!-- BEGIN GENERATED: rerun-multimap -->

| Flag | Default | Help |
|------|---------|------|
| `--run-dir RUN_DIR` | *required* | Completed source run. It is never modified. |
| `--output, -o OUTPUT` | *required* | New result directory; must not already be non-empty. |
| `--multimap-method {conservative,em,equal}` | *required* | Multimapping resolution method to apply: conservative, em or equal (see 'viralscan run --help' for what each does). |
| `--cores, -c CORES` | `6` | Number of cores for snakemake workers. Default: 6. |
| `--multimap-em-max-iter N` | *(none)* | (em only) Maximum EM iterations. Default: preserves existing config value. |
| `--multimap-em-tol TOL` | *(none)* | (em only) EM convergence tolerance. Default: preserves existing config value. |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->

### viralscan rerun-programs

<!-- BEGIN GENERATED: rerun-programs -->

| Flag | Default | Help |
|------|---------|------|
| `--run-dir DIR` | *required* | A completed viralscan run directory (the one holding kb-python/ and results/). |
| `--programme-min-breadth N` | *(none)* | Distinct non-overlapping overlap groups required before a programme is called. Must be >= 1. Default: the run's own config (2 if it set none). |
| `--programme-latent-min-breadth N` | *(none)* | Distinct latent overlap groups required before 'latent' is called. Must be >= 1. Default: the run's own config (1 if it set none). |
| `--programme-min-umi X` | *(none)* | A marker counts as evidence only at >= X molecules; 0 means any nonzero value. Must be >= 0. Default: the run's own config (0.0 if it set none). |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->

### viralscan hostresponse

<!-- BEGIN GENERATED: hostresponse -->

| Flag | Default | Help |
|------|---------|------|
| `--output, -o OUTPUT` | *required* | Existing viralscan sample output directory (contains config.yaml). |
| `--host-h5ad PATH` | *required* | Host gene-expression h5ad (cells × genes, matched to the viralscan run). |
| `--n-seeds N` | *(none)* | Random seeds for multi-seed L2 regression, 1-6 (default: from config or 6). |
| `--n-stab-iter N` | *(none)* | Selection-frequency (randomized-Lasso) iterations (default: from config or 100). |
| `--use-hvg, --no-use-hvg` | `True` | Use highly variable genes as features (default: on). --no-use-hvg uses all genes. |
| `--stab-min-prob P` | *(none)* | Min uncalibrated selection frequency to flag a gene; a ranking threshold, not an error rate (default: 0.6). |
| `--top-n-genes N` | *(none)* | Top N above-threshold genes to pass to pathway enrichment (default: 50). |
| `--min-umi N` | *(none)* | Minimum selected-method viral molecule estimate for the per-cell positive label (default: from config or 10). |
| `--label LABEL` | `raw` | Candidate-support label: 'raw' (default, depth-confounded estimate>=threshold) or depth-normalized 'cpm'/'fraction' (prevalence-matched burden per host molecule). |
| `--depth-match` | `False` | Restrict analysis to a coarsened-exact depth-matched cohort; this reduces measured total-depth imbalance but does not establish independence from depth confounding. |
| `--mito-control, --no-mito-control` | `True` | Add %mito as a covariate to the per-gene E-values (default: on; --no-mito-control to disable). |
| `--gene-symbols` | `False` | Annotate output CSVs with HGNC symbols from Ensembl IDs via mygene.info (network). |
| `--differential` | `False` | Write a genome-wide depth/%mito-adjusted differential table (<virus>_differential.csv). |
| `--enrichment` | `False` | Run pathway enrichment via gget (requires viralscan[enrichment]). |
| `--enrichment-db DB` | *(none)* | gget.enrichr database (default: GO_Biological_Process_2023). |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |
| `--cv CV` | *(none)* | CV scheme (default: auto). auto = group-disjoint CV on --groups, or on the one donor/sample-like obs column if there is exactly one, else cell-level CV with a leakage warning (also if group CV cannot run on that auto-selected column); cell = legacy cell-level CV (a donor can span train and test); group = group-disjoint CV, requires --groups, fails closed. |
| `--groups GROUPS` | *(none)* | obs column declaring donor/sample groups; used by --cv group and --cv auto (ignored, with a warning, under --cv cell). |
| `--cv-folds CV_FOLDS` | *(none)* | Number of group-disjoint folds (default: 5). |
| `--cell-type-column CELL_TYPE_COLUMN` | *(none)* | obs column used for cell-type baselines and strata. |
| `--cell-types CELL_TYPES` | *(none)* | Declared cell types to analyze separately (space- or comma-separated). |
| `--panel-in-fold, --no-panel-in-fold` | *(none)* | Fit the selection-frequency panel within each training fold. |
| `--permutations PERMUTATIONS` | *(none)* | Structured-null replicates (default: 0). |
| `--permutation-unit PERMUTATION_UNIT` | *(none)* | Explicit exchangeability unit for the null. |
| `--permutation-block PERMUTATION_BLOCK` | *(none)* | obs column restricting label exchangeability. |
| `--min-negative-cells MIN_NEGATIVE_CELLS` | *(none)* | Minimum negative cells per stratum (default: 10). |
| `--min-groups MIN_GROUPS` | *(none)* | Minimum groups per stratum (default: 2). |

<!-- END GENERATED -->

### viralscan check-whitelist

<!-- BEGIN GENERATED: check-whitelist -->

| Flag | Default | Help |
|------|---------|------|
| `--r1 R1` | *required* | R1 FASTQ (barcode read). |
| `--whitelist, -w PATH` | *required* | Barcode whitelist (optionally .gz). |
| `--technology, -x TECHNOLOGY` | `10xv3` | Single-cell technology (default: 10xv3). |
| `--min-match-rate F` | `0.5` | Minimum acceptable barcode match rate (default: 0.5). |
| `--n-sample N` | `100000` | Number of R1 reads to sample (default: 100000). |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->

### viralscan check-chemistry

<!-- BEGIN GENERATED: check-chemistry -->

| Flag | Default | Help |
|------|---------|------|
| `--r1 R1` | *required* | R1 FASTQ (barcode read). |
| `--r2 R2` | *(none)* | R2 FASTQ (for the pilot). |
| `--index, -i INDEX` | *(none)* | kallisto index (enables the strand pilot). |
| `--transcripts, -t TRANSCRIPTS` | *(none)* | t2g file (enables the strand pilot). |
| `--whitelist, -w PATH` | *(none)* | Barcode on-list. |
| `--technology, -x TECHNOLOGY` | *(none)* | Technology you plan to use; checked. |
| `--cores, -c CORES` | `4` | Threads for the pilot. |
| `--json PATH` | *(none)* | Also write the report as JSON. |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->

### viralscan doctor

<!-- BEGIN GENERATED: doctor -->

| Flag | Default | Help |
|------|---------|------|
| `--profile PROFILE` | `full` |  |
| `--json PATH` | *(none)* | Also write the report as JSON. |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->

### viralscan validate-run

<!-- BEGIN GENERATED: validate-run -->

| Flag | Default | Help |
|------|---------|------|
| `--run-dir RUN_DIR` | *required* | ViralScan output directory containing run_manifest.json. |
| `--no-verify-inputs` | `False` |  |
| `--json PATH` | *(none)* | Also write the report as JSON. |
| `--verbose` | `False` | Enable DEBUG-level logging. |
| `--quiet` | `False` | Suppress INFO messages. |

<!-- END GENERATED -->
