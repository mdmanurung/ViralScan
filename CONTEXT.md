# ViralScan

Domain and architecture vocabulary for ViralScan — a Snakemake-driven CLI that
quantifies viral load from paired-end single-cell FASTQ samples using
`kb-python` (kallisto + bustools). This file pins terms so architecture reviews
and design conversations stay consistent; see `LANGUAGE.md` in the
`improve-codebase-architecture` skill for the architecture vocabulary
(module, interface, seam, depth).

## Language

### A run and its context

**Run**:
A single invocation of viral quantification for one sample pair — from FASTQ
input through the six Snakemake rules to the per-sample output directory.
_Avoid_: job, pipeline run, execution.

**Run Config**:
The validated, typed description of a Run — every parameter the rules need
(paths, thresholds, flags, multimap policy). Coerced and validated exactly once,
at the CLI→workflow edge (`from_cli`); read back by downstream rules as a
trusted typed value (`from_yaml`).
_Avoid_: config dict, settings, params, options.

**Run Context**:
A passive value passed to each worker script's `run(ctx)`, bundling the
**Run Config** with the **Kb Count Outputs** for that Run. Answers "what" and
"where"; it does not perform I/O.
_Avoid_: session, environment, state, handle.

**Kb Count Outputs**:
The module that owns the on-disk layout `kb-python` writes under a Run's output
directory — the named file paths (adata, bus, ec, barcodes, genes…) and the
resolution of which adata is current for a Run.
_Avoid_: paths, file map, output dir.

### Viral evidence

**Viral Gene**:
One gene of the index, identified by its `gene_id` (the GTF `gene_id`, column 2
of the kb `t2g`). The unit of counting: every count matrix column is a gene, and
a gene is viral or host. Whether it is viral is decided once, by the **Virus
Identity Table**.
_Avoid_: viral accession (that now means the genome), virus id.

**Viral Accession**:
The genome accession a **Viral Gene** sits on (the GTF seqname, column 5 of the
kb `t2g`, for example `NC_007605.1`). It is the join key to the virus catalogue
and so to the virus. Before MECH-A this word meant the gene_id; that is now
**Viral Gene**.
_Avoid_: viral gene id, genome id.

**Virus**:
One reported organism: all **Viral Genes** sharing one `virus_key`. The key is the
NCBI taxid for non-segmented viruses (segments of one isolate share it; only the
species-level taxid 11320 "Influenza A virus" splits, by strain), `genus:<Genus>`
for anelloviruses, `accession:<acc>` for a genome the catalogue does not list,
and `name:<legacy>` for a pre-v3 three-column `t2g`. Results are grouped by key
and named by `virus_name`; names are unique across keys.
_Avoid_: species, virus name (as an identifier), virus group.

**Virus Identity Table**:
The per-**Run** table `results/virus_identity.tsv`, one row per indexed **Viral
Gene** or host gene, written once by the `analysis` rule (`virus_identity.py`).
It carries `viral`, `virus_key`, `virus_name`, taxonomy, `sibling_group` and
`risk_class`. `multimap`, `detection`, `umap`, `hostresponse`, `evidence` and
gene programmes read it; none re-derives a gene's virus from its ID string. A
Run with no table (made before it existed) falls back to `log/analysis.txt` and
the legacy prefix rules in `virus_grouping.py`.
_Avoid_: name map, gene-to-virus map, analysis.txt (the old viral-gene list).

**Index Kind**:
Whether the kb index is `combined` (host and viral genes) or `virus_only`, read
off the **Virus Identity Table** (any host row means combined). It decides which
multimap policies apply.
_Avoid_: reference type, index mode.

**Detection**:
Calling a Viral Accession present in a sample because its total UMI count meets
the `detection_threshold` (inclusive `>=`).
_Avoid_: hit, call, positive.

**Multimapper Correction**:
Redistribution of UMI mass from equivalence classes that span multiple genes,
per the configured `multimap_method`. Owned by `multimapping.py`.
_Avoid_: ambiguity resolution, dedup.

## Relationships

- A **Run** has exactly one **Run Config**.
- A **Run Context** bundles one **Run Config** and one **Kb Count Outputs**.
- **Run Config** is validated once at `from_cli`; every later read is `from_yaml`
  against a file only `menu._write_run_config` writes.
- **Kb Count Outputs** resolves the current adata for a Run from the
  **Run Config**'s `multimapping` flag — not from file existence. If the flag is
  set the multimap adata must exist; a missing file is an error, not a fallback.
- **Detection** reads counts for each **Viral Gene** from the current adata,
  groups them by **Virus** through the **Virus Identity Table**, and sums per
  virus.
- A **Viral Gene** belongs to exactly one **Virus** (one `virus_key`); a **Virus**
  has one or more **Viral Accessions** (segments, isolates, genomes).

## Example dialogue

> **Dev:** "When `detection` needs the count matrix, does it decide which adata
> to read, or does the **Run Context** hand it the current one?"
> **Maintainer:** "The **Kb Count Outputs** resolves the current adata — there's
> one rule for it, so `detection` and `umap` can't drift. `detection` just asks
> for `ctx.outputs.current_adata`."

## Flagged ambiguities

- "config" was used for three different things: the argparse Namespace, the
  Snakemake `--config` dict, and the YAML on disk. Resolved: the typed value is
  the **Run Config**; the others are representations it is built from, not the
  thing itself.
- "viral gene" and **Viral Accession** were used interchangeably. Resolved
  (MECH-A): **Viral Gene** is the `gene_id`; **Viral Accession** is the genome
  accession (t2g column 5) it sits on. They are different things.
- Mapping a **Viral Accession** to a virus name was done two ways — substring
  match (`detection`) and prefix match (`umap`) — which disagreed on 151/2692
  bundled gene IDs and were both buggy. Resolved: one **boundary-aware** rule in
  `virus_grouping.py` (key matches when the gene ID equals it, or starts with it
  and the next char is `_` or a digit; longest key wins).
  Since MECH-A step 4 that rule is the legacy fallback, used only when a Run has
  no **Virus Identity Table**; the table assigns the virus from the genome
  accession instead.
