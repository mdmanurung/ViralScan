# Translated search of RNA-seq against viral proteins

This exploratory workflow ports the indexing and counting method in
[Laura Luebbert's kallisto tutorial](https://kallisto.readthedocs.io/en/latest/translated/notebooks/virus_detection_sc.html).
`kb ref --aa --workflow custom` builds the protein index, and
`kb count --aa --parity single` translates the sequencing reads and performs
pseudoalignment. Kallisto and kb-python implement the method; ViralScan does
not implement a new aligner. Snakemake tracks preparation, host masking,
indexing, and optional counting.

## Reference and adaptation

The tutorial uses PalmDB RdRP proteins. Here the input is **all annotated CDS
proteins of the RNA and DNA viruses represented in the local GTF corpus**.
This includes retroviral, hepatitis-delta, and DNA-virus proteins. It is an
adaptation of the reference scope, not a reproduction of the PalmDB benchmark.
No proteins are sequence-clustered. All proteins for each reference genome
are included, rather than restricting them to individual GTF feature IDs.

`scripts/prepare_viral_proteins.py` uses the already cached GenBank records.
Each record must match both the accession version and the catalogue's
`gb_sha256`. Biopython reads the complete CDS `/translation` annotations,
including multiline translations; the package's minimal GenBank parser only
reads their first line and is deliberately not used for proteins. Existing
nonempty reference directories are preserved; rebuilds require a new directory.

Outputs:

- `proteins.faa`: annotated amino acid sequences, identified by genome
  accession version plus protein accession version.
- `viruses.t2g`: protein targets mapped to catalogue taxids, for per-virus UMI
  counting; related viruses remain ambiguous when their proteins share matches.
- `proteins.tsv`: protein, product, gene, source genome, taxid, and annotation
  location, for interpreting each match.
- `coverage.tsv`: one row for every GTF, including records without annotated
  translations. A protein-free reference record is a coverage gap, not a
  negative biological result.
- `reference.json`: catalogue, script, GenBank, and output hash provenance,
  method attribution, counts, and annotation gaps.

The 2026-10-14 local corpus has 393 GTFs and 299 distinct catalogue labels
(labels can represent genotypes/serotypes, not just taxonomic species).
It yields 4,124 proteins from 389 reference accessions, covering 296 labels:
492 proteins cover 111/113 RNA-virus labels, and 3,632 cover 185/186 DNA-virus
labels. `JN860200.1` (Alkhumra), `D90457.1` (Coxsackievirus A24), `M12737.1`
(HPV8), and `NC_005225.1` (Puumala L) have no annotated translations. Puumala
has proteins from other segments. Modelled nucleotide GTF features are not
silently translated into protein references; notably the Alkhumra record
describes a nonfunctional, mutated polyprotein. These gaps remain explicit.

## Run

Use an existing environment with Python >=3.10, Biopython, kb-python,
kallisto with translated-search support, and Snakemake. No package dependency
or production ViralScan CLI was changed. The checked local environment uses
kallisto 0.52.0, kb-python 0.30.2, and Snakemake 9.23.1. The index uses
`-k 31`, matching the linked tutorial. Counting reads k from the index;
kb-python 0.30.2 no longer accepts the tutorial's `kb count -k` option.
Both steps use the same explicit kallisto executable.

From the repository root:

```bash
cp analysis/translated_search/config.example.yaml analysis/translated_search/config.local.yaml
# Edit config.local.yaml: executable paths, host files, and optionally samples.
# Use a writable runtime cache if the environment's default cache is read-only.
export XDG_CACHE_HOME="$PWD/analysis/translated_search/work/cache"
PYTHONPATH=src snakemake --snakefile analysis/translated_search/Snakefile \
  --configfile analysis/translated_search/config.local.yaml --cores 2 --dry-run
PYTHONPATH=src snakemake --snakefile analysis/translated_search/Snakefile \
  --configfile analysis/translated_search/config.local.yaml --cores 2
```

`samples: {}` builds only the reference index. Supplying sample R1/R2 paths
also produces per-sample unfiltered virus-by-barcode H5AD counts. R1 must
contain barcodes/UMIs and R2 cDNA according to `technology`. `--parity single`
is required for translated search; the two single-cell input files retain
their technology-specific roles. This is not a paired-fragment protein
alignment. For bulk paired cDNA reads, the tutorial requires separate read
searches; this single-cell workflow does not expose a bulk paired mode.

Supply a **host-only nucleotide genome and cDNA reference**, never a combined
host-plus-virus reference as the host mask. The genome must be uncompressed
FASTA and cDNA gzip-compressed FASTA; the workflow concatenates them with a
record separator. The index receipt records their combined checksum, protein
reference checksum, kallisto/kb versions, and the exact build command. Outputs,
local configuration, and Snakemake state are ignored by Git. Run larger builds
through the institution's existing scheduler setup.

## Running arm P on many libraries

`make_configs.py` reads a DSR manifest plus one `chem/<sample>.json` per library and writes one
Snakemake config per group of identical technology, strand and whitelist (the workflow takes one
`technology` and at most one `whitelist` per config):

```bash
python analysis/translated_search/make_configs.py --manifest MANIFEST.tsv --chem-dir CHEM_DIR \
  --base BASE.yaml --out-dir analysis/translated_search/work/configs \
  --taxa analysis/translated_search/expected_taxa.tsv
```

- Technology is `chemistry.chemistry` of the chem JSON; the manifest's free-text `hint` is never read.
  (GSM4682311 is therefore 10xv2, the 2026-10-06 user decision, although its manifest hint says 10xv3.)
- R1 must be the barcode+UMI read. Each pair is opened and R1's modal read length must equal the chem
  JSON's `r1_length`; samples that fail are excluded and the exit status is 1. For ERR13027027 the barcode read is
  the `_2` file and the manifest already lists it as R1, so no swap is applied (`--swap SAMPLE` exists for
  manifests that get it wrong).
- A library whose chem JSON says `user on-list` needs the manifest whitelist; it is passed as the optional
  Snakefile key `whitelist`, which adds `kb count -w <file>`. Without the key, kb's bundled on-list is used
  (unchanged behaviour).
- `strand` is recorded in the config but not passed: kallisto's `--aa` mode searches all six frames of
  each read (`ProcessReads.cpp`, kallisto 0.50.0), so read orientation does not change the match.
- `expected_virus.tsv` maps each library to its expected taxon ids from `expected_taxa.tsv` (data, not code;
  taxon ids checked against `viruses.t2g` and `virus_catalog.tsv`). Add rows as datasets are discovered.

`summarize_translated.py` writes a per-library summary (UMIs and barcodes for the expected virus, top viruses,
top proteins, ambiguous multi-virus records) from a `kb count` directory; it needs `bustools` on PATH.
`--cells FILE` restricts it to a frozen called-cell list. `check_barcodes.py OBSERVED FROZEN` reports
overlap, suffix/length/case differences and the reverse-complement overlap, and exits 2 below `--abort-below`
(default 0.90) so a barcode-format mismatch cannot silently shrink a cell-restricted summary.
On the EBV 1M-pair pilot output the summary reproduces the AA-02 numbers (2,312 UMIs in 798 barcodes).

Only the Snakefile `whitelist` key changes the workflow. Keep outputs of real libraries under the ignored
`work/` directory.

## Interpretation

RNA-seq detects expressed protein-coding regions of DNA viruses as well as
RNA viruses. It does not measure viral DNA or physical particles. Noncoding
RNAs, including EBV EBERs and latency transcripts without protein products,
cannot match this protein reference. A protein match does not alone establish
infection, productive replication, or reactivation. Host masking, raw-read
verification, negative controls, and comparison to the nucleotide workflow
remain necessary before making biological claims. Per-virus counts here are
an exploratory output, separate from ViralScan's calibrated detection calls.

## Index provenance

Host mask of the shared kallisto 0.50.0 protein index. The shared protein index (`index/receipt.json`, kallisto 0.50.0, kb-python 0.30.2) records the plain
nucleotide mask `host_cdna_dna.fa` (sha256 `8c788453...601b`: GRCh38 genome plus Ensembl cDNA) as its D-list.
A retry note says an "encoded 0.52 copy" of the mask was reused. These two statements were checked
on 2026-10-14 and the receipt is the one that holds:

- The `index` rule of the retry ran `kb ref ... --d-list <plain mask>`; its Snakemake log, receipt and
  command all name the plain file, whose hash equals the receipt's. No recorded command passes the encoded file.
- The 0.52 attempt left two files in its `tmp/`: a 5.5 MB encoded copy of the viral proteins (the
  "`encoded_target_sha256`" in the note is this file, not the mask) and a 26 GB comma-free-code (CFC)
  encoded mask, sha256 `90e93fc7...4965`, whose hash no receipt records.
- kallisto 0.50.0 reads the D-list FASTA as is (`KmerIndex::DListFlankingKmers`); it does not translate
  it. A nucleotide mask is therefore compared with a graph built from CFC-encoded protein sequence.
  Controlled test (same proteins, `-k 31`, same kallisto/kb): the first 300 MB of the plain mask (chr1) gives
  108 distinguishing flanking k-mers; the first 1.75 GB of the encoded mask (the same chr1 region, six frames)
  gives 2,208. The flanking set is a union over D-list records, so it can only grow with more input: the shared
  index has 483, fewer than the 2,208 of this prefix of the encoded mask alone, so it cannot have used the
  encoded mask.

Consequence: the shared v0500 index used the **plain** mask, and in kallisto 0.50.0 that masks far fewer host
k-mers than the translated mask (about 20-fold fewer flanking k-mers on the chr1 comparison). Host masking in
that index is therefore largely inactive and host-derived off-target viral matches are not expected to be
suppressed; this is consistent with, not shown by, the single papillomavirus protein that ranks among the top
off-target hits in the HSV-1 pilot. The intended mask is the encoded one; building the index
with it under 0.50.0 (or giving 0.52 more than the one hour it was allowed) is not done here. Until then, treat
per-virus calls from this index as unmasked.

## Verification (2026-10-14)

Four focused extraction tests pass: complete multiline translations, DNA-virus
and retrovirus inclusion with annotation gaps, cache-corruption rejection,
accession-version rejection, and preservation of nonempty existing outputs.
Changed Python files pass lint/format, and the repository data-governance check
passes. These checks do not qualify the full ViralScan package for release.

A real-tool Snakemake run rebuilt the reference and protein index and counted
40 synthetic 10xv3 reads encoding two real reference proteins. The resulting
H5AD has one barcode and 296 virus columns, with exactly 20 UMIs for Aichi
virus (`taxon_1313215`) and 20 for HSV-2 (`taxon_10310`), and no other positive
counts. This confirms that an expressed protein from a DNA virus can be
detected by this workflow. The mask for this test contains synthetic sequences;
it is not a host-specificity experiment. The generated smoke receipt and logs
are under `work/smoke_2026-10-14/`.

The populated, ignored `config.local.yaml` points to the available human
GRCh38-2024-A genome and Ensembl 116 cDNA and selects no real samples. Its DAG
passes dry-run. The full human-mask index build and real-sample counts have
not been executed. The reference preparation receipt is
`reference_2026-10-14_v3/reference.json`.

## Attribution

- Tutorial: Laura Luebbert, *Viral proteins in single-cell RNA seq data*,
  updated January 7, 2025; linked above, consulted October 14, 2026.
- Method: Luebbert L, Sullivan DK, Carilli M, Hjorleifsson KE, Winnett AV,
  Chari T, Pachter L. *Detection of viral sequences at single-cell resolution
  identifies novel viruses associated with host gene expression changes*.
  Nature Biotechnology (2025).
  [doi:10.1038/s41587-025-02614-y](https://doi.org/10.1038/s41587-025-02614-y).
- Upstream method/code: [pachterlab/LSCHWCP_2023](https://github.com/pachterlab/LSCHWCP_2023)
  (BSD-2-Clause). No upstream source files are vendored here; the workflow
  invokes the documented kb/kallisto commands, with the adaptation described
  above. Existing dependency licenses continue to apply.
- Protein annotations: NCBI GenBank/RefSeq, pinned by accession version and
  source-record checksum in the catalogue and generated tables. PalmDB
  sequences are not used by this targeted reference. If a separate PalmDB
  discovery arm is added, cite Edgar et al., *Petabase-scale sequence alignment
  catalyses viral discovery*, Nature (2022),
  [doi:10.1038/s41586-021-04332-2](https://doi.org/10.1038/s41586-021-04332-2).
