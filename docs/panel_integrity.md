# Panel integrity checks

The bundled reference builder checks the local GTF corpus hash manifest before
downloads. The default is
`analysis/panel_expansion/gtf_corpus_manifest.tsv`; use `--gtf-manifest` to select
another reviewed corpus. The manifest records hashes and sizes, with unknown
source dates. It identifies local files and does not imply GTF bundling in the
wheel or validate their biological annotations.

Run the offline structure and catalogue matrix check with:

```bash
PYTHONPATH=src python scripts/panel_integrity.py --output /tmp/panel_integrity.tsv
PYTHONPATH=src python scripts/write_gtf_manifest.py \
  --check analysis/panel_expansion/gtf_corpus_manifest.tsv
```

The matrix includes every shipped accession, deliberate index exclusions, and
reverse GTF accessions absent from the catalogue. Real anellovirus gene-table
exons are included because the builder emits them during assembly. Missing gene
annotations remain visible even when the builder can create a whole-genome
fallback from its FASTA. Duplicate records, duplicate gene/transcript ownership,
malformed rows and transcripts without exons fail the structure check. These
checks are separate from hashing: a matching hash does not repair annotation.

By default the check inspects the normalised view, which is what the bundled
builder gives `kb ref` (REF-13, `src/viralscan/gtf_normalise.py`; the shipped GTF
files are never rewritten). `--raw` gates on the files as they are, and the
summary line always prints their count. The normaliser gives every CDS-only gene
a gene-scoped transcript (ID = `gene_id`) with an exon, retags its CDS rows, and
removes introns only for single-protein spliced CDS; the rules and their limits are
in the module docstring.

The current local corpus has 2,626 reported structure problems in the shipped
files (legacy CDS transcripts without exons and reused transcript IDs) and **7** in
the normalised view, plus 47 anellovirus annotation gaps. The 7 are CDS rows whose
transcript has no exon inside a gene that already has exons (Ebola 2, EBV, HCMV, HHV-8,
measles, West Nile); they are left as shipped because no rule derives their exons. Its reverse catalogue check has no unlisted GTF accessions.
The four stale gene rows for the excluded `AB303562.1` duplicate have been
removed; the canonical `NC_038359.1` annotations remain. (2,628 became 2,629 on
2026-10-10: Human astrovirus `NC_001943.1` wrote `Non structural gene` as its
column-3 feature type, so two CDS rows were invisible to this check; as `CDS`
they now show the same exon-less shape as the other legacy records. 2,629 became 2,626 when the duplicate `NC_002076.2` GTF
`Torque_teno_virus_NC_002076.gtf` was removed in favour of `Alphatorquevirus homin1.gtf`.)

Only an explicit `gene_biotype` or `transcript_biotype` of `whole_genome` labels a
feature as such; gene names and long spans do not classify it. As of 2026-10-10
the local corpus tags no gene as `whole_genome`: HPV-8 `M12737.1` and HPV-71
`NC_039089.1` were given modelled CDS (`scripts/model_hpv_cds_gtfs.py`, QC in
`analysis/panel_expansion/hpv_modelled_cds_qc.tsv`), and the genome-spanning
`misc_RNA` / `prim_transcript` of HAV `NC_001489.1` and HTLV-2 `NC_001488.1` were
removed (`scripts/fix_thin_gtfs.py`), keeping their real CDS. `whole_genome`
genes still arise for FASTA-only user references; excluding them from
gene-level outputs and programmes (`DEF-04`) is deferred because the
whole-genome gene is such a virus's only unit. This check does not silently change counts or gene membership.

The EC budget helper accepts a normalized JSON object with integer
`max_ec_size` and `discarded_ec_count` values:

```bash
PYTHONPATH=src python scripts/panel_integrity.py --output /tmp/panel_integrity.tsv \
  --ec-summary /tmp/normalized_inspect.json --max-ec-size 20 --max-discarded-ecs 0
```

The numbers above illustrate explicit budgets, not calibrated defaults. The
helper rejects missing/invalid metrics and exceeded budgets. It does not run
or parse `kallisto inspect`; CAT-40 needs the real index experiment and an
audited normalization of that result. `--max-ec-size`/`--max-discarded-ecs`
without `--ec-summary` is a usage error (the budgets would otherwise be ignored),
and a catalogue with no `shipped` rows or an unknown `panel` value is rejected
rather than checked vacuously.

`viral_summary.tsv` and the HTML report share `reference_risk_flags` from the
run Virus Identity table. Endogenous and contaminant roles produce
`endogenous_reference` and `vector_reagent_reference`; the curated `HHV-6`
sibling group carries `iciHHV6_possible`. An EVE-risk annotation additionally
carries `endogenous_overlap_possible`. These labels require interpretation;
they do not identify the source, diagnose integration, suppress counts or
establish infection. Automatic HERV exclusion and a calibrated CAR-T/vector
filter remain separate unfinished decisions.

SV40 (`NC_001669.1`) and AAV1/7/8 (`NC_002077.1`, `NC_006260.1`, `NC_006261.1`) were fetched from NCBI,
catalogued as `shipped` + `contaminant` and given generated GTFs on 2026-10-10 (PLAN `PANEL-01` WP2). The GTFs
stay local and are listed in `analysis/panel_expansion/gtf_corpus_manifest.tsv`.
`analysis/panel_expansion/promotions_pending_records.tsv` is the overlay that recorded the decisions;
`scripts/panel_promote.py ... --check` passes now that all four accessions are in the catalogue. The catalogue
generator writes the requested set rather than merging it, so never run it against the live catalogue without
`--out`.

## REF-13 normaliser QC

`scripts/ref13_normaliser_qc.py` generates the `kb ref` cDNA with the installed
`ngs_tools` before and after normalisation, on the cached genomes, and writes
`analysis/panel_expansion/ref13_normaliser_qc.tsv` (one row per gene: rule, transcript
count and cDNA length before and after, `bp_removed`, `identical`). Over the 393 local GTFs
(2026-10-10): 4,147 genes, **4,040 identical** and 107 changed.

- 101 genes are `spliced`: about 99 kb of intron removed from 100 single-copy genes, and
  UL111A (`HUM_CYTO_HHV5wtgp098`), which keeps both its spliced and its span transcript, so it
  gains 531 bp. Thirteen collapsed joins in 11 HHV-6B genes (DR1 and DR6 once per copy) need the cached
  flatfile, because the shipped HHV-6B GTF collapsed them into one CDS row each.
- 6 HHV-6B genes (DR1, DR3, DR6, B1, B2, B3) have one transcript per terminal-repeat copy
  (`G`, `G-c2`), so `bp_removed` is negative for them: the first copy was lost before.
- HBV P and S (origin-wrapping, two gene rows with `part`) are unchanged: `ngs_tools`
  concatenates exons in genomic order, so the wrap cannot be represented and the last gene
  row still wins (`wrap_last_row`). Their first-part sequence (2309-3182 for P, 2850-3182 for S)
  is still absent from the index.
- B19V `HUM_PARVO_unassigned_gene_1`: its CDS rows are re-keyed onto the gene-only gene
  `HUM_PARVO_B19V_gp4` at the same span (4890-5174), which is already indexed, rather than
  given a second gene row that would duplicate the locus.
- Adenovirus gp01/04/07/08/09/10/13/14/16 and other multi-protein genes keep the span.
- UL111A: Jenkins et al. 2004 (J Virol 78:1440) reports a latency transcript with a splicing
  pattern different from the lytic one; its abstract does not say that an intron is retained.
  The span transcript is therefore kept as a conservative second transcript, unverified.

A viral-only `kallisto index -k 31` of all 393 GTFs, before and after (`kallisto inspect`):
targets 4,206 to 4,213; unitigs 32,246 to 32,285; k-mers 4,653,849 to 4,638,071; max EC
size 17 and 0 discarded ECs in both.

Regenerating `gene_programs.tsv` with the normalised exons changes the overlap groups
(EBV independent groups 5 to 2, HHV-6B 5 to 3; 63 of 75 rows differ in some column). The committed
TSV has not been regenerated.
