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

The current local corpus has 2,628 reported structure problems, mostly legacy
CDS transcripts without exons and reused transcript IDs, plus 47 anellovirus
annotation gaps. Its reverse catalogue check has no unlisted GTF accessions.
The four stale gene rows for the excluded `AB303562.1` duplicate have been
removed; the canonical `NC_038359.1` annotations remain.

Only an explicit `gene_biotype` or `transcript_biotype` of `whole_genome` labels a
feature as such; gene names and long spans do not classify it. The local corpus
explicitly tags `M12737_gene1`. HAV and HTLV-2 genome-spanning RNA annotations
lack this tag and remain unresolved. Tagging an entire mixed accession would
also exclude its genuine CDS. Whole-genome gene/programme exclusions and their
catalogue integration therefore require the verified feature classification;
this check does not silently change counts or gene membership.

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

SV40 and AAV1/7/8 have no cached NCBI flatfiles in this checkout. The accepted
overlay decisions are ready in
`analysis/panel_expansion/promotions_pending_records.tsv`. To fetch the records
when network work is authorized:

```bash
PYTHONPATH=src python -c 'import os; from viralscan.scripts.ncbi_fetch import fetch_reference; fetch_reference(["NC_001669", "NC_002077", "NC_006260", "NC_006261"], "/tmp/viralscan_pending_records", email=os.environ["NCBI_EMAIL"])'
```

Import those cached records into the existing catalogue while preserving its
rows before applying the overlay:

```bash
python scripts/panel_promote.py analysis/panel_expansion/promotions_pending_records.tsv --check
```

The check fails until all four accessions exist in the catalogue. The catalogue
generator writes the requested set rather than merging it, so passing only
these four accessions with the live catalogue as its output would destroy the
remaining rows. GTF promotion, catalogue merge and broad-discovery re-sync are
still pending their records; no speculative entries were added.
