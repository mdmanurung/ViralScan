# F-015 — The catalogue is not the index, and VIRTUS2's real value is HPV types + EBV type 2

**Date:** 2026-09-27
**Status:** confirmed (all numbers measured, reproducible offline)
**Tags:** reference-panel, catalogue, index, EBV, HPV, VIRTUS2, VirDetect, provenance

## Claim

Three separate findings. The first is the largest actionable defect in the
reference and has nothing to do with any external tool.

### 1. Most of the catalogue is not in the index (biggest gap)

`src/viralscan/data/virus_catalog.tsv` has 2,249 accessions, but the built
reference `covid_viralscan/viralscan_ref/viral_genome.dedup.fa` (2,312 records)
does **not** contain many of them. A catalogued virus that is not indexed is
**undetectable** — no amount of threshold work can recover it.

| Family | Catalogued | Not indexed |
|---|---:|---:|
| Papillomaviridae | 16 | **13** |
| Retroviridae | 5 | **5** (all of them) |
| Orthomyxoviridae | 30 | 7 |
| Coronaviridae | 7 | 4 |
| Polyomaviridae | 8 | 3 |
| Flaviviridae / Hepadnaviridae / Adenoviridae | 10 / 1 / 1 | 0 |

The 13 unindexed HPV are types **18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59,
66, 69** — including HPV18 and HPV31, the two highest-risk oncogenic types.

**All five Retroviridae are unindexed**: HTLV-1 (`NC_001436.1`), HTLV-2
(`NC_001488.1`), HIV-1 (`NC_001802.1`), HIV-2 (`NC_001722.1`), and simian
foamy virus (`NC_001364.1`). Both HIV genomes *and* both HTLV genomes are
therefore undetectable today. This bears directly on the tonsil work, since
HTLV-1/2 are T-cell-tropic and a plausible confounder for a T-cell compartment
result.

Herpes coverage, by contrast, is complete: all 8 catalogued Orthoherpesviridae
(HSV-1/2, VZV, HHV-5/6A/7, KSHV) are indexed, alongside EBV type 1.

**Implication:** the highest-value reference work is *reconciling the catalogue
against the built FASTA*, not importing a new genome collection. An import that
does not also fix the build reproduces the same invisibility at larger scale.

### 2. EBV type 2 is missing, and it is type-discriminating

The index carries **only** `NC_007605.1` (EBV type 1 / B95-8, 171,823 bp).
`NC_009334.1` (EBV type 2 / AG876, 172,764 bp) is in neither our catalogue nor
the index. Measured over distinct unambiguous 31-mers:

```
EBV1 distinct 31-mers  144,283
EBV2 distinct 31-mers  144,736
shared               115,215   79.9% of EBV1 / 79.6% of EBV2
EBV1-unique           29,068   20.1%
EBV2-unique           29,521   20.4%
```

So ~20% of the k-mer space is type-specific. A type-2 infection is currently
**uncallable**; its reads either drop out or smear as multimaps against type 1.
Adding `NC_009334.1` is a one-accession change with a real, measurable
consequence, and it makes EBV typing possible rather than merely EBV detection.

**Correction to an earlier subagent report.** That report stated "77% of EBV1
31-mers are type-1-unique." Measured, the ~77–80% figure is the **shared**
fraction; the unique fraction is **20.1%**. The direction of the conclusion is
unchanged (adding type 2 gains resolvable evidence) but the magnitude is ~20%,
not ~77%. Recorded here so the number is not re-cited from the report.

### 3. The host reference carries no viral contigs (negative result)

`VirDetect`'s documentation warns that some hg38 builds include a `chrEBV`
contig, which would make EBV silently invisible to a host-subtraction pipeline.
Checked both host inputs actually used for the D-list:

- `refdata-gex-GRCh38-2024-A/fasta/genome.fa` — 194 contigs, **zero** contigs
  matching `EBV|chrEB|herpes|virus|viral|HHV|KSHV|TTV|torque|HBV|HCV`.
- `Homo_sapiens.GRCh38.116.gtf.gz` — **zero** viral seqnames.

The `chrEBV` hazard **does not apply** to this build. Recorded as a closed
negative so it is not re-investigated, and so any future host-reference swap
re-runs the check rather than assuming it.

## CORRECTION (2026-09-28) — the HPV type counts above were wrong

An earlier version of this finding reported **92** VIRTUS2 HPV types and **76**
new. Both are wrong. The type-extraction regex was `HPV\d{1,2}`, which matches
only one- or two-digit type numbers and silently missed every three-digit type
(HPV100–HPV171). Re-counted with `HPV\d+`:

```
total distinct HPV types in VIRTUS2 : 178   (my own recount)
VIRTUS2 HPV types                   : 182   (independent recount, 3 methods)
already catalogued by type          :  16
genuinely new                       : 166
```

So the HPV gap is roughly **twice** the size first reported, not modestly larger.
The direction of every conclusion here is unchanged — the gap is real, and it is
not an accession-matching artifact — but the magnitude was understated and any
plan sized on "76 types" was sized wrong. See F-016 for the built panel, which
resolved 182/182 types and measured the multimapping risk directly.

## What VIRTUS2 actually adds (and does not)

`200830_viruses.txt` = 762 lines; `viruses.fasta` = 762 records. The list uses
**space-separated** accessions (`NC 000883.2`) while the FASTA uses underscores
(`NC_001499.1_...`) — a naive accession regex silently drops ~96 of 762 records,
which is how an earlier pass of this analysis undercounted the collection.

After normalizing both forms, against our 2,249-accession catalogue:

| | count |
|---|---:|
| already in our catalogue | 169 |
| new to us | 497 |
| — HPV types (**by type, not accession**) | **166** |
| — influenza (legacy lab strains) | 37 |
| — TTV / anellovirus | 27 |
| — animal retroviruses | 5 |
| — hepatitis | 5 |
| — coronaviruses | 2 |
| — herpesvirus (**= `NC_009334.1`, EBV type 2**) | 1 |
| — other (largely segmented / animal) | 325 |

HPV type comparison is the honest one, because VIRTUS2 carries HPV as
`gi|…|lcl|HPV##REF.1` records whose accessions will never match our `NC_`
accessions even for the same type:

```
our catalogued HPV types   16   (HPV1,4,16,18,31,33,35,39,45,51,52,56,58,59,66,69)
VIRTUS2 HPV types          182
genuinely new            166
```

So the HPV expansion is real, not an accession artifact: **16 → 182 types**.

## Recommendations

1. **Reconcile catalogue vs index before importing anything** (largest win).
   Add a build-time assertion that every catalogued accession intended for
   detection is present in the emitted FASTA. This is the defect that made 13
   HPV and all HIV invisible.
2. **Add `NC_009334.1`** (EBV type 2). Cheap, high value, enables typing.
3. **HPV: prefer RefSeq `NC_` accessions over VIRTUS2 `lcl|` records** for the
   same 76 types. The types are what we need; the `lcl|` provenance is a
   liability (see caveats), and RefSeq accessions give us GTFs for free.
4. **Skip VIRTUS2's influenza.** It is PR8 / H9N2 HK/97 / H5N1 goose 1996 /
   H3N2 NY/2004 / B/Lee/1940 / C/Ann Arbor — legacy lab strains, no better than
   what we already hold, and the 37 "new" ones are mostly segmented negatives.
5. **Skip the 27 anelloviruses.** Already covered with far better provenance by
   the upstream 2,023 hardmasked representatives.
6. **Take the method, not the genomes, from VIRTUS2/VirDetect:** host
   subtraction before viral quantification, per-strand counting, and explicit
   multimap handling. These are pipeline properties, not reference properties.

## Caveats

- `lcl|` records are RefSeq *local* submissions with weaker curation than `NC_`
  curated records. Verify each before indexing; do not bulk-import.
- The 325 "other" records were bucketed by name only; several are animal or
  segmented viruses of little value for human scRNA-seq. Not audited
  individually.
- EBV type 1/2 k-mer sharing was computed on whole genomes including the
  internal repeat regions, which are more conserved than the unique-sequence
  space. 20% is therefore a conservative floor for type-discriminating evidence.
