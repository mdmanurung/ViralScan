# ANELLO-15: draft request to the clareaulab authors (NOT SENT)

Status: draft for the maintainer to review and send. Nothing here has been sent to anyone.

## What changed upstream since the row was written

`clareaulab/human_anellovirus_pangenome` gained a `LICENSE` on 2026-10-08 (commit `517c264`, "Create LICENSE"):
the MIT licence, "Copyright (c) 2026 Lareau Lab, Memorial Sloan Kettering". The row said the repo had no licence;
that part of open item (1) is answered. The README still gives no paper, DOI or citation text.
Still unanswered: a citable release, and permission to ship the `orf1_genus` column inside ViralScan.
The MIT text probably covers redistribution of `simple_anello_metadata_V2.csv`; ask anyway, because the
column is analysis output with unstated authorship.

## Draft message

Subject: Citing and reusing the human anellovirus pangenome in ViralScan

Dear Lareau lab,

ViralScan is an open-source viral-load tool (kallisto + bustools) with an anellovirus panel. Its accession list
follows your non-redundant set of 2,023 representatives (`clareaulab/human_anellovirus_pangenome`, commit
`3ed77e1`, later `517c264`). We fetch sequence and CDS from NCBI ourselves, so we ship none of your FASTA.
Thank you for the resource, and for adding the MIT licence on 2026-10-08.

We have three questions.

1. Is there a paper, preprint, DOI or Zenodo release we should cite? Today we cite the URL and commit hash.
2. May we ship the `orf1_genus` column of `simple_anello_metadata_V2.csv` (keyed by accession) in our package
   data, with attribution? It resolves the genus for 576 of the 702 anellovirus rows our NCBI-derived catalogue
   leaves generic.
3. Fifteen of your `orf1_genus` calls differ from the genus we hold from NCBI (table attached). Which method and
   tree version produced `orf1_genus`, and do you consider the nine Beta-to-Heto calls final?

Best regards,
(maintainer name)

## Disagreement table (measured 2026-11-04)

Inputs: our `src/viralscan/data/virus_catalog.tsv` (family Anelloviridae, genus set) against the lab CSV at
commit `517c264`. 3,422 of our 3,452 anellovirus rows join on `accession_version` (1:1). The 15 differences
below are the same on all 3,422 rows and on the 2,023 `cdhit_representative` rows. This table is
read-only: no catalogue, FASTA, GTF, t2g or membership file changed, so the frozen cat42d panel inputs
(G4) are untouched.

| Our genus (NCBI) | Lab `orf1_genus` | n | Accessions |
|---|---|---:|---|
| Betatorquevirus | Hetorquevirus | 9 | MN769141.1, MN769193.1, MN769310.1, MN769649.1, MN770248.1, MN771065.1, MN771204.1, MN771334.1, MN771373.1 |
| Sadetorquevirus | Hetorquevirus | 2 | MN771248.1, MN769236.1 |
| Yodtorquevirus | Hetorquevirus | 1 | MN770091.1 |
| Gammatorquevirus | Betatorquevirus | 1 | MH017566.1 |
| "Small anellovirus" (no genus) | Gammatorquevirus | 1 | AY622909.1 |
| Lamedtorquevirus | Unclassified | 1 | MN774989.1 |

Reading. The "9 Beta" in the row reproduces exactly (all nine go to Hetorquevirus). The "3 Gamma" does not
reproduce on this snapshot: only one of our Gammatorquevirus rows differs (MH017566.1, to Beta), and one more
Gamma call sits on a row we hold under a non-genus label (AY622909.1). The row's figure probably came from a
different join or an earlier CSV; the table above is the current measurement. The "171 of 185" generic-row
figure also drifts: 576 of 702 generic catalogue rows (576 of 581 among the 2,023 representatives) carry an
`orf1_genus`.

## Decisions left for the maintainer

- Send the request, and from which address.
- Whether to adopt `orf1_genus` at all. It would change genus labels for catalogue rows and therefore needs a
  new panel version; cat42d stays frozen until then (G4).
- Which genus wins on the 15 conflicts (NCBI taxonomy vs. lab phylogeny). Not decided here.
