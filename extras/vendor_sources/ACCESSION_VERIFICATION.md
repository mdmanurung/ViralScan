# ViralScan reference panel — live NCBI accession verification

**Date of verification:** 2026-09-27
**Sources:** NCBI E-utilities (`eutils.ncbi.nlm.nih.gov/entrez/eutils`) and NCBI Datasets
REST API v2/v2alpha (`api.ncbi.nlm.nih.gov/datasets`). All requests carried
`tool=viralscan&email=…` and were rate-limited to ≥0.4 s intervals.
**Databases queried:** `nuccore` (GenBank/RefSeq INSDC), `biosample`, `sra`, `taxonomy`.

**Status vocabulary**

| Status | Meaning |
|---|---|
| VERIFIED | Accession exists; organism, length, title and requested qualifiers read live from the record |
| VERIFIED (partial) | Accession exists and the requested qualifiers check out, but the record has a caveat noted in the row |
| NOT-FOUND | No such record found in `nuccore` within the searches actually run |
| UNVERIFIABLE | Cannot be established from INSDC/NCBI records; the label is not carried by NCBI |

> Nothing in this file is filled in from memory. Where a value could not be read
> from a live NCBI record it is written as **NOT-FOUND** or **UNVERIFIABLE** and
> the searches that produced the negative are stated.

---

## Two facts that shape the whole exercise

**1. Pango lineage is *not* a Pango attribute in INSDC, and Entrez cannot filter on it.**
`nuccore` exposes only 34 Entrez fields and none of them is `pango_lineage`. Two
traps were found and ruled out by direct experiment:

| Attempted query | Result | Why it fails |
|---|---|---|
| `pango_lineage[BA.3.2]` in `nuccore` | 5 hits | Entrez **silently drops the bracketed value**; it translates to bare `pango_lineage[All Fields]` and returns the same 5 records for *every* lineage queried |
| `pango_lineage[BA.3.2]` in `biosample` | 0 | The BioSample attribute is named **`lineage`**, not `pango_lineage` |
| `lineage[BA.3.2]` in `biosample` | 129 440 (constant) | Same bracketed-value drop; constant for every lineage queried |
| `filters.pangolin_classification=…` on `/virus/taxon/2697049/dataset_report` | ignored | Returns `total_count=9215231` for every value, including deliberately bogus filter names. (The `genome` endpoint *does* honour filters — `filters.assembly_level` cuts 12 413 → 2 733 — so the convention is right, this endpoint just does not implement it.) |

The only machine-readable Pango label available per INSDC accession is
`virus.pangolin_classification` from `GET /datasets/v2alpha/virus/accession/{acc}/dataset_report`.
This is **computed by NCBI**, not deposited: for `PP848088.1` (labelled `XBB.1.5.18`)
the GenBank flatfile contains no occurrence of "pango" or "lineage" anywhere.
The `completeness` field of that endpoint is also unreliable as a length proxy —
`PP848088.1` is 29 749 bp (a full genome) yet reports `PARTIAL`, while the 3 810 bp
record `PZ120670.1` also reports `PARTIAL`. Length was therefore used as the gate,
not `completeness`.

**2. Fine-grained Influenza clade labels are absent from INSDC records.**
Occurrence of the clade string in the whole `nuccore` text index, restricted to
`"Influenza A virus"[Organism] AND hemagglutinin[Title]`:

| Clade label | Records matching |
|---|---|
| `2.3.4.4b` | 1471 |
| `2.3.4.4` | 547 |
| `2.3.4` | 2288 |
| `6B` | 10 (all coincidental — see below) |
| `3C.3a` | 4 (all coincidental) |
| `B3.13` | 12 (genuine) |
| `D1.1` | 16 (genuine) |
| `6B.1`, `2a.3a.1`, `3C.2a1b.1a` | **0** |

The `6B` and `3C.3a` "hits" are false positives — they are isolate names, not clade
labels (e.g. `MF807096.1` is `A/Canada/6B/2014(H3N2)`, whose `6B` is a place token;
the `3C.3a` hits are `A/Czech Republic/56/2015(H3N2)` etc.). Consequently the
H1N1/H3N2 sub-clade rows below are marked **UNVERIFIABLE** in the *clade* column
while the accession/isolate/subtype/segment facts are fully VERIFIED. The H5N1
genotype rows (B3.13, D1.1) *are* supportable, and are backed by sequence identity
computed here rather than by any label in the record.

---

## Task 1 — SARS-CoV-2 representative genomes (INSDC / open access)

### Corpus counts

| Quantity | Value | Source |
|---|---|---|
| `"SARS-CoV-2"[Organism]` records in `nuccore` | **9 218 691** | esearch, live |
| … of which RefSeq (`srcdb_refseq[PROP]`) | **1** | esearch, live |
| The single RefSeq record | `NC_045512.2` (UID 1798174254) | esummary, live |
| `virus` records for taxon 2697049 (Datasets) | 9 215 231 | dataset_report `total_count` |
| SARS-CoV-2 `biosample` records | not separately counted; `JN.1` alone = 21 693 | esearch, live |

The entire RefSeq representation of SARS-CoV-2 in `nuccore` is **one** record.

### Verified representatives

`Pango label` is `virus.pangolin_classification` (NCBI-computed), read live.
`Deposited pango attr` records whether a depositor-supplied Pango label was
actually retrievable.

| Lineage | Accession | Pango label (NCBI) | Len (bp) | Isolate / strain | Collection date | Geo | Deposited pango attr | Status |
|---|---|---|---|---|---|---|---|---|
| JN.1 | `PP832909.1` | JN.1.45 | 29 776 | SARS-CoV-2/human/MAR/IPM404-2024/2024 | 2024-01-04 | Morocco | no | VERIFIED |
| JN.1 | `PP357841.1` | JN.1.55 | 29 767 | SARS-CoV-2/human/FRA/045-0011-D-P_2401V120213/2024 | 2024-01-05 | France | no | VERIFIED |
| XFG | `PZ323188.1` | **XFG** (exact) | 29 764 | SARS-CoV-2/human/USA/PHL2-M-A-A07-20260417/2026 | 2026-03-30 | USA: Pennsylvania, Philadelphia | yes — `EPI_ISL_20434955; Lineage: XFG.2` in DBLINK | VERIFIED |
| XFG | `PZ056993.1` | XFG.23.1.3 | 29 763 | SARS-CoV-2/human/USA/PHL2-M-B-H08-20260218/2026 | 2026-02-02 | USA: Pennsylvania, Philadelphia | no | VERIFIED |
| **BA.3.2** | `OZ536352.1` | **BA.3.2.2** | 29 872 | *(none — "genome assembly, complete genome: monopartite")* | *(none)* | United Kingdom | no | VERIFIED (partial) |
| NB.1.8.1 | `PZ201453.1` | **NB.1.8.1** (exact) | 29 598 | SARS-CoV-2/human/USA/FL-LC1141514/2025 | 2025-08-20 | USA: Florida | **yes** — BioSample `lineage` = NB.1.8.1 | VERIFIED |
| NB.1.8.1 | `PZ201531.1` | **NB.1.8.1** (exact) | 29 690 | SARS-CoV-2/human/USA/PA-LC1141738/2025 | 2025-08-20 | USA: Pennsylvania | **yes** — BioSample `lineage` = NB.1.8.1 | VERIFIED |
| BA.4.6 | `OR325407.1` | **BA.4.6** (exact) | 29 660 | SARS-CoV-2/human/FRA/IHUCOVID-095741/2022 | 2022-10 | France | no | VERIFIED |
| XBB.1.5 | `PP848088.1` | XBB.1.5.18 | 29 749 | SARS-CoV-2/human/FRA/IHUCOVID-102576/2023 | 2023-01 | France | no | VERIFIED (partial) |
| XBB.1.5 | `PP846642.1` | **XBB.1.5** (exact) | 29 755 | SARS-CoV-2/human/FRA/IHUCOVID-S0584/2023 | 2023-07 | France | no | VERIFIED (partial) |
| Wuhan-Hu-1 | `NC_045512.2` | *(not exposed)* | 29 903 | Wuhan-Hu-1 | — | China (Wuhan) | no | VERIFIED |
| **PQ.16.1.1** | — | — | — | — | — | — | no | **NOT-FOUND** |

Caveats on the "(partial)" rows, stated precisely:

- **`OZ536352.1`** carries no `/isolate` and no `/collection_date`; it is a
  *genome assembly* rather than a single-segment clinical isolate. Its DEFINITION
  does contain the literal string "complete genome", and Datasets reports
  `completeness=COMPLETE`. It is the **only** BA.3.2-labelled record found.
- **XBB.1.5** rows: `PP848088.1` DEFINITION is the long
  `…ORF1ab polyprotein (ORF1ab), … complete cds; …` enumeration form and therefore
  does **not** contain the literal phrase "complete genome", although it is a
  29 749 bp whole genome. `"XBB.1.5"[All Fields] AND "complete genome"[Title]`
  returns **0** — *no* XBB.1.5 record in `nuccore` carries that literal phrase.

### PQ.16.1.1 — documented negative

Four independent searches, all zero:

| Search | Result |
|---|---|
| `"PQ.16.1.1"[All Fields]` in `nuccore` (no organism restriction) | 0 |
| `"PQ.16.1.1"[All Fields]` / `[ATTR]` / `[WORD]` / `[TITL]` in `biosample` | 0 / 0 / 0 / 0 |
| `"PQ.16.1.1"[All Fields]` in `sra` | 0 |
| `pangolin_classification` in 1 000-record Datasets pages | 0 of ~901 000 records scanned (release dates 2023-07 → 2026-09; 4 523 distinct Pango labels observed in that window, including `XFG.*`, `PQ.17`, `PY.1.1.1`, `KP.3.1.1`, `XEC`, `XFZ`, `RK.1`) |

The Datasets listing is ordered by accession descending (not by date), so the scan
covered a contiguous, most-recent-first ~10 % of the corpus. This is a **bounded**
negative, not proof of global absence: **`PQ.16.1.1` is NOT-FOUND in the window
searched**, and no Pango label was inferred for it.

### BA.3.2 — how it was found, given the negative search

`"BA.3.2"[All Fields]` in `nuccore` = 0, and the 3 `biosample` hits are
unrelated (a *Biomphalaria* snail and two non-coronavirus samples). The record was
recovered by paging the Datasets `/virus/taxon/2697049/dataset_report` listing
(page cap 1 000 records/request, 429/backoff-aware) and matching
`pangolin_classification` client-side. This is slow — ~2.3 s per 1 000 records —
and is the reason BA.3.2 yields only one representative.

### Deposited-`lineage` vs NCBI-computed disagreement (1 observed)

`PZ201583.1` was reached through BioSample `SAMN56708780` whose depositor attribute
reads `lineage = NB.1.8.1`, but NCBI's own `pangolin_classification` for it is
**`PQ.17`**. For a lineage-labelled reference panel this is a live hazard: the
depositor-declared and NCBI-computed labels disagree. `PZ201583.1` was **excluded**
from the panel above for that reason.

---

## Task 2 — Influenza A/B clade representatives

### 2a. The critical check: `OP212288` — user hypothesis REFUTED

The user believed `OP212288` is *A/Texas/61/2022* HA, clade 2.3.4.4b genotype
B3.13 (H5N1), "deposited with (H3N2) in its organism field".

**Finding: the record is genuinely H3N2, and it is not an organism-field error.**
It is a name collision. Verified three independent ways:

| Line of evidence | Result |
|---|---|
| ORGANISM field (source feature) | `Influenza A virus`, `db_xref=taxon:11320` — **no subtype in the organism field at all**. The `(H3N2)` appears in the *DEFINITION/isolate string*, not in ORGANISM. |
| `/serotype` qualifier | `H3N2` |
| CDC/WHO `##FluData##` block | `NAME :: A/Texas/61/2022`, **`TYPE :: H3`**, `LOCATION :: United States / Texas`, `COLLECT_DATE :: 20-Apr-2022`, `SENDER_LAB :: DC Public Health Lab`, `EPI_ISL_13284199` — the depositing authority itself assigned type **H3**. |
| CDS translation length | **567 aa** = canonical H3 HA (H5N1 2.3.4.4b HA is 579–580 aa) |
| Protein identity vs references | **87.7 % / 82.2 %** to H3N2 HA (`AF251419.1`, `AY531037.1`); **42.4 % / 42.7 %** to H5N1 HA (`AF216729.1`, `AF216737.1`) |
| HA0 cleavage motif (aa 312–333) | query `NVNRITYGACPRYVKQSTLKLATGMRNVPEKQTR` — **monobasic**, matching H3N2 ref `NVNRITYGACPKYVKQKTLKLATGMRNVPEKQTR`; H5N1 ref is polybasic `PLTIGECPKYVKSNRLVLATGLRNTPQRERRRKK` |

Whole-isolate check — **all 8 segments** of `A/Texas/61/2022` are `H3N2`:

| Segment | Accession | `/serotype` |
|---|---|---|
| PB2 | `OP212285.1` | H3N2 |
| PB1 | `OP212286.1` | H3N2 |
| PA | `OP212287.1` | H3N2 |
| **HA** | **`OP212288.1`** | **H3N2** |
| NP | `OP212289.1` | H3N2 |
| NA | `OP212290.1` | H3N2 |
| M | `OP212291.1` | H3N2 |
| NS | `OP212292.1` | H3N2 |

So: there is **no organism-field/subtype mismatch** on `OP212288`; the `(H3N2)`
label is correct and independently corroborated by sequence. The name
*A/Texas/61/2022* is shared with the well-known H5N1 clade 2.3.4.4b human case,
which is the trap.

### 2b. Second name collision found: A/Brisbane/59/2007

A/Brisbane/59/2007 is a widely cited **H3N2** reference, but the `nuccore` records
are **H1N1** (`/serotype="H1N1"`, 8/8 segments, `08-Jul-2007`/`01-Jan-2008`).
Recorded here because a naive subtype→isolate table would have mis-assigned it.

### 2c. Verified complete segment sets

All rows are VERIFIED: organism, serotype, host, geo, collection date and segment
list were read live from the source features of the cited accessions.

| Group | Isolate | Subtype (record) | Host | Geo | Collection date | Segments | Clade in record? |
|---|---|---|---|---|---|---|---|
| H1N1pdm09 2009 pandemic | A/California/07/2009 | H1N1 | Homo sapiens (M, 54) | USA: California | 2009-04-09 | **8/8** | no |
| H3N2 | A/Perth/16/2009 | H3N2 | Homo sapiens (M, 6) | Australia: Western Australia | 2009-04-07 | **8/8** | no |
| H1N1 | A/Indiana/09/2016 | H1N1 | Homo sapiens | USA | 2016-01-23 | **8/8** | no |
| H5N1 human (B3.13) | A/California/147/2024 | H5N1 | Homo sapiens | USA: California | 2024-10-07 | **8/8** | no |
| H5N1 human (D1.1) | A/Washington/239/2024 | H5N1 | Homo sapiens | USA: Washington | 2024-10-18 | **7/8 — PB2 absent** | no |
| H7N9 human | A/Changsha/58/2017 | H7N9 | Homo sapiens (male) | China | 2017-03-17 | **8/8** | no |
| H9N2 human | A/Guangdong/GZ179/2016 | H9N2 | Homo sapiens | China | 2016-07-29 | **8/8** | no |
| H9N2 swine | A/swine/Jiangsu/C1/2008 | H9N2 | swine | China | 2008-10-02 | **8/8** | no |
| H1N2 swine | A/swine/Mexico/AVX61/2013 | H1N2 | Sus scrofa scrofa | Mexico | 2013-09-30 | **8/8** | no |
| H2N2 | A/Mallard Duck/MI/23-034361-019-original/2023 | H2N2 | Mallard Duck | USA: MI | 2023-10-31 | **8/8** | no |
| Flu B | B/Alaska/12/2015 | Influenza B virus | Homo sapiens | USA | 2015-04-30 | **8/8** | no |
| Flu B | B/Yamagata/16/1988 | Influenza B virus | Homo sapiens | Japan | 1988 | **8/8** | no |

**Segment accessions (one row per accession, as requested):**

| Isolate | Subtype | Seg | Accession |
|---|---|---|---|
| A/California/07/2009 | H1N1 | PB2 | `NC_026438.1` |
| A/California/07/2009 | H1N1 | PB1 | `NC_026435.1` |
| A/California/07/2009 | H1N1 | PA | `PX715226.1` |
| A/California/07/2009 | H1N1 | HA | `NC_026433.1` |
| A/California/07/2009 | H1N1 | NP | `NC_026436.1` |
| A/California/07/2009 | H1N1 | NA | `NC_026434.1` |
| A/California/07/2009 | H1N1 | M | `NC_026431.1` |
| A/California/07/2009 | H1N1 | NS | `NC_026432.1` |
| A/Perth/16/2009 | H3N2 | PB2 | `KJ609203.1` |
| A/Perth/16/2009 | H3N2 | PB1 | `KJ609204.1` |
| A/Perth/16/2009 | H3N2 | PA | `KJ609205.1` |
| A/Perth/16/2009 | H3N2 | HA | `KP457178.1` |
| A/Perth/16/2009 | H3N2 | NP | `KJ609207.1` |
| A/Perth/16/2009 | H3N2 | NA | `KJ609208.1` |
| A/Perth/16/2009 | H3N2 | M | `KJ609209.1` |
| A/Perth/16/2009 | H3N2 | NS | `KJ609210.1` |
| A/Indiana/09/2016 | H1N1 | PB2 | `KX005726.1` |
| A/Indiana/09/2016 | H1N1 | PB1 | `KX005727.1` |
| A/Indiana/09/2016 | H1N1 | PA | `KX005728.1` |
| A/Indiana/09/2016 | H1N1 | HA | `KX005729.1` |
| A/Indiana/09/2016 | H1N1 | NP | `KX005730.1` |
| A/Indiana/09/2016 | H1N1 | NA | `KX005731.1` |
| A/Indiana/09/2016 | H1N1 | M | `KX005732.1` |
| A/Indiana/09/2016 | H1N1 | NS | `KX005733.1` |
| A/California/147/2024 | H5N1 | PB2 | `PQ468764.1` |
| A/California/147/2024 | H5N1 | PB1 | `PQ468760.1` |
| A/California/147/2024 | H5N1 | PA | `PQ468762.1` |
| A/California/147/2024 | H5N1 | HA | `PQ468757.1` |
| A/California/147/2024 | H5N1 | NP | `PQ468759.1` |
| A/California/147/2024 | H5N1 | NA | `PQ468763.1` |
| A/California/147/2024 | H5N1 | M | `PQ468758.1` |
| A/California/147/2024 | H5N1 | NS | `PQ468761.1` |
| A/Washington/239/2024 | H5N1 | PB1 | `PQ573552.1` |
| A/Washington/239/2024 | H5N1 | PA | `PQ573557.1` |
| A/Washington/239/2024 | H5N1 | HA | `PQ573554.1` |
| A/Washington/239/2024 | H5N1 | NP | `PQ573556.1` |
| A/Washington/239/2024 | H5N1 | NA | `PQ573555.1` |
| A/Washington/239/2024 | H5N1 | M | `PQ573553.1` |
| A/Washington/239/2024 | H5N1 | NS | `PQ573551.1` |
| **A/Washington/239/2024** | H5N1 | **PB2** | **NOT-FOUND** (no segment-1 record deposited for this isolate) |
| A/Changsha/58/2017 | H7N9 | PB2 | `MG572501.1` |
| A/Changsha/58/2017 | H7N9 | PB1 | `MG572502.1` |
| A/Changsha/58/2017 | H7N9 | PA | `MG572503.1` |
| A/Changsha/58/2017 | H7N9 | HA | `MF370259.1` |
| A/Changsha/58/2017 | H7N9 | NP | `MG572504.1` |
| A/Changsha/58/2017 | H7N9 | NA | `MF370260.1` |
| A/Changsha/58/2017 | H7N9 | M | `MF370261.1` |
| A/Changsha/58/2017 | H7N9 | NS | `MF370262.1` |
| A/Guangdong/GZ179/2016 | H9N2 | PB2 | `KX808593.1` |
| A/Guangdong/GZ179/2016 | H9N2 | PB1 | `KX808594.1` |
| A/Guangdong/GZ179/2016 | H9N2 | PA | `KX808595.1` |
| A/Guangdong/GZ179/2016 | H9N2 | HA | `KX808596.1` |
| A/Guangdong/GZ179/2016 | H9N2 | NP | `KX808597.1` |
| A/Guangdong/GZ179/2016 | H9N2 | NA | `KX808598.1` |
| A/Guangdong/GZ179/2016 | H9N2 | M | `KX808599.1` |
| A/Guangdong/GZ179/2016 | H9N2 | NS | `KX808600.1` |
| A/swine/Jiangsu/C1/2008 | H9N2 | PB2 | `KX867822.1` |
| A/swine/Jiangsu/C1/2008 | H9N2 | PB1 | `KX867823.1` |
| A/swine/Jiangsu/C1/2008 | H9N2 | PA | `KX867824.1` |
| A/swine/Jiangsu/C1/2008 | H9N2 | HA | `KX867825.1` |
| A/swine/Jiangsu/C1/2008 | H9N2 | NP | `KX867826.1` |
| A/swine/Jiangsu/C1/2008 | H9N2 | NA | `KX867827.1` |
| A/swine/Jiangsu/C1/2008 | H9N2 | M | `KX867828.1` |
| A/swine/Jiangsu/C1/2008 | H9N2 | NS | `KX867829.1` |
| A/swine/Mexico/AVX61/2013 | H1N2 | PB2 | `KU976668.1` |
| A/swine/Mexico/AVX61/2013 | H1N2 | PB1 | `KU976613.1` |
| A/swine/Mexico/AVX61/2013 | H1N2 | PA | `KU976799.1` |
| A/swine/Mexico/AVX61/2013 | H1N2 | HA | `KU976596.1` |
| A/swine/Mexico/AVX61/2013 | H1N2 | NP | `KU976879.1` |
| A/swine/Mexico/AVX61/2013 | H1N2 | NA | `KU976577.1` |
| A/swine/Mexico/AVX61/2013 | H1N2 | M | `KU976501.1` |
| A/swine/Mexico/AVX61/2013 | H1N2 | NS | `KU976522.1` |
| A/Mallard Duck/MI/23-034361-019-original/2023 | H2N2 | PB2 | `QB056508.1` |
| A/Mallard Duck/MI/23-034361-019-original/2023 | H2N2 | PB1 | `QB056509.1` |
| A/Mallard Duck/MI/23-034361-019-original/2023 | H2N2 | PA | `QB056510.1` |
| A/Mallard Duck/MI/23-034361-019-original/2023 | H2N2 | HA | `QB056511.1` |
| A/Mallard Duck/MI/23-034361-019-original/2023 | H2N2 | NP | `QB056512.1` |
| A/Mallard Duck/MI/23-034361-019-original/2023 | H2N2 | NA | `QB056513.1` |
| A/Mallard Duck/MI/23-034361-019-original/2023 | H2N2 | M | `QB056514.1` |
| A/Mallard Duck/MI/23-034361-019-original/2023 | H2N2 | NS | `QB056515.1` |
| B/Alaska/12/2015 | Influenza B | PB1 | `KT866400.1` |
| B/Alaska/12/2015 | Influenza B | PB2 | `KT866401.1` |
| B/Alaska/12/2015 | Influenza B | PA | `KT866402.1` |
| B/Alaska/12/2015 | Influenza B | HA | `KT866403.1` |
| B/Alaska/12/2015 | Influenza B | NP | `KT866404.1` |
| B/Alaska/12/2015 | Influenza B | NB/NA | `KT866405.1` |
| B/Alaska/12/2015 | Influenza B | M1/BM2 | `KT866406.1` |
| B/Alaska/12/2015 | Influenza B | NS/NEP | `KT866407.1` |
| B/Yamagata/16/1988 | Influenza B | PB1 | `OQ034469.1` |
| B/Yamagata/16/1988 | Influenza B | PB2 | `OQ034470.1` |
| B/Yamagata/16/1988 | Influenza B | PA | `OQ034471.1` |
| B/Yamagata/16/1988 | Influenza B | HA | `OQ034472.1` |
| B/Yamagata/16/1988 | Influenza B | NP | `OQ034473.1` |
| B/Yamagata/16/1988 | Influenza B | NB/NA | `OQ034474.1` |
| B/Yamagata/16/1988 | Influenza B | M1/BM2 | `OQ034475.1` |
| B/Yamagata/16/1988 | Influenza B | NS/NEP | `OQ034476.1` |
| B/Victoria/2/1987 | Influenza B | 6/8 — HA (seg 4) and NB/NA (seg 6) absent | — | — | — | 6/8 | — |

> **Note on Influenza B segment count.** The brief asked for 2 segments for type B.
> In `nuccore`, Influenza B is deposited as **8 segments** (PB1, PB2, PA, HA, NP,
> NB+NA, M1+BM2, NS+NEP) — as shown for both B rows above. No 2-segment-only
> convention was found.

### 2d. Genotype confirmation (B3.13 and D1.1) — sequence-based, not label-based

Neither genotype is written into the `PQ…` records, so HA protein identity was
computed here against records that *do* carry the genotype string.

| Comparison (HA CDS, global alignment) | Identity |
|---|---|
| `PQ468757.1` (A/California/147/2024, B3.13 human set) vs `PV602015.1` (A/swine/Kansas/ExpPig61-NS5DPI/2024, record contains `B3.13`) | **99.5 %** |
| `PQ468757.1` vs `OP212288.1` (H3N2) | 41.3 % |
| `PQ468757.1` vs `AF251419.1` (H3N2 ref) | 41.7 % |
| `PQ573554.1` (A/Washington/239/2024, D1.1 human set) vs `PZ350417.1` (A/Bovine/Wisconsin/P3b/2025, record contains `D1.1`) | **99.6 %** |
| `PQ573554.1` vs `PQ468757.1` (B3.13) | 97.9 % |

Conclusions: the user's **B3.13** attribution for `PQ468757`–`PQ468764` is strongly
supported (99.5 %), and the **D1.1** attribution for `PQ573551`–`PQ573557` is
strongly supported (99.6 %). The two genotypes are 97.9 % identical to each other,
so they are close but distinguishable — a single HA sequence comparison at ~98 %
identity is **not** sufficient to separate them; the 99.5/99.6 % figures are what
carries the call, and a fuller genotype marker panel would be needed for routine
discrimination.

### 2e. Requested accession ranges — both exist

| Range | Exists? | What they are |
|---|---|---|
| `PQ468757`–`PQ468764` (8 accessions) | **YES, all 8** | Complete 8-segment set of *A/California/147/2024* (H5N1), Homo sapiens, USA: California, 2024-10-07. HA = `PQ468757.1`. Matches the B3.13 hypothesis. |
| `PQ573551`–`PQ573557` (7 accessions) | **YES, all 7** | 7 of the 8 segments of *A/Washington/239/2024* (H5N1), Homo sapiens, USA: Washington, 2024-10-18. Matches the D1.1 hypothesis. **PB2 is missing from both the deposited range and `nuccore` for this isolate.** |

### 2f. Clade labels the records cannot support

| Requested clade | Status | Evidence |
|---|---|---|
| 3C.3a, 6B, 6B.1, 6B.2 (H1N1pdm09) | **UNVERIFIABLE** in INSDC | 0 in-record occurrences for 6B.1; the 6B / 3C.3a "hits" are isolate place-name tokens. A/California/07/2009 (8/8) and A/Indiana/09/2016 (8/8) are supplied as era-appropriate complete sets; their clade labels are **external designations, not read from NCBI**. |
| 3C.2a, 3C.2a1, 3C.2a1b.1a, 2a, 2a.3, 2a.3a.1, 2a.3a.1(J.2.4.1) (H3N2) | **UNVERIFIABLE** in INSDC | 0 in-record occurrences for `3C.2a1b.1a`, `2a.3a.1`. A/Perth/16/2009 (8/8) supplied as an era-appropriate H3N2 set. No 2a-era human H3N2 set was asserted. |
| 2.3.4, 2.3.4.4, 2.3.4.4b (H5N1) | Clade strings present (2288 / 547 / 1471 HA hits) but **not on the chosen isolates** | 2.3.4.4b/D1.1 and B3.13 assignments above rest on computed sequence identity. |
| Influenza B Victoria / Yamagata | **UNVERIFIABLE** in INSDC | `ORGANISM` is `Influenza B virus` (taxid 11520) for all B records; NCBI has **no separate taxid** for B/Victoria or B/Yamagata (taxonomy search returns no such node). A `Victoria` string match inside *A/Perth/16/2009* records proved to be free-text noise (0 literal tokens on re-check) and must not be used. `B/Yamagata/16/1988` (8/8) and `B/Alaska/12/2015` (8/8) are supplied with lineage as an **external designation**. |

### 2g. Isolate-designation probes that returned nothing

Tested and confirmed **absent from `nuccore`** (not guessed at):
`A/California/07/2009`'s companions aside, the following isolate designations return
0 records: `A/Shanghai/46682/2013`, `A/Hangzhou/469/2013`, `A/Hong Kong/3399/2015`,
`A/Vietnam/119/2004`, `A/Vietnam/06/2004`, `A/Cambodia/2/2013`, `A/Canberra/332/2019`,
`A/Singapore/21/1957`, `B/Melbourne/12/2015`, `B/Yamagata/12/2011`,
`B/Victoria/4/1989`, `B/Victoria/2570/2019`, `B/Shanghai/10/1993`, `B/Hong Kong/8/1995`.

Consequences: the **famous 2013 human H7N9 isolates are not in INSDC** under those
names (A/Changsha/58/2017 was substituted, 8/8, human), and **historical human
H2N2 (1957/1968) is not retrievable** — the only H2N2 records in `nuccore` are avian
and a 2023 mallard duck. `B/Victoria/2/1987` exists but has only 6/8 segments.

### 2h. Organism-field inconsistency affecting subtype searches

Non-human-adapted subtypes are **not** filed under `ORGANISM = Influenza A virus`
(taxid 11320). They carry `ORGANISM = "Influenza A virus subtype H7N9"` etc., under
separate serotype taxids (`H7N9 subtype` 333278, `H9N2 subtype` 102796,
`H1N2 subtype` 114728, `H5N1 subtype` 102793; **no H2N2 node exists**). This is why
`"Influenza A virus"[Organism] AND H7N9[Title] AND segment[Title]` returns **0** while
`txid333278[Organism:exp] AND segment[Title]` returns 6 315. Any panel builder that
filters on `Organism="Influenza A virus"` will silently miss H7N9/H9N2/H1N2 entirely.
Separately, `H7N9[Title]` alone is swamped by patent/synthetic-construct noise
(18/20 sampled hits were `synthetic construct`).

---

## Task 3 — the four endemic coronaviruses + companions

| Accession | Taxon | Organism (as filed) | Len (bp) | Title | Status |
|---|---|---|---|---|---|
| `NC_045512.2` | 2697049 | Severe acute respiratory syndrome coronavirus 2 | 29 903 | Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome | VERIFIED |
| `NC_004718.3` | 227984 | **SARS coronavirus Tor2** | 29 751 | SARS coronavirus Tor2, complete genome | VERIFIED |
| `NC_019843.3` | 1335626 | Middle East respiratory syndrome-related coronavirus | 30 119 | Middle East respiratory syndrome-related coronavirus isolate HCoV-EMC/2012, complete genome | VERIFIED |
| `NC_038294.1` | 1263720 | **Betacoronavirus England 1** | 30 111 | Betacoronavirus England 1 isolate H123990006, complete genome | VERIFIED |
| `NC_002645.1` | 11137 | Human coronavirus 229E | 27 317 | Human coronavirus 229E, complete genome | VERIFIED |
| `NC_005831.2` | 277944 | Human coronavirus NL63 | 27 553 | Human Coronavirus NL63, complete genome | VERIFIED |
| `NC_006213.1` | 31631 | Human coronavirus OC43 | 30 741 | Human coronavirus OC43 strain ATCC VR-759, complete genome | VERIFIED |
| `NC_006577.2` | 290028 | Human coronavirus HKU1 | 29 926 | Human coronavirus HKU1, complete genome | VERIFIED |
| `NC_028752.1` | **1699095** | **Camel alphacoronavirus** | 27 395 | Camel alphacoronavirus isolate camel/Riyadh/Ry141/2015, complete genome | VERIFIED |

All 9 VERIFIED — 9/9, no failures.

### The `NC_028752.1` / taxid 11137 hypothesis: REFUTED

The brief asked to confirm that `NC_028752.1` is filed under taxid 11137 alongside
HCoV-229E. It is **not**.

| Claim | Verified value |
|---|---|
| `NC_028752.1` taxid (esummary and flatfile agree) | **1699095** — Camel alphacoronavirus |
| `NC_002645.1` taxid | **11137** — Human coronavirus 229E |
| Records under `txid11137[Organism:exp]` | 1 409, all *Human coronavirus 229E* |

They are distinct species under distinct taxids (`Camel alphacoronavirus` 1699095 in
genus *Alphacoronavirus*; `Human coronavirus 229E` 11137). If the panel groups
`NC_028752.1` under taxid 11137, any taxid-based reference lookup will mis-bucket it.

Two further naming points worth carrying into the panel:

- SARS-CoV-1 (`NC_004718.3`) is filed as **"SARS coronavirus Tor2"** (taxid 227984),
  not as a *Severe acute respiratory syndrome coronavirus* name matching SARS-CoV-2.
- MERS: `NC_019843.3` is *Middle East respiratory syndrome-related coronavirus*
  (1335626), but `NC_038294.1` is filed as **"Betacoronavirus England 1"** (1263720),
  not under a MERS name. Both were listed in the brief as "MERS-CoV"; the second one
  carries no MERS label in NCBI at all.

---

## Task 4 — adjacent human ssDNA viruses

| Accession | Taxon | Organism (as filed) | Len (bp) | Host | Geo | Collection date | Title | Status |
|---|---|---|---|---|---|---|---|---|
| `KF764702.1` | 93678 | **TTV-like mini virus** | 2 853 | Homo sapiens | Netherlands | 1991 | TTV-like mini virus isolate D50, complete genome | VERIFIED |
| `NC_076339.1` | 93678 | TTV-like mini virus | 2 907 | Homo sapiens | Viet Nam | 2015-05-18 | TTV-like mini virus strain vzttmv4, complete genome | VERIFIED (RefSeq) |
| `NC_076166.1` | 93678 | TTV-like mini virus | 2 943 | Homo sapiens | China | 2016-01-01 | TTV-like mini virus isolate zhenjiang, complete genome | VERIFIED (RefSeq) |
| `OP549860.1` | 432261 | Torque teno midi virus | 3 213 | Homo sapiens | South Africa | 2017-04-11 | Torque teno midi virus isolate UC156V1_C016, complete genome | VERIFIED |
| `MN774988.1` | 432261 | Torque teno midi virus | 3 301 | Homo sapiens | Tanzania | 2015 | MAG: Torque teno midi virus isolate SAfiA-353-59, complete genome | VERIFIED |
| `KM593803.2` | 136966 | SEN virus | 3 840 | Homo sapiens | France | May-1997 | SEN virus strain HDMU-97, complete genome | VERIFIED |
| `KM593802.2` | 136966 | SEN virus | 3 840 | Homo sapiens | France | Feb-2013 | SEN virus strain HDMU-13, complete genome | VERIFIED |
| `NC_135454.1` | **3025750** | **Human circovirus 1** | 2 021 | Homo sapiens | France | 2022-01-12 | Human circovirus 1 strain Paris, complete genome | VERIFIED (partial) |
| `PP968832.1` | 3298505 | Human circovirus | 2 021 | Homo sapiens | Hong Kong | 2022-01-13 | circovirus strain HK P1, complete genome | VERIFIED |
| `OR905605.1` | 3298505 | Human circovirus | 2 024 | Homo sapiens | Switzerland | 2022-10-05 | Human circovirus strain Liestal, complete genome | VERIFIED |
| `OZ282268.1` | 3298505 | Human circovirus | 2 024 | *(no /host)* | Switzerland: BL | 2022-10-05 | Human circovirus isolate Circovirus/Switzerland/BL-KSBL-Path-M2022.1259/2022 genome assembly, **chromosome: 1** | VERIFIED (partial) |

All 11 VERIFIED (2 with caveats below). Search yields: 2 482 TTMV complete genomes
with `Homo sapiens` host, 63 TTMDV, and exactly **2** SEN-virus complete genomes.

Caveats, stated precisely:

- **TTMV organism name.** The organism in `nuccore` is **"TTV-like mini virus"**
  (taxid 93678), not "Torque teno mini virus". Searching the organism string
  `Torque teno mini virus` still resolves (Entrez matches the lineage), but any
  exact-string organism filter must use `TTV-like mini virus`. TTMDV, by contrast,
  *is* filed as "Torque teno midi virus" (432261).
- **`NC_135454.1` is a different species name from the other three.** It is
  **"Human circovirus 1"** (taxid **3025750**), whereas `PP968832.1`,
  `OR905605.1` and `OZ282268.1` are **"Human circovirus"** (taxid **3298505**).
  The brief listed all four together; they are not one taxon in NCBI.
- **`OZ282268.1` does not say "complete genome".** Its DEFINITION is
  `…genome assembly, chromosome: 1.` It is a single-chromosome assembly at
  2 024 bp and carries no `/host` qualifier, so it fails a strict
  "complete genome + human host" filter while `OR905605.1` (same country, same
  date, 2 024 bp, complete genome, human host) is the cleaner panel entry.
- **SEN virus (taxid 136966)** has 303 records total, but only **2** are complete
  genomes (`KM593803.2`, `KM593802.2`). The other 301 are partial
  `ORF1 gene, partial cds` clones. The brief asked for 2–3 accessions; only 2
  qualifying ones exist.

---

## Summary for panel design

**Verified:** 11 SARS-CoV-2 accessions (5 of 7 requested lineages, plus
Wuhan-Hu-1), ~95 Influenza A/B segment accessions across 12 complete isolates,
9/9 endemic-coronavirus accessions, 11/11 ssDNA accessions.

**Failed / not found:**

| Item | Status |
|---|---|
| `PQ.16.1.1` representative | **NOT-FOUND** (0 in nuccore/biosample/sra; 0 in ~901 k Datasets records spanning 2023-07→2026-09) |
| `PQ5735xx` PB2 | **NOT-FOUND** (A/Washington/239/2024 is 7/8) |
| `6B.1`, `3C.2a1b.1a`, `2a.3a.1` clade labels | **UNVERIFIABLE** — absent from INSDC |
| Influenza B Victoria/Yamagata lineage | **UNVERIFIABLE** — no NCBI taxid, no in-record label |
| Historical human H2N2 (1957/1968) | **NOT-FOUND** — only avian/2023-duck H2N2 in INSDC |
| 2013 human H7N9 isolates (`A/Shanghai/46682/2013` etc.) | **NOT-FOUND** — substituted A/Changsha/58/2017 |
| BA.3.2, second representative | Only 1 record found (`OZ536352.1`); no isolate/date qualifiers on it |

**Surprises worth acting on:**

1. **`OP212288` is H3N2, not H5N1 — REFUTED.** There is no organism-field error:
   `ORGANISM` is bare "Influenza A virus" (taxid 11320); the `(H3N2)` sits in the
   DEFINITION isolate string, and the depositing CDC/WHO centre's own `##FluData##`
   block says `TYPE :: H3`. Confirmed by 567 aa CDS length, 87.7 % identity to H3N2
   vs 42.4 % to H5N1, and a monobasic (H3-like) HA0 cleavage site. All 8 segments
   `OP212285`–`OP212292` are `/serotype="H3N2"`. The isolate name collides with the
   H5N1 clade 2.3.4.4b Texas case.
2. **A/Brisbane/59/2007 is H1N1 in `nuccore`,** not the H3N2 its citation implies.
3. **`NC_028752.1` is taxid 1699095 (Camel alphacoronavirus),** not 11137; taxid
   11137 is HCoV-229E alone.
4. **`NC_135454.1` is "Human circovirus 1" (3025750),** a different species name
   from the other three circovirus records (3298505).
5. **Depositor `lineage` and NCBI `pangolin_classification` can disagree**
   (`PZ201583.1`: depositor NB.1.8.1 vs NCBI PQ.17). Pin which one the panel trusts.
6. **Entire SARS-CoV-2 RefSeq representation in `nuccore` is one record**
   (`NC_045512.2`) out of 9 218 691.
7. **Subtype searches keyed on `Organism="Influenza A virus"` silently miss
   H7N9/H9N2/H1N2** — those use distinct serotype taxids.
8. **Entrez bracketed-field syntax silently drops its value**
   (`pango_lineage[X]`, `lineage[X]`), returning the same result for every `X`.
   This is the single easiest way to draw a false conclusion from NCBI.
9. `XBB.1.5` has **no** `nuccore` record with the literal phrase "complete genome"
   in the title; whole-genome XBB.1.5 entries use the "…complete cds…" enumeration.
