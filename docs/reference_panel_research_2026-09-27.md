---
title: Viral reference panel design — influenza and human respiratory viruses
subtitle: Verified literature and NCBI-resource research for ViralScan
date: 2026-09-27
status: research findings, not an implemented design
---

# Viral reference panel design — influenza and human respiratory viruses

Research note for ViralScan. Everything below is either (a) fetched from a
primary source, (b) queried live against NCBI E-utilities on **2026-09-27**, or
(c) explicitly flagged as unverified. No accession list here was written from
memory.

**Verification legend**

| Tag | Meaning |
|-----|---------|
| `[FETCH]` | I downloaded and read the primary page/PDF/XML myself |
| `[EBI]` | Queried live against NCBI E-utilities on 2026-09-27; reproducible with the query given |
| `[SNIP]` | I could only retrieve the page through the search index (page blocked direct fetch); treat details as second-hand |
| `[UNVERIFIED]` | Claimed by a secondary source but not confirmed here — must be checked before it enters a panel manifest |

---

## 0. Why this matters for the current panel

`src/viralscan/data/` today contains, for influenza, exactly the RefSeq set:

* `Influenza_A_virus_NC_002016…NC_002023` = **A/Puerto Rico/8/1934(H1N1)**, 8 segments
* `Influenza_B_virus_NC_002204…NC_002211` = **B/Lee/1940**, 8 segments
* `Influenza_C_virus_NC_006306…NC_006312` = **C/Ann Arbor/1/50**, 7 records

That is **one IAV strain from 1934**, one IB strain from 1940, and one ICV
strain from 1950. There is no A(H1N1)pdm09, no A(H3N2) after 2004, no B/Victoria,
no B/Yamagata, and no avian subtype beyond what §4 shows RefSeq holds. The panel
is functionally a genome-completeness fixture, not a detection panel.

---

## 1. IAV subtype situation (current epidemiology)

### 1.1 Human seasonal subtypes — unchanged since 1977

Human seasonal influenza today is **A(H1N1)pdm09, A(H3N2), and B/Victoria**.
A(H2N2) caused the 1957–1968 pandemics and has not circulated in humans since
the 1968–69 season `[SNIP: CDC trivalent-vaccine page, cdc.gov/flu/vaccine-types/trivalent.html]`.
A(H1N1)pdm09 has been the sole H1N1 in humans since 2009, replacing the
pre-2009 seasonal H1N1 lineage and displacing classical swine H1N1 from the
human slot `[FETCH: WHO fact sheet 2026-01-21; CDC archive]`.

**Recommended panel treatment:** H1N1pdm09 and H3N2 are *mandatory and
time-varying*; H2N2 is *historic* (include for archival deep-sequencing projects
only, e.g. paleogenomics of 1957/1968 pandemics — otherwise omit).

### 1.2 Avian/swine subtypes with verified human infections

`[FETCH: WHO "Influenza (avian and other zoonotic)" fact sheet, updated 21 January 2026]`
lists: H5N1, H5N6, H5N2, H5N8, H6N1, H7N9, H7N2, H7N3, H7N4, H7N7, H9N2,
H10 (and, for swine, H1 and H3).

`[SNIP: CDC "Reported Human Infections with Avian Influenza A Viruses", updated
2026-03-20]` — direct fetch returned HTTP 403, so the following is from the
search index of that page:

* LPAI subtypes with virologically confirmed human infections: **A(H3), A(H5),
  A(H6), A(H7), A(H9), A(H10)**.
* A(H3N8) — small numbers of people in China since 2022 (mild to critical).
* A(H5N2) — one person, hospitalised and died, Mexico 2024.
* A(H5N6) — sporadic, mostly China, since 2014; **93 cumulative human cases
  as of 17 Jan 2026, 92 in the Chinese mainland** `[SNIP: CHP Hong Kong Avian
  Influenza Report vol.22 wk03]`.
* A(H5N8) — one asymptomatic poultry worker, Russia 2020.
* A(H6N1) — one person, moderate LRTD, Taiwan 2013.
* A(H7N9) — >1,500 in China 2013–2017; **1,568 confirmed, 616 fatal (CFR 39%),
  last case 5 Apr 2019** `[SNIP: WHO WPRO ai_20250530.pdf; CHP AIR vol.21 wk51;
  WHO fact sheet 2026-01-21]`.
* A(H9N2) — >100 people since 1998 across China, HK, Bangladesh, Cambodia,
  Egypt, India, Oman, Pakistan, Senegal, Vietnam. **117 cumulative by 9 May 2025**
  `[SNIP: EID, Changsha H9N2 report, PMC12620571]`.
* A(H10N3) — two people, severe pneumonia, China 2021 and 2022; **7 cumulative
  since 2021 as of Jan 2026** `[FETCH: WHO human-animal interface, 23 Jan–31 Mar 2026 PDF]`.
* A(H10N5) — one person, died, China late 2023 `[SNIP: CDC page]`.
* A(H10N7) — small numbers, Egypt 2004 and Australia 2010 `[SNIP: CDC page]`.

**ECDC, most recent overview (28 Feb – 4 Jun 2026)** `[FETCH]`:
19 human cases in six countries/territories — Bangladesh 2×A(H5N1) (1 fatal),
Cambodia 3×A(H5N1) (1 fatal), India 1×A(H5N1), Italy 1 imported A(H9N2)
(*first H9N2 human case in the European Region*), China 10×A(H9N2) +
1 fatal A(H5N6), Taiwan 1×A(H7N7). First detection of A(H9N2) **clade G5.5** in
European poultry. A(H5N5) detected in a polar bear and a walrus in Norway. No
sustained human-to-human transmission documented.

### 1.3 IAV H1N2 — has it "arisen" in humans?

**Yes, three separate times, but never as a sustained human seasonal subtype.**

1. **A(H1N2) seasonal, 1988–89** — sporadic cases over one winter in China
   `[EBI/PMC5176240: Emerging Infect Dis, whole-genome characterisation of a
   novel human A(H1N2) variant, Brazil]`.
2. **A(H1N2) seasonal, 2000–2003** — emerged in the human population, became
   widespread in Europe, sporadic cases through Europe/Asia/Africa/Americas
   2001–2003. Reassortant of human seasonal H1N1 and H3N2 `[EBI/PMC5176240]`.
3. **A(H1N2)v (swine-origin), ongoing sporadic human cases** — US Michigan
   agricultural fair Aug 2023; UK first case Nov 2023 (clade 1b.1.1); Brazil
   2015 (H1N2v); ~50 human H1N2v cases globally since 2005
   `[SNIP: cdc.gov/swine-flu; BMJ 2023;383:p2819; PMC5176240]`.
   **Most recent: China notified WHO of a laboratory-confirmed human
   A(H1N2)v case on 3 February 2026** `[FETCH: WHO human-animal interface
   summary and assessment, 23 January – 31 March 2026, PDF, line 260]`.

The same WHO report also records **A(H1N1)v** (China, 20 Mar 2026) and
**A(H3N2)v** (Brazil, 26 Jan 2026) `[FETCH, same PDF, lines 254–292]`.

**Recommended panel treatment:** A(H1N2) is not a human-seasonal subtype today.
Include H1N2v references only if the tool's remit includes swine-variant
detection. Document explicitly that a H1N2 hit is *not* evidence of
pandemic-risk reassortment.

---

## 2. IAV lineages that matter for a reference panel

### 2.1 A(H1N1)pdm09 — HA clade/subclade

The 3C/6B lineage is dead for pdm09; operational naming is now:

| Level | Name | Notes |
|---|---|---|
| Clade | `5a.2a` (alias `C.1`) | H1N1pdm09 since ~2014 |
| Clade | `5a.2a.1` (alias `D`) | dominant since 2019 |
| Subclade | `D.2` | +R113K; 19% of 5a.2a.1 in 2024–25 EU/EEA `[SNIP: ECDC wk40 2024–wk33 2025]` |
| Subclade | `D.3` | 69% of 5a.2a.1; carries T120A; **no reference strain** `[SNIP: same]` |
| Subclade | `D.3.1` | NH 2025-26 and **2026-27** CVV A/Missouri/11/2025 `[SNIP: FDA VRBPAC 12 Mar 2026 deck]` |
| Subclade | `D.3.1.1` | co-dominant in 2024–25 `[SNIP: FDA VRBPAC 12 Mar 2026 deck]` |

`[SNIP: FDA VRBPAC March 12 2026]` — "the vast majority of HA genes … belonged
to clade 5a.2a1 subclades D.3.1 and D.3.1.1." Serology used
5a.2a (C.1.9.3), 5a.2a.1 (D.3.1, D.3.1.1).

### 2.2 A(H1N1) pdm09 vs pre-2009 seasonal vs classical swine

* **pdm09** — triple reassortant: Eurasian avian PB2/PB1/PA, North American
  classical-swine H1, N1, M, NS `[SNIP: cdc.gov archive 2009 Q&A]`. GenBank
  reference set: `NC_026431…NC_026438` (A/California/07/2009) — RefSeq, all
  eight segments `[EBI]`.
* **Pre-2009 seasonal H1N1** — the 6B/6B.1 lineage; extinct as a human
  subtype but historically important. GenBank segment sets exist
  (e.g. A/Brisbane/59/2007 = `CY064971`–`CY064978`; A/Beijing/502/2009 =
  `GQ290105`–`GQ290112`) `[EBI]`.
* **Classical swine H1N1** — still enzootic; `A/Sw/1976/1976` is the prototypic
  strain. `[UNVERIFIED]` — I could not confirm a canonical classical-swine
  H1N1 full 8-segment accession set by name search; resolve from GISAID/IVSR
  before including.

### 2.3 A(H3N2) clade evolution — **the naming has changed and is actively changing**

This is the single most error-prone part of a panel manifest. Verified history:

| Older name | Current short name | Status |
|---|---|---|
| `3C.2a1b.2a` | — | pre-2020 |
| `3C.2a1b.2a.1` | — | 2021–22 |
| `3C.2a1b.2a.2` | `2a.3` (renamed **Feb 2023**) | `[SNIP: ECDC/WHO influenza virus characterisation, Mar 2024]` |
| `3C.2a1b.2a.2a` | `2a.3a` (E50K vs A/Norway/24873/2021) | `[SNIP: ECDC wk40 2024–wk33 2025]` |
| `3C.2a1b.2a.2a.3a.1` | `2a.3a.1` (alias **J**) | >99% of EU/EEA 2024–25 `[SNIP: ECDC wk40 2024–wk33 2025]` |
| `3C.2a1b.2a.2a.2a.3a.1` | `2a.3a.1` | `[SNIP: Melbourne WHO CC report 2023]` |
| `3C.2a1b.2a.2b` | `2a.3b` / subclade `G.1.3.1` | 15 viruses, minor `[SNIP: ECDC wk40 2024–wk33 2025]` |
| — | `2a.3a.1 (J.2)` | NH 2025-26: A/Croatia/10136RV/2023 (egg), A/District of Columbia/27/2023 (cell) `[SNIP: FDA VRBPAC 12 Mar 2026]` |
| — | `2a.3a.1 (J.2.3)`, `(J.2.4)`, **`(J.2.4.1)` = "K"** | K emerged **Aug 2025**, dominated 2025-26; 47% of EU/EEA GISAID submissions 1 May–17 Nov 2025; **NH 2026-27 CVV = A/Darwin/1454/2025 (egg), HA subclade K** `[SNIP: FDA VRBPAC 12 Mar 2026; ECDC TA Brief 2025; MDPI Pathogens 15(2):37]` |

**There is a *proposed* fourth naming system.** `[FETCH: Europe PMC full text,
PMC12904685, Influenza Other Respir Viruses, doi 10.1111/irv.70230]` — a GISRS
paper "Nomenclature for Tracking of Genetic Variation of Seasonal Influenza
Viruses" proposes algorithmic, dynamically-assigned short labels (e.g. `G.1.3`,
`C.1`, `J.2`, `K`) supported by Nextclade, with the subclade definitions in
supplementary Tables S1–S6. The main text confirms the operational GISRS
practice in 2026 is *both* the legacy `2a.3a.1 (J.2.4.1)` style and the short
letter/digit alias.

**Panel implication:** do **not** hard-code `3C.2a1b.2a.2a.3a.1` strings. Store
both a stable internal `h3n2_subclade` label and the as-of-date GISRS label, and
re-derive the manifest at each release.

### 2.4 A(H5N1) clade 2.3.4.4 and 2.3.4.4b

`[SNIP: J Virol 99(6) e00424-25, "Clade 2.3.4.4b … knowns, unknowns, and
challenges"]`:

* Clade 2.3.4.4 (Gs/Gd/96 lineage) has **eight subclades, 2.3.4.4a–2.3.4.4h**;
  **2.3.4.4b is the dominant subclade**.
* By end of 2015, clade 2.3.4.4 H5Nx viruses (N1–N9) replaced all earlier H5N1.
* 2.3.4.4b H5N8 was prevalent in wild birds until reassortment produced a novel
  H5N1 in 2020.
* >100 reassortant H5N1 genotypes in North America, >50 in Europe.

**2024–2025 genotype replacement inside 2.3.4.4b** `[SNIP: FAO/WHO/WOAH joint
public health assessment, Jul 2025; EID 31(12); EID 29(6)]`:

| Genotype | Description | Verified GenBank accessions `[EBI]` |
|---|---|---|
| **B3.13** | Eurasian 2.3.4.4b backbone with North American LPAI PB2/PB1/NP/NS. Caused the 2024 US dairy-cattle panzootic (1,074 herds / 17 states, Mar 2024–1 Jul 2025) and the US human cases. | A/California/147/2024 = `PQ468757`–`PQ468764` (8 seg); A/Texas/61/2022 = `OP212285`–`OP212292` (8 seg) |
| **D1.1** | Separate reassortant, primarily poultry; two severe human cases (BC Canada; fatal, Louisiana) | A/Washington/239/2024 = `PQ573551`–`PQ573557` (7 seg) + `PQ525415` |
| **B3.2, B3.6** | North America / South America (marine mammals) | `[UNVERIFIED]` accessions not resolved here |
| **EA-2021-AB, EA-2022-BB, EA-2023-DG** | European genotypes (Fusaro et al., *Virus Evolution* 2024, 1,956 European genomes) | `[UNVERIFIED]` |

**Other H5 clades that still matter for human surveillance:**
`[FETCH: WHO human-animal interface 23 Jan–31 Mar 2026]` and `[SNIP: Lam et al.,
bioRxiv 2025.11.23.690055]`:

* **2.3.2.1a** (Bangladesh 2026 case; A/Victoria/149/2024-like triple
  reassortant) — a live human clade.
* **2.3.2.1c** — human cases in China, Nepal, Cambodia, Indonesia, Vietnam.
* **2.3.2.1d**, **2.3.2.1e** — 2025 Cambodia cases identified as 2.3.2.1e.
* **2.3.4.4e, 2.3.4.4g, 2.3.4.4h** — East/SE Asia 2013–2020; human infections
  reported for most 2.3.4.4 subclades **except 2.3.4.4c and 2.3.4.4f**.

⚠️ **NCBI organism-string trap (verified).** `OP212288` — the segment-4 HA of
**A/Texas/61/2022**, the first human infection from a US dairy cow — is
deposited as `Influenza A virus (A/Texas/61/2022(**H3N2**))` `[EBI: efetch]`.
Its HA is clade 2.3.4.4b H5. Any panel builder that types by the `[Organism]`
field will file this critical reference under H3N2. Type by **sequence**, not by
the organism string. `[EBI]` also shows `"A/Texas/61/2022"[Organism]` returns
**0 hits** — the searchable organism string is the parenthesised form only.

### 2.5 H7N9, H9N2, H10N7, H3N8

* **H7N9** — `NC_026422`–`NC_026429` (A/Shanghai/02/2013, RefSeq, all 8
  segments) `[EBI]`. 1,568 human cases 2013–2019, CFR 39%, **no human case since
  5 Apr 2019**. Two waves (Shanghai 2013, Guangdong/HK 2014–2017) plus
  intermediate/pre-2013 precursors. `[UNVERIFIED]` — accessions for the second
  wave (e.g. A/Hong Kong/4801/2014 is H3N2, not H7N9; H7N9 wave-2 GenBank
  accessions not resolved here).
* **H9N2** — two things matter, and they are phylogenetically distinct:
  * the **G1 lineage** (prototype `A/quail/Hong Kong/G1/1997`), split
    G1-W ("Western") and G1-E ("Eastern") — most human isolates;
  * the **BJ94 / Y280 / G9 lineage** (prototype `A/chicken/Beijing/1/94`,
    also called `A/Duck/Hong Kong/Y280/97`) — a second stable poultry lineage.
  All human H9N2 isolates carry HA from G1-W, G1-E, or BJ94 `[SNIP: Viruses
  2019, "A Global Perspective on H9N2 AIV", PMC6669617]`. 2025 Changsha cases
  were Y280-like `[SNIP: PMC12620571]`.
  RefSeq: `NC_004905`–`NC_004912` (A/Hong Kong/1073/99, G1) `[EBI]`.
  GenBank: A/Duck/Hong Kong/Y280/97 = `AF156394`, `AF156405`, `AF156419`,
  `AF156433`, `AF156447`, `AF156461`, `AF156475` (partial CDS for several) `[EBI]`.
  `[SNIP: ECDC]` also reports a clade designation **`G5.5`** for A(H9N2) in
  European poultry (first detection there) — a newer Gx.y scheme. `[UNVERIFIED]`
  — the Gx.y H9N2 nomenclature source was not retrieved; do not use it in a
  manifest without locating the primary naming document.
* **H10N7** — human cases Egypt 2004 and Australia 2010 `[SNIP: CDC page]`.
  NCBI holds 7,445 records mentioning "H10N7" `[EBI]`, but these are dominated
  by wild-bird surveillance submissions (e.g. QB056286–QB056291,
  A/Green-winged Teal/UT/23-038097-001-original-repeat/2023(H10N7)) `[EBI]`.
  `[UNVERIFIED]` — the specific 2004 Egypt / 2010 Australia human-case genomes
  were not located.
* **H3N8** — LPAI, human cases in China since 2022 `[SNIP: CDC page]`.
  `[UNVERIFIED]` — human-case accessions not located; 4,011 records mention
  "H3N8" with "segment 4" in the title `[EBI]`, but they are wild-bird
  surveillance isolates.

---

## 3. Influenza B — is Yamagata gone?

**Verified: no confirmed B/Yamagata circulation since March 2020.**

* `[SNIP: CDC trivalent vaccine page, updated 2026-05-11]` — "There have been no
  confirmed detections of influenza B/Yamagata lineage viruses after March 2020.
  Reports of B/Yamagata detections after March 2020 for which samples were
  available were confirmed as naturally occurring B/Victoria lineage viruses or
  were identified as the B/Yamagata lineage component of live attenuated
  vaccines." US vaccines went trivalent from 2024-25. **"As of the 2026 Southern
  Hemisphere VCM there will no longer be updated recommendations for the
  B/Yamagata lineage component."** Explicitly: "it is not known whether
  influenza B/Yamagata lineage viruses are extinct at this time."
* `[FETCH: WHO candidate-vaccine-virus index page]` — the **February 2026-2027
  NH** recommendation lists only A(H1N1)pdm09, A(H3N2) and **B/Victoria**
  (egg + cell) — **no B/Yamagata**. B/Yamagata is still listed for the
  **September 2026 Southern Hemisphere** season (and for 2025, 2024, 2023,
  2022, 2021, 2020 back-compat entries).
* `[SNIP: CDC ACIP Influenza B/Yamagata Update, Feb 2024; FDA VRBPAC Feb 2024]`
  — "no confirmed influenza B/Yamagata viruses in global surveillance since
  March 2020"; WHO and FDA VRBPAC recommended excluding B/Yamagata.
* `[SNIP: Caini et al., Lancet Microbe 2024, "Probable extinction of influenza
  B/Yamagata and its public health implications"]` — the peer-reviewed statement
  of "probable", not "confirmed", extinction.
* `[SNIP: Euro Surveill 2022]` — post-2020 FluNet B/Yamagata reports trace to
  LAIV.

**B/Victoria HA subclades:** `V1A` (clade 1) → `V1A.1`, `V1A.2`, `V1A.3`
(clade 1A) → `V1A.3a` (1A) → `V1A.3a.2` (3a.2) / `V1A.3a.3` (3a.3). NH 2025-26
and 2026-27 CVVs are `V1A.3a.2` and `C.3.1` respectively
`[SNIP: FDA VRBPAC 12 Mar 2026 deck]`. B/Yamagata clades `Y1`, `Y2`, `Y3`.

**Recommended panel treatment:** B/Victoria **mandatory**, with ≥3 subclade
representatives (1A.3 / 3a.2 / 3a.3). **B/Yamagata optional-but-recommended** as
a *retained-informational* reference: it costs one genome, it is still the
component of every LAIV product used in any geography, and its absence would be
an unflagged behavioural change. Label it `retired_from_surveillance` in the
manifest rather than deleting it. This is a judgement call, not a literature
finding — state it as such in the methods.

RefSeq for IB is one strain: **B/Lee/1940**, `NC_002204`–`NC_002211`
`[EBI]`.

---

## 4. What NCBI actually holds (verified counts, 2026-09-27)

### 4.1 RefSeq is **not** a diversity resource

`[EBI]` `esearch db=nuccore term=txid11320[Organism:exp] AND srcdb_refseq[PROP]`:

> **Count = 56.** All 56 are individual **segments**, resolving to **7 strains**:

| Strain | Subtype | RefSeq accessions | Complete? |
|---|---|---|---|
| A/Puerto Rico/8/1934 | H1N1 | `NC_002016`–`NC_002023` | 8/8 segs |
| A/Hong Kong/1073/99 | H9N2 | `NC_004905`–`NC_004912` | segs 5,7,8 full; 1,2,3,4,6 CDS-only |
| A/Goose/Guangdong/1/1996 | H5N1 | `NC_007357`–`NC_007364` | 8/8 (segs 7,8 as "complete sequence"; 1–6 CDS) |
| A/New York/392/2004 | H3N2 | `NC_007366`–`NC_007373` | 8/8 segs |
| A/Korea/426/1968 | H2N2 | `NC_007374`–`NC_007382` | 8/8 segs |
| A/Shanghai/02/2013 | H7N9 | `NC_026422`–`NC_026429` | 8/8 segs |
| A/California/07/2009 | H1N1 (pdm09) | `NC_026431`–`NC_026438` | 8/8 segs |

Six subtypes, two of them H1N1, **no H1N2, no H5N6, no H5N8, no H10N7, no H6N1,
no H3N8, and no post-2013 pdm09 or 2.3.4.4b diversity.** This is the same set
ViralScan already ships.

Other RefSeq counts `[EBI]`: Influenza B = **8** (B/Lee/1940 segments
`NC_002204`–`NC_002211`); Influenza C = **7** (C/Ann Arbor/1/50
`NC_006306`–`NC_006312`); SARS-CoV-2 = **1** (`NC_045512.2`, Wuhan-Hu-1);
**Human bocavirus 1 = 0**.

RefSeq overall context `[SNIP: ncbi.nlm.nih.gov/refseq]`: Release 237, 184,752
organisms. RefSeq viral FASTA release file
`ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz`, 170 MB,
dated 2026-09-03 `[EBI: HTTP HEAD]`.

**Design consequence:** build the panel from **GenBank/IVSR**, not RefSeq. Use
RefSeq only where you want NCBI's curated gene models (which matters for the
GTF side of a kallisto reference).

### 4.2 GenBank volumes (verified, 2026-09-27)

Total records: IAV = **1,666,387** `[EBI]`. Counts of records whose full text
mentions the subtype string `[EBI]`, `"HxNy"[AllFields]`:

| Subtype | Records mentioning subtype | Subtype | Records |
|---|---|---|---|
| H1N1 | 452,423 | H9N2 | 70,987 |
| H3N2 | 594,895 | H3N8 | 34,197 |
| H5N1 | 228,507 | H5N2 | 16,041 |
| H7N9 | 18,048 | H5N8 | 12,791 |
| H1N2 | 49,717 | H6N1 | 7,025 |
| H5N6 | 7,504 | H10N7 | 7,443 |
| H2N2 | 3,149 | H3N2-adjacent minor types | — |

⚠️ **These counts overlap and are not genome counts.** A single reassortant
record can mention two subtype strings; `H1N1` also matches H1N10–H1N19;
`H5N2` counts exclude many H5N2 that appear only in feature tables. The
subtype sums exceed the 1,666,387 total precisely because of this overlap. Use
them only to order subtypes by data availability, never as a denominator.

### 4.3 How IVSR maps to a detection panel

`[FETCH: ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=genomeset]`
plus the linked database page. The page carries this banner:

> "Important Update — NCBI plans to redirect the Influenza Virus Resource to
> NCBI Virus. For most up-to-date and accurate virus data, go to NCBI Virus."

The **Influenza Virus Genome Set** tool is *by construction* a "one genome per
strain, ordered by segment" resource. Its controls map onto panel design as
follows:

| IVSR control | Panel-design meaning |
|---|---|
| Genome sets: **Complete only** (default) | keep only strains with all 8 segments present — this is your assembly unit |
| **Full-length only** (vs "Full-length plus") | drop partial-CDS segment records; mandatory, or your FASTA will have segment-level holes |
| Required segments 1–8 checkboxes | per-segment completeness control |
| Include/Exclude **Only Pandemic (H1N1) viruses** | lets you partition pdm09 from pre-2009 seasonal H1N1 — exactly the §2.2 split |
| Include/Exclude **Lab strains**, **Vaccine strains**, **Lineage defining strains**, **Mixed subtype** | curation filters; a good default is *exclude lab strains, include vaccine and lineage-defining* |
| **Collapse identical sequences** | NCBI's own redundancy filter — a check on any CD-HIT threshold you choose |
| Host = Avian / Human | splits the avian-reservoir block from the human block |

**Practical workflow:** use IVSR Genome Set to *enumerate* strains per subtype
and per host, per year; then fetch the eight segment accessions per strain via
E-utilities; concatenate in segment order 1→8 with segment namespacing; then
cluster with CD-HIT at your chosen threshold. Plan for the IVSR→NCBI
Virus/Datasets migration now; the CGI is session-based and will break.

Reproducible recipes `[EBI]`:

```bash
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
# RefSeq IAV records (segments), as counted above
curl -s "$E/esearch.fcgi?db=nuccore&rettype=count&retmode=json&term=txid11320%5BOrganism%3Aexp%5D+AND+srcdb_refseq%5BPROP%5D"
# Total IAV GenBank records
curl -s "$E/esearch.fcgi?db=nuccore&rettype=count&retmode=json&term=txid11320%5BOrganism%3Aexp%5D"
# A complete 8-segment strain: find segment 4, then siblings by strain string
curl -s "$E/esearch.fcgi?db=nuccore&retmode=json&retmax=20&term=%22A%2FShanghai%2F02%2F2013%22%5BAllFields%5D"
```

---

## 5. Other human respiratory viruses — verified RefSeq holdings and lineage lists

`[EBI]` RefSeq record counts and accessions (2026-09-27):

| Virus | RefSeq n | RefSeq accessions (complete genomes) | GenBank records |
|---|---|---|---|
| RSV (A + B) | **3** | `NC_001781` HRSV-B; `NC_001803` HRSV-A/UK/S2.ts1C/1995; `NC_038235` HRSV-A/USA/001/198X | 74,218 |
| Human metapneumovirus | **1** | `NC_039199` isolate 00-1 | 15,587 |
| HPIV-1 | 1 | `NC_003461` | 1,485 |
| HPIV-2 | 1 | `NC_003443` ("Human rubulavirus 2") | 1,168 |
| HPIV-3 | 3 | `NC_001796`, `NC_075446` (ZHYMgz01), `NC_038270` (= Simian Agent 10) | 4,449 |
| HPIV-4a | 1 | `NC_021928` (strain M-25) | 671 |
| Rhinovirus A | 2 | `NC_001617` (RV89), `NC_038311` (RV1, ATCC VR-1559) | 12,526 |
| Rhinovirus B | 2 | `NC_001490` (RV-B14), `NC_038312` (RV3) | 2,483 |
| Rhinovirus C | 2 | `NC_009996`, `NC_038878` (polyprotein CDS only) | 8,852 |
| **HBoV-1** | **0** | — | 1,168 |
| SARS-CoV-2 | 1 | `NC_045512.2` | 9,218,691 |
| HAdV-A | 2 | `AC_000005`, `NC_001460` (both "Human mastadenovirus A") | 587 |
| HAdV-B | 5 | `AC_000010` (SAdV-21), `AC_000018` (HAdV-7), `AC_000019` (HAdV-35), `NC_011202` (HAdV-B2), `NC_011203` (HAdV-B1) | 7,772 |
| HAdV-C | 4 | `AC_000007` (Ad2), `AC_000008` (Ad5), `AC_000017` (Ad1), `NC_001405` | 6,403 |
| HAdV-D | 3 | `AC_000006`, `NC_010956`, `NC_012959` (HAdV-54) | 3,155 |
| HAdV-E | 3 | `AC_000011` (SAdV-25), `NC_003266`, `NC_017825` (chimp AdV Y25) | 1,031 |
| HAdV-F | 1 | `NC_001454` | 3,949 |
| HAdV-G | 1 | `NC_006879` (**SAdV-1**, *not* a human genome) | 187 |

Note the NCBI NCBI-species placeholder records: `NC_001460`, `NC_001405`,
`NC_010956`, `NC_003266`, `NC_001454` are species-level
`Human mastadenovirus {A,C,D,E,F}` entries `[EBI: efetch]` — verify which
serotype each actually is before use; the `DEFINITION` gives no serotype.

### 5.1 RSV

Two genotypes, **A and B**, antigenically and genomically distinct `[SNIP: WHO
RSV fact sheet, 19 Dec 2025]`. The AG-2 lineage notation (GA-1…GA-3,
GB-1…GB-4) is the older scheme `[UNVERIFIED — not confirmed here]`. Two
sublineages are often used operationally: RSV-A ON1 and RSV-A ON2, and
RSV-A(BA)/RSV-B(BA) for the 2023-era global BA lineage `[UNVERIFIED]`.
**For a detection panel, A + B with 2–3 genomes each is sufficient** — RSV is
genetically far more conserved across the season than influenza.

### 5.2 Human metapneumovirus

Two antigenic subgroups, **A and B**, each split into two genetic sublineages
**A1, A2, B1, B2** (van den Hoogen 2004; Biacchesi 2003) `[SNIP: ICTV
naming proposal 2012.012V.N.v1; PMC446134]`. ICTV's naming proposal lists
**A1, A2, B1, B2**. A later, more granular scheme used in the Netherlands
surveillance literature is **A1, A2.1, A2.2.1, A2.2.2, B1, B2** (six
lineages) `[SNIP: PMC9973309, "Emergence and Potential Extinction of Genetic
Lineages of HMPV between 2005 and 2021"]` — A1 not detected in the Netherlands
after 2006; A2.2 viruses with 180-nt and 111-nt G-gene duplications became
dominant. A competing **A2a/A2b/A2c** scheme also exists `[SNIP: PMC12385779]`.
**Three naming schemes coexist.** Recommended: A1, A2.2.2 (111-nt dup),
A2.2.2 (180-nt dup), B1 — four genomes covers the currently circulating space.

### 5.3 Parainfluenza 1–4

HPIV-3 is the genetically most diverse: lineages **A1, A1.1, A1.2, A2, A2.1,
A2.2, A3, B, C, D, E, F, G** circulate; the EID study of SARS-CoV-2-negative
rapid-antigen-positive samples resolved A1, C, E, F, G lineages from
GISAID reference sets `[SNIP: EID 31(5):24-1191]`. HPIV-1 and HPIV-2 are
comparatively conserved. **Recommended: HPIV-1 ×2, HPIV-2 ×2, HPIV-3 ×4
(A1, A2, B, C/E), HPIV-4a ×2** (`NC_021928` is 4a; 4b also exists —
`[UNVERIFIED]`, accessions not resolved here).

### 5.4 Rhinovirus — the genotype problem is the panel problem

Three species, **A, B, C** `[SNIP: ICTV; PMC4441521]`. Genotype assignment
requires **>86–87% aligned nucleic acid identity in VP4/VP2 or VP1**
`[SNIP: PMC4441521, Study Group on Rhinovirus nomenclature]`. Genotype counts in
circulation differ by source and are inconsistent:

* 167 serotypes total `[SNIP: ScienceDirect topic page]`
* 169 serologic subtypes across 3 species `[SNIP: Johns Hopkins ABX Guide]`
* **RV-A 80, RV-B 32, RV-C 55 serotypes** `[SNIP: J Virol, "Human Rhinovirus
  Diversity and Evolution"]`
* 75 / 25 / 3+ `[SNIP: ScienceDirect topic page]`
* "over 200 distinct rhinoviruses may exist" `[SNIP: Yale PDF]`

**Do not attempt one genome per genotype** — that is 100–170 genomes and it
still would not cover reassortant genotypes. RV-A vs RV-B genome identity is
~51–58% `[SNIP: PMC1892812]`, i.e. they are essentially unrelated at the
genome level; within a species, 1A capsid region identity is ~70–96% and VP1
~80–96%. **Recommended: 4–6 genomes per species selected for maximum VP1/VP4
phylogenetic spread within the species, not for serotype coverage.** NCBI holds
only 2 RefSeq genomes per RV species `[EBI]`, so the panel must come from
GenBank (`[EBI]` e.g. `NC_009996`, and 173 "complete genome"-titled RV-C
records exist under txid463676).

### 5.5 Adenovirus

* **Seven human species (A–G)**, ~100 types `[SNIP: CDC MMWR 73(50), 18 Dec 2024]`;
  ICTV's Mastadenovirus report confirms species separated phylogenetically
  `[FETCH-via-index: elliot1.ictv.global]`.
* The **51 classic serotypes** are all sequenced; types **52–90+** were defined
  genomically, and by 2018 there were **90 distinct genotypes in 7 species**,
  >half of them species D `[SNIP: PMC6141750, "Adenoviromics"]`.
* The subset that matters for a **respiratory** panel: **B (3, 7, 11, 14, 16,
  21, 34, 35, 50), C (1, 2, 5, 6), E (4), and F (40, 41)** `[SNIP: PMC3955829
  clinical table]`. Species A (12, 18, 31) and D (8–10, 13, 15, 17, 19, 20,
  22–30, 32, 33, 36–39, 42–49, 51, 53) are predominantly
  keratoconjunctivitis/gastroenteritis; G is just HAdV-52.
* US surveillance 2017–2023, top respiratory types by share of respiratory
  specimens: **3 (23.7%), 4 (18.3%), 2 (13.9%), 7 (13.4%), 1 (11.3%)** `[SNIP:
  CDC MMWR 73(50)]`. That is a defensible, citable basis for a 5-genotype core.
* Note: Recombination within a species is common and generates new types —
  HAdV-3, 7, 11, 14 are all recombinants in species B `[SNIP: ICTV report]`.
* **Recommended: 1, 2, 3, 4, 5, 7, 11, 14, 55/56** ≈ 9–12 genomes, with a
  species-level rep each for D and E. Verified GenBank examples `[EBI]`:
  HAdV-14 complete genome `PQ657852`–`PQ657854`; HAdV-B1 `NC_011203`;
  HAdV-B2 `NC_011202`; HAdV-54 `NC_012959`.

### 5.6 Human bocavirus 1

**Zero RefSeq genomes** `[EBI]`. 1,168 GenBank records. Complete genomes exist,
e.g. `LC651171` (H565, complete) and `LC832987`
(`HBoV1/human/Japan/NIID24-146/2024`, complete) `[EBI]`. HBoV1 is single
serotype, so **1–2 genomes is enough**; the interesting HBoV axis is co-detection
with RSV, not diversity.

### 5.7 SARS-CoV-2

`NC_045512.2` is the only RefSeq genome `[EBI]`, and 9.2 M GenBank records exist.
For a *detection* panel, one genome (Wuhan-Hu-1) is generally sufficient and is
the RVDB reference `[SNIP: mSphere 2025, "Refinement of the RVDB"]`. Note
that lineage-defining sets (BA.1, BA.2, XBB, JN.1, LP.8.1) are a *phylogenetic*
need, not a detection need — a kallisto k-mer index against Wuhan-Hu-1
recalls any SARS-CoV-2 within ~0.1–0.5% divergence. Add 2–3 lineage
representatives only if the tool reports lineage.

---

## 6. Published guidance on reference count and clustering thresholds

**Honest answer: there is no citable published rule of thumb for "N reference
genomes per respiratory virus."** What exists:

1. **Reference count, qualitative.** `[SNIP: Virus Evolution 2(2):vew022, "Challenges
   in the analysis of viral metagenomes"]` — the field's explicit bottleneck:
   "the relatively few documented viral reference genomes compared to the
   estimated number of distinct viral taxa renders classification
   problematic," and "observed viral genomes often deviate considerably from
   reference genomes demanding use of exhaustive alignment approaches." It
   recommends using **multiple viral reference sequences** to mitigate
   single-reference failure. This is a *direction*, not a number.

2. **Database size, measured.** `[SNIP: PMC8953373, "Performance of Five
   Metagenomic Classifiers for Virus Pathogen Detection Using Respiratory
   Samples"]` — the single most directly relevant paper found. They used **one
   genome per human respiratory virus** (including HPIV-4, HCoV-NL63,
   HCoV-229E, influenza A, influenza B) as a common reference database across
   five classifiers, and obtained **83–100% sensitivity, 90–99% specificity,
   91–98% AUC**. That is the empirical anchor: one genome per virus is a
   *working* panel for clinical metagenomics.

3. **Clustering thresholds actually used in the literature:**

| Threshold | Tool | Source | Domain |
|---|---|---|---|
| **98%** nucleotide | CD-HIT-EST | `[SNIP: Goodacre et al., mSphere 3(2):e00069-18, 2018]` — RVDB, 561,676 creps | broad viral HTS, adventitious-agent detection |
| **98%** nucleotide | CD-HIT-EST | `[SNIP: Chin et al., mSphere, 2025, "Refinement of the RVDB"]` — SARS-CoV-2 handled separately at 98% then MMseqs2 | same |
| **95%** nucleotide | CD-HIT-EST (`-c 0.95 -aS 0.85`) | `[SNIP: PMC11694666, coral reef virome]` — viral contig dereplication | environmental virome |
| **95% ANI + 80% shared genes** | Virathon | `[SNIP: same]` — viral *population* dereplication | vMAG-level |
| **99%** amino acid | MMseqs2 (v15-6f452) | `[FETCH: PalmDB paper, Nat Biotechnol, Luebbert et al.]` — clustering RdRP "palmprints" into sOTUs | protein, not nucleotide |

**PalmDB is the most important single result for ViralScan specifically.**
`[FETCH: Nature Biotechnology, "Detection of viral sequences at single-cell
resolution identifies novel viruses associated with host gene expression
changes", Luebbert, Sullivan, Carilli, Hjörleifsson, Winnett, Chari, Pachter]`:

* kallisto was extended with **translated search** (nucleotide → amino-acid
  reference), retaining barcode/UMI single-cell resolution.
* **PalmDB = 296,623 unique RdRP-containing amino-acid sequences representing
  146,973 virus species**, clustered at 99% protein identity with MMseqs2.
* 36 MB reference; detects **16× more viruses than the 8,694 Riboviria
  reference genomes** in RefSeq.
* Recall up to 27.5% more viral RdRP sequences than Kraken2 translated search.
* Paper explicitly notes existing methods "are limited to (NCBI) reference
  genomes" and that cross-species contamination of reference databases
  misclassifies host reads — arguing for **host read removal before viral
  alignment**.

### 7. Recommended panel design

**Retention philosophy (ViralScan-specific judgement, not literature-derived).**
ViralScan's own SENTINEL negative-control study found a held-out plant
recovered only **8–40% of reads from TTV strains ≤84% identical** to any
reference. The dominant failure mode in a small-index kallisto detector is
**divergent low-level homologues**, so: retain completeness, drop only
provable redundancy. **Inclusion is broad; exclusion must be provable.**

**Deduplication: within-subtype only, 95% for viruses, 98% for slow-evolving
dsDNA viruses.**

Rationale: collapsing across subtypes conflates H3N8 human-infecting
(α2,3) with H3N2 (α2,6); within a subtype the HA/NA surface genes cluster far
below whole-genome 95% anyway, so CD-HIT retains the clade representatives that
actually matter. 98% for adenovirus and papillomavirus (slow, compact dsDNA;
near-duplicate serotypes are redundant). Run **at the segment level and then at
the concatenated level** — a "complete strain" is a labelled unit, so collapsing
strains loses the reassortment structure a segmented genome exists to express.

**Per-virus genome recommendations** (all segments complete, verified complete
genomes):

| Virus / group | Genomes | Composition rule | Indexed bp |
|---|---|---|---|
| **IAV — A(H1N1)pdm09** | **5–6** | 1 emergent (D.3.1.1), 1 recent NH CVV (D.3.1), 1 D.3, 1 founding clade (CA/07/2009), 1 early 6B.1A.5a | ~90 kb |
| **IAV — A(H3N2)** | **6–8** | **one per subclade**: 2a.3 (older), 2a.3a, 2a.3a.1 J.2, **J.2.4/K (emergent 2025-26)**, 2a.3b, 1a/1b residual | ~110 kb |
| **IAV — A(H1N2)v** | 2 | 1 North American, 1 Eurasian. Document as `variant_only` | ~27 kb |
| **IAV — A(H2N2)** | 1–2 | optional; `NC_007374`–`NC_007382` | ~27 kb |
| **IAV — A(H5N1) avian** | **10–14** | 2.3.2.1a/c/d/e ×1 each + **2.3.4.4b B3.13, D1.1, B3.2, B3.6** + 2.3.4.4 (a,e,g,h) — all H5Nx | ~190 kb |
| **IAV — A(H7N9)** | **2** | wave 1 (`NC_026422`–`NC_026429`) + wave 2 / intermediate | ~27 kb |
| **IAV — A(H9N2)** | **3** | G1-W, G1-E, BJ94/Y280 | ~40 kb |
| **IAV — A(H5N6)** | 2 | clade 2.3.4.4h (also 2.3.4.4g/e); 93 human cases | ~27 kb |
| **IAV — minor poultry** | 1 each | H3N8 (equine), H10N7, H10N3, H10N5, H6N1, H5N2, H5N8, H7N7, H3N6 | ~270 kb |
| **IAV total** | **35–42** | | **~550–600 kb** |
| **IB — Victoria** | 3 | V1A.3a (1A), V1A.3a.2 (3a.2), V1A.3a.3 (3a.3) | ~36 kb |
| **IB — Yamagata** | 1 | retained, flagged `retired_from_surveillance` | ~12 kb |
| **ICV** | 1 | `NC_006306`–`NC_006312`; C/SC — no human public-health concern | ~13 kb |
| **RSV A / RSV B** | 3 / 3 | ON1/ON2 + BA-era; A2 + long; B-WaDC-18537-1962 | ~91 / ~91 kb |
| **hMPV** | 4 | A1, A2.2.2-111dup, A2.2.2-180dup, B1 | ~53 kb |
| **HPIV-1 / 2 / 3 / 4** | 2 / 2 / 4 / 2 | HPIV-3: A1, A2, B, C | ~94 kb |
| **Rhinovirus A / B / C** | 5 / 4 / 5 | max VP1 phylogenetic spread per species | ~100 kb |
| **Adenovirus B / C / E** | 4 / 4 / 2 | B: 3, 7, 11, 14; C: 1, 2, 5, 6; E: 4 | ~255 / ~255 / ~65 kb |
| **Adenovirus D** | 1 | species rep only | ~35 kb |
| **HBoV-1** | 2 | `LC651171`, `LC832987` | ~11 kb |
| **SARS-CoV-2** | 1 (+3 optional lineage reps) | `NC_045512.2` | 30 / 120 kb |
| **Seasonal coronaviruses** | 1 each | 229E, NL63, HKU1, OC43 | ~110 kb |
| **Respiratory total** | **~45** | | **~1.1 MB** |

**Grand total ≈ 85 genomes ≈ 1.7 MB** (plus host contigs). Index build at
`k=31` is a few minutes. That is a **55% genome increase** over the current
195-entry panel and roughly **2.5× the total bases** while fixing every lineage
gap in §0.

**Per-segment counting must be per-genome, not per-accession.** RefSeq IAV = 56
accessions = 7 genomes; non-segmented viruses = 1 accession = 1 genome. Report
both numbers and never mix them. Your `reference_manifest.json` should carry
`genome_id`, `segments[]` with per-segment accession/length/completeness,
`lineage` (internal label), `gisrs_label_asof`, and `narrow_mappability_flags`.

**Mandatory engineering constraints, each traceable to a failure mode above:**

1. **Segment-namespaced targets.** IAV is 8 segments; concatenating without
   namespacing creates artifactual junction k-mers. Use per-segment FASTA
   records with a segment index in the transcript name.
2. **Never derive subtype from the NCBI organism field.** `OP212288` (A/Texas/
   61/2022 HA) is typed `(H3N2)` in the organism string but is clade 2.3.4.4b
   H5 `[EBI]`. `"A/Texas/61/2022"[Organism]` returns 0 hits.
3. **Model crosses/subtype spans.** For a respiratory panel most CI genomes are
   single-segment; the only cross-genome CI risk is IAV's 8-segment
   concatenation, handled by namespacing rather than by CI.
4. **Keep any viral CI and label the ambiguity** ("multimapping: IAV HA
   H1/H3/H5/H7/H9/H10" must be a first-class output, not a discarded bucket —
   the class exists for exactly this virus).
5. **Cap HA/NA/NA-adjacent segment identity loss.** Because HA/NA divergence is
   where a 95% whole-genome collapse is most likely to hide a clade difference,
   cluster HA and NA at a stricter threshold than internal genes.
6. **Host subtraction before viral quantification**, per the PalmDB paper.
7. **Negative-control expectation must be published.** The TTV 8–40%
   held-out-divergent result sets the realistic ceiling; a panel manifest that
   claims ≥90% recall on divergent homologues is not credible.

**Unresolved and explicitly not answered by this research:**
`[UNVERIFIED]` classical swine H1N1 canonical accession set; `[UNVERIFIED]`
H7N9 wave-2 accessions; `[UNVERIFIED]` H10N7 and H3N8 human-case accessions;
`[UNVERIFIED]` H9N2 `Gx.y` nomenclature primary source; `[UNVERIFIED]` B3.2/B3.6
and EA-* European H5 genotype accessions; `[UNVERIFIED]` HAdV species-level
RefSeq serotype identity; `[UNVERIFIED]` HPIV-4b; `[UNVERIFIED]` RSV ON/BA
sublineage naming primary source.

---

## 8. Citations

**WHO**
1. Influenza (avian and other zoonotic), fact sheet, 21 Jan 2026. https://www.who.int/news-room/fact-sheets/detail/influenza-(avian-and-other-zoonotic) `[FETCH]`
2. Influenza at the human-animal interface — Summary and risk assessment, 23 January to 31 March 2026. https://cdn.who.int/media/docs/default-source/influenza/human-animal-interface-risk-assessments/influenza-at-the-human-animal-interface-summary-and-assessment--from-23-january-to-31-march-2026.pdf `[FETCH]`
3. Recommended composition of influenza virus vaccines for use in the 2026-2027 northern hemisphere influenza season, 27 Feb 2026. https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2026-2027-northern-hemisphere-influenza-season
4. Candidate vaccine viruses and potency testing reagents. https://www.who.int/teams/global-influenza-programme/vaccines/who-recommendations/candidate-vaccine-viruses `[FETCH]`
5. Recommendations for influenza vaccine composition. https://www.who.int/teams/global-influenza-programme/vaccines/who-recommendations
6. Influenza Update No. 463, 22 Jan 2024. https://cdn.who.int/media/docs/default-source/influenza/influenza-updates/2023/2024_01_22_surveillance_update_463.pdf
7. Recommended composition of influenza virus vaccines, 2023-2024 NH. https://cdn.who.int/media/docs/default-source/influenza/who-influenza-recommendations/vcm-northern-hemisphere-recommendation-2023-2024/202302_seasonal_recommendation_a.pdf
8. Human infection with avian influenza A(H5) viruses (WPRO), 30 May 2025. https://cdn.who.int/media/docs/default-source/wpro---documents/emergency/surveillance/avian-influenza/ai_20250530.pdf
9. Respiratory syncytial virus (RSV) fact sheet, 19 Dec 2025. https://www.who.int/news-room/fact-sheets/detail/respiratory-syncytial-virus-(rsv)
10. Standardization of terminology for influenza virus variants infecting humans: Update (Feb 2021). https://cdn.who.int/media/docs/default-source/influenza/global-influenza-surveillance-and-response-system/nomenclature/standardization_of_terminology_influenza_virus_variants_update.pdf
11. Updated unified nomenclature system for the HPAI H5N1 avian influenza viruses (WHO/OIE/FAO H5N1 Evolution Working Group, Oct 2011). https://cdn.who.int/media/docs/default-source/influenza/global-influenza-surveillance-and-response-system/nomenclature/updated_nomenclature_system_h5n1_avian_influenza_viruses.pdf

**CDC**
12. Reported Human Infections with Avian Influenza A Viruses, updated 2026-03-20. https://www.cdc.gov/bird-flu/php/surveillance/reported-human-infections.html `[SNIP — direct fetch 403]`
13. Global Summary of Recent Human Cases of H5N1 Bird Flu, 4 Aug 2025. https://www.cdc.gov/bird-flu/spotlights/h5n1-summary-08042025.html
14. Types of Influenza Viruses, 26 Sep 2025. https://www.cdc.gov/flu/about/viruses-types.html
15. Trivalent Influenza Vaccines, updated 2026-05-11. https://www.cdc.gov/flu/vaccine-types/trivalent.html
16. Influenza Activity in the United States during the 2024-25 Season and Composition of the 2025-26 Influenza Vaccine. https://www.cdc.gov/flu/whats-new/2025-2026-influenza-activity.html
17. Influenza Activity in the United States during the 2023-24 Season and Composition of the 2024-25 Vaccine. https://www.cdc.gov/flu/whats-new/flu-summary-2023-2024.html
18. Influenza B/Yamagata Update (ACIP, 28 Feb 2024). https://stacks.cdc.gov/view/cdc/148509/cdc_148509_DS1.pdf
19. Surveillance of Human Adenovirus Types and the Impact of the COVID-19 Pandemic on Reporting — US, 2017-2023. MMWR 73(50), 18 Dec 2024. https://www.cdc.gov/mmwr/volumes/73/wr/mm7350a1.htm
20. Swine flu: CDC Confirms Another Human Infection with Flu Virus from Pigs. https://www.cdc.gov/swine-flu/comm-resources/swineflu-infection.html
21. Flu Viruses of Special Concern (archived). https://archive.cdc.gov/www_cdc_gov/flu/pandemic-resources/monitoring/viruses-concern.html
22. Influenza Risk Assessment Tool (IRAT) — clade 2.3.4.4b A(H5N1), March 2025. https://www.cdc.gov/pandemic-flu/media/pdfs/2025/IRATA-California-Washington.pdf
23. About Human Metapneumovirus, 10 Feb 2026. https://www.cdc.gov/human-metapneumovirus/about/index.html

**FDA**
24. Influenza Vaccine Composition for the 2026-2027 U.S. Influenza Season (VRBPAC, 12 Mar 2026). https://www.fda.gov/vaccines-blood-biologics/vaccines/influenza-vaccine-composition-2026-2027-us-influenza-season
25. VRBPAC March 12, 2026 — Global Influenza Virus Surveillance and Characterization (CDC slide deck). https://www.fda.gov/media/191518/download
26. VRBPAC March 5, 2024 — Global Influenza Virus Surveillance and Characterization. https://www.fda.gov/media/176782/download

**ECDC / EFSA / EU**
27. Avian influenza overview March–May 2026, 26 Jun 2026. https://www.ecdc.europa.eu/en/publications-data/avian-influenza-overview-march-may-2026 `[FETCH]`
28. Influenza virus characteristics, week 40 2024 to week 33 2025, EU/EEA, 30 Sep 2025. https://www.ecdc.europa.eu/en/publications-data/influenza-virus-characteristics-week-40-2024-week-33-2025
29. Influenza virus characterisation — ECDC/WHO, 31 Mar 2024. https://www.ecdc.europa.eu/sites/default/files/documents/influenza-ECDC-WHO-Report-March-2024.pdf
30. Assessing the risk of influenza for the EU/EEA (TA Brief, 2025) — "Subclade K accounts for 47% of A(H3N2) sequences … 1 May to 17 November". https://www.ecdc.europa.eu/
31. Avian influenza overview December 2024 – March 2025. https://health.ec.europa.eu/document/download/e04590a4-6b71-4620-b528-37578b1a7a16_en

**NCBI**
32. Influenza Virus Database. https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database
33. Influenza Virus Genome Set. https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=genomeset `[FETCH]`
34. RefSeq. https://www.ncbi.nlm.nih.gov/refseq/ — Release 237
35. RefSeq viral release. https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
36. E-utilities (esearch/esummary/efetch) — all counts in §4/§5 `[EBI]`

**Literature**
37. Luebbert L, Sullivan DK, Carilli M, Hjörleifsson KE, Winnett AV, Chari T, Pachter L. Detection of viral sequences at single-cell resolution identifies novel viruses associated with host gene expression changes. *Nature Biotechnology* 44:100–109 (2026). https://www.nature.com/articles/s41587-025-02614-y `[FETCH]` — PalmDB: 296,623 RdRP amino-acid sequences, 146,973 species, 99% MMseqs2 clustering, 36 MB, 16× more viruses than 8,694 Riboviria RefSeq genomes
38. Nomenclature for Tracking of Genetic Variation of Seasonal Influenza Viruses. *Influenza Other Respir Viruses*, doi 10.1111/irv.70230. https://pmc.ncbi.nlm.nih.gov/articles/PMC12904685 `[FETCH, Europe PMC]`
39. Goodacre N, Aljanahi A, Nandakumar S, Mikailov M, Khan AS. A Reference Viral Database (RVDB) to Enhance Bioinformatics Analysis of HTS for Novel Virus Detection. *mSphere* 3(2):e00069-18 (2018). https://pmc.ncbi.nlm.nih.gov/articles/PMC5853486 — CD-HIT-EST at 98% nt
40. Chin P-J, Bhavsar JD, Bosma TJ, MacDonald ML, Polson SW, Khan AS. Refinement of the Reference Viral Database (RVDB) for improving bioinformatics analysis of virus detection by HTS. *mSphere* (2025). https://pmc.ncbi.nlm.nih.gov/articles/PMC12306153
41. Performance of Five Metagenomic Classifiers for Virus Pathogen Detection Using Respiratory Samples from a Clinical Cohort. https://pmc.ncbi.nlm.nih.gov/articles/PMC8953373 — one genome per respiratory virus, 83–100% sens
42. Challenges in the analysis of viral metagenomes. *Virus Evolution* 2(2):vew022. https://academic.oup.com/ve/article/2/2/vew022/2797616
43. Laboratory validation of a clinical metagenomic next-generation sequencing assay for respiratory virus detection and discovery. *Nat Commun* (2024). https://pubmed.ncbi.nlm.nih.gov/39532844
44. Clade 2.3.4.4b highly pathogenic avian influenza H5N1 viruses: knowns, unknowns, and challenges. *J Virol* 99(6):e00424-25 (2025). https://journals.asm.org/doi/10.1128/jvi.00424-25
45. H5N1 Clade 2.3.4.4b: Evolution, Global Spread, and Host Range Expansion (2025). https://pmc.ncbi.nlm.nih.gov/articles/PMC12472894
46. H5N1 2.3.4.4b: a review of mammalian adaptations and risk of pandemic emergence. *J Gen Virol* 106(6) (2025). https://www.microbiologyresearch.org/content/journal/jgv/10.1099/jgv.0.002109
47. Fusaro A, et al. High pathogenic avian influenza A(H5) viruses of clade 2.3.4.4b in Europe. *Virus Evolution* 10:veae027 (2024). https://pubmed.ncbi.nlm.nih.gov/38699215 — 1,956 European genomes
48. Lam TT, et al. Nomenclature updates to the hemagglutinin gene clade designations … clades 2.3.2.1c and 2.3.4.4. bioRxiv 2025.11.23.690055. https://www.biorxiv.org/content/10.1101/2025.11.23.690055v1
49. Highly Pathogenic Avian Influenza A(H5N1) Clade 2.3.4.4b Virus Infection in Poultry Farm Workers, Washington, USA, 2024. *EID* 31(12). https://wwwnc.cdc.gov/eid/article/31/12/pdfs/25-1118-combined.pdf
50. Chen H et al. Whole-Genome Characterization of a Novel Human Influenza A(H1N2) Virus Variant, Brazil. *EID* (PMC5176240). https://pmc.ncbi.nlm.nih.gov/articles/PMC5176240
51. Cogdale J, et al. A case of swine influenza A(H1N2)v in England, November 2023. (PMC10797662). https://pmc.ncbi.nlm.nih.gov/articles/PMC10797662
52. Zhao D, et al. First report of human infection caused by swine-origin influenza A(H1N2)v virus in mainland China, Kunming (Yunnan), 2025. *Front Public Health* 14:1893536 (2026). https://pmc.ncbi.nlm.nih.gov/articles/PMC13461623
53. Caini S, et al. Probable extinction of influenza B/Yamagata and its public health implications. *Lancet Microbe* (2024). https://www.sciencedirect.com/science/article/pii/S2666524724000661
54. Paget J, Caini S, Del Riccio M. Has influenza B/Yamagata become extinct? *Euro Surveill* 27(39):2200753 (2022). https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2022.27.39.2200753
55. Kitamura N, et al. Proposals for the classification of human rhinovirus species A, B and C into genotypically assigned types. *J Gen Virol* (2013); and the Study Group on Rhinovirus nomenclature site (PMC4441521). https://pmc.ncbi.nlm.nih.gov/articles/PMC4441521
56. Human Rhinovirus Diversity and Evolution. *J Virol* (PMC / doi 10.1128/jvi.01659-16). https://journals.asm.org/doi/10.1128/jvi.01659-16
57. Genome-wide diversity and selective pressure in the human rhinovirus. (PMC1892812). https://pmc.ncbi.nlm.nih.gov/articles/PMC1892812
58. ICTV taxonomy report chapter: Family Adenoviridae. https://elliot1.ictv.global/report/chapter/adenoviridae/adenoviridae ; Genus Mastadenovirus. https://elliot1.ictv.global/report/chapter/adenoviridae/adenoviridae/mastadenovirus
59. Human adenovirus: Viral pathogen with increasing importance. (PMC3955829) — clinical serotype table. https://pmc.ncbi.nlm.nih.gov/articles/PMC3955829
60. Adenoviromics: Mining the Human Adenovirus Species D Genome. (PMC6141750) — 90 genotypes in 7 species. https://pmc.ncbi.nlm.nih.gov/articles/PMC6141750
61. van den Hoogen BG, et al. Recovery of human metapneumovirus genetic lineages A and B from cloned cDNA. (PMC446134). https://pmc.ncbi.nlm.nih.gov/articles/PMC446134
62. Emergence and Potential Extinction of Genetic Lineages of Human Metapneumovirus between 2005 and 2021. (PMC9973309). https://pmc.ncbi.nlm.nih.gov/articles/PMC9973309
63. ICTV proposal 2012.012V.N.v1 — naming convention for hMPV strains (A1, A2, B1, B2). https://ictv.global/sites/default/files/web-files/General_Information/2012.012V.N.v1.metapneumovirus_names.pdf
64. Epidemiological and Genetic Characterization of Three H9N2 Viruses Causing Human Infections — Changsha, April 2025. (PMC12620571). https://pmc.ncbi.nlm.nih.gov/articles/PMC12620571
65. A Global Perspective on H9N2 Avian Influenza Virus. *Viruses* (PMC6669617) — G1 / BJ94 / Y439 lineages. https://pmc.ncbi.nlm.nih.gov/articles/PMC6669617
66. Updated joint FAO/WHO/WOAH public health assessment of recent influenza A(H5) events in animals and people, Jul 2025. https://www.woah.org/app/uploads/2025/07/25728-fao-woah-who-h5-assessment.pdf
67. Updated FASTA/GenBank assembly study — Influenza A(H3N2) Subclade K (J.2.4.1). *Pathogens* 15(2):37. https://www.mdpi.com/2036-7449/18/2/37
68. Sullivan DL, et al. kallisto, bustools and kb-python for quantifying bulk, single-cell and single-nucleus RNA-seq. *Nat Protoc* (2024). https://pmc.ncbi.nlm.nih.gov/articles/PMC10690192

**Other**
69. Centre for Health Protection, Hong Kong. Avian Influenza Report, Vol 22 Wk 03 / Vol 21 Wk 51. https://www.chp.gov.hk/files/pdf/2026_avian_influenza_report_vol22_wk03.pdf
70. ATCC Respiratory Pathogen Panels (MP-36 inclusivity, MP-37 exclusivity) — useful for positive-control design. https://www.atcc.org/-/media/product-assets/documents/panels/microbiology/respiratory-pathogen-panels.pdf
71. Falk S, et al. INSaFLU: an automated open web-based bioinformatics suite "from-reads" for influenza WGS-based surveillance. *Genome Medicine* 10:49 (2018). https://link.springer.com/article/10.1186/s13073-018-0555-0
72. MCRL: using a reference library to compress a metagenome … *Bioinformatics* 38(3):631. https://academic.oup.com/bioinformatics/article/38/3/631/6390794
73. CD-HIT User's Guide. https://home.cc.umanitoba.ca/~psgendb/birchhomedir/doc/cd-hit/cdhit-user-guide.pdf
