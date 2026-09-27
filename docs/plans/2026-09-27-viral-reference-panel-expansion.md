# ViralScan viral reference panel — breadth and strain-diversity design

> **SUPERSEDED IN PART — read `PLAN.md` `WP4H` and `.living/findings/F-014` first.**
> This document is kept for the evidence it records, but four of its
> recommendations were **retracted** on 2026-09-27 after adversarial review, and
> the decision that motivated them was **overturned**. `PLAN.md` is the tracker;
> where the two disagree, `PLAN.md` wins.
>
> **Retracted**
> 1. *§4.3 / `CAT-20`: "ICTV is now 37 genera / 243 species, not the 8 the shipped
>    table encodes", implying a taxonomy refresh is needed.* It is not. All 16
>    ICTV genera absent from the panel have **zero** complete genomes with
>    `"Homo sapiens"[Host]` (Lambda 408/0, Eta 256/0, Iota 139/0, Kappa 24/0) —
>    swine, feline, canine, tupaia and pinniped viruses. Adding them buys no
>    sensitivity and adds false-positive surface. The refresh is **deleted**.
> 2. *§4.3: "our genus labels are unreliable predictions".* **Backwards.** The
>    582 panel genomes with no GenBank genus are the ~584 that the upstream ORF1
>    phylogeny **resolved**. Our labels are better than GenBank's. Do **not**
>    reconcile them against GenBank lineages.
> 3. *§2.3: "re-pull by family taxid, est. 20–40 k new records / ~60 Mbp".* Wrong
>    twice over: `clareaulab/human_anellovirus_pangenome` already performed that
>    reduction (3,545 → **2,023** representatives, 4.80 % hardmasked, **0**
>    pure-homopolymer k-mers), and the 42,755 raw `Anelloviridae[Organism]`
>    records include non-human hosts.
> 4. *§3.1 / `CAT-18`: "only 1 accession carries both real genes and a
>    placeholder", from grouping `t2g` by accession prefix.* Wrong key. There is
>    a second duplicate — `AB303562.1` / `NC_038359.1`, byte-identical, the
>    GenBank and RefSeq copies of one Gammatorquevirus genome, **present in the
>    upstream panel too**. `CAT-05` must key on **sequence**.
>
> **Overturned decisions**
> - *§5.1 `CAT-18` (two-index architecture)* — not adopted. Measured on the built
>   `t2g`: 470,468 rows = 468,425 real + 2,043 placeholders (0.4 %), and the
>   99.8 % single-bucket failure is *within*-anellovirus, not gene-vs-genome. One
>   index, one feature set per accession.
> - *§4.2 / `CAT-20` (translated ORF1 `--aa` index)* — not adopted; a single
>   nucleotide index was chosen. Consequence accepted: **anellovirus genus is a
>   label from the catalogue, never a measurement**, and is not recoverable from a
>   nucleotide panel (measured: no 31-mer is shared by ≥1,000 of the genomes).
> - *§8 Q1 (SARS-CoV-2 lineage genomes)* — resolved: ENA/INSDC, 4–6 genomes.
> - *§8 Q3 (two indexes)* — resolved: one. *Q4 (translated ORF1)*: declined.
>   *Q5 (bundle)*: bundle, and the size crisis in §7.3 is void — the whole
>   hardmasked anellovirus panel is **6.03 MB**.
>
> **Superseded by measurement.** The document's §4 framing — "the anellovirus
> panel is probably structurally incomplete" — is the right instinct, landed on
> by the wrong evidence. The real, measured defects are in F-014 (the panel
> carries **170 pure-homopolymer 31-mers**; the 1.44 % false-read floor is a
> **k-mer** property, not an N-masking one) and in the `CAT-05` duplicate above.
>
> **Also corrected here:** §2.2's Serratus anchor is **17** viruses, not 18 —
> `sources/viral_panel/Serratus_hits_all_viruses.tsv` has exactly 17 rows with
> `reactivation_candidate == TRUE`, and HHV-6B is not among them.

Status: **design proposal (partly superseded)**
Date: 2026-09-27
Branch: `codex/viralscan-v3`
Operational tracker: [`PLAN.md`](../../PLAN.md) → `WP4H`, rows `CAT-09`…`CAT-30`
(this document supplies the evidence and the rows that were missing from it; it does not
supersede `WP4H`)

---

## 0. What this document adds to WP4H

`WP4H` already asks the right question and has the right two-tier decision. Six of its
open rows are under-specified in ways that would produce a wrong panel, and five
failure modes found in this session are not in it at all. Summary:

| # | Gap in WP4H | Where |
|---|---|---|
| G1 | `CAT-12` says "SARS-CoV-2: RefSeq plus one per WHO variant lineage". **RefSeq holds exactly one SARS-CoV-2 genome and no lineage labels anywhere.** The instruction is not executable. | §3.1 |
| G2 | `CAT-12` says "Influenza A: one strain per relevant subtype". The axis that carries the diversity is the **clade within a subtype**, and the current 8 IAV RefSeq records are all segments of *one 1934 lab strain* — so "one per subtype" still leaves H3N2 with no circulating isolate. | §3.2 |
| G3 | **Not in the plan at all:** a naive whole-genome anellovirus panel produces **1.44 % of R2 reads as false hits**, 100 % attributable to homopolymer / tandem-repeat k-mers. Measured. | §4.3 |
| G4 | **Not in the plan:** `kb count` **discards** multimapping UMIs by default in the scRNA-seq path, so adding near-identical genomes converts viral reads into *dropped* UMI rather than split counts. Expansion is not monotone-good. | §5.2 |
| G5 | **Not in the plan:** anelloviruses have **no conserved nucleotide k-mer space at all** (0 of 2,042 genomes share any 31-mer with ≥1,000 others), so a nucleotide panel cannot be genus-sensitive. The only route is translated search on ORF1. | §4.2 |
| G6 | **Not in the plan:** the anellovirus panel is probably **structurally incomplete**. `Anelloviridae[Organism]` = 42,755 records but genus-name queries reach only 21,283; SENV (taxid 136966) sits outside every genus and would never be retrieved. ICTV is now 37 genera / 243 species, not the 8 genera the shipped table encodes. | §4.4 |
| G7 | **Not in the plan:** the **Serratus reactivation panel as an explicit, checkable acceptance criterion** — the 18 latent-DNA-virus accessions from Lareau et al. 2023. This is the "Serratus activation screen" anchor and belongs in `CAT-10` as a testable list, not as prose. | §2.2 |
| G8 | **Not in the plan:** `kb --workflow=custom` / `--aa` **silently disables host masking** unless `--d-list` is passed explicitly, and CoV 3′ poly-A tails are clipped at index build. | §5.3 |

---

## 1. Where the panel stands

Measured 2026-09-27 against
`references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.fa`
and `src/viralscan/data/`.

| | |
|---|---:|
| Records in the widest build | 2,313 |
| Anelloviridae + Gyrovirus | 2,041 (88 %) |
| Everything else | 272 records → 107 distinct species |
| …of which HHV-6B pseudocontig fragments | 97 (not genomes) |
| Bundled GTFs (`src/viralscan/data/*.gtf`) | 195 files, 2,692 gene IDs |
| Anellovirus accessions | 2,042 (1,995 annotated, 2,515 genes) |
| Anellovirus genomes, total bp | 5,994,773 (median 2,896; min 1,141; max 3,996) |

Effective non-anellovirus content: **~175 genomes over ~106 species, almost all
single-genome.**

Absent entirely: **SARS-CoV-2**, HIV-1/2, HTLV-1/2, HCoV-OC43/NL63/HKU1, hMPV,
bocavirus, influenza D, TSPyV, HPyV6/7, HPV18 and every high-risk HPV type except 16.
Influenza A is a single 1934 lab strain (PR8); influenza B a 1940 strain; influenza C a
1950 strain.

**Scale is not the constraint.** The full target panel is ≈ **30 Mbp** of viral
sequence against a 3.15 Gb host genome — under 1 % by mass. The binding constraints
are *semantic*: annotation quality, equivalence-class behaviour, and host
cross-talk. `CAT-13`'s size gate should be reframed accordingly (see §7.3).

---

## 2. Target scope

### 2.1 Tier structure (adopt WP4H's decision; refine the rules)

| Tier | Rule | Purpose | Expected genomes |
|---|---|---|---:|
| **T1 — breadth floor** | one representative per human-host RefSeq virus species | close the absence list in one step | ~350–450 |
| **T2 — persistence set** | literature-curated, **deep strain/clade diversity**, cited per entry | latent/reactivation calling, gene programmes | ~600–900 |
| **T3 — anellovirus cloud** | full family panel + ICTV exemplars | presence-only, genus-resolution via translated search | 2,042 + 243 |

The Serratus evidence (§2.2) defines **T2**, and the ICTV/MSL41 evidence (§4.4) defines
**T3**. T1 is the floor that makes "cover as many viruses as possible" true.

### 2.2 T2 anchor: the Serratus reactivation panel

The "Serratus activation screen" is **Lareau CA, …, Satpathy AT. "Latent human
herpesvirus 6 is reactivated in CAR T cells." *Nature* 623(7987):608–615 (2023).**
doi:[10.1038/s41586-023-06704-2](https://doi.org/10.1038/s41586-023-06704-2) ·
PMID 37938768 · code `github.com/caleblareau/serratus-reactivation-screen`.

Design facts worth copying: input panel = **129 curated human viruses** (ViralZone /
Hulo 2011); reactivation criterion = **DNA-genome-only, known latent replication
cycles that can reactivate *in vivo*** (their ref: Traylen et al., *Future Virology*
6:451, 2011), narrowed to **17 DNA viruses from Herpesviridae, Polyomaviridae,
Adenoviridae and Parvoviridae**; positive call = **≥100 reads AND ≥50 % mean mapped
identity** per sample–virus pair.

> Every reactivation-capable agent in this literature is a **DNA** virus. That is
> precisely why an RdRP-only platform cannot find them, and the strongest argument for
> ViralScan's separate reference.

**All 18 listed accessions are already in the panel** (verified against the GTF
filenames) — T2's *latent DNA* core needs no new accessions, only strain diversity
and real gene structure. It becomes the acceptance criterion for the panel, not new
breadth:

| Family | Accessions (verified in panel) |
|---|---|
| Herpesviridae (10) | `NC_006273` HCMV · `NC_001664` HHV-6A · `NC_000898` HHV-6B · `NC_001716` HHV-7 · `NC_009333` KSHV · `NC_001806` HSV-1 · `NC_001798` HSV-2 · `NC_001348` VZV · `NC_007605` EBV · `NC_006560` Cercopithecine herpesvirus |
| Polyomaviridae (5) | `NC_001538` BK · `NC_001699` JC · `NC_009238` KI · `NC_009539` WU · `NC_010277` Merkel cell |
| Adenoviridae (1) | `NC_001405` |
| Parvoviridae (2) | `NC_000883` B19 · `NC_001401` AAV |

(Upstream labels this "17 DNA viruses" but enumerates 18; the discrepancy is in the
source and is flagged rather than resolved.)

T2 beyond the Serratus anchor, per WP4H `CAT-10`: polyomaviruses TSPyV/HPyV6/HPyV7;
HIV-1/2, HTLV-1/2; HBV, HDV; HPV; parvovirus B19; measles (SSPE); HCV; pegivirus.
Each needs a `persistence_class` (`latent-episomal`, `latent-integrated`,
`chronic-productive`, `persistent-commensal`, `recurrent-lytic`).

**Known Serratus weaknesses to inherit deliberately:** the input panel contains **no
SARS-CoV-2** and only **3 HPV genotypes**; HHV-6 was filed as HHV-6A only and needed a
separate HHV-6B query; TTV appears with 206 incidental hits and is never flagged
reactivation-relevant. T2 must be a *superset* of Serratus, never a copy.

---

## 3. Strain and lineage diversity

The governing principle, from the only benchmark that measured it
(van Bemmelen, Nika & Baaijens, *BMC Genomics* 2026, doi:10.1186/s12864-026-12874-w):

> For genomes >99 % identical, **accuracy falls as the reference set grows** — Spearman
> ρ = **−0.73** (F1) and **−0.57** (abundance) against set size. The largest single
> gain came from **restricting to geographically plausible diversity**
> (**+109 %** abundance, **+240 %** F1), not from finer clustering.

So: **add genomes only where they add k-mer space.** Concretely —

### 3.1 SARS-CoV-2 — the hard blocker (G1)

**RefSeq contains exactly one SARS-CoV-2 genome: `NC_045512.2`** (Wuhan-Hu-1, 29,903 bp).
Verified three ways: `txid2697049[Organism] AND srcdb_refseq[PROP]` → 1;
`"Severe acute respiratory syndrome coronavirus 2"[Organism] AND refseq[filter]` → 1.
**Zero RefSeq records database-wide carry a `pango_lineage` attribute** — it is a
BioSample attribute, not a nuccore one. `NC_045512.2`'s flatfile contains no occurrence
of `pango`, `gisaid` or `lineage`.

The other 9.2 M SARS-CoV-2 sequences exist only in INSDC, **unlabelled**.

Consequences for `CAT-12`:

1. The row as written ("RefSeq plus one per WHO variant lineage") **cannot be
   executed** from RefSeq. Either accept a 1-genome panel, or source from
   INSDC/ENA with lineage assignment done locally (pangolin/UShER), which introduces
   a GISAID data-use-agreement dependency and a CI reproducibility problem.
2. **Decision required from the user** (see §8 Q1).
3. **Recommended cap: 4–6 genomes, not 12–15.** At k=31, genomes differing by ~30–70 nt
   are indistinguishable in conserved regions; extra genomes add essentially zero
   k-mer space while multiplying equivalence-class size. Proposed set:
   `NC_045512.2` (ancestral, carries UTR coverage later lineages lack) + one
   JN.1-descendant + one **BA.3.2** (the only non-JN.1 branch in circulation — a
   saltation with ~40 spike substitutions, genuinely different) + one XFG + one
   pre-Omicron VOC (BA.4.6 or XBB.1.5) for re-analysis of 2021–2025 datasets.

Current WHO designations (page updated 2026-07-28): VOI **JN.1**; VUMs **XFG**,
**NB.1.8.1**, **PQ.16.1.1**, **BA.3.2**. **Pango is depth-capped at 4 nodes**
(measured: 0 lineages deeper than 3 dots) while WHO counts >3,700 JN.1 descendants —
so "one genome per lineage" is not a well-posed target at any depth.

Verified companion accessions: SARS-CoV-1 `NC_004718.3`; MERS-CoV `NC_019843.3`
**and** `NC_038294.1` (2 RefSeq records — include both); HCoV-229E `NC_002645.1`;
NL63 `NC_005831.2`; OC43 `NC_006213.1`; HKU1 `NC_006577.2`. Note
`NC_028752.1` (camel alphacoronavirus) is filed by NCBI under taxid 11137 alongside
HCoV-229E — keep it deliberately as a **cross-reactivity control**.

### 3.2 Influenza — clade, not subtype (G2)

Current IAV RefSeq content is **7 strains**, all segments of lab/legacy isolates, with
**no H3N2 and no circulating isolate**. "One strain per relevant subtype" is therefore
insufficient: H3N2's diversity lives in its clades.

⚠ **Clade names were renamed in Feb 2023** — `3C.2a1b.2a.2a.3a.1` → `2a.3a.1`; the
current emergent clade is `2a.3a.1 (J.2.4.1)` = "**K**" (Aug 2025; 47 % of EU
submissions; NH 2026-27 CVV). H1N1pdm09 is `5a.2a.1` D.3.1 / D.3.1.1. **Do not hard-code
`3C` strings**; store both and read the current name from a dated source.

| Agent | Budget | Axis that carries diversity |
|---|---:|---|
| IAV H1N1 pdm09 | 6 | pre-2009 seasonal vs 2009 pandemic; clades 3C.3a, 6B, 6B.1, 6B.2 |
| IAV H1N1 (non-pdm09) | 2 | classical swine; 2009-era seasonal |
| IAV H3N2 | 8 | 3C.2a, 3C.2a1, 3C.2a1b.1a, 2a, 2a.1, 2a.3, **2a.3a.1**, **2a.3a.1 (J.2.4.1)** |
| IAV H5N1 | 6 | 2.3.4, 2.3.4.4, **2.3.4.4b**; human-case genotypes B3.13 and D1.1 |
| IAV H7N9 | 3 | 2013 wave + post-2014 |
| IAV H9N2 | 4 | human + avian lineages (G.x.y) |
| IAV H1N2 | 2 | seasonal 2000-03; sporadic H1N2v (China notified WHO of a human H1N2v case **3 Feb 2026**) |
| IAV H2N2 | 1 | historical only (gone since 1968-69) |
| IB Victoria | 4 | lineage/ subclade |
| IB Yamagata | 1 | **flag `retired_from_surveillance`** — no confirmed detection since March 2020; dropped from the NH 2026-27 CVVs but still in the Sept 2026 Southern Hemisphere list, and it is the LAIV component. Keep 1, flagged. |
| IC / ID | 1 each | breadth only |

**8 segments per strain, always.** `CAT-03` now groups segments into one species, but
each strain still needs all 8 (or 7 for IB, 2 for ID/PB1-F2) or the k-mer space is
incomplete.

⚠ **Type by sequence, not organism string.** `OP212288` — A/Texas/61/2022, the first
dairy-cattle human H5N1 case, clade 2.3.4.4b B3.13 — is deposited with `(H3N2)` in its
organism field. Verified H5N1 human-case accessions: B3.13 = `PQ468757`–`PQ468764`;
D1.1 = `PQ573551`–`PQ573557`.

### 3.3 Other respiratory and enteric breadth (T1/T2)

Verified RefSeq anchors: RSV `NC_001781` (add A and B, 6 total); hMPV (4: A/B, 2
clades each); parainfluenza 1–4 (`NC_003461` present, expand to 10); rhinovirus
**by species not genotype** (A/B/C, 14 total — pan-RV primers are species-level);
adenovirus **by species, all 7** (11–15 total; species C/B/E are the respiratory ones;
98 % lower-respiratory AdV in children is C); bocavirus 1 (2).

**Clustering rule: within-subtype/species only, 95 % nt identity, 85 % coverage**
(`cd-hit-est -c 0.95 -aS 0.85`) for RNA viruses; 98 % for adenovirus.
Justification: CD-HIT was *excluded* from the van Bemmelen benchmark for runtime, so
there is no benchmark-derived CD-HIT threshold to cite — 95/85 is the closest documented
viral-derep default (`votuderep`, ANI ≥95 % + target coverage ≥85 %, the CheckV method)
cross-checked against RVDB's 98 % CD-HIT-EST and a virome paper's 95 %. The
directionally relevant cautionary result is PalmDB (*Nat Biotechnol* 2025): clustering
viral references at 99 % made **67.4 % of true taxa undetectable** vs 3.3 % for an
unclustered index, so they abandoned clustering and grouped by taxonomy instead —
which is what `CAT-03` already does. That experiment was in amino-acid space, so it is
directional, not quantitative, for us.

**Inclusion broad, exclusion only where provable.** Nothing in the literature supports
excluding an animal-reservoir virus that has caused a documented human case.

---

## 4. Anelloviridae (T3) — the deepest part of the design

### 4.1 What we measured this session

Built the exact-31-mer panel from the 2,042 shipped genomes and screened the EBV LCL
sample `SRR12682296` (EBV-immortalised LCLs, 10x v2; R1 = 26 bp barcode, R2 = 98 bp
UMI+insert), on 2,000,000 R2 reads:

| Probe | Result |
|---|---|
| EBV `NC_007605.1` (288,566 unique 31-mers) | **29,207 reads = 1.46 % of R2** — strong positive, file and method validated |
| TTV prototype `NC_002076.2` (6,392 unique 31-mers) | **0 reads** |
| All 2,042 anelloviruses, no complexity filter (9,690,411 31-mers) | 7,175 / 500,000 = **1.44 %** — apparently a strong positive |
| …same panel, low-complexity-masked (maxrun ≤ 11, ≥3 distinct bases, no tandem repeat) | **0 / 6,437 captured hit reads** |

Probe validation: exact recovery on synthetic controls at 1 %, 0.1 % and 0.02 %
abundance (1000/1000, 100/100, 20/20) with zero false positives among ~100,000 random
reads.

**The 1.44 % was 100 % artifact.** 89.5 % of the captured hit reads carry a
homopolymer run ≥ 31 bp (median longest run 49, max 97) — the poly-A/poly-T tails that
dominate 10x R2. 5,320 low-complexity k-mers in the panel (all C-rich) were matching
them. After masking, not one hit read had a single genuine anellovirus k-mer.

### 4.2 There is no conserved nucleotide k-mer space (G5)

| k-mers found in N of the 2,042 genomes | count |
|---:|---:|
| 1 | 8,789,070 |
| 2–4 | 809,350 |
| 5–19 | 79,974 |
| 20–49 | 5,866 |
| 50–99 | 1,930 |
| 100–199 | 918 |
| 200–399 | 396 |
| 400–999 | 54 |
| **≥1,000** | **0** |

**No 31-mer is shared by half the panel.** Genera share at most ~44 % *ORF1 amino-acid*
identity. A nucleotide k-mer panel therefore cannot be genus-sensitive: a read is
detectable only if it matches the specific genotype it came from. This is why
`CAT-11` measured `P(90 bp fragment captured)` per genus at 0.788 (Alpha) down to
**0.064 (Gyrovirus)**, and `CAT-12`'s per-genus bars (Beta 0.498, Gamma 0.332,
Samek 0.255, Mem 0.247).

**The only genus-sensitive route is translated search on ORF1 protein.** kallisto
supports `--aa` (PalmDB used it for exactly this: 4.4 GB RAM, 593 MB index with a
D-list). Constraints to respect: `--parity single` is mandatory for `--aa`; `--aa`
implies `--workflow=custom`, which **disables host masking unless `--d-list` is passed
explicitly** (§5.3).

### 4.3 The panel is probably structurally incomplete (G6)

| Query (`db=nuccore`, 2026-09-27) | Count |
|---|---:|
| `Anelloviridae[Organism]` | **42,755** |
| …organism string naming an ICTV genus | 21,283 |
| …organism string naming **no** genus | **21,472** |
| `Anelloviridae[Organism] AND srcdb_refseq` | 178 |
| **RefSeqViral** | **0** |
| `"unclassified Anelloviridae"` (taxid 363628) | 21,463 |

Entrez `[Organism]` matches the taxonomic lineage, so a query built on **genus names
structurally cannot reach the ~21,000 pre-2005 records** whose organism string is
`Torque teno virus` (12,382), `Torque teno mini virus` (6,718), `Torque teno midi
virus` (5,735) or `SEN virus` (303). If the 2,042 panel was genus-derived, it missed
them. **Re-pull by family taxid 687329 / `Anelloviridae[Organism]`**, and filter for
human host at the **flatfile `/host` qualifier** level — `"Homo sapiens"[Organism]`
returns 0 for anelloviruses and will mislead.

ICTV **MSL41** now has **37 genera / 243 species** (was 14 genera pre-2021, 34 genera /
173 species in MSL39.3, +4 genera ratified Mar 2025). Demarcation: genus from ORF1
*aa* phylogeny; species from ORF1 identity **< 69 %**. NCBI Taxonomy is ahead of ICTV
(660 species, incl. `Betatorquevirus sp. 'homini939'`). Every one of the 243 species has
a RefSeq exemplar — **that is the curated-provenance backbone to seed with.**

Adjacent human ssDNA viruses to add explicitly (they will not arrive from a genus
query): **TTMV** (Betatorquevirus, 6,718 records), **TTMDV** (Gammatorquevirus, 5,735),
**SEN virus** (taxid 136966, outside every genus, 303 records), **human circovirus**
(4 genomes: `NC_135454.1`, `PP968832.1`, `OR905605.1`, `OZ282268.1` — Circoviridae,
formally not Anelloviridae).

⚠ **Licence/provenance:** RefSeqViral is empty for this family; the 178 RefSeq records
are exactly the ICTV exemplars, and everything else is INSDC bulk-submitted
(public domain). **Panel composition is currently dictated by whoever bulk-submitted
last, not by ICTV or prevalence.** Seeding with the 243 ICTV exemplars and adding bulk
data on top fixes the provenance story; per-accession SHA-256 + retrieval date
(`REF-13`) is mandatory for anything redistributed.

**Genogroup is not recoverable.** ICTV abolished the concept; 0 of 2,042 records carry
`/genogroup`; the 2 that carry `/genotype` contradict their own organism. Legacy
"genogroup 1" is a many-to-one collapse of what are now several ICTV species. Assign
species by **ORF1 identity against the 243 exemplars** (the SCANellome V2 method:
CD-HIT 90 % nt on complete ORF1 → representatives), and **ignore `/genotype` entirely.**

### 4.4 What to report, and why the null is the expected answer

TTV in an immunocompetent adult: **~10²–10³ copies/mL plasma** (*Viruses* 17(2):140,
2025, PMID 40143262), prevalence **>90 %** (PMID 12721794); "probably the entire human
population is AV infected" (Kaczorowska & van der Hoek, *FEMS Microbiol Rev* 44:305,
2020). TTMV prevalence 40.1 % (PMID 30016934).

But scRNA-seq is 3′-biased and polyA-selected, TTV is ssDNA, and the panel has no
conserved k-mer space. The measured EBV LCL result agrees: **EBV 1.46 % of R2 reads,
zero genuine anellovirus reads.** So:

- **Report anellovirus as presence-only, binary, `screening_only`** (already `REF-10`).
- **Attach the expected-value calculation to every negative**, do not publish a bare
  zero. This is the same discipline already required for the tonsil EBV screen
  (`SENS-06`).
- **Never publish per-genotype anellovirus abundance.** The 87 % of annotated genomes
  carrying a single CDS spanning the whole genome are still competition buckets
  (`ANELLO-12`).
- Add a `detection_bound` column so a zero is quantitative, not silent.
- Genus-level resolution, if wanted, comes from the **translated ORF1 index** (§4.2) —
  never from the nucleotide panel.

### 4.5 Proposed T3 build

| Set | n | Source | Purpose |
|---|---:|---|---|
| Shipped anellovirus panel | 2,042 | `anellovirus_accessions.tsv` | continuity |
| ICTV MSL41 exemplars | 243 | MSL41 OLS `has narrower match` | curated provenance backbone |
| Re-pull by family taxid 687329, human host, complete genome | **TBD** (est. 20–40 k) | E-utilities + flatfile `/host` | close the 21k-record gap |
| SCANellome V2 representatives | 3,864 | Zenodo `10.5281/zenodo.7937276` | cross-check against a second community set |
| TTMV / TTMDV / SENV / human circovirus | ~13 k records, 4 in-panel | NCBI | completeness |

**Dedup: do not cluster.** `CAT-11`/`REF-01` already established that 20 genomes retain
1.1 % of k-mer space, 204 retain 11.3 %, 2,042 retain 100 % ⇒ keep all. Adding the
21k-record re-pull is ~60 Mbp — still ~2 % of the host genome. **Cheap.**

---

## 5. Index and quantification architecture

### 5.1 Two index tiers, not one equivalence-class space

The 99.8 % single-bucket failure (`MW455439.1_gene1`, 1,167,103 / 1,169,272 anellovirus
UMI reported as "Alphatorquevirus") is a *competition-bucket* failure: whole-genome
pseudo-transcripts and real CDS transcripts for the same virus compete in one EC space.
Fix: **two `kb ref` builds and two `kb count` passes on the same reads.**

| Index | Contents | Used for | Never used for |
|---|---|---|---|
| `viral_gene` | real CDS transcripts only (GenBank path, genome-scoped IDs) | gene programmes, latent/lytic, gene-level breadth | genome-level abundance |
| `viral_genome` | one whole-genome pseudo-transcript per accession | presence/absence, `detection_bound`, `accession_breadth` | per-genome abundance, gene programmes |

Both write into the same v3 output schema with a `tier` column. This is the structural
version of `CAT-01`, which fixed annotation but not the EC collision.

### 5.2 Expansion is not monotone-good (G4)

- kallisto k is hard-capped at **31** (`k > 31` needs `-DMAX_KMER_SIZE=64`). Genomes
  differing by ~30–70 nt share k-mers in conserved regions.
- In the scRNA-seq path kallisto/bustools **discard UMIs that map to multiple genes**;
  `--multimapping` distributes them uniformly instead. So growing a near-identical panel
  converts viral reads into **dropped** UMI. Assert `max EC size` and
  `number of ECs discarded` from `kallisto inspect` on every build.
- `--distinguish` (custom workflow, zero-indexed numeric target names) is the
  purpose-built escape hatch for exactly this — separates k-mers by target instead of
  collapsing them. **Worth one experiment** before committing to a per-virus cap.
- `kb count --multimapping` (not `--em`) is the correct setting; kallisto's authors
  advise against `--em` in scRNA-seq.

### 5.3 Build mechanics that will silently break (G8)

- `--workflow=custom` / `--aa` **sets `dlist = None`** → **no host masking at all**
  unless `--d-list` is passed explicitly. Must pass host genome **and** transcriptome.
- A D-list FASTA **with** a header filters only distinguishing flanking k-mers;
  **without** a header it filters every k-mer in the file.
- kallisto **clips poly-A tails longer than 10 nt at index build** — harmless, but it
  makes 5′/3′ UTR k-mers asymmetric. Coronaviruses are the pathological case, and
  3′-biased scRNA-seq reads depend on the 3′ end.
- Default `--workflow=standard` extracts only cDNA from GTF-annotated CDS, losing all
  UTR and intergenic sequence. **Use `--workflow=custom` with full viral genomes.**
- `kb count -m` "may also cause it to crash due to memory."
- Duplicate FASTA IDs or duplicate sequences break `kallisto index` outright (it
  happened once already with `NC_002076.2`). `CAT-05` must run on **every** build path,
  including `scripts/build_bundled_panel_ref.py`, which currently has no guard.

### 5.4 Low-complexity masking — new build step (G3)

Required at index-build time, based on §4.1:

- Drop panel k-mers with homopolymer run > 11, fewer than 3 distinct bases, or a
  perfect tandem repeat of a unit ≤ 5 bp.
- Calibrate on the host transcriptome: kallisto's own guard is `max 10` nt, so our
  threshold is deliberately slightly looser.
- **Add a regression test:** build a small panel containing a C-rich/homopolymer
  anellovirus region, assert synthetic poly-A reads yield **0** hits. Without this, the
  next panel expansion silently reintroduces a 1.4 % false-signal floor.
- Report the masked-k-mer count per accession in the reference manifest — a genome that
  loses a large fraction of its k-mers is a degenerate annotation and should be flagged,
  not silently shipped.

---

## 6. Reporting contract changes

1. `tier` column (`gene` / `genome`) on every viral row — the reader must be able to
   tell which equivalence-class space produced the number.
2. `detection_bound` for every negative (see §4.4).
3. Genus-level anellovirus rollup **only** from the translated index; the nucleotide
   panel reports family/genus presence, never abundance.
4. `expected_value_note` for anellovirus and for any low-prevalence latent agent in a
   tissue where E[infected cells] < 1 — the discipline already required for tonsil EBV
   (`SENS-06`).
5. `accession_breadth` recomputed over **all** index genes of the virus, not detected
   ones (`ANDET-01`).
6. Corpus-level negative-control reporting stays mandatory: a host-only sample must
   produce no call (`CAT-14`).

---

## 7. Validation

### 7.1 Acceptance criteria

| # | Criterion |
|---|---|
| A1 | All 18 Serratus latent-DNA accessions present, gene-level annotated, ≥1 real multi-exon gene each |
| A2 | Zero unnamed gene IDs (already 0 after `CAT-03`; must stay 0) |
| A3 | Synthetic poly-A/tandem-repeat reads → **0** viral calls (§5.4) |
| A4 | Host-only negative → **0** calls (`CAT-14`) |
| A5 | `max EC size` and discarded-EC count recorded and within budget on every build |
| A6 | `measure_kmer_capture.py` per-genus `p_fragment` ≥ 0.95 and zero-coverage ≤ 5 % (`CAT-12`) |
| A7 | EBV LCL `SRR12682296` re-run reproduces `PROG-07` (895 latent / 236 productive / 311 mixed / 3,094 indeterminate) (`CAT-15`) |
| A8 | Anellovirus reported as expected-negative **with** the §4.4 calculation attached |

### 7.2 Truth panel (extends `WP5`)

Needs ≥1 orthogonally confirmed positive per persistence class, and must include a
**known TTV-positive specimen** — the field's own method is rolling-circle amplification
+ NGS with ORF1 phylogenetic assignment (up to 15 TTV species resolved in a single
plasma sample, PMID 38543797). Without one, every anellovirus result stays
`screening_only` (`REF-10`) and no code change lifts that.

### 7.3 Re-scope `CAT-13`

The size gate should assert on the **viral payload** (gzipped) and keep the 60 MB PyPI
limit in mind, but should record that the *real* build cost is the **host genome
D-list** (3.15 Gb, ~64 GB RAM, ~8 h, `scripts/build_genome_panel_ref.sh`), not the virus.
Decide explicitly: genome-D-list (best artifact control, expensive) vs
cDNA-D-list (cheap, but `F-005` showed it manufactures ~90 % spurious Alphatorquevirus
signal).

---

## 8. Open questions for the user

| # | Question | Default if unanswered |
|---|---|---|
| **Q1** | **SARS-CoV-2 lineage genomes — in or not?** RefSeq gives exactly 1. Getting 4–6 means INSDC/ENA + local lineage assignment, a GISAID DUA dependency, and a non-reproducible CI. | Ship `NC_045512.2` only; document the limit. Detection is unaffected (lineage assignment is not a ViralScan feature); only re-analysis of historical datasets benefits. |
| **Q2** | **Genus-name anellovirus queries → re-pull by family taxid 687329?** Est. 20–40 k new records / ~60 Mbp. Changes the panel from "clareaulab collection" to "family census". | Yes, re-pull; store every accession with SHA-256 + retrieval date. |
| **Q3** | **Two `kb ref` builds / two `kb count` passes** (double counting cost, structural fix for the 99.8 % bucket) vs one build with a tier column in t2g (cheaper, leaves the EC collision)? | Two builds. The single-bucket failure is the most damaging defect in the panel's history. |
| **Q4** | **Translated ORF1 index for anellovirus genus resolution** — extra index, `--parity single`, extra D-list handling. In or defer? | Defer to a separate subcommand; ship the nucleotide panel as presence-only first. |
| **Q5** | **Sequence bundling** — `WP4H` decision 3 bundles sequences in the package (self-contained), against PR 8's move-data-off-package direction. Viral payload ≈ 30 Mbp gzipped. Confirm. | Bundle Tier 2 + anellovirus catalogues; Tier 1 sequences via the fetch path. |

---

## 9. Sequencing

Unchanged blocking order from `WP4H` (`CAT-01`, `CAT-03` first), with these inserted:

| Order | Row | Depends on |
|---|---|---|
| 1 | `CAT-01`, `CAT-03` *(done)* | — |
| 2 | **`CAT-17`** low-complexity masking + poly-A regression test | nothing — **do this first**, it is cheap and it is a live correctness bug |
| 3 | `CAT-05` duplicate guard on every build path | — |
| 4 | `CAT-09` T1 breadth floor (adds SARS-CoV-2, HIV, HTLV, OC43/NL63/HKU1, hMPV, bocavirus, flu D, TSPyV, HPyV6/7) | `CAT-01` |
| 5 | `CAT-10` T2 persistence set, with the §2.2 Serratus list as an explicit row | `CAT-09` |
| 6 | **`CAT-18`** clade-level influenza budget (§3.2) replacing the "one per subtype" wording | `CAT-09` |
| 7 | **`CAT-19`** two-index architecture (§5.1) | `CAT-05` |
| 8 | **`CAT-20`** anellovirus family-taxid re-pull + ICTV MSL41 exemplars + SENV/TTMV/circovirus (§4.5) | `CAT-09` |
| 9 | `CAT-14` host cross-talk gate | `CAT-09`, `CAT-20` |
| 10 | `CAT-12` measured capture bars, **now against the two-index design** | `CAT-18`, `CAT-19`, `CAT-20` |
| 11 | `CAT-15`, `CAT-16`, `CAT-13` | all |

`CAT-11` is done and its numbers are the input to `CAT-12`; note the correction it
already made (analytic `p_fragment` 1.0000 vs **measured 0.5097**) — `CAT-12` must be
gated on the measured value.

---

## 10. Sources

**Repo:** `PLAN.md` (`WP4G`, `WP4H` incl. expansion, `WP4I`, `WP4J`; `REF-01`, `REF-11`,
`REF-13`, `ANELLO-12`, `ANDET-01`, `CAT-01`…`CAT-16`, `SENS-01`, `SENS-06`, `F-005`),
`docs/reference_panel.md`, `docs/faq.md`, `docs/cli_reference.md`,
`src/viralscan/{anellovirus,constants,virus_grouping,virus_catalog}.py`,
`src/viralscan/scripts/{ncbi_fetch,build_reference,analysis}.py`, `scripts/build_bundled_panel_ref.py`,
`extras/build_*.py`, `.living/findings/`, `FINDINGS_REGISTRY.md`.

**Measured in this session** (`~/ttv_probe/`, not committed): 31-mer panel
construction over 2,312 genomes; EBV `SRR12682296` R2 screens; synthetic positive
controls at 1 %/0.1 %/0.02 %; low-complexity re-scoring of 6,437 hit reads.
Scripts: `ttv_kmer_probe.py`, `probe.py`, `anello_full.py`, `anello_filtered.py`,
`build_complex_panel.py`, `rescore_hits.py`, `final_scan.py`.

**Literature:**
Lareau CA et al. *Nature* 623:608–615 (2023) doi:10.1038/s41586-023-06704-2 ·
Edgar RC et al. *Nature* 602:142 (2022) doi:10.1038/s41586-021-04332-2 ·
Traylen A et al. *Future Virology* 6:451 (2011) doi:10.2217/fvl.11.21 ·
Laubscher F, Kaiser L, Cordey S. *Viruses* 16:1349 (2024) PMID 39339826 ·
Kraberger S et al. ICTV Anelloviridae chapter, *J Gen Virol* 107(3) (2026) doi:10.1099/jgv.0.002222 ·
ICTV MSL41 · varsani S et al. *Arch Virol* 166:2943 (2021) PMID 34383165 ·
van Bemmelen M, Nika M, Baaijens SA. *BMC Genomics* (2026) doi:10.1186/s12864-026-12874-w ·
PalmDB, *Nat Biotechnol* (2025) doi:10.1038/s41587-025-02614-y ·
Brani et al. *Viruses* 17(2):140 (2025) PMID 40143262 ·
Kaczorowska A, van der Hoek L. *FEMS Microbiol Rev* 44:305 (2020) PMID 32188999 ·
PMID 12721794 (prevalence) · PMID 30016934 (TTMV prevalence) · PMID 38543797 (RCA-NGS,
15 species/sample) · PMID 18094586 (genogroup 5-assay survey) ·
WHO SARS-CoV-2 variant tracking (page updated 2026-07-28) · WHO influenza fact sheet
(2026-01-21) · Nextstrain ncov clades · kallisto 0.50 docs (k≤31, EC size, `--distinguish`).
