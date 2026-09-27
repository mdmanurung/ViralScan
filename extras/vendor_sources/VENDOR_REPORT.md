# Vendor Sources — Provenance & Measurement Report

Vendored external reference assets for the ViralScan reference-panel work.
Everything in this directory is a **byte-exact upstream copy** pinned to a specific
commit. Nothing here was edited, reformatted, or "fixed"; the one known data defect
found upstream is documented below and deliberately left intact.

- Generated: 2026-09-27
- Manifest: [`VENDOR_MANIFEST.tsv`](VENDOR_MANIFEST.tsv) (41 rows + header, 7 columns)
- Fetch method: `curl -sSL --fail --retry 3` against `raw.githubusercontent.com`
- All file digests are SHA-256, computed after download from the bytes on disk.

## Integrity verification

| Check | Result |
|---|---|
| Files downloaded | **41** |
| Download failures (HTTP 4xx/5xx) | **0** |
| Manifest rows | 41, all 7 columns populated |
| SHA-256 recomputed from disk vs. manifest | 41/41 match |
| Byte size on disk vs. manifest | 41/41 match |
| Byte size vs. GitHub tree-API blob `size` | **41/41 match** — all files are byte-complete, none truncated |
| Pinned-commit assertion (`trees` API `sha` == pinned ref) | 5/5 match |

The last row is the load-bearing one: it re-derives every expected file size
independently from the GitHub tree API and confirms the downloaded bytes are complete.

## Totals

| Source repo | Pinned commit | Files | Bytes | MiB |
|---|---|---:|---:|---:|
| `clareaulab/human_anellovirus_pangenome` | `3ed77e19` (supplied) | 22 | 12,135,246 | 11.57 |
| `caleblareau/pan-viral-reactivation` | `74136de5` (supplied) | 8 | 3,447,879 | 3.29 |
| `clareaulab/ad-hsv-mapping` | `53801369` (HEAD, 2026-08-26) | 4 | 349,459 | 0.33 |
| `yyoshiaki/VIRTUS3` | `b7873791` (HEAD, 2026-03-16) | 2 | 164,405 | 0.16 |
| `huangyh09/ViralScan` | `d8279c37` (HEAD, 2025-08-06) | 5 | 15,571,256 | 14.85 |
| **TOTAL** | | **41** | **31,668,245** | **30.20** |

### How the three dynamic pins were resolved

For sources 3–5 the pin is `GET /repos/{repo}/git/trees/HEAD?recursive=1` → `.sha`.
Worth recording: that `sha` field is **byte-identical to the HEAD commit sha**
(`/repos/{repo}/commits?per_page=1` → `.sha` returns the same 40 hex chars), and
`raw.githubusercontent.com` accepts it as a ref. So the tree endpoint's `sha` is usable
as a commit pin — no extra resolution step was needed. Each resolved tree was also
non-truncated, so no pagination follow-up was required.

## Missing files / 404s / substitutes

**None.** All 42 requested paths were requested minus one consolidation (below); 41
files were downloaded and **0 returned 404**. No substitute path was needed anywhere.

- Source 4 (`yyoshiaki/VIRTUS3`): the task warned `data/NC_007605.1_CDS_EBER12.fa`
  "may not exist at that exact path." It **does** exist, at exactly that path, at the
  pinned commit. The tree listing of `data/` was fetched anyway and shows the two
  relevant EBV files: `NC_007605.1_CDS_EBER12.fa` (160,472 B) and
  `NC_007605.1_CDS_EBER12.tgMap.tsv` (3,933 B). There is also a sibling
  `NC_007605.1_CDS.fa` (CDS-only, 94 records) and a `*_salmon_index/` directory; the
  CDS-only variant was **not** downloaded as it was not requested, but note the EBER12
  file is precisely that file plus the 2 EBER records (see below).
- Source 5 (`huangyh09/ViralScan`): `viral_reference/` contains exactly 4 blobs, all
  far under the 50 MB cap (7.58 MB, 60.5 KB, 7.60 MB, 329 KB), so **all 4 were
  downloaded and nothing was skipped**.
- Requested-path arithmetic: 22 + 8 + 4 + 2 + (1 script + 4 directory) = 41 files.

---

# Measurements

## 1. `anello/ref/hardmasked_cdhit_rep_anelloviridae_061725_genomic.fa`

The hardmasked Anelloviridae pangenome (2023 CD-HIT representative genomes).

| Metric | Value |
|---|---|
| Contig count (`grep -c '^>'`) | **2023** (independent Python FASTA parse also 2023 — agrees) |
| Total bp | **5,927,006** |
| Min contig length | **1,141** bp |
| Median contig length | **2,895** bp (exact; n=2023 is odd so this is a single observed value) |
| Max contig length | **3,996** bp |
| Total N count | **284,762** |
| N as % of bp | 4.80 % |
| Contigs containing ≥1 N | 2,005 of 2,023 (18 are N-free) |
| GC | 38.17 % |

First 5 contig IDs: `PQ438050.1`, `PP728782.1`, `PP728781.1`, `PP728779.1`, `PP728776.1`.

**Cross-check:** `ref/anello_sizes.txt` has exactly 2023 entries and every one matches
the FASTA length for the same contig ID — the two shipped reference files are
self-consistent.

## 2. `anello/ref/anello_t2g.txt`

kb-style target-to-gene file. **No header row**; 8 tab-separated columns, consistently
8 on all 2023 lines.

| Metric | Value |
|---|---|
| Line count (non-empty) | **2023** |
| Distinct `target_id` (col 1) | 2023 |
| Distinct `gene_id` (col 2) | 2023 |
| **Distinct accession-ish prefixes** | **2023** |

"Prefix" is resolution-dependent, so all three readings are given:

- Accession **with** version (e.g. `PQ438050.1`): **2023** distinct.
- Accession **minus** version (e.g. `PQ438050`): **2023** distinct.
- Also collapsing the RefSeq `NC_` underscore (e.g. `NC`): **2018** distinct.

The 2023 → 2018 gap is fully explained: exactly 6 target IDs carry an `NC_`-style
accession — `NC_038359.1`, `NC_038351.1`, `NC_038350.1`, `NC_038346.1`, `NC_038345.1`,
`NC_038337.1`. Every other accession is a letter-prefixed 6-digit form with no
underscore, so it is invariant under the split. Bottom line: **there is one target per
contig, and no accession collisions** — the file maps 1:1 onto the 2023 pangenome
contigs. (Contrast with the panviral file below, which *does* have a collision.)

Sample row: `PQ438050.1_transcript  PQ438050.1  PQ438050.1  (blank)  PQ438050.1  1  3062  +`

## 3. `anello/simple_anello_metadata_V2.csv`

| Metric | Value |
|---|---|
| Data rows (excluding header) | **3,545** |
| Column count | **18** (all 3,545 rows have exactly 18 fields — no ragged rows) |

Columns, in order:

1. `Accession`
2. `Species`
3. `Genus`
4. `Family`
5. `Isolate`
6. `Virus Taxonomic ID`
7. `Virus Name`
8. `Length`
9. `Isolate Lineage source`
10. `leiden`
11. `genus`
12. `species`
13. `infer_genus`
14. `cdhit_representative`
15. `orf_id`
16. `aa_length`
17. `phylo_cluster`
18. `orf1_genus`

Note 3,545 rows ≫ 2,023 contigs: this table is at **ORF granularity** (it carries
`orf_id` and `aa_length`), not contig granularity. It also mixes taxonomic columns
(`genus`/`species`) with assigned-cluster columns (`leiden`/`phylo_cluster`) and
`infer_genus` — worth reading `docs/plans/2026-09-27-viral-reference-panel-expansion.md`
before treating `genus` as authoritative.

## 4. `anello/example/expected_output/full_SRR32170409/anello_detected.tsv` — KNOWN POSITIVE CONTROL

Full contents (5 lines: 1 header + 4 detections, 311 bytes):

```tsv
target_id	genus	Species	est_counts	tpm	phylo_cluster
MH649122.1	Alphatorquevirus	Anelloviridae sp.	2308.89	850801	29
MZ286128.1	Alphatorquevirus	Anelloviridae sp.	234.83	94884.2	29
MW455390.1	Alphatorquevirus	Anelloviridae sp.	96.2165	45304.2	29
KP343840.1	Alphatorquevirus	Torque teno virus	16.0608	7727.26	29
```

All 4 positives are genus **Alphatorquevirus** (torque teno virus) in phylo cluster
**29**, spanning ~2.7 orders of magnitude in `est_counts` (2308.89 → 16.06). Three of
the four are annotated only to `Anelloviridae sp.` at species level, so **genus-level
annotation, not species-level, is what this control actually tests.**

Companion `run_info.json` for the same run (SRR32170409): 36,109,130 reads processed,
2,659 pseudoaligned (**0.0074 %**), 2,084 unique, kallisto 0.52.0, k=31, index v13.
This is a very low pseudoalignment rate, so it is a good null-ish control as well as a
positive one: a pipeline that emits the 4 Alphatorquevirus calls *and* stays quiet on
the other 2,019 contigs is behaving correctly.

## 5. `panviral/reference/pan_virus_annotation_plain.tsv` — CORRUPTION CONFIRMED

| Metric | Value |
|---|---|
| Rows (excluding header) | **724** |
| Columns | `EC`, `Gene`, `Nuccore`, `Virus` |
| Distinct `Virus` values | **9** |
| Distinct `Nuccore` values **as shipped** | **10** |
| Distinct `Nuccore` **after stripping the stray bracket** | **9** |
| Malformed lines | **1** |

`EC` is 0-based and strictly complete: 0 … 723, 724 unique values, no gaps, no
duplicates.

### The malformed line, verbatim

The single defective row is **file line 532** (the `EC == 530` record):

```
530	Jvgp6	[NC_001699.1	Jcpolyomavirus
```

Split on tab, the four fields are:

| # | Field | Value |
|---|---|---|
| 1 | `EC` | `530` |
| 2 | `Gene` | `Jvgp6` |
| 3 | `Nuccore` | **`[NC_001699.1`** ← stray leading `[` |
| 4 | `Virus` | `Jcpolyomavirus` |

**Your suspicion is correct, and it is exactly one row — no others.** Scanning all
724 data rows for `[`, `]`, or leading/trailing commas returns precisely this one
line. No other field in the file is malformed.

### Why it is a corruption, and the likely mechanism

The `Virus → Nuccore` mapping is otherwise strictly 1:1, and this one row is the
violation:

| Virus | Accession | Gene records |
|---|---|---:|
| AAV2 | `NC_001401.2` | 9 |
| HHV1 | `NC_001806.2` | 77 |
| HHV3 | `NC_001348.1` | 73 |
| HHV4 | `NC_007605.1` | 94 |
| HHV5 | `NC_006273.2` | 169 |
| HHV6B | `AF157706.1` | 103 |
| HHV7 | `U43400.1` | 107 |
| HHV8 | `MK733606.1` | 86 |
| Jcpolyomavirus | `NC_001699.1` | 5 |
| Jcpolyomavirus | **`[NC_001699.1`** | **1** ← same virus, 6th gene |

So `Jcpolyomavirus` has 6 gene records, of which one carries the corrupted
accession. That is why the as-shipped distinct-accession count is 10 instead of 9:
`[NC_001699.1` and `NC_001699.1` are two *strings*, so they count as two accessions
for one virus. Stripping the bracket restores the clean 9 viruses ↔ 9 accessions
bijection.

The sibling file `reference/pan_virus_annotation.txt` (the GenBank feature-table
variant, 724 lines) contains the same record at the same `EC` and shows the correct
value — record `EC 530` is:

```
>lcl|NC_001699.1_cds_NP_043513.1_6 [locus_tag=Jvgp6] [db_xref=GeneID:1489521] [protein=hypothetical protein] [protein_id=NP_043513.1] [location=complement(4495..5013)] [gbkey=CDS]
```

Note `NC_001699.1` is clean there, and the record's distinguishing feature is that its
accession is immediately followed by `_cds_...` after the standard `lcl|` prefix. That
strongly suggests a one-off regex/parse slip when the plain TSV was derived from the
bracket-annotated feature table: the converter was stripping `[...]` attribute blocks,
and on this single row a `[` leaked into the extracted `Nuccore` field.

**Recommendation:** do **not** ship this file as-is to anything that joins on
accession. A stray `[` in a join key will silently produce a *missing* join rather
than an error, so the JCV `Jvgp6` gene will vanish from any annotation-join downstream
with no warning. Repair is a one-character strip of `[` on line 532. This is recorded
here and in the manifest `note`; the vendored copy is left byte-identical to upstream
so the defect stays auditable.

## 6. `hsv/reference-genomes/HSV1-LATonly.fasta` — HSV-1 LATENCY-ONLY REFERENCE

| Metric | Value |
|---|---|
| Sequence count | **3** |
| Total bp | **8,347** |
| GC | — |

| # | Header | Length | Locus span (NC_001806.2) |
|---|---|---:|---|
| 1 | `NC_001806.2:118805-119464:+:exon-HHV1gp00s02-1` | **660** | exon 1 of ICP0/LAT |
| 2 | `NC_001806.2:119465-121420:+:intron-HHV1gp00s02-1` | **1,956** | intron 1 |
| 3 | `NC_001806.2:121421-127151:+:exon-HHV1gp00s02-2` | **5,731** | exon 2 |

All three records are the same locus, gene `HHV1gp00s02` = **ICP0 / LAT**, forward
strand, spanning 118,805–127,151 contiguously.

**Important caveat before this is used as a reference panel entry:** the file is *not*
a contiguous genomic segment and *not* a spliced transcript — it is **3 separate FASTA
records, one of which is an intron**. Two consequences:

1. The 1,956 bp intron is present as its own k-mer-countable record. Reads from
   unspliced/latent HSV-1 transcripts will map to it, so this "LAT-only" reference is
   not restricted to the spliced LAT mRNA despite the name.
2. If a downstream tool assumes one-record-one-locus, it will treat 3 distinct
   targets. The S2G/t2G companion (`viral_reference/HSV1-LATonly.idx`, also in the
   upstream repo) encodes whatever convention upstream used — check it before
   assuming.

For contrast, the two coding references in the same directory: `HSV1-coding.fasta`
= 77 records / 123,288 bp, and `VZV-coding.fasta` = 73 records / 113,370 bp. The 77
and 73 gene counts match the HHV1 and HHV3 rows of the panviral annotation table
exactly (see §5) — a useful independent cross-check that these repos were built
against the same RefSeq annotation release.

## 7. `virtus3/data/NC_007605.1_CDS_EBER12.fa` — EBV (HHV-4), EBER-AUGMENTED

| Metric | Value |
|---|---|
| Total sequences | **96** |
| **EBER / non-coding** | **2** (2.08 %) |
| Protein-coding CDS | 94 (97.92 %) |
| Distinct gene symbols | 86 |

The two EBER records, in full:

| Header | Length | First 40 nt |
|---|---:|---|
| `rna-HHV4_EBER-1 [gene=EBER-1]` | **167** | `AGGACCTACGCTGCCCTAGAGGTTTTGCTAGGGAGGAGAC` |
| `rna-HHV4_EBER-2 [gene=EBER-2]` | **173** | `AGGACAGCCGTTGCCCTAGTGGTTTCGGACACACCGCCAA` |

Both are named by the `rna-` prefix (vs. `lcl|..._cds_...` for the CDS records), carry
`[gene=EBER-1]` / `[gene=EBER-2]`, and are appended at the **end** of the file after
the 94 CDS records. They are small RNAs, so they have no `gbkey=CDS` and no
`protein_id`. This is the file's whole point: it is the stock `NC_007605.1_CDS.fa`
(94 CDS) **plus** EBER1 and EBER2, so that an EBV quantifier catches EBER expression —
which is diagnostically relevant and invisible to a CDS-only reference.

**Cross-check:** 94 CDS here exactly equals the 94 HHV4 gene records in the panviral
`Virus` table (§5). Both repos agree on the EBV CDS set.

Duplicate gene symbols (gene symbols are **not** unique — 96 records, 86 symbols):

| Symbol | Records | Why |
|---|---:|---|
| `BWRF1` | 8 | EBV BWRF1 has multiple predicted isoforms (several `partial=5'` records) |
| `LMP-1` | 3 | BNLF2a / BNLF2b / BNLF2c are three differently-spliced LMP-1 products |
| `EBNA-3B/EBNA-3C` | 2 | 3B and 3C come from one bicistronic transcript |

So use the FASTA record ID (or the tgMap col-1 IDs) as the join key, **not** the gene
symbol — joining on `[gene=]` will fan out and double-count.

`NC_007605.1_CDS_EBER12.tgMap.tsv` (3,933 B) matches the FASTA exactly: **96 lines**,
no header (line 1 is `lcl|NC_007605.1_cds_YP_401631.1_1  LMP-2A`), of which **2** are
EBER rows. Target IDs align 1:1 with the FASTA record headers.

---

## Provenance caveat

These are third-party research repositories vendored as **reference inputs**, not as
code we run. Note that `anello/ref/kallisto_idx.sh` hardcodes a `module load
lareauc/kallisto/0.48.0` and the `run_info.json` call string embeds an absolute path
under `/Users/gutierj6/...` — these are records of how upstream built the panel, kept
for reproducibility context, not runnable scripts. Reproduce the index with your own
kallisto version (upstream used 0.52.0 / k=31 / index v13) if you need bit-comparable
numbers to the `anello_detected.tsv` control.
