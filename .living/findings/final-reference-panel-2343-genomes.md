# F-016 — The final reference panel, built: 2,343 genomes, and what it cost to get

**Date:** 2026-09-28
**Status:** confirmed (panel built; index in progress at time of writing)
**Tags:** reference-panel, build, CAT-17, CAT-31, CAT-05, EBV, HPV, kb, dustmasker, EBER

## What was built

`scripts/build_bundled_panel_ref.py` with the pinned toolchain, producing a
**2,343-genome** viral panel:

| Component | Count | Source |
|---|---:|---|
| GTF-backed non-anellovirus genomes | 323 | 195 bundled + 35 verified additions + 94 HPV |
| Anellovirus (clareaulab) | 2,020 | CAT-05 applied, 2,041-row table |
| **Total** | **2,343** | |

The 35 additions are the 34 catalogued-but-unindexed accessions of F-015 plus
EBV type 2. All 15 IARC high-risk HPV types are now present.

## CAT-31: the reconciliation guard works, and it found the CAT-05 pair

Running the new guard during the build:

```
catalogued accessions NOT in the panel : 2   (was 34)
  AB303562   UNEXPLAINED  Torque teno midi virus 12
  NC_038359  UNEXPLAINED  Torque teno midi virus 12
indexed but absent from the catalogue  : 96
```

The 35 additions closed **32 of the 34** gaps in one build. The 2 survivors are
precisely the CAT-05 duplicate pair, which is the correct outcome: a duplicate
genome must not be indexed twice.

**CAT-05 is settled.** Live NCBI (both 3,187 bp, taxid 2065053, isolate MDJN51,
identical md5 `5e78d2d32ea911d61962cd8ad6255a4f`; `NC_038359.1` carries an
`assemblyacc` of `AB303562`; the RefSeq COMMENT reads *"The reference sequence is
identical to AB303562"*). **Keep `NC_038359.1`, drop `AB303562.1`** — the RefSeq
copy has 4 curated `gene` features and the GenBank record has zero. Dropping the
GenBank twin costs no gene coverage. Table is now 2,041 rows.

The 96 "indexed but uncatalogued" are mostly the newly added INSDC-only HPV plus
`AF157706.1` (HHV-6B, the source of the 97 `HUM_HERP6B_*` pseudo-contigs). They
have no catalogue row, so no provenance or tier. The catalogue needs extending to
match the panel.

## Three latent bugs found on the way

1. **20 duplicate `transcript_id`s across the 195 bundled GTFs** — a real `kb ref`
   crash (`SegmentCollection.copy`). `unassigned_transcript_1` appeared in 20
   files, `ROTA_B_`/`ROTA_C_`/`ROTA_A_`/`BANNA_` in 10–11 each. Namespacing by
   genome accession fixed 301 occurrences across 111 GTFs; a re-scan confirms
   **0** remaining. This was a live crash waiting for a build that used them.

2. **EBER1/EBER2 were silently dropped.** `ncbi_fetch.py` filtered features with
   `if key == "CDS"`, discarding the `misc_RNA` EBER features that are the
   highest-abundance latent EBV transcripts. Live evidence: `NC_007605.1` carries
   EBER-1 at 6629..6795 (167 nt) and EBER-2 at 6956..7128 (173 nt). Now emitted;
   the regenerated EBV GTF has 96 genes (94 CDS + 2 EBER), matching VIRTUS3's
   independently built 96-record reference exactly.

3. **kallisto ≥0.50 writes a single-file index** — no `.td`/`.nd`. Any
   artifact check asserting those files exist is wrong, not the build.

## dustmasker cannot do the job the gate asks of it

The CAT-17 gate failed the build twice even with masking on. The cause is that
**dustmasker masks essentially nothing here**: it removed **785 of 11,048,660
bases (0.01 %)** and left all 85 pure-homopolymer 31-mers in `NC_001479.1`
intact. DUST scores compositional complexity over a window; it is not built to
strip a homopolymer tract, so the poly-A/poly-T 10x tail that manufactured
1.44 % of R2 reads (F-014) passes straight through it.

Two changes, in the builder, both measured:

- **Targeted homopolymer-run masking** after dustmasker: every run of ≥31
  identical bases becomes N. A "pure-homopolymer 31-mer" is by definition such a
  run, so this removes exactly the class the gate fails on. Panel result:
  `NC_001479.1` **85 pure → 0**; **panel-wide pure-homopolymer 31-mers = 0**
  across all 2,343 records, for 473 masked bases (0.004 %).
- **The fraction default `0.0 → 0.05`.** A tandem repeat is legitimate sequence;
  measured across the panel the worst record is at **0.0205** (median 0.0035,
  p90 0.0100). A 0.0 ceiling is unachievable by any masking and only guaranteed a
  red build. The primary control remains the *absolute* pure-homopolymer count,
  which is exactly what F-014 recommended.

The fraction gate is now a backstop on `long_run`/`few_bases`/`tandem`
(4,902 / 2,391 / 2,466 k-mers panel-wide), not the thing that gates the build.

## The HPV multimapping fear was wrong — measured

A subagent warned that adding ~160 HPV genomes would cause pervasive
cross-type multimapping, because alpha-papillomavirus L1 is conserved and
`kb count` discards multimapping UMIs. **Measured at k=31, this does not hold:**

```
unique-31-mer fraction   mean 0.9711   median 0.9806   min 0.8670 (HPV94)
k-mer space type-discriminating      98.68%
max pairwise 31-mer sharing          12.0%     max types per 31-mer  15
L1 vs non-L1 unique fraction         0.9838 vs 0.9683   (L1 is MORE specific)
100-nt read assignability            99.96% (worst type 99.15%)
types with <10% unique k-mers        0
```

HPV16/HPV18 are 51.5% identical at nucleotide level and **share zero 31-mers**;
54 of 78 high-risk pairs share zero. L1 is a conserved *protein*, not a
conserved *nucleotide* sequence. The closest pair of distinct types is 90.4%
identical, well clear of the >99 % regime where CAT-19 measured collapse.

Adopted: **109 human-pathogen types** (alpha + beta genera), dropping 73
gamma/mu/nu, which costs ≤0.36 pp of per-type assignability. Caveat: the
kallisto pseudo-inverse *threshold* effect is invisible to this metric, so the
build is gated on `kallisto inspect` (max EC size) rather than on k-mer overlap
alone.

## Provenance reality check on HPV

Only **68 of 182** HPV types have a curated RefSeq genome. **114 are INSDC-only**,
including **8 of the 15 IARC group-1 high-risk types** (39, 51, 56, 58, 59, 66,
68, 69). Those 8 are INSDC because RefSeq curation never reached them, not
because they are low quality. The panel includes them, with `source` and
`refseq` recorded per row in `hpv_types.tsv`, and `virus_catalog.tsv` needs the
same treatment before the HPV rows can be trusted as curated.

## Toolchain facts worth keeping

- `kb` **exits 0 after kallisto fails.** `build_reference.py`'s
  `CalledProcessError` handler never fires. Artifact checks are the only gate —
  never `$?`.
- **kallisto and bustools are vendored inside the `kb_python` wheel.** Pinning
  the conda `kallisto`/`bustools` specs is cosmetic; pin `kb-python` itself. The
  drifted env was running conda `bustools` 0.45.1 on PATH while `kb` actually
  executed the bundled 0.51.1.
- The pinned env resolved to **Python 3.8**, but the repo requires 3.9+
  (`dict[str, str]`). The working invocation is the repo's Python 3.12 with the
  pinned env's `bin` first on `PATH`, so `kb`/`kallisto`/`dustmasker` are pinned
  while the builder runs on a supported interpreter.
- The shipped covid `index.idx` is **healthy** (`kallisto inspect` OK, 470,468
  targets, positive control 900→842 pseudoaligned, host reads 0). The crash was
  in `kb ref` **rebuilding**: `combined.gtf` is *chromosomal* (seqnames 1..22,
  4.66 GB) paired with `combined.fa`, a *cDNA* FASTA keyed by ENST. Point the
  build at `combined_cdna.gtf` instead, and preflight that every GTF seqname
  exists as a FASTA record.
- RefSeq taxonomy drift: `NC_009334.1` is placed in **Orthoherpesviridae** by
  current NCBI, not Herpesviridae (post-split). The panel's older family labels
  are stale in the same way.

## Method note — an agent hallucinated accessions

A subagent was told not to guess accessions and did anyway: it reported
`NC_001856.1` as a polyomavirus (it is *Borrelia burgdorferi*), `NC_014968.1` as
HCoV-OC43 (a sweet-potato leaf-curl plant virus), and `NC_007455.1` as EBV (a
bocavirus). Every accession in the final panel is instead derived from
NCBI-sourced tables in this repo and re-verified live. The lesson generalises:
on a curation task, verify the returned `DEFINITION` string against the intended
species, because a plausible accession resolves to a real but wrong organism.

## Open items this build did not close

- `virus_catalog.tsv` is behind the panel by 96 rows and mislabels family for the
  herpesviruses.
- HPV rows need `source`/`refseq`/`oncogenic_class` propagated into the
  catalogue.
- `M12737` (HPV-8) has no CDS upstream and falls back to a whole-genome
  pseudo-transcript — the exact shape that collapsed 99.8 % of anellovirus UMI
  into one bucket in `CAT-18`. `NC_039089` (HPV-71) is similarly thin (RefSeq
  annotates only E1).
- HAV `NC_001489.1` and HTLV-2 `NC_001488.1` now emit genome-scale `misc_RNA`
  features (100 % and 94 % of their genomes) because of the EBER fix. They
  belong to the two-index tier split (`CAT-18`), not to this panel.
- `kallisto inspect` (max EC size, discarded EC count) has not yet been run on
  the new index.
