# Discovery versus mg2sc/Cellsnake: parameters and per-UMI aggregation

Date: 2026-10-21, rerun after plan step R5 (per-UMI rule changed to deepest-consistent) on top of DISC-09 (`ad1a534`).
Script: `scripts/compare.py`.
Outputs: `outputs/comparison.tsv` (sha256 `881a98e0…`), `outputs/umi_detail.tsv.gz` (sha256 `7d03ef77…`).
A second run from scratch reproduced both files byte for byte. (Earlier version at `0f219eb`: `1ff9fc46…` / `37872351…`, columns `viralscan_lca` and `viralscan_deepest_consistent_audit`; the `viralscan_deepest_consistent` numbers are identical to that audit column.)

**Scope.** This is a synthetic benchmark with known truth. It compares *aggregation rules and Kraken2 settings* under
otherwise identical inputs. It is not an estimate of real-data sensitivity or specificity (plan §12.2, gate G13).

## Setup

**Reference.** A Kraken2 2.17.2 database built from six 30 kb random genomes on a 16-node taxonomy:
- E. coli and S. flexneri share 40 % of their sequence (same family, different genera);
- the EBV genome carries a 2 kb segment identical to the human genome (an EVE-like homology);
- S. aureus and human are unrelated;
- a sixth genome ("novel") is **not** in the database.

**Sample.** 60 called cells (real 10x v3 on-list barcodes) with 20 UMIs each: 1,200 UMIs and 3,641 R2 reads of
90 bp, with 1 % substitutions and an `N` run in 5 % of reads.

| UMI kind | n | 1–5 reads from |
|---|---|---|
| ebv | 305 | anywhere in EBV |
| ebv_homolog | 71 | only the human-identical 2 kb |
| ecoli | 138 | anywhere in E. coli |
| ecoli_shared | 81 | only the region shared with Shigella |
| shigella | 82 | anywhere in Shigella |
| saureus | 137 | S. aureus |
| human | 170 | human |
| novel | 95 | the genome absent from the database |
| collision | 121 | two different species alternating (two transcripts sharing a CB/UMI) |

**Compared, everything else held equal.**

Kraken2 settings:
- ViralScan default `0.0 / 2` (Kraken2's own defaults);
- Cellsnake's documented `0.05 / 4` (`scrna-workflow` `config.yaml:50-51`, `cellsnake/command_line.py:142-143`);
- the Snakefile fallback `0.01 / 3` (`workflow/Snakefile:167-168`).

Aggregation per (CB, UMI):
1. ViralScan `resolved_taxid` (`viralscan discover` run end to end), now the **deepest-consistent** rule (user
   decision 2026-10-21): the deepest taxon when every classified read lies on one root-to-leaf path, else the LCA;
   root-level results stay ambiguous. Column `viralscan_deepest_consistent`.
2. The **previous LCA-only rule**, recomputed from the `taxids` column of `umi_groups.tsv.gz` (no longer what discovery
   outputs): the LCA of classified reads. Column `viralscan_lca_only_previous`. Kept as a baseline.
3. mg2sc's vote: `max(set(labels), key=labels.count)` over all reads, unclassified `"0"` included (`k2sc.py:119-121`,
   identical in Cellsnake's fork `e744c2c`). It is applied to discovery's own classifications and corrected
   barcodes, so only the aggregation differs. It is run under 10 `PYTHONHASHSEED` values.

Not emulated:
- mg2sc's BAM front end (Cell Ranger `CB`/`UB` tags, `samtools view -f 4`, the positional BAM/Kraken2 `zip`);
- the fork's optional fastp, bowtie2 and second-Kraken2 host filters;
- Cellsnake's cell-level filters (`min.features = 3`, `Homo*` name removal).

## Results (setting `0.0 / 2`; the other two settings are within ±8 UMIs on every row)

| Outcome (of 1,200 UMIs) | ViralScan deepest-consistent (current) | previous LCA-only | mg2sc vote |
|---|---:|---:|---:|
| Single-source UMIs (1,004): exact species | **799** | 606 | 752 |
| … correct but coarser (an ancestor of the truth) | 108 | 231 | 144 |
| … no taxon (unclassified / ambiguous) | 6 | 76 | 1 |
| … **root-level taxon reported as a feature** | 0 | 0 | 16 |
| … wrong lineage | 0 | 0 | 0 |
| EBV-homolog-only UMIs (71): no taxon | 71 | 71 | 0 |
| … **root-level taxon reported as a feature** | 0 | 0 | 71 |
| Collision UMIs (121): **one species claimed for both transcripts** | 8 | 2 | **118** |
| … resolved to a common ancestor | 25 | 25 | 3 |
| … no taxon | 88 | 94 | 0 |
| Novel-genome UMIs (95): any taxon assigned | 0 | 0 | 0 |
| UMIs whose label changes with the hash seed | 0 | 0 | **132** (11 %) |

Full grid in `outputs/comparison.tsv`; one row per UMI in `outputs/umi_detail.tsv.gz`.

## Findings

1. **Kraken2 settings barely matter here.**
   - Across `0.0/2`, `0.05/4` and `0.01/3`, exact calls moved by at most 3 UMIs and "no taxon" by at most 8, for
     every aggregation.
   - Random genomes with 1 % error are not where these thresholds bite. This benchmark gives no reason to adopt
     Cellsnake's `0.05/4`, nor to reject it.
   - The real-data grid (plan §12.2, M5) is still needed.
   - On real data they do matter: on EBV SRR12682296, 0.05/4 kept 246,233 EBV groups in called cells against
     623,096 at 0.0/2 (PLAN `DISC-R10`, `docs/discovery.md`).
2. **mg2sc's vote is non-deterministic.**
   - 132–139 of 1,200 UMIs (11 %) changed label across 10 `PYTHONHASHSEED` values.
   - Ties follow `set` iteration order, which Python randomises per process for strings.
   - ViralScan output was byte-identical across reruns.
3. **mg2sc claims one species for colliding transcripts.**
   - It did so for 118 of 121 two-species UMIs.
   - ViralScan's deepest-consistent rule did so for 8 to 12 (the previous LCA-only rule: 2 to 5). These happen when
     one transcript's reads were all unclassified or ancestral, so nothing contradicts the other. The new rule trades
     about 6 more collision overcalls for 193 more exact calls.
4. **mg2sc reports root as a taxon.**
   - 87 UMIs became a "root" feature: all 71 human-homolog EBV UMIs, plus 16 others.
   - Cellsnake's `Homo*` name filter would not remove them. ViralScan marks them ambiguous and keeps them out of
     the matrix.
5. **The deepest-consistent rule recovers the exact calls the LCA-only rule gave up.**
   - LCA-only made 606 exact calls against mg2sc's 752. Most of the gap was lineage-*consistent* UMIs: one read is
     root- or family-level, from the human-identical EBV segment or the E. coli/Shigella shared region, and another is
     species-level. The LCA climbed to that ancestor, although no read contradicts the species.
   - Deepest-consistent gives 799 exact calls (more than mg2sc), 6 with no taxon, 0 root features and 8 collision
     overcalls (mg2sc: 118).
   - Truly conflicting reads still go to their LCA.
   - Limit: a species unit can rest on a single species-specific read when the other reads are ancestral.
6. **The novel-genome test is uninformative.** A random genome shares no k-mers with the database, so nothing
   assigned it. A related-but-absent organism (for example a species withheld from its genus) is needed to test
   false discovery.

## Decision

Decision 4 was changed by the user on 2026-10-21: the per-UMI taxon is the deepest classified taxon when all
classified reads lie on one root-to-leaf path, else their LCA; root-level results stay ambiguous. Gate G7 is now "no
group is counted below the LCA of its *conflicting* classified reads". Implemented in `_group_shard` (plan step R5);
the numbers above are from that code. This benchmark is synthetic; it does not estimate real-data accuracy.
