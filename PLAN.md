# ViralScan 3.0 execution tracker

Status: **active**

Branch: `codex/viralscan-v3`

Last reconciled: 2026-09-27

Release target: `3.0.0rc1`, then `3.0.0`

Scientific target: methods-grade manuscript generated only from frozen v3 artifacts

This is the single operational tracker. The accepted scientific and product
contracts are in
[`docs/plans/2026-07-22-viralscan-v3-correctness-and-publication.md`](docs/plans/2026-07-22-viralscan-v3-correctness-and-publication.md).
The former pre-v3 tracker is retained only in the local governance archive and
excluded from Git and source distributions because it contains institutional
paths and private-analysis notes. It must not be used to establish v3
completion.

## Next action

**2026-10-04 (later): body census done, and the read-side measures leak.**
Over all 19,785 covid anellovirus-aligned reads, no read body longer than 26 nt
places on any anellovirus genome. The genuine component is ≤ 13 reads
(CP95 0.11 %), all at chance level. But `has_tso` misses edge-truncated TSO,
and `read_body` scores a 5′ TSO/TruSeq fragment as a complex body. Next: the
reagent-leak fix under `ANELLO-PRIOR.3`, which needs one user decision on
column semantics. Then `REL-16` and `ANDET-09f`.

**2026-10-04 (latest): `ANDET-09` closed; the covid TTV signal is an artefact.**
The STARsolo anellovirus alignment branch is implemented and tested, and ships
**off** (`ANDET-09e` criterion 1 failed 1/8 genomes: `SRR2037085_NODE_7436`,
kallisto 14, alignment 0; criteria 2 and 3 passed). The commissioned adversarial
review returned and *upheld the artefact conclusion while correcting its stated
reasoning* — 0/30 read bodies align to any anellovirus genome, and the
pileup/poly-A/`NM = 0` triad is non-diagnostic
(`.living/findings/covid-anellovirus-signal-is-polyg-reads.md`, update
2026-10-04). Next: `ANDET-09f` (filter sweep on the existing plant set),
`REL-16` (the kallisto segfault), and the open user decision on scope for the
covid and SFL-tonsil re-quantification.

**2026-10-03 (night): `ANDET-09` started.** A STARsolo anellovirus
alignment branch that does not depend on kallisto, from a grill of an external
bulk/rustar plan.

**2026-10-03 (latest, evening): DEF-02 chemistry module landed.** `-x` now
defaults to read-based detection and a contradicted `-x` stops the run
(`--force-technology` overrides); SW-21 closed (Drop-seq gets `-w None`); the
3′ strand pilot picks forward. Waiting on HSV-1 job 25695959 for the
end-to-end check. Per-chemistry default *values* stay with DEF-00.

**2026-10-03 (latest): CAT-42 closed; cat42b is the current panel.**
Every swap gate passed (see `CAT-42`). The covid 5′ strand reruns failed in
emptyDrops, because the `test_viralscan` env's `Rscript` has no Matrix. They
were resubmitted with `--resume` as 25695488. Still open in this pass:
GOV-06 closed 2026-10-03. The 5′ reruns finished and `--strand auto` landed on
2026-10-03 (see DEF-02 and F-020). TTMDV reverse-strand read check: planned, see F-020.
PROG-07 was re-measured but reopened (`PROG-17`: layer 2 ignored the
called-cell set; both closed later the same day). The REF-13 manifest landed 2026-10-03.
PROG-08 closed 2026-10-03: KSHV is `complete` (LANA cluster, K12, LANA2);
real-data validation is `PROG-18`. PROG-17 and PROG-07 closed 2026-10-03:
layer 2 scores called cells only; EBV LCL = 932 cells, 526 latent, 0 inversions.

**2026-10-02 (latest): finish-pass over every started (`[~]`) row.**
The reviewed plan is at `~/.claude/plans/read-last-handoff-document-moonlit-noodle.md`.
- W0 (housekeeping and blocker notes) is done.
- W1 runs three implementers: `v3/sw24-mask` (SW-24 plus the CAT-42 composition mask), `v3/cli` (SW-06, DOC-02, later `--strand auto`) and `v3/gov06` (LVC-11/12, protocol 1.2.0).
- User decisions: CAT-42 replaces the shipped panel if every check passes; GOV-06 v3 gets an explicit `-gtf` (protocol 1.2.0); no PyPI installs; no push.

**2026-10-01 (latest): CAT-42 panel rebuild running (job 25691642,
homopolymer mask ≥20).** Checks 2 and 3 under `CAT-42` are queued behind it
(25691838 GSE189670, 25691839 EBV/HSV-1). When the build ends, read
`kallisto inspect` (check 1), then the two runs, then return to the
5′ `--strand` reruns and `DEF-00`. `MECH-A2` landed while it ran. Handoff:
[`.living/HANDOFF_2026-10-01.md`](.living/HANDOFF_2026-10-01.md).

**2026-09-29 (latest): every run without `-w` skipped barcode correction
(`SW-13`, fixed).** Snakemake turned an empty whitelist into `kb count -w None`,
which means "bypass". Rerun any result that needs absolute numbers.

The mechanism review (`WP1D`) ranks the next work:

1. ~~`SW-14`…`SW-20` defect fixes~~, all landed 2026-09-29.
0. **Grill decisions confirmed 2026-09-29** — see `WP1E`. They bind all work
   below. On 2026-09-30, MECH-A step 4, `DEF-03` and `--strand` landed.
   Steps 4a–4c (the guard false positive, the NC_000898.1 catalogue row,
   old-run backfill) also landed, which closes MECH-A. Next: the
   5′ reruns with `--strand`,
   then `DEF-00`, the protocol amendment.
2. `MECH-A`, the Virus Identity table, is done (2026-09-30). Step 1 (the catalogue
   merge, 4,128 accessions with taxid and `panel` scope) and step 2
   (`virus_identity.py`, golden-tested on 3 real indexes) landed on 2026-09-29.
   Step 3 (the `analysis` rule writes `results/virus_identity.tsv`) also
   landed. Step 4, re-pointing the consumers, landed 2026-09-30, and so did
   steps 4a–4c. Both guard rules are decided (grill Q8 = B, adding an
   index manifest; Q8b = A2, genus plus "Anelloviridae (genus unassigned)").

The 4,127-genome max panel (`viral_panel_max_2026-09-28/`) built and passed
every gate, but must not become the default before `MECH-A`/`MECH-B`. EBV
type 2 took 16 % of EBV molecules through the equal split of shared k-mers
(F-017).

**2026-09-28: the final reference is BUILT — 2,343 genomes.** `CAT-31`
through `CAT-35` shipped together; the panel is 323 GTF-backed genomes plus
2,020 anelloviruses, and the CAT-31 guard closed 32 of 34 catalogued gaps inside
that build. Full measurements in
[`.living/findings/final-reference-panel-2343-genomes.md`](.living/findings/final-reference-panel-2343-genomes.md)
(F-016). Four things now matter, none of them "add more viruses":

1. **`CAT-40` — gate the index on `kallisto inspect`** (max EC size, discarded
   ECs) before trusting it. The HPV k-mer metrics are computed outside kallisto
   and cannot see the pseudo-inverse threshold effect, which is `CAT-19`'s real
   failure mode. This is the one open item that could still invalidate the build.
2. **`CAT-36` — the catalogue is behind the panel by 96 rows**, so 96 indexed
   genomes have no provenance, family or tier. `CAT-31`'s reverse check will keep
   reporting them.
3. **`CAT-38`/`CAT-39` — three whole-genome pseudo-transcripts** entered the
   panel (`NC_001489.1` HAV, `NC_001488.1` HTLV-2, `M12737` HPV-8). This is the
   exact shape that collapsed 99.8 % of anellovirus UMI into one bucket, so they
   need a deliberate decision, not a silent default.
4. **`CAT-30` is still the only row blocked on you** — retract the published
   Alphatorquevirus number, or leave the published results frozen? The science
   is closed; the retraction never happened.

**2026-09-27 (later): F-015 — the catalogue is not the index, and it outranks
every reference-import question above.** Measured against the built
`viral_genome.dedup.fa`: **13 of 16 catalogued HPV are absent** (including
HPV18 and HPV31), and **all 5 Retroviridae — both HIV-1 and HIV-2 — are
absent**, plus 7/30 influenza, 4/7 coronaviruses, 3/8 polyomaviruses. A
catalogued virus that is not indexed is *undetectable*, so this silently caps
sensitivity today and would do the same to any import.

1. **`CAT-31` shipped 2026-09-28 (guard only — the 34 misses are still open).**
   The build now reconciles the catalogue against the panel it actually emits
   and writes `catalogued_not_indexed.tsv`, so the gap is reported on every
   build instead of being invisible. The guard is deliberately **red**: it
   reports all 34 catalogued-but-unindexed accessions as `unexplained` and fails
   under `--strict-reconciliation`, because none of them is a documented
   exclusion — they are F-015's measured sensitivity loss. Each one now needs a
   genome or a named decision, not a silent default. Next: work the report.
2. **`CAT-32` is the cheapest real EBV win.** EBV type 2 (`NC_009334.1`) is
   missing; 31-mer overlap with type 1 is 79.9% shared, leaving ~20%
   type-discriminating k-mers. One accession turns EBV detection into EBV
   *typing*. The `CAT-31` report will confirm it lands.
3. **`CAT-33` reframes the HPV question.** VIRTUS2 offers 92 HPV *types* vs our
   16 → 76 genuinely new, but as weak-provenance `lcl|` records; source the same
   types from curated RefSeq `NC_` accessions. Start with the **13 already
   catalogued** types in the `CAT-31` report — that is 13 real genomes, not 76
   speculative ones.
4. **`CAT-34` closes a worry.** No `chrEBV`/viral contig exists in our actual
   host build (GRCh38-2024-A, 194 contigs; GRCh38.116 GTF), so the VirDetect
   silent-invisibility hazard does not apply here.
5. **VIRTUS2's other 400 records are not worth importing** (`CAT-35`, done):
   legacy influenza lab strains and 27 anelloviruses we already cover better.
   Take the *method* (host subtraction, strand-aware counting), not the genomes.

**2026-09-27: `CAT-17` shipped and the panel expansion was adversarially
reviewed — five new rows (`CAT-25`…`CAT-30`), four retractions, 41 external
files vendored.** Three things now need a decision rather than more work:

1. **`CAT-30` is the most urgent and is blocked on you.** `main` @ `4fcd748`
   already closed F-005 as an accession-level artifact, but
   `covid_viralscan/results/*/results/viral_summary.tsv` **still publishes**
   `Alphatorquevirus` 1,167,103 and 1,605,631 UMI at 17.6–19.0 % of cells. The
   published number and the closed finding contradict each other today.
2. **`CAT-05` has a live blocker.** The new gate rejected *both* the upstream
   2,023-rep panel and the shipped panel on one duplicate —
   `NC_038359.1` / `AB303562.1`, byte-identical. It is an upstream
   dereplication failure we inherited, and a duplicate sequence has broken
   `kallisto index` before.
3. **`CAT-20` stays gated** on the F-013/F-011 reconciliation, per the
   previous session's own instruction. Adopting the 2,023 reps is the same
   unmeasured question as restoring the 1,522, approached from the other side.

An advisory note on sequencing: the panel question was framed as "cover more
viruses", but for anelloviruses — 88 % of the panel — more genomes is provably
*not* the answer, because the ceiling is set by upstream dereplication. The
half of the request with a clear answer is breadth outside Anelloviridae
(`CAT-21` SARS-CoV-2 lineages, `CAT-22` influenza clades, `CAT-27` EBER,
`CAT-29` HSV-1 latency), and that half is cheap.

**Earlier, superseded:** WP4H gained a design document and 8 rows
(`CAT-17`…`CAT-24`) — read
[`docs/plans/2026-09-27-viral-reference-panel-expansion.md`](docs/plans/2026-09-27-viral-reference-panel-expansion.md)
before touching `CAT-09`…`CAT-16`.** It corrects three rows as
under-specified and adds six failure modes that were not tracked:

1. **`CAT-17` is the next thing to build.** A naive whole-genome anellovirus
   panel manufactures **1.44 % of R2 reads as false hits** in the EBV LCL
   `SRR12682296` — 6,437/6,437 captured hit reads had *zero* genuine
   anellovirus k-mer once low-complexity k-mers are masked. It is cheap, has
   no dependencies, and every panel expansion silently reintroduces it.
2. **`CAT-12` is not executable as written.** RefSeq holds **exactly one**
   SARS-CoV-2 genome and **zero** lineage-labelled RefSeq records
   database-wide, so "one per WHO variant lineage" cannot be built from
   RefSeq. → `CAT-21`, blocked on a user decision (document §8 Q1).
3. **"One influenza strain per subtype" is the wrong axis** — the diversity is
   in the clade, and the 8 current IAV RefSeq records are all segments of *one
   1934 lab strain*. → `CAT-22`.
4. **Panel expansion is not monotone-good** (`CAT-19`): `kb count` *discards*
   multimapping UMIs by default in the scRNA-seq path, and the one benchmark
   that measured this found F1 falls ρ = −0.73 with reference-set size for
   >99 %-identical genomes.
5. **The anellovirus panel is structurally incomplete** (`CAT-20`): genus-name
   queries reach 21,283 of 42,755 family records; SENV sits outside every
   genus; ICTV is now 37 genera / 243 species, not the 8 the shipped table
   encodes.
6. **The Serratus anchor becomes checkable** (`CAT-23`): the 18 latent-DNA-virus
   accessions of Lareau et al. *Nature* 2023 are already in the panel, so they
   are an acceptance criterion rather than new breadth.

`CAT-18` (two-index architecture) is the structural fix for the 99.8 %
single-bucket failure that `CAT-01` did not address: whole-genome and real-CDS
transcripts for the same virus still share one equivalence-class space.

**New, 2026-09-27: WP4H gains a breadth-and-diversity expansion
(`CAT-09`…`CAT-16`).** A census of the current reference found 2,313 records but
only ~106 non-anellovirus species, nearly all single-genome, with SARS-CoV-2,
HIV-1/2, HTLV-1/2, 3 of 4 seasonal coronaviruses, hMPV, bocavirus, influenza D
and every high-risk HPV except 16 **absent entirely**, and influenza A present
only as one 1934 lab strain. `CAT-01` (build-ref discards the real GTF) and
`CAT-03` (nothing groups segments) block all of it and come first. Diversity
targets are set by measured 31-mer capture, not by quota.

**New, 2026-09-27: `WP4J` — the SFL tonsil pool shows no TTV and no HPV in any
of the 24 donors (`TONSIL-01` done).** The screen used the host-subtracted
cellranger BAM: 118.9 M unmapped GEX reads, and 23 reads survived to the end,
none of them anellovirus or HPV. A positive-control plant recovered 99 %. A
held-out plant recovered only 8–40 % of reads from TTV strains ≤ 84 % identical
to any reference, so divergent low-level TTV is not ruled out. Only
`x223` is gene expression (`x225` is ADT), and the library is 5′ v3 R2-only.
None of these zeros is an informative negative (`SENS-06`). A native `viralscan`
run waits for `TONSIL-02` (strand and whitelist support). The WP4I → WP4G → WP4H
order below is unchanged.

**New, 2026-09-27: three work packages opened, in this order — `WP4I`
(latent/lytic), `WP4G` (anellovirus detection), `WP4H` (human-virus catalogue).**
`PROG-10` is done; next is `PROG-11` (catalogue
biology), because the `complete` labels and the `PROG-07` numbers are unreliable
until they land. `CAT-01` (build-ref discards the real GTF) blocks every natively
built index, so it precedes any `CAT-08` build and `ANELLO-13`. Do not flip the
anellovirus default (`ANDET-07`) before `ANDET-01`–`ANDET-04`. Housekeeping that
preceded this: the WP4B2–WP4F work was squashed into one commit so every commit
passes its own suite, and the schema 1.1.0 governance migration was finished
(see the evidence log).

**New, 2026-09-26: WP4F (`ANELLO-01`–`ANELLO-13`) landed** — the Anelloviridae
panel now carries real NCBI gene structure instead of one placeholder gene per
genome. The root cause was a *discarded* GenBank GTF, not a missing one:
`_fetch_one` already wrote a real annotation for 1,995 of 2,042 accessions and
both panel builders threw it away and rebuilt a placeholder from the FASTA. A
second, independent defect would have survived that fix — gene IDs were the bare
`/gene=` value, so `ORF1` (150 genomes) would have collapsed 1,995 genomes onto
2,316 columns. `extras/build_anellovirus_genes.py` generates
`src/viralscan/data/anellovirus_genes.tsv`; measured CDS coverage is 205/206
(99.5 %) across all eight genera, and the low CDS *count* is biology, not a gap:
75 % of the panel is Betatorquevirus TT-mini genomes that genuinely carry one
ORF. NCBI carries **no** genogroup for Anelloviridae, so that column ships empty
rather than inferred. **WP4F also fixes the two things `HPV-09` was blocked on**
(genome-scoped `_genbank_to_gtf` IDs, public `fetch_genbank()`); that row is left
for its owner to close against their own tests. `ANELLO-12`/`ANELLO-13` record
the limit: the 99.8 %-in-one-bucket covid artifact is now explainable and
testable, but the test needs `kb` and a re-run, and per-genotype anellovirus
quantification is still not defensible for the one-ORF majority of the panel.

**New, 2026-09-26: WP4E (`HPV-01`–`HPV-10`) landed** — HPV ORFs are now named.
All 14 high-risk genotypes plus HPV1/HPV2, 124 ORFs in
`src/viralscan/data/hpv_genes.tsv`, names taken from each record's own `/gene` or
`/product` qualifier rather than from a coordinate table, because papillomavirus
genomes are linearised circles cut at the submitter's choice of point and a
coordinate table is wrong for 15 of the 16 genotypes here. Two follow-ups are
recorded rather than done: `HPV-09` needs a `ncbi_fetch.py` fix (out of scope
for that change) and `HPV-11` needs the real index rebuilt to *measure* the L1
cross-mapping this row currently only predicts. **Do not publish a per-type HPV
number before `HPV-11` closes.**

**New, 2026-09-26: WP4D (`PROG-01`–`PROG-07`) landed** — layer 2 gene-programme
inference, opt-in via `--gene-programs`, with the EBV LCL regression test
pinning the design. `PROG-08`/`PROG-09` are the open scope questions. WP4E does
not change `PROG-09`: HPV still has no latency/lytic dichotomy, so it stays out
of the gene-programme catalogue. Also in this change: WP4C (`SENS-01`–`SENS-05`)
and WP4B2 (`HOST-01`–`HOST-04`) landed alongside the `REF-01`/`REF-13`
reference-visibility fixes.** An audit of detection sensitivity found that
`limit_of_detection` existed only as an unrun endpoint in the protocol schema,
that the bundled 20-genome TTV panel captured
**1.36 %** of the real anellovirus 31-mer space (median per-genome coverage
**0.00 %**, 85.8 % of genomes sharing zero 31-mers), and that 49.2 % of gene IDs
in the covid panel resolved to no virus name — so 9 of 17 published
`viral_summary.tsv` rows were bare gene IDs rather than viruses.

The immediate consequence: **`informative_negative` is `false` on essentially
every run, by design**, because certifying a negative needs a measured k-mer
capture term and no shipped workflow plants a control. `SENS-06` is the row that
fixes that, and it should be built together with `VAL-01` — a generator that
plants a target at known abundance *is* the positive control `SENS-04` consumes.
Until then, no negative result from this package may be reported as an absence.

**Governance follow-up required at commit time (not done here, deliberately).**
`analysis/v3_artifact_inventory.tsv` has no rows for the artifacts this work
introduces, and adding them before the commit would break
`test_artifact_inventory`: `governance_utils` resolves each row's `path` from
the **git tree at the recorded `git_sha`**, so a row pointing at an untracked
file can never validate. At commit, add rows for
`src/viralscan/sensitivity.py`, `src/viralscan/gene_programs.py`,
`src/viralscan/data/gene_programs.tsv`, `src/viralscan/scripts/gene_programs.py`,
`results/sensitivity.tsv`, `results/positive_control.json`,
`results/gene_program_summary.tsv` and `results/gene_program_cells.tsv` (and
re-derive the five already-failing rows:
`artifact-inventory-schema`, `claim-registry`, `claim-registry-schema`,
`output-reference-doc`, `v3-counting-contract-doc`, which fail today only
because those files are modified-but-uncommitted in the working tree).

Pre-existing and still first in the queue: **`G0` archive-build evidence, then
`G1`.** `GOV-03`,
`GOV-04`, and `GOV-05` now have fail-closed local evidence: the sanitized
artifact inventory, claim graph and public ship-scope validators pass, and a
wheel plus sdist built directly through the installed setuptools backend match
the positive allowlist. `G0` remains partial because the locked Python lacks the
PyPA `build` frontend required by the frozen validation command. No dependency
was downloaded or installed without approval.

### Do now

1. Supply a compatible PyPA `build` frontend to
   `benchmark_runs/legacy_v2_v3/env_full` through an explicitly approved network
   install or a user-provided offline artifact.
2. Run the prescribed `python -m build --no-isolation` command and the exact
   wheel/sdist/context membership check; rerun the full local gate and frozen
   dirty-path hashes before promoting `G0`.
3. Treat the separate integration-marker failures as `G1` evidence: the real CLI
   currently requires a `conda` executable for `--use-conda`, and this runtime
   does not provide one; the installed Snakemake also emits no rule listing with
   the tests' `--quiet` dry run.

## How to use this tracker

1. Take the first unchecked item whose dependencies are all `[x]`.
2. Add or update its test before changing production behavior.
3. Run the stated acceptance command and save the named artifact.
4. Mark `[x]` only with dated evidence; use `[~]` for partial and `[!]` for an
   external blocker.
5. Update **Next action**, the work-package row, and the claim registry in the
   same commit.

Status legend: `[x]` verified; `[~]` partial; `[ ]` ready/not started; `[!]`
blocked on a person, credential, private datum, external service, or unavailable
compute.

## Critical path

| Order | Work package | Status | Depends on | Exit gate |
|---:|---|:---:|---|---|
| 0 | Governance and legacy freeze | `[~]` | none | `G0` |
| 1 | Core software contracts | `[~]` | `G0` | `G1` |
| 2 | Install, lock, and artifact parity | `[~]` | `G1` | `G2` |
| 3 | Preregister scientific validation | `[~]` | `G0` | `G3` |
| 4 | Freeze production references | `[~]` | `G3` | `G4` |
| 5 | Build and validate the truth panel | `[ ]` | `G3`, `G4` | `G5a` |
| 6 | Run public positives and comparators | `[ ]` | `G2`, `G4`, `G5a` | `G5b` |
| 7 | Score, audit, and freeze results | `[ ]` | `G5b` | `G5` |
| 8 | Regenerate docs and claims | `[~]` | `G1`; numbers require `G5` | `G6a` |
| 9 | Regenerate and review manuscript | `[ ]` | `G5`, `G6a` | `G6` |
| 10 | Publish and test `3.0.0rc1` | `[ ]` | `G2`, `G5`, `G6a` | `G7` |
| 11 | Publish stable `3.0.0` and submit | `[ ]` | `G6`, `G7` | `G8` |
| 12 | Observe and maintain | `[ ]` | `G8` | ongoing |

Work packages 1-4 may proceed in parallel. Do not start truth-panel outcome
analysis before `G3`, regenerate quantitative documentation before `G5`, or tag
a release before its gate.

## Verified baseline

- [x] `BASE-01` — v3 branch and locked correctness roadmap exist.
- [x] `BASE-02` — molecule-safe streaming allocation, distinct-gene EC
  projection, five explicit methods, contracted H5AD layers, and exact mass
  audits pass synthetic/property tests.
- [x] `BASE-03` — retained EBV input completed with 103,145,071 BUS records,
  57,957,364 input molecules, and exact partition/matrix conservation. Evidence:
  [`analysis/v3_ebv_baseline/README.md`](analysis/v3_ebv_baseline/README.md).
- [x] `BASE-04` — exact-target evidence completes extraction, competitive
  alignment/BLAST, deduplication, QC, coverage, BAM indexing, and IGV on the tiny
  real-tool fixture.
- [x] `BASE-05` — the 2026-07-22 local gate passed 720 tests with 21 deselected;
  Ruff check/format, data governance, and protocol validation also passed. Rerun
  the full gate after each remaining software slice rather than treating this
  development-branch result as release evidence.

## WP0 — Governance and legacy freeze

Objective: prevent unsafe pre-v3 results, private data, or unsupported claims
from entering the release. Estimated remaining effort: 1-2 days.

- [x] `GOV-01` — mark pre-v3 corrected/combined counts scientifically
  incompatible and ineligible for v3 claims.
- [x] `GOV-02` — add a ship-scope data-governance check and exclude private
  clinical outputs from package/container/source-distribution scope.
- [x] `GOV-03` — inventory every claim-bearing input, reference, intermediate,
  result, and scheduler record with software version and SHA-256. Write
  `analysis/v3_artifact_inventory.tsv`.
- [x] `GOV-04` — define the public ship-scope allowlist, then remove or quarantine
  institutional defaults from every included script and manifest.
- [x] `GOV-05` — expand `claims/registry.json` into a validated claim graph with
  source location, artifact digest, Git SHA, input/reference hashes, schema,
  layer, denominator, generation command, scope, and status.
- [x] `GOV-06` — **closed 2026-10-03 (diagnostic-only, outcome-ineligible).**
  LVC-11 to LVC-14 are `[x]` in `analysis/legacy_v2_v3/TRACKER.md`. The
  interpretation is in `FRESH_CONTROLS.md` (merged from `v3/gov06-close`,
  17f59f3).
  - **v2:** all five fresh rows reproduce the archive on every record within
    1e-6, including the attempt-3 highmem row SRR6825024.
  - **v3:** all five fresh rows completed with `validate-run` 0. They used
    the frozen wheel and an explicit `-gtf` (protocol 1.2.0, packet
    `66918c50…`). Exact-read evidence succeeded for every expected target.
  - **Fresh v3 vs archived v3:** 836 of 862 records match. The 26
    mismatches:
    - 20 HIV reference differences (the arms differ on HIV);
    - 3 Cercopithecine herpesvirus 1→0 molecules (unresolved);
    - 3 EBV: 2 float noise, and SRR12682298 at −2 molecules (unresolved).
  - **HIV:** v2's HIV result is structurally negative, not a sensitivity
    comparison.
  - **EBV:** v3 is 0.46–0.59× v2. These are different quantities and are
    not rescaled. EBV recovery is qualitative positive-control evidence
    only.
  - **Caveats:** the frozen v3 wheel predates SW-13 and carries SW-24. The
    GTF/t2g parity residual was accepted fail-closed. No truth, comparator,
    calibration, release or publication gate closes on this.
  - **Retained failure:** v3 evidence attempt 1 failed because `kallisto`
    and `bustools` were not on PATH; it is kept.
  - **Open:** for the HIV controls, the evidence helper's "largest
    non-target" is another HIV gene, so no non-HIV candidate was traced.
  - Previous text follows as history. Execute the outcome-ineligible
    ViralScan 2.2.0 versus v3
  diagnostic in `analysis/legacy_v2_v3/`. The identical-BUS arm is complete
  (44/44 valid rows) and the five-control input gate is closed. Fresh
  matched-FASTQ attempt-2 arrays `25331035` and `25331037` are **terminal with
  all ten tasks failed** across three independent causes: a nested-output
  validator defect that manufactured exit 65 on four genuinely successful v2
  rows, an unpopulated Zenodo viral-data cache that killed all five v3 rows
  before quantification, and one genuine v2 out-of-memory on `SRR6825024`.
  Attempt 3 is fully wired and tested but **not frozen**, because the cache pin
  it requires is blocked on `REF-11`; `LVC-13`–`LVC-14` remain unstarted and have
  no comparison tooling yet. A same-input head-to-head with evonk's released
  2.2.0 is tracked separately as `CMP-06` (WP6B).
  - 2026-10-02 (`v3/gov06`, e3a08b9).
    - **Tooling:**
      - highmem tier wired in, via `--highmem-task`;
      - explicit `-gtf` path, so the cache is optional;
      - GTF/t2g parity check;
      - LVC-12 exact-read evidence step;
      - `compare_legacy_v2_v3.py fresh-vs-archive`. The four attempt-2 v2 rows match the archive on all 300 records, within 1e-6.
    - **LVC-11:** the attempt-3 v2 highmem packet is submitted as job 25694919.
    - **LVC-12:** blocked on a user decision.
      - The parity check of the v2-arm GTF (195 files, 2,692 genes) against the frozen index t2g fails:
        - 89 GTF-only genes;
        - 200 non-Ensembl t2g genes not in the GTF, including all of HIV-1.
      - The v2 panel has no HIV GTF, so the archive could never call HIV.
      - The protocol 1.2.0 text is drafted, not applied, until that decision.
    - **LVC-12 unblocked (2026-10-02, `v3/gov06-hiv`, 325cc58).** The user decided to add an HIV GTF.
      - Source: evonk `Serratus_v2/gtf/Human_immunodeficiency_virus_NC_001802.gtf`. Its gene IDs and coordinates equal the index's IMMUNO_HIV1gp1–10.
      - Parity after the addition: 2,613 of 2,702 genes in the t2g. The residual (89 GTF-only, 190 t2g-only) is accepted fail-closed (`--accept-parity-residual`).
      - Protocol **1.2.0** applied, `outcome_triggered: true`; the arms now differ on HIV.
      - Packet `fresh_control_packet_attempt3b_v3/` (5 v3 rows) frozen and submitted.

`G0` passes when the governance scan is green, every public quantitative claim
is registered, all pre-v3 quantitative claims are rejected or historical, and
the inventory identifies enough retained BUS/reference/FASTQ material to rebuild
each eligible result.

Evidence to record:

```text
python3 scripts/check_data_governance.py
python3 scripts/validate_claim_registry.py --coverage
sha256sum analysis/v3_artifact_inventory.tsv claims/registry.json
```

## WP1 — Core software contracts

Objective: close remaining correctness and workflow-consistency gaps before
scientific-scale execution. Estimated remaining effort: 4-7 engineering days.

### WP1A — Schemas, assignments, and reruns

- [x] `SW-01` — move/package all v3 schemas inside the installed `viralscan`
  distribution, load them with `importlib.resources`, and fail closed when a
  required schema is missing. Add wheel and sdist tests.
- [x] `SW-02` — enforce every public JSON/TSV/H5AD v3 schema at write and
  `validate-run` boundaries; remove generic silent skips. Three of the six
  shipped schemas had no reader at any boundary: `count_audit`,
  `reference_manifest`, and `evidence_manifest` shipped without ever validating
  a document. `h5ad_contract.json` is not a JSON Schema — it declares no
  keywords, so handing it to a validator would have accepted everything while
  looking like enforcement; `validate_json_schema` now refuses it by code
  (`not_a_json_schema`) and `_matrix_issues` reads `required_layers` and
  `required_uns` from it instead of from literals, closing a two-sources-of-truth
  gap that had left `quantification_unit` and `multimap_method` unchecked.

  Write boundaries raise (`SchemaContractError`) rather than returning issues;
  `validate-run` reports. The split is deliberate and is about authorship: at
  `validate-run` the artifact is input and a violation is a finding, while at a
  write boundary ViralScan is the author and a violation is a defect in this
  code, so publishing the file anyway would ship it under a schema it does not
  meet. Do not re-litigate this into a uniform policy.

  Note the h5ad count invariants were already fail-closed at *construction*
  (`multimapping.py:698` raises); `_matrix_issues` re-derives them from bytes on
  disk, which is a different guarantee — it catches a truncated write, a
  hand-edited file, or an artifact from another version rather than a compute
  bug. Two audit fields (`resolved_molecules`, `ignored_read_multiplicity`) and
  the contract's third invariant (unique mass equals audited unique molecules)
  had no reconstruction-side check at all; they do now.

  Silent skips removed: the `if manifest:` guard that let an empty run manifest
  pass schema validation and fingerprint checks, and the `ImportError` branch in
  `validate_json_schema` that turned a broken install of a hard dependency into
  a soft finding.
- [ ] `SW-03` — add optional compressed molecule-assignment evidence containing
  CB, UMI, ECs, distinct genes, ambiguity class, method, weights, and exclusion
  reason without changing default matrix mass.
- [x] `SW-04` — make `rerun-multimap` regenerate every method-dependent artifact
  in a new result tree: matrix/layers, count audit, summaries, evidence tiers,
  UMAPs, and host-response inputs.

  Audited 2026-07-28, one artifact at a time, because "regenerate everything" is
  not checkable without knowing which artifacts are actually method-dependent:

  - **matrix/layers** — regenerated. `_swap_multimap_layer` rewrites `X`, both
    compositional layers, and `uns["multimap_method"]`.
  - **count audit** — *not* stale, and does not need regenerating. Every field in
    `molecule_audit` is method-invariant: the molecule counts describe resolution
    rather than allocation, and `allocated_ambiguous_mass` is invariant because
    `host_conservative` divides `count / sum(cons_eligible)` across eligible genes
    rather than dropping mass — total allocated mass equals the ambiguous molecule
    count under every non-EM method, which `MoleculeAudit.validate` requires
    anyway. An earlier reading of this row assumed the swap corrupted the audit;
    it does not, and the fixture used to "reproduce" it was a state the pipeline
    cannot produce.
  - **summaries** — was broken, and not only on rerun. `multimap.py` opened
    `summary.txt` with mode `"w"` and wrote three totals; `detection.py` runs
    later in the DAG and opened the same path the same way, so those totals were
    destroyed on **every** run and never published. Nothing consumed them — no
    test, doc, or notebook. detection is now the sole writer and recomputes them
    from the H5AD, so they are both published and correct after a layer swap.
  - **evidence tiers, UMAPs, host-response inputs** — audited 2026-07-28 by a
    multi-agent review. Three further defects, one shared root cause: the
    command's invalidation list was incomplete, so `shutil.copytree` left old
    artifacts in a tree labelled with the new method.
    - `run_manifest.json`: **the reported defect was not real, and the fix for it
      was a regression.** The review claimed the manifest lives beside each
      sample's `config.yaml`, so the rewrite at `output_dir/` (the tree root) was
      a no-op. Running the pipeline end to end for `SW-10` showed the opposite —
      a completed run of `viralscan -o out` produces `out/run_manifest.json`
      alongside `out/<sample>/config.yaml`, exactly one manifest, at the root.
      The original code was correct. The per-sample rewrite shipped in `e6315cb`
      moved it to a path that never exists, so the manifest stopped being updated
      at all. Reverted 2026-07-29, with the layout now pinned by
      `test_run_manifest_is_at_the_tree_root`. The finding was accepted without
      being run; five of the review's findings were spot-checked and this was not
      one of them.
    - `plots/` was never cleared. Which viruses clear `detection_threshold` is
      method-dependent, and `generate_html_report` globs the directory, so a
      demoted virus's figure was re-embedded into a report whose own table no
      longer listed it. Only detection-owned patterns are cleared; `umap.py`
      writes into the same directory and its output is left alone.
    - `hostresponse/` was never cleared and `log/hostresponse.done` was never
      dropped, so it re-ran only if snakemake happened to judge it stale by
      mtime, and a virus falling below `MIN_VIRUS_CELLS` kept its old CSVs.

  Closed by `SW-05`.
- [x] `SW-05` — add an integration test proving no stale artifact survives a
  method change and the source result remains untouched.
  `tests/integration/test_rerun_no_stale_artifacts.py`, eight cases over the
  scenario that actually orphans files: a virus that clears the detection
  threshold under the source method and falls below it under the new one, so it
  is skipped rather than rewritten. Asserts the demoted virus leaves no plot or
  CSV, the surviving virus is regenerated rather than merely kept, provenance
  names the new method, artifacts owned by other rules survive the cleanup, every
  method-dependent sentinel is dropped, and the source tree is byte-unchanged.

  Both fixes were mutation-tested before the row was flipped: making
  `clear_stale_virus_outputs` a no-op fails
  `test_the_demoted_virus_leaves_no_hostresponse_csv_behind`, and restoring the
  old root-level manifest path fails `test_provenance_names_the_new_method`.

  Scope limit, stated in the module docstring: Snakemake is not invoked, so this
  proves nothing stale survives *when the rules re-run*, not that the DAG
  re-executes. That needs `SW-10`.

### WP1B — Workflow safety and architecture

- [x] `SW-06` — `doctor`, `validate-run`, fingerprints, resume/overwrite safety,
  and atomic manifests exist; finish whole-workflow staging and an atomic
  completion marker.
  - Done 2026-10-02 (`v3/cli`, deaf5c9). The user accepted that the
    completion marker serves as the whole-workflow staging.
    - `run_complete.json` holds the fingerprint, version and samples, plus the
      sha256 of `results/{viral_summary,virus_identity,multimap_evidence}.tsv`
      and `adata_multimap.h5ad` when present. It is written atomically after
      the sample loop and after rerun-multimap.
    - `prepare_output_directory` clears it on resume, overwrite and the
      rerun-multimap copy. rerun-programs and hostresponse re-stamp it.
    - `completion_marker: true` is added to the manifest after the
      fingerprint hash, so old `--resume` still matches.
      `validate_run` errors only for manifests that declare the marker; for
      old runs a missing marker is only a warning.
    - Schema `run_complete.schema.json` is in `REQUIRED_V3_SCHEMAS`.
    - Tests: 7 unit tests, plus an e2e check that deleting the marker fails
      validate-run.
- [~] `SW-07` — exact-fragment STAR filtering and mate synchronization exist;
  emit a reason for every retained/removed fragment plus lost-truth and
  host-virus-ambiguous boundary counts.
  - 2026-10-02 (`v3/tests-rel`, cbd5274):
    - `fragment_lineage.tsv.gz` now has one row per input fragment:
      `retained/host_unmapped` or `removed/host_mapped`. Memory is
      O(retained), by streaming input IDs against the retained set.
    - `lost_truth_counts()` is a Python hook (no CLI flag) that gives
      D15/D16 counts. Tests include a real STAR integration test.
    - Documented in `output_reference.md`.
    - Still open: the host–virus-ambiguous boundary count (deferred until
      defined), and the real lost-truth numbers, which need VAL-01.
- [x] `SW-12` — range-check the EM parameters. `multimap_pseudocount` was
  guarded; `multimap_em_max_iter` and `multimap_em_tol` were not. A budget of
  zero makes `range(1, max_iter + 1)` empty, so `em_gene_abundances` and
  `em_cell_abundances` return their pre-loop seed weights while the H5AD records
  `multimap_method` as an EM method. Mass conservation still holds, so
  `MoleculeAudit.validate` cannot see it, and the only trace is a
  `converged: false` diagnostic that nothing reads. Both now raise in
  `RunConfig.from_snakemake_config`.
- [x] `SW-08` — remove unsafe kallisto CB-UMI-wide host filtering from the stable
  CLI because exact fragment identifiers are unavailable.
- [ ] `SW-09` — split the oversized CLI into thin parsers plus importable service
  functions; convert Snakemake scripts to minimal wrappers without changing
  outputs.
- [x] `SW-10` — run one tiny paired-end fixture through documented CLI commands:
  preflight, reference, quantification, molecule allocation, cell calling,
  summaries, evidence, BAM/BLAST/plots/IGV, and `validate-run`.
  Executed for real on 2026-07-29 against `tests/data/evidence_tiny` (two read
  pairs, two reference sequences) using the tool binaries in
  `benchmark_runs/legacy_v2_v3/env_full/bin` (kallisto 0.50.1, bustools 0.43.2,
  kb_python 0.28.2, snakemake 8.20.5). `create_config` → `kb_count` → `analysis`
  → `multimap` → `detection` completes with exit 0, publishes all thirteen
  expected artifacts, and `validate-run` returns `ok: true` with zero issues —
  which exercises every `SW-02` schema check against a real artifact rather than
  a fixture. Codified as `tests/integration/test_tiny_end_to_end.py`, 17 cases,
  skipped when the binaries are absent.

  Two things the run established that no unit test could:
  - **`REF-11` is a hard blocker on the default path, demonstrated rather than
    inferred.** Without `-gtf` the run dies in `analysis`: the bundled panel is
    fetched from the unregistered Zenodo DOI. The fixture now ships its own
    minimal `viral.gtf`. When `REF-11` resolves, add a variant that drops `-gtf`.
  - **`run_manifest.json` lives at the tree root**, beside the per-sample
    directories rather than inside one. See `SW-04` below.

  Not covered: the evidence/BAM/BLAST/IGV leg, which needs `blastn`,
  `makeblastdb`, and `minimap2` — none present in this environment. Those legs
  are exercised by `test_exact_lineage.py` and `test_evidence_chain.py`. The row
  stays `[~]` until they run in one sequence.
  - Done 2026-10-02 (`v3/tests-rel`, 8149423). `TestDocumentedSequence` runs
    the tiny fixture end to end:
    1. `doctor --profile full`.
    2. The documented `-ref -fasta -gtf` build (`kb ref`). `build-ref` needs
       a network download, so it isn't used. The fixture needs a distinct
       viral record name and `transcript` rows.
    3. A run with visuals (plots and HTML).
    4. `--cell-calling external`.
    5. `viralscan evidence --blast` on the real run (all 15 artifacts).
    6. `validate-run` before and after evidence.

    Correction to the text above: env_full now has blastn, makeblastdb,
    minimap2, samtools, Rscript, STAR and cd-hit-est. The no-`-gtf` variant
    still waits on REF-11.
- [x] `SW-11` — make production cell calling fail closed: caller exceptions,
  zero-match external lists, invalid barcode geometry, and canonical collisions
  must never silently turn every barcode into a cell; `none` remains explicit.
  `detection.py` caught every exception and continued with `called_mask=None`,
  which `compute_stats` expands to all-ones, so a caller failure silently
  replaced the called-cell denominator with every barcode and still labelled the
  result a called-cell rate. All four paths now raise `CellCallingError`;
  `--cell-calling none` remains the explicit way to report over all barcodes.
  Completed 2026-07-27 by the review follow-ups: `emptydrops_seed` and
  `emptydrops_niters` are declared configuration with CLI flags rather than
  function-signature defaults, so the protocol's frozen `seeds.cell_calling`
  actually reaches `set.seed()` in `emptydrops.R`; and `Rscript` is preflighted
  whenever the resolved caller is `emptydrops`, because failing closed at the end
  of a multi-hour run for a knowable reason is the wrong place to fail.

- [x] `SW-13` — resolve `snakemake --config` wire values in Python, not in shell
  templates (2026-09-29).
  - Snakemake parses an empty `whitelist=` as `None`, and `{config[whitelist]:q}`
    rendered it as the literal `None`. **Every run without `-w` therefore called
    `kb count … -w None`, which kb reads as "bypass barcode error correction".**
    Seen in the SRR12682296 `kb_info.json`: bus → sort → inspect → count, no
    `bustools correct`, 793,308 barcodes. All reference-strategy and max-panel
    runs since then are uncorrected and need rerunning for absolute numbers.
  - Booleans: Snakemake converts only `True`/`False`, so `gene_programs=false`
    arrived as the non-empty string `"false"` and the rule was always planned.
  - Now a rule `params` lambda resolves the whitelist, and the gene_programs gate
    compares `str(...).lower() == "true"`.
  - Regression tests in `tests/test_snakefile_dag.py`: a template check, plus
    `snakemake -n -p` renders `WL=` when unset and `WL=/fake/wl.txt` when set.
    Both fail on the old Snakefile.
  - Root-cause removal is `MECH-C` (single Run Config writer).
  - Follow-up: `tests/data/evidence_tiny/R1.fastq` carried off-list barcodes
    (`AAAA…`/`GGGG…`), and the tiny e2e run had passed only because correction
    was bypassed. With correction on, kb drops them and cell calling fails
    closed on 0 barcodes. The fixture now uses two v3 on-list barcodes.
- [x] `SW-20` — (2026-09-29) the other half of `SW-13`. ViralScan's own matrix still skipped
  barcode correction when no `-w` is given.
  - `multimap.prepare_resolved_bus` sorts the **raw** `output.bus` and runs
    `bustools correct` only for a user whitelist. The multimap-derived X, which
    feeds detection, is therefore still uncorrected, while kb's
    `counts_unfiltered` is now corrected.
  - Use kb's corrected and sorted `output.unfiltered.bus`, or its copied on-list
    (e.g. `10x_version3_whitelist.txt`), when `-w` is absent.
  - Count the dropped off-list records in the audit.
  - Test with a fixture barcode one mismatch from the on-list.
  - Fixed. `multimap.select_bus_input` picks kb's `output.unfiltered.bus`
    (`KbCountOutputs.kb_corrected_bus`) when no user on-list is given, and warns
    when neither kb nor the user corrected. Checked end to end on the tiny
    fixture plus one read with barcode `GAACCCAAGAAACACT`, which is one
    mismatch from on-list `AAACCCAAGAAACACT` and has a unique neighbour:
    - fixed code: the read is corrected, 3 input molecules, 2 viral molecules;
    - previous code: the raw barcode is dropped **before** the audit counts it,
      2 input molecules, 1 viral molecule, and the loss is invisible.
  - Still open, under `MECH-F`: count off-list drops in the audit rather than
    before it.
- [x] `SW-21` — **closed 2026-10-03** (user decision: bypass). See DEF-02 (c).
  Original text: (found 2026-09-29, F-018) since `SW-13`, a Drop-seq run
  without `-w` is pre-filtered to a data-derived cell list.
  - Cause: Drop-seq has no official on-list. `kb count` therefore runs
    `bustools allowlist`, a knee on the data, and corrects against it. Every
    other barcode is discarded, both in kb's matrix and in the corrected BUS
    that `SW-20` feeds to multimap.
  - Measured on HSV-1 SRR8315713 (job 25666447, v1 index):
    - kept barcodes: 1,809,992 → 5,468;
    - HSV-1 molecules: 32,404 → 23,170 (−28 %);
    - infected barcodes: 10,227 → 1,127.
  - The 5,468 barcodes become the "all barcodes" denominator, and ViralScan's
    own cell calling runs after this implicit one.
  - 10x runs behave as intended. EBV SRR12682296 (10xv2, official on-list)
    keeps 885,933 of 906,202 molecules (−2.2 %), and the `summary.txt` headline
    now matches, at 885,939.
  - Decision needed, under `MECH-D`, the chemistry module: a technology with no
    official on-list either keeps no correction (pass kb its bypass value) or
    keeps kb's allowlist and documents the pre-filter. This is not changed in
    the Snakefile until then.
- [ ] `SW-22` — (found 2026-09-29) `rerun-multimap` does not skip `kb_count`.
  - Cause: `_run_rerun_multimap` rewrites the copy's `config.yaml`, which is an
    input of `kb_count`. Under snakemake's default mtime trigger, `kb_count`
    (kallisto on the FASTQs), `analysis` and `multimap` all rerun, which also
    undoes the fast layer swap.
  - Shown by dry-run on a skeleton copy of `combined_corrected/SRR8315713`.
    Without the config rewrite, only detection, umap and all are scheduled.
  - This predates MECH-A. A fix (for example keeping the config's mtime, or
    moving the multimap parameters out of `kb_count`'s inputs) belongs with
    `MECH-C`, the single Run Config writer.
- [x] `SW-14` — (2026-09-29) the Snakemake invocation in `menu.main` and `_run_rerun_multimap`, now one helper, `menu._snakemake_run_command`, covered by `tests/test_cli.py::TestSnakemakeRunCommand`:
  - put the `all` target before `--quiet`, because snakemake 9 lets
    `--quiet [...]` consume the target;
  - drop the unconditional `--use-conda`, since conda is not required at run
    time.
- [x] `SW-15` — (2026-09-29) STAR host-filter geometry. Done in `host_filter.starsolo_barcode_args` (pure) and `_plain_whitelist`, with tests in `TestStarsoloBarcodeArgs`. Merging with the benchmark builder is left to `MECH-D`:
  - drop the invalid `--outSAMflag None`;
  - add `--soloBarcodeReadLength 0`, because 10x 5′ with a 150 bp R1 aborted
    with "barcode length 150 ≠ 28" on SRR20710647's library;
  - decompress a gzipped whitelist;
  - share `reference_strategy.starsolo_geometry_args`.
- [x] `SW-16` — (2026-09-29) the `summary.txt` headline reported "Viral molecules … 0; Cells
  with viral reads 0/793308" beside 906,202 EBV molecules in
  `viral_summary.tsv`. `_headline_totals` filters virus *names* against
  `var_names`. Fixed: the call site passes the detected gene IDs, and
  `_headline_totals` raises if handed IDs that match no `var_name`. A new tiny-e2e
  assertion (`test_summary_headline_matches_viral_summary_total`) compares the
  headline with the `viral_summary.tsv` total; it fails on the old code.
- [x] `SW-17` — (2026-09-29) the `evidence --virus hhv6a/hhv6b/hhv8/kshv` selectors raised. `VIRUS_ALIASES` now lists every label each virus resolves to (legacy, catalogue, RefSeq) and matches case-insensitively, until `MECH-A` replaces it with taxid lookup. The cause was:
  `VIRUS_ALIASES` overwrites `hhv8`, and there is a `6B`/`6b` case mismatch.
- [x] `SW-18` — (2026-09-29) `anellovirus.gtf_text_for` returned `''` for uncatalogued
  accessions when `fasta_texts` is absent, so genomes silently lose their GTF.
  Its docstring promises a 1 bp exon.
  - Fixed: it now raises `ValueError` naming the unsized accessions.
  - It also iterated `accessions` twice, so a generator argument silently lost
    every placeholder. The argument is now materialised once.
  - Tests are in `tests/test_anellovirus_reference.py`.
- [x] `SW-19` — (2026-09-29, reproduced then fixed) `rerun-multimap` rewrote the *source* run's h5ad, because
  `from_yaml` runs before the `output` rewrite (menu.py ~595 vs 616). Write an
  e2e test through `_run_rerun_multimap` first, and fix only if it reproduces.
  - **Reproduced.** `TestRerunLeavesSourceUntouched` drives the real function
    with snakemake stubbed, and the source h5ad's sha256 changed.
  - Fix, part 1: the fast path resolves the h5ad under the copy's sample
    directory.
  - Fix, part 2: `cfg["output"]` keeps its trailing separator. Without it the
    Snakefile's `f"{config['output']}log/…"` paths would have been
    `…SAMPLElog/`.
  - Fix, part 3: the `--config` args come from
    `RunConfig.to_snakemake_config_args()` instead of a private `_arg_val`.

`G1` passes when all count invariants, schema checks, safety scenarios, rerun
consistency, and the full tiny workflow are green. No known correctness or data-
loss defect may remain.

Required gate:

```text
NUMBA_CACHE_DIR=/tmp/viralscan-numba-cache PYTHONPATH=src python3 -m pytest tests/ -q
python3 -m ruff check .
python3 -m ruff format --check .
PYTHONPATH=src python3 -m pytest -m "integration and not network" -q
```

### WP1D — Mechanism consolidation (architecture review, new 2026-09-29)

A read-only review found three walks of the code whose scientific decisions are
re-derived at every call site from string formats:

- virus identity and detection;
- the counting path and Run Config;
- reference construction.

The helpers pass their tests while the call sites produce wrong outputs.
`SW-13`…`SW-19` are the defects it observed. The user chose `MECH-A` first,
and settled three decisions in review:

- **A "virus" row** is one NCBI organism with rollup columns. The key is the
  taxid; segmented viruses are keyed by (species, strain or isolate).
- **An uncatalogued indexed gene** counts as viral if it is in `--gtf`, with a
  warning.
- **The catalogue** is merged into the packaged `virus_catalog.tsv`.

The review report (candidates A–F, the defect table, the builder × gate
matrix and MECH-A progress) is published as a private artifact:
https://claude.ai/artifact/QXWXwBk3BiJNSBnioYiUKH (2026-09-29).

- [x] `MECH-A` — (2026-09-30) per-Run **Virus Identity table**, `src/viralscan/virus_identity.py`,
  built once by the `analysis` rule. It maps gene_id → genome accession (t2g
  column 5) → catalogue row → viral status, virus key, name, family, sibling
  group and risk class.
  - [x] Step 1 (2026-09-29): catalogue merge. `virus_catalog.tsv` now holds
    4,128 accessions (504 species, every row with a taxid), up from 2,249.
    `build_virus_catalog.py` adds `taxid`, `organism`, `strain` and
    `serotype`, and carries forward the overlay columns `common_name`,
    `sibling_group`, `role` and `panel`. `extras/seed_catalogue_overlays.py`
    seeded them:
    - `common_name` on 303 rows, from the legacy prefix map;
    - `sibling_group` by taxid: HSV, HHV-6, and HHV-4 (EBV-1/2, F-017);
    - `risk_class=eve` on 3,448 rows;
    - `role=decoy` on 9 rows;
    - `panel`: `shipped` on 2,345 rows, `max` on 1,783.
  - Decision: CAT-31 now reconciles only `panel=shipped` rows.
    `catalogue_detection_targets(panel="shipped")` is the default, and
    `panel=None` checks every row. The max-panel rows are there so the identity
    table can name them; they are not claims of the shipped panel.
  - Fix: the builder keeps the newest version of each accession. The pre-merge
    catalogue had `NC_006312.1`, but the bundled GTF uses `.2`.
  - Finding for step 2: strain text is not a usable key. Two of the eight PR8
    segments spell the strain "A/Puerto Rico/8/1934(H1N1)". The isolate-level
    taxid (211044) is the same on all eight, but 24 Influenza A rows carry only
    the species taxid 11320. The segmented key therefore needs a normalised
    strain, with the `(HxNy)` suffix stripped.
  - [x] Step 2 (2026-09-29): `src/viralscan/virus_identity.py` and
    `tests/test_virus_identity.py` (32 tests).
    - `build_identity_table(t2g, gtf_gene_ids)` gives one row per indexed gene,
      with status catalogued / uncatalogued / host / legacy_prefix.
    - It fails when no gene is viral, and warns with a count for uncatalogued
      genes.
    - Golden tests pass on `t2g_v2`, the final `panel.t2g` and `t2g_max`, and
      skip where those files are absent:
      - every viral gene is catalogued, and every host gene is ENSG (41,145);
      - EBV-1 and EBV-2 are two viruses in sibling group HHV-4;
      - HPV45's two records form one virus;
      - the 18 SARS-CoV-2 genomes form one virus;
      - PR8 is one virus with 8 segments;
      - no virus holds a segment twice;
      - no two viruses share a display name.
    - The virus counts are 113 (v2), 236 (final) and 421 (max).
  - Decision, revising the "segmented = (species, strain)" key: segments are
    keyed by **taxid**. A taxid is split by normalised strain only when it holds
    the same segment twice.
    - Reason: strain text is unreliable inside one RefSeq set. It is empty on
      some segments ("Hantavirus Z10" M, "Pichinde" S) and spelled "…/1/96" or
      "…/1/1996" on others (goose/Guangdong, taxid 93838).
    - Only the species-level taxid 11320 ("Influenza A virus", 4 isolates)
      splits.
  - Decision (user, grill Q8b = A2, 2026-09-29): anelloviruses are keyed by
    **genus** (`ANDET-05`) and not by taxid. Genomes with no assigned genus are
    named "Anelloviridae (genus unassigned)", because a bare family row beside
    genus rows reads as their total.
    - 4 of 120 anellovirus taxids span several genera. They are the catch-all
      bins that hold most rows:
      - 2055263 "Anelloviridae sp." (Alpha, Beta, Gamma);
      - 68887 "Torque teno virus" (Alpha, Beta, Gamma, Samek);
      - 93678 "TTV-like mini virus" (Beta, Gamma, Het);
      - 432261 "Torque teno midi virus" (Gamma, Mem, Samek).
    - Catalogued anellovirus rows carry no `common_name`. (An earlier count of
      "15" treated the table's genus label and the same genus from the
      catalogue as different labels.)
    - The genus comes from the anellovirus accession table. Otherwise it falls
      back to the catalogue genus, when that is an ICTV genus name (one word
      ending in "virus"), then to "Anelloviridae". This keeps the non-genus
      lineage token "Small anellovirus" out of the keys.
  - Step 3 caveat for step 4: a run directory made before 35940ec has no
    `results/virus_identity.tsv`. Snakemake does not rebuild it, because no rule
    consumes it yet (checked by dry-run on a skeleton copy of
    `combined_corrected/SRR8315713`: only detection, umap and all are
    scheduled). Once step 4 makes it an input, `rerun-multimap` and resume must
    build it in process for old runs, or the `analysis` rule reruns.
  - Guard (user, grill Q8 = B): keep the structural guard below and add an
    index build manifest that records the host and viral gene sets. The
    manifest overrides `--gtf`, and a run is refused when `--gtf` would mark a
    manifest-host gene as viral. The manifest, fixing `reference_strategy.py`'s
    `-gtf`, and detection reading the table are all step 4 work.
    Structural guard: a `--gtf` gene whose t2g
    column 5 is a transcript of the index is host. This is the host-cDNA row
    shape of every combined index measured: 465,769 host rows, 0 viral rows.
    - Reason: `reference_strategy.py:704` passes the combined GRCh38+viral GTF
      as `-gtf`. Today's `analysis.py` makes every host gene in that path viral.
      Our SLURM runs passed viral-only GTFs and are unaffected.
  - [x] Step 3 (2026-09-29): the `analysis` rule writes
    `results/virus_identity.tsv` as a declared output, next to
    `log/analysis.txt`, which is unchanged.
    - `analysis.write_identity_table` builds the table from `config.transcripts`.
    - A run whose index has no viral gene now stops at this step.
    - The column reference is in `docs/output_reference.md`.
    - The tiny e2e (20/20) checks that the fixture's 3-column t2g takes the
      legacy fallback.
    - No consumer reads the table yet.
  - [x] Step 4 (2026-09-30, branch `v3/mech-a-consumers`: 63c3de1 phase 1,
    160a4db phase 2): the consumers read the table.
    - `multimap.py` takes its viral/host partition from the `viral` column.
    - `detection.py` groups by `virus_key`, names by `virus_name`, and takes
      the sibling note from `sibling_group`. The weaker members of a group
      are flagged against its dominant member, which generalises the old
      pairs. `eve_risk` comes from `risk_class`.
    - `umap`, `hostresponse`, `evidence` (`--virus` resolves by key, taxid,
      handle, name, organism or family) and `gene_programs` (by `virus_key`)
      also read the table.
    - The consumer rules declare `virus_identity.tsv` as an input.
    - `SIBLING_VIRUS_PAIRS`, `EVE_RISK_GENERA` and `merged_name_map` remain
      only as the fallback for runs with no table. CONTEXT.md has the new
      terms.
    - Tests: `test_mech_a_consumers.py` and `test_mech_a_presentation.py`.
      They cover every retired sibling pair and EVE genus, per-gene →
      per-group conservation including `genus:` keys, the partition,
      selectors and programme lookup. Main tree: 1,519 passed, 0 skipped.
    - Decision (EVE, user-confirmed 2026-09-30): a virus is flagged when any
      of its genes has `risk_class=eve`. An empty `risk_class` on a catalogued virus means no
      EVE risk (680 rows are empty; flagging them would flag EBV).
      Uncatalogued and legacy viruses fall back to the genus-name test, with
      a warning.
    - Count parity was checked against 815838a on copies of the EBV, HSV-1,
      HHV-6B and covid x213 runs (`viralscan_work/parity_mechA/`).
      `count_audit`, `found_genes`, `analysis.txt` and `positive_control` are
      identical in all four; every multimap layer is identical for EBV,
      HSV-1 and HHV-6B.
    - Explained naming diffs:
      - `Anelloviridae` → "Anelloviridae (genus unassigned)", 274 genes;
      - "Torque teno virus" and raw `D1P6x_gpN` IDs → `Alphatorquevirus`;
      - covid `VARVgp184`: VZV → Variola virus, which is correct (1 UMI);
      - `n_viral_accessions_in_reference` 9193 → 5093, because only indexed
        genes are counted.
    - Name checks pass: EBV "Epstein-Barr virus"; HHV-2 `possible_em_bleed`
      on HSV-1 at 665:1; the EBV-1/EBV-2 split is unchanged.
  - [x] Step 4a (2026-09-30) — **structural-guard false positive** (covid
    x213): gene `HUM_HERP6B_DR1` has t2g column 5 equal to its own
    transcript ID. The guard therefore marked a viral gene as host: 1 UMI in
    1 cell, `counts_unique_viral` −1.
    - t2g shape, measured: the covid index has 97 **self-named** rows
      (transcript = gene = column 5), all `HUM_HERP6B_*`. They are VIRTUS-style
      records, one per gene, so DR1 was only the one with a UMI. Host cDNA
      rows that point at a transcript never have transcript = gene (`ENST` ≠
      `ENSG`): 0 of 465,769 host rows in each of the five stored combined
      t2g files (covid, v1 full panel ×2, final 2,343, max).
    - Fix (`virus_identity.read_t2g`): a self-named row names no genome. Its
      accession is `""`, and it is not structural-host. `_resolve_genes`
      names a GTF gene with no accession by the legacy prefix map (status
      `legacy_prefix`, key `name:<virus>`), so the 97 genes are one virus,
      "Human herpesvirus 6b", not 97 accession-keyed ones. The exemption
      does not use column 3: 3,037 viral rows of the final panel have it empty.
    - Covid parity (`viralscan_work/parity_mechA/covid/final4a/`, job
      25684236): the 97 genes form one viral row, and HHV-6B gets its 1 UMI
      back. Every multimap layer, including `counts_unique_viral`, is now
      identical to 815838a. The HHV-6 note does not fire: HHV-6 has 38 UMI against HHV-6B's
      1, which is under the 50:1 threshold, and legacy-prefix rows carry no
      `sibling_group`. Detection warns that `eve_risk` comes from the genus
      fallback, as designed for uncatalogued viruses.
    - Tests: two in `tests/test_virus_identity.py` (a self-named row, and
      one prefix-named virus next to an `ENST`/`ENSG` self row that stays
      host). Full suite: 1,527 passed.
  - [x] Step 4b (2026-09-30) — **catalogue row for NC_000898.1** (HHV-6B,
    taxid 32604, `common_name` "Human herpesvirus 6b", `sibling_group`
    HHV-6). The catalogue held HHV-6B only as AF157706.1. On the stored
    SRR20710641 index HHV-6B therefore read `NC_000898.1`, and the HHV-6A/6B
    note did not fire.
    - `panel=legacy`, a new value, not `shipped` as first written here.
      `NC_000898.1` is in neither built panel. The shipped `viral.fa` has
      HHV-6B only as AF157706.1, and max-panel dedup dropped it as an
      `exact_hash` duplicate of AF157706.1 (`dedup.tsv`), the CAT-05 pattern.
      Only pre-v3 stored indexes (VIRTUS2-sourced) carry it. A `shipped` or
      `max` value would claim a panel that does not index it, and a
      reconciliation of that panel would report it as catalogued but not
      indexed. Reconciliation skips any other `panel` value. The identity
      table does not filter by `panel`, so naming is unaffected. It is the
      only catalogue row missing from its panel's FASTA.
    - Checked on a copy of the stored run with the working tree
      (`viralscan_work/parity_mechA/hhv6b/final/`, job 25683797). It now
      reads "Human herpesvirus 6b" (5,596 UMI), and the note fires on
      "Human herpesvirus 6" at 323:1.
  - [x] Step 4c (2026-09-30): **in-process backfill for old run dirs**.
    - `virus_identity.backfill_identity_table(config)` builds the table from
      the run's own `log/analysis.txt`, not a fresh GTF glob, so catalogue
      or GTF changes since the run cannot move its partition. The file takes
      `analysis.txt`'s mtime.
    - `rerun-multimap` calls it after rewriting the config, and so does
      `--resume` before snakemake starts. `write_identity_table` moved from
      `scripts/analysis.py` into `virus_identity.py` (re-exported), so
      `menu.py` can call it without the script's logging setup.
    - Tests in `tests/test_identity_backfill.py`, including a byte-identical
      match with the table the `analysis` rule writes.
    - Dry run on a timestamp skeleton of the stored SRR8315713 run, with
      `detection.done` and `umap.done` removed (what `rerun-multimap`
      does): without backfill, analysis + multimap + detection + umap are
      scheduled; with it, only detection + umap.
    - A complete old run with every target present schedules nothing either
      way, because snakemake does not rebuild a missing intermediate for
      up-to-date targets.
  - MECH-A stays `[~]` until 4a and 4b land and the HHV-6B name check passes
    on the stored index. All three were met on 2026-09-30, and MECH-A is
    closed.
  - It replaces 7 prefix-matching call sites, the GTF-only viral/host
    partition, the name-keyed `SIBLING_VIRUS_PAIRS`, and the substring EVE test.
  - Observed failures it fixes:
    - the HHV-6 sibling note is empty at 297:1 (SRR20710641);
    - one genome appears under two names (`TTVgp1` vs `NC_002076.2_gene1`);
    - the 1,912 max-panel accessions would report as raw IDs;
    - EBV-2 bleed shows as 80 unflagged per-gene rows (F-017).
- [x] `EXPL-R2.10`: exploratory reruns on the current reference (grill
  R2.10 = D, 2026-09-29). These do not tune anything (R2.2, R2.4). Submitted
  on 2026-09-29:
  - 25672273 `scripts/slurm_quant_max_corrected.sh`: F-017 redone on the max
    panel with barcode correction on, for EBV and HSV-1. Output goes to
    `viral_panel_max_2026-09-28/runs/combined_max_corrected/`.
  - 25672274 `scripts/slurm_quant_hhv6b_5p.sh`: HHV-6B SRR20710641 as 5′
    `-x 10xv2` on the v1 index. Its R1 has the TSO at base 27, so the barcode
    is 16 bp and the UMI 10 bp. Output goes to
    `ebv_latest_ref_2026-09-27/runs/combined_corrected_5p/`.
  - 25672275 `scripts/slurm_quant_covid_exploratory.sh`: covid x213 and x216
    with current code, the published index, and the CellRanger whitelist. R1 is
    28 bp with no TSO, so `-x 10xv3` geometry is correct. Barcode correction
    was already on in the published run, and there is no read filter yet, so
    the F-019 poly-G signal is expected to come back. Output goes to
    `covid_viralscan/results_v3_exploratory/`.
  - 25672276 `scripts/slurm_strand_test.sh`: 4M read pairs of each 5′ library
    (HHV-6B, x213, x216), run through `kb count --strand`
    forward/reverse/unstranded. Output goes to
    `viralscan_work/strand_test/*/strand_summary.tsv`.
  - Results (2026-09-30):
    - **F-017 is confirmed with barcode correction on:** EBV-2 takes 16.1 %.
    - **HSV-1 max panel:** +0.05 %.
    - **HHV-6B as 10xv2:** 5,596 molecules in 3,556 cells. It still reports as
      `NC_000898.1_gene1`, the MECH-A naming defect.
    - **Covid on current code:** 1.08M / 1.50M Alphatorquevirus, which
      reproduces F-019.
    - **F-020, the strand default:** kallisto's forward default keeps only
      6.5–8.9 % of reads on every 5′ library. `--strand` (DEF-02) is now a
      priority.
    - **Tonsil read origin:** reproduces the F-019 artefact, with 3.02M reads
      and 0 clean full-length viral reads.
- [x] `EXPL-TTMDV` — **Result 2026-10-03: artefact** (F-020 "EXPL-TTMDV
  result"). The reproduction gate passed exactly (x213 55,046 / 24,658;
  x216 7,391 / 1,829).
  - 0 clean, full-length, viral-best reads in either sample. 100 % carry
    homopolymers and 96 % carry the 10x TSO: they are TSO–oligo-dT–poly-G
    concatemers.
  - They hit a 26-nt poly-A tail after `AATAAA` in the covid index's
    AB303552.1/AB303557.1, which is the F-021 sink class. cat42b masks it.
  - gget blast: 0 of 14 reads hit an anellovirus.
  - It says nothing about donor TTV carriage, which is commensal and
    expected. My earlier "no homopolymer over 12 nt" was measured on the
    masked cat42b FASTA, and is corrected.
  - Original row: (2026-10-03, user-requested) read-level check of the
  reverse-strand Gammatorquevirus signal in covid x213/x216 (F-020 update):
  55,046 + 24,658 unique molecules on AB303552.1/AB303557.1 under
  `--strand reverse`. Exploratory; it tunes nothing. Anelloviruses are
  commensal (F-022), so the check is neutral on the prior.
  - Method, reused from F-019 (`viralscan_work/f005_readorigin/`), outputs in
    `viralscan_work/ttmdv_readcheck/<sample>/`:
    1. `kallisto bus -n --rf-stranded` (0.51.1, same index, `-x 10xv3`).
    2. Capture the anellovirus ECs, then correct.
    3. Extract the exact R2 reads.
    4. Competitive `minimap2 -ax sr` against GRCh38 + the covid panel.
    5. `classify`, then a spliced alignment to the TTMDV genomes.
    6. gget blast (login node) of the per-sample consensus and about 20
       position-stratified clean reads.
  - **Reproduction gate (before classifying):** the captured molecules on
    AB303552.1/AB303557.1 must match the reverse run (x213: 55,046 /
    24,658). Stop if they do not.
  - **Pre-registered reading:**
    - **Real TTMDV:** most molecules clean (no homopolymer ≥15, no TSO, not
      low entropy), full-length and viral-best; reads spread across the
      transcribed region (not one hotspot); blast top hits are anelloviruses.
      Variants recurring across independent CB-UMIs point to a donor strain.
    - **Artefact:** most molecules low-complexity, TSO or host-best; a single
      hotspot, especially in the GC-rich non-coding region; blast hits human
      or vector.
    - A consensus identical to the reference at every covered site is a
      contamination/reference signature.
    - Not evidence: strand (guaranteed sense by `--rf-stranded`), and
      identity ≥ 0.9 to AB303552 (donor strains diverge).
  - **Gate run 1 failed (job 25695673, 2026-10-03).** Re-counted vs run,
    in unique molecules:
    - x213: 57,889 vs 55,046 (+5.2 %) and 25,109 vs 24,658 (+1.8 %);
    - x216: 7,947 vs 7,391 (+7.5 %) and 1,956 vs 1,829 (+6.9 %).

    Always high, which fits a gate defect: the capture kept only
    target-EC records, so a molecule with other-EC reads looked
    target-unique. The tolerance is unchanged (2 %). Run 2 (job 25695740)
    gates on every record of each target-touching UMI (`bustools capture
    -u`). If it still fails, stop.
    - Descriptive only: identity distribution; called vs empty-droplet
      barcodes (ambient); shared CB-UMIs between x213 and x216. Donor
      identity of x213/x216 is unknown (batch 1/2), so a cross-sample
      consensus comparison is descriptive.
- [x] `EXPL-HPV16` — (2026-10-01) HPV16 positive control on GSE189670 (Bedard et al.,
  Nat Commun 2023, PMID 37031202), user-requested 2026-09-30. Isogenic NIKS
  keratinocyte rafts, 10x 3′ v3: SRR19537341 (GSM5705760, "HPV16 infected
  keratinocytes") vs SRR19537339 (GSM5705759, "Normal keratinocytes"), one
  run each. Shipped 2,343-genome panel (`viral_ref_final/build/panel.idx`),
  `-x 10xv3`. Exploratory; it tunes nothing.
  - Truth is qualitative only: the authors see HPV16 early and late genes in
    the HPV16 rafts and none in the parental rafts, with no per-cell % given.
    Pass = HPV16 present in the HPV16 sample, and ~0 in the normal one.
  - ENA serves R2 only for these runs, so reads come from SRA with
    `fasterq-dump --include-technical --split-files`.
  - Gene identity: the bundled GTF's `gene` attribute maps gp1=E6, gp2=E7,
    gp3=E1, gp4=E2, gp5=E1^E4, gp6=E5, gp7=L2, gp8=L1. The genes overlap,
    and 3′ reads pile up at the early and late polyA sites, so per-gene
    counts reflect polyA position, not ORF identity.
  - Scoping of the other user-supplied datasets (GSE164690, GSE208653,
    CELLxGENE HPV/CMV, E-MTAB-13687 tonsil):
    `viralscan_work/dataset_scoping_2026-09-30/`.
  - Result (2026-10-01, job 25684772, HEAD f40f72a): **HPV16 PASS on
    presence.** SRR19537341 (HPV16) has 7,725 molecules in 2,442 called
    cells. SRR19537339 (normal) has 8 molecules in 6 cells, about 970:1.
    - Per gene, HPV16 sample, as unique + allocated:
      E5 2,374; E1^E4 1,501 (1,478 allocated); E7 1,440; E1 1,100 (all
      allocated); E2 696; E6 451; L1 112; L2 51. Early genes dominate, which
      matches the authors' "mostly early genes". E5 sits next to the early
      polyA site, as expected for 3′ reads.
    - Both array tasks were marked FAILED only because the script's final
      `test -s $OUT/results/...` was wrong. viralscan nests results under
      `$OUT/<sample>/`. Snakemake completed every rule.
    - The cell-calling knee is not usable: total ≥10 gives 60,475 /
      1,569,733 (HPV16) and 86,135 / 1,292,165 (normal) barcodes. This is
      `SW-23`.
    - Indicative per-cell result, cut at the barcode-rank steepest descent
      (total ≥952 / ≥739, computed outside viralscan): 1,869 / 8,023 HPV16
      cells are HPV16+ (23.3 %), holding 7,034 of 7,725 molecules. The
      normal raft has 0 / 11,045, and all 8 of its molecules sit in empty
      droplets.
    - **emptyDrops rerun** (2026-10-01, job 25689583, detection only on
      copies in `runs_emptydrops/`; FDR 0.01, lower 100, seed 100):
      | | HPV16 raft | normal raft |
      |---|---:|---:|
      | called cells | 16,079 (knee 9,895, inflection 1,052) | 11,884 |
      | HPV16+ called | **2,012 (12.5 %)** | **0** |
      | HPV16+ comparable (≥200 host UMI) | 1,927 / 9,086 (21.2 %) | 0 / 9,793 |
      | Gammatorquevirus+ called | 2,552 (15.9 %) | 3,314 (27.9 %) |
      These are the reportable numbers. The rerun reused the run's stored
      `virus_identity.tsv`, so it still shows the old label "16,18"; new runs
      get the fixed name. The anellovirus row is the F-021 artefact, now
      27.9 % of the normal raft's called cells.
    - **Specificity failure, anellovirus.** Gammatorquevirus has 5,688
      molecules (HPV16) and 10,037 (normal) in an anellovirus-free cell
      line. In the normal sample all 10,037 sit on one placeholder gene,
      `KP343824.1_gene1` (an "UNVERIFIED" TTV isolate S57 record whose
      first 29 nt are poly-T; the builder masks only runs ≥31). This is the
      same single-bucket sink as F-005/F-019. Read-level check: of the
      first 20M R2 reads, 27,098 carry a 31-mer from that poly-T head. They
      are host mRNA ends running into **poly-A tails**, not the poly-G reads
      of F-019. See F-021 and `CAT-42`. KP343824.1 is also in the max
      panel, yet the 3′ EBV (10x v2) and HSV-1 (Drop-seq) runs put 0 on it,
      so the capture depends on the library, not on 3′ chemistry alone.
    - The label was wrong: HPV16 reported as "Human papillomavirus 16,18".
      Fixed in this commit (catalogue `common_name` and both legacy prefix
      maps). `HUM_PAP_1618_*` genes are HPV16 only; HPV18 is NC_001357
      under `HPV18_*`.
    - Low-level calls present in both samples at similar levels (HPV118 15–18,
      HPV29 17, HHV-6 11–12) are background, not HPV16-specific.
  - **Read-level validation (2026-10-02, job 25694035; scripts in
    `hpv16_gse189670/readcheck/`): HPV16 is real.**
    - The scan found 13,119 R2 reads with a non-homopolymer HPV16 31-mer
      in the HPV16 raft, against 10 in the normal raft.
    - Competitive minimap2 against GRCh38 + the viral panel:
      - 11,349 (86.5 %) are clean, full-length and HPV16-best;
      - 38 (0.3 %) are host-best;
      - 33 are low-complexity.
    - 99.2 % read in the mRNA-sense direction (11,254 vs 95), so this is
      RNA, not DNA.
    - The pile-up sits just upstream of the early polyA site: 4,173 reads
      start at nt 3,000–3,250 of the indexed layout, which numbers from the
      E1 ATG, so E6 = 7,125.
    - **Canonical HPV16 splice junctions** (K02718 numbering, on a genome
      rotated to start at nt 7,500):
      - 880^3358 (E1^E4): 81 reads;
      - 880^2709: 29;
      - 226^409 (E6*I): 25;
      - 880^3391/3361: 22;
      - 226^526: 1.
    - gget BLAST of 12 reads sampled along the genome: 11 are full-length
      HPV16 at e = 7e-38, and 1 is a 71-nt HPV16 chimera.
    - The normal raft's molecules are index hopping: 4 of its 9 clean
      HPV16 CB-UMI pairs also occur in the HPV16 raft.
    - Caveat: ~5,000 reads sit in E6/E7/E1, far from the early polyA. They
      are still sense-strand and clean, so they are most likely internal
      priming on viral RNA. Not investigated.
- [x] `SW-24` — `MoleculeAudit.validate` (`multimapping.py:72`) compares
  allocated ambiguous mass with an absolute `atol=1e-9`. On EBV SRR12682296
  with the CAT-42 index, 10,250,998 ambiguous molecules sum to
  10,250,998.000000002 (diff 1.9e-9, relative 1.8e-16, float rounding), and
  the run fails (job 25691839, probe 25694067). Fix: a relative tolerance,
  plus a test with ~10M fractional shares. Blocks CAT-42 check 3.
  - Done 2026-10-02 (`v3/sw24-mask`, 08ee19a): `rtol=1e-9, atol=1e-9`.
    `tests/test_multimapping_tolerance.py` covers 3.3M thirds in a csr matrix,
    a +0.5 mismatch that still raises, and the observed 10,250,998 + 2e-9
    case, which `rtol=0` would reject.
- [x] `SW-23` — `--cell-calling knee` puts the knee in the empty-droplet
  tail (found 2026-10-01). `cellcalling.knee_cells` takes `argmin` of the
  signed distance (the point furthest *below* the chord). On a barcode-rank
  curve the cell plateau lies above the chord and the empty tail below, so
  the knee lands at about `knee_min_umi` (10).
  - Every knee run this week logged "knee at total>=10": HHV-6B 15,800,
    covid 118,061, HPV16 60,475, normal 86,135 cells.
  - On GSE189670 the steepest descent sits at rank ~8,000 / ~11,000, and
    `argmax` (above the chord) at rank 3,061 / 4,690.
  - Every `infected_called` / `pct_infected_called` from a knee run is
    inflated in its denominator.
  - **Decision (user, 2026-10-01): emptyDrops for every reported number.**
    `knee` stays in the code only because the frozen protocol lists it
    under `sensitivity_only_callers`. Its estimator is not fixed. Every knee
    run now logs a WARNING (sensitivity-only, SW-23), and the help text and
    docs say so (`cli_reference.md`, `faq.md`, `output_reference.md`).
  - The three tracked SLURM scripts (`slurm_quant_covid_exploratory.sh`,
    `slurm_quant_hhv6b_5p.sh`, `slurm_quant_max_corrected.sh`) now use
    `--cell-calling emptydrops`. They put
    `conda/envs/R4_51/bin` (R 4.5.1, DropletUtils 1.30.0) at the end of
    `PATH`, because the bench env has no R. The quickstart vignette is also
    switched.
  - `TestKneeCells` passes with the bug because its data is a trivial
    bimodal split. It was left unchanged, and the new test only asserts the
    warning.
  - Check before quoting: the covid `pct_infected_called` figures under
    `HOST-03` (143,243 → 28,921 called cells) may come from knee runs.
- [ ] `ANELLO-PRIOR` — anelloviruses are a commensal virome, detectable in
  most people without disease (user, 2026-10-01). So an anellovirus call in
  human tissue or blood is biologically expected and must not be filtered as
  noise by family. Any gate has to be read-level: are the assigned reads
  low-complexity (homopolymer/poly-A/poly-G), and do they land on viral
  sequence along the genome, near the viral polyA, or only on a homopolymer
  tract? This applies to F-005, F-019, F-021, CAT-42 and DEF-01 wording.
  - First application (2026-10-01): on GSE189670 normal raft, 0 of 27,098
    KP343824.1-hitting reads (20M scanned) carry any non-homopolymer viral
    31-mer, so that call is the poly-A sink (F-021). The scan script
    (whole-genome k-mer position profile) is the template for the gate.
  - Detection audit (2026-10-01, F-022), against Kane et al. 2026 (plasma
    DNA prevalence 79 % young / 100 % older). No human anellovirus call
    has been read-validated yet. The default call is 1 UMI on one gene,
    with no read-level gate. The CDS-only models leave out the 5′ region
    (median 513 nt) and, for ~10 % of genomes, ≥200 nt before polyA. The
    tonsil whole-genome null (F-010) still says cellular mRNA is rare.
    Sub-steps, in order:
    - [ ] `ANELLO-PRIOR.1` plant positive control: held-out genomes,
      3′/5′-end reads into a real host library. Includes a TTV-high
      public dataset if one exists.
    - [ ] `ANELLO-PRIOR.2` rerun `measure_kmer_capture.py` on CDS exons
      vs exons + UTR. Its F-011/F-012 numbers were measured on whole
      genomes. Extend models to the mRNA extent only (never whole genome),
      judged on the negative controls too.
    - [~] `ANELLO-PRIOR.3` default read-level **diagnostic** (not a gate).
      - **Re-specified 2026-10-04.** The original spec — "homopolymer/entropy
        fraction and position along the genome" — names two measures the
        adversarial review ruled non-diagnostic. 10x 3' chemistry sees only
        transcript ends, so a genuine anellovirus read is poly-A rich by
        definition; and an artefactual poly-A read lands on the genome's
        longest templated A-tract, which in an anellovirus *is* the polyA
        site. Position and homopolymer fraction cannot separate the
        hypotheses. See F-019, update 2026-10-04.
      - **User decision, 2026-10-04: tune for not missing a real infection.**
        Flag, never filter; the confirmatory checks come afterwards. This
        follows `ANELLO-PRIOR` — anelloviruses are commensal, a call is
        expected, and every homopolymer-based filter removes precisely the
        genuine 3'-end reads needed to prove one. Such evidence *bounds* a real
        infection; it never excludes one.
      - **Done 2026-10-04.** Three measured columns, all labels:
        - `alignment_complex_body_fraction` — the decisive measure. The body
          5' of the first >=15-nt run must be >=20 nt at >=2.0 bits
          dinucleotide entropy. Thresholds measured, not assumed: over 5,000
          random ACGT draws the genuine minimum is 2.21 bits at 20 nt, while
          the artefact classes top out at 1.70 (poly-A 0.00, AC 1.00, CAG
          1.58, A-rich 1.70). 2.0 sits in that gap. The review's 3.5 was
          measured on 90-nt reads and would discard half of all genuine 25-nt
          bodies. `MIN_BODY_LEN` is 20, not the review's 25, because a 3'UTR
          ending in A's merges them into the tail (<=6 nt lost, measured).
        - `alignment_tso_fraction` — TSO in R2 in either orientation at <=2
          mismatches. It cannot occur in a genuine molecule: 10x 3' chemistry
          sequences only fragments carrying the bead-oligo end. Entropy does
          not catch the TSO (H = 3.56), which is why both columns exist.
        - `alignment_median_query_coverage` — identity without the length it
          was measured over is not evidence. **Needed before `ANDET-09f`,**
          which sweeps the very coverage filter that currently bounds it.
      - Also fixed: `scripts/plant_anello_10x.py` gave its `3p` plants no
        untemplated poly-A tail, so any poly-A-sensitive measure scored on that
        set would have looked free. `truth.tsv` now records `body_len` and
        `tail_len` per read, so recovery can be scored against body length.
        Pinned by `tests/test_plant_anello_tails.py`.
      - Caught in review, fixed before the columns were trusted: SAM stores
        SEQ reverse-complemented on a reverse-strand record, which **inverts**
        both read-side measures — a genuine `[body][poly-A]` arrives as
        `[poly-T][rc body]` and scores as artefact, while a `[poly-A][TSO]`
        chimera arrives TSO-first and scores as complex. With `--soloStrand
        Unstranded` that is about half of all records. `Alignment.read_seq()`
        now restores sequencing orientation; pinned by a test that fails
        against the unfixed code.
      - [x] Remaining, and **the four arms as specified will not give a
        reading**: covid x213 and the synthetic negative yield ~0 aligned reads
        through the branch (its 0.80 coverage filter already rejects the
        chimeras), so their columns come out *empty*, which is "not measured",
        not "artefact"; and HPV16 is not an anellovirus, so the branch never
        touches it. Validate instead by running `is_complex_body` / `has_tso`
        directly over the 19,785 covid anellovirus-aligned reads already
        extracted — which is also the full-set body census that turns the
        review's 0/30 into a bound. Un-revcomp any read taken from BAM SEQ.
        - **Done 2026-10-04** (`scripts/anello_body_census.py`; F-019 update
          "full body census"). Bound: no body > 26 nt places on any
          anellovirus genome; ≤ 13/19,785 reads (CP95 0.11 %) / ≤ 13/2,791
          CB+UMI (CP95 0.80 %) are even candidates, all 20–26 nt, chance-level.
          **The measures failed validation**: `has_tso` misses TSO truncated
          at the read edge, and `read_body` takes a 5′ TSO/TruSeq fragment as
          the body (H ≈ 3.5). 6,171 reads kept, 98.9 % of them with a body that
          places nowhere. Fix tracked in the row below.
      - [ ] Fix the reagent leak found by the census (2026-10-04). `has_tso`
        should match a TSO overlapping the read edge (flag-only, fits the
        locked decision). Pin it with the census's real reads as fixtures, and
        re-run `scripts/anello_body_census.py` to show the leak closed.
        **User decision:** do TruSeq R1 and other reagent hits go into
        `tso_fraction`, or into a renamed `reagent_fraction`? And should
        `read_body` strip a leading reagent before measuring the body?
      - [ ] The columns are reachable only with `--anello-align`, which ships
        off, so a default run's anellovirus call still carries nothing but
        `artifact_risk`. Apply the same pure functions in `viralscan evidence`,
        which `output_reference.md` already names as the route for checking a
        kallisto call's reads. **This puts `REL-16` (the segfault that breaks
        `viralscan evidence`) ahead of `ANDET-09f`.**
      - [ ] Expect the regenerated planted arm to lose roughly a quarter of its
        3' reads: a read starting >228 nt into the 300-nt window carries a tail
        of >=19 nt, dropping coverage below 0.80. That is the branch filter, not
        the new columns, and it is exactly the sensitivity cost the permissive
        decision rejects — so score 3' recovery by `body_len` bin from
        `truth.tsv` and make that the axis of `ANDET-09f`. Re-check criterion 3
        as well: A-rich tails may now be lost to the host filter.
    - [x] `ANELLO-PRIOR.4` (user decision) relabel `eve_risk` for
      Anelloviridae; its basis is the F-005 mechanism F-019 revised.
      - **Done 2026-10-03** (user chose "clear + interim note" in a grill).
        Anellovirus `eve_risk` is now False.
        `viral_summary.tsv` gains `artifact_risk`, set to `low_complexity`
        (F-019/F-021). It is a diagnostic label, never a filter, until
        `ANELLO-PRIOR.3`'s measured metric replaces it.
      - Implementation:
        - All 3,448 catalogue `risk_class=eve` rows became
          `low_complexity`. They are Anelloviridae only, and no other column
          changed.
        - `EVE_RISK_GENERA` is now empty. The mechanism is kept for a real
          EVE family; inherited ciHHV-6 is parked.
        - The new `LOW_COMPLEXITY_RISK_GENERA` seeds the label, and the
          genus-name fallback covers uncatalogued names.
      - Basis: a quick Europe PMC search found no germline human anellovirus
        EVE; the one integration is somatic, in the SKNO-1 cell line (PMID
        42671192).
      - cat42b probe: 8 anellovirus groups are labelled, the other 222
        viruses are not, and `eve_risk` is False everywhere.
    - MECH-B (group → sum → threshold) must close before any anellovirus
      sensitivity claim.
- [x] `CAT-42` — homopolymer mask misses runs under 31 nt (F-021,
  2026-10-01). KP343824.1 starts with 29 T and captured poly-A tail reads
  in GSE189670 (10,037 molecules in an anellovirus-free cell line). The 3′
  EBV and HSV-1 max-panel runs put 0 on it, so the capture is
  library-dependent.
  - Panel scan (2026-10-01): 38 runs of ≥20 identical bases in 38 of 2,343
    genomes (32 poly-A, 3 T, 2 G, 1 C; 9 are 25–30 nt). Every
    Gamma/Betatorquevirus molecule in both GSE189670 rafts (5,688 + 941 and
    10,037 + 727) landed on one of those genomes. The ~250 "genus
    unassigned" molecules did not.
  - `build_bundled_panel_ref.py --homopolymer-run-length` default is now
    20 (was 31). Rebuild job 25691632 → `viral_ref_cat42/build/`. After it
    finishes:
    1. check `kallisto inspect` against CAT-40 (max EC size, ECs
       discarded);
    2. rerun GSE189670 (both rafts) for anellovirus ~0 and HPV16
       unchanged;
    3. rerun the EBV/HSV-1 regression (molecules within ±5 %).
  - Queued 2026-10-01 with `afterok:25691642`:
    `viral_ref_cat42/quant_gse189670.sbatch` (25691838, emptyDrops, output
    in `hpv16_gse189670/runs_cat42/`) and `viral_ref_cat42/regression.sbatch`
    (25691839, output in `viral_ref_cat42/runs/`). Both use
    `viral_ref_cat42/viral_panel.gtf`, cut from the new `combined.gtf`. In
    it, anellovirus gene IDs carry their accession prefix, so 124 gene IDs
    differ from the old panel GTF; HPV16 is unchanged.
  - The rename hits 124 RefSeq anellovirus genes (`TTV3_gp1` →
    `NC_014081.1_TTV3_gp1`, `D1P64_gp1` → `NC_038336.1_…`). Identity is not
    affected, because t2g column 5 is still the genome accession, so those
    genes stay catalogued by accession. The output gene IDs do change.
    Before this index ships, find the commit that changed the builder's
    gene-ID naming (not yet traced).
  - Pass criterion for check 2, judged by accession, not by genus label:
    molecules on KP343824.1 and on the other 37 genomes with a ≥20-nt run
    are ≈0, and HPV16 stays 7,725 molecules / 2,012 cells. *Corrected
    2026-10-02:* the ~250 "genus unassigned" molecules are **not** a
    positive signal. They appear equally in the normal raft, an
    anellovirus-free cell line, so they are artefacts too.
  - **Composition mask (2026-10-02, `v3/sw24-mask`, 665b274):**
    - `build_bundled_panel_ref.py --lowcomplexity-kmer-mask`, off by default.
      It N-masks each 31-nt window with ≤2 distinct bases and a base ≥20, or
      any base ≥28, in **Anelloviridae records only**. A review simulation of
      the panel-wide rule masked EBNA-2 and HSV-1 s-gene sequence. It also
      writes `lowcomplexity_mask.tsv`.
    - Dry run on `viral_ref_cat42/build/viral.fa`: 1,586 bases masked across
      46 records. HPV16, EBV-1/2, HSV-1 and HSV-2 lose 0 bases. KP343822.1
      loses 34 bases, all inside nt 2405–2460. KP343842.1 loses 35 of the 38
      nt at 3306–3343.
    - The rebuild goes into `viral_ref_cat42b/`. It is checked against cat42,
      so the mask is the only variable.
  - **Results (2026-10-02):**
    - Build 25691642: `BUILD_OK`, 1,346 bases masked. `kallisto inspect`:
      max EC 9921, 0 ECs discarded, 471,944 targets. **Check 1 passes.**
    - GSE189670 (25691838): Gamma/Betatorquevirus 10,037 → **0**. HPV16
      is still 7,725 molecules (2,007 / 15,502 called cells). "Genus
      unassigned" stays at 246 / 273 in both rafts (KP343822.1,
      KP343842.1), so the artefact persists below the 20-nt mask.
      **Check 2 passes for KP343824.1, with one residual artefact.**
    - Regression (25691839): HSV-1 completed. EBV failed in multimap, on
      the conservation check, with `SW-24`. **Check 3 is blocked.**
  - **cat42b results (2026-10-03; build 25694852, checks 25694853/25694854,
    code d01043c). Every gate passes, so cat42b is now the current panel**
    (user decision 2026-10-02: replace if every check passes). Compared with
    cat42, so the composition mask is the only variable:
    - `kallisto inspect`: 0 ECs discarded, 471,944 targets, 84,860,192 k-mers
      (cat42: 84,860,149; kallisto fills N with pseudorandom bases).
    - GSE189670: "genus unassigned" 246 → **0** (HPV16 raft) and 273 → **0**
      (normal raft). KP343822.1 220 → 0, KP343842.1 26 → 0. HPV16 stays at
      **7,725** molecules (2,006 / 15,615 called cells; cat42 2,007 / 15,502.
      emptyDrops runs on the changed matrix, so the called set moved). The only anellovirus left is 1 Samektorquevirus molecule in the
      normal raft, as on cat42.
    - EBV SRR12682296: EBV-1 747,531 → 747,532, EBV-2 142,960 → 142,960.
      Gene-level changes: KP343822.1 65 → 0, KP343821.1 1 → 0, BGLF3/3.5/4
      +0.33 each. EBNA-2 unchanged. Against v1 (885,933): EBV-1 + EBV-2 =
      890,492.
    - HSV-1 SRR8315713: 23,187 → 23,187 (v1 23,170). KP343822.1 19 → 0. HHV-2
      is still `possible_em_bleed` (666:1). HSV-1 s-genes unchanged.
  - **Swap:** the golden identity tests now include `cat42b`
    (`tests/test_virus_identity.py`, `final` kept). CMP-06's "latest
    reference" is now `viral_ref_cat42b/build/` (see the note there).
  - **Gene-ID rename, traced:** the `{accession}_{token}` gene-ID rule
    (`ncbi_fetch.py:350`) came in with 3379b7c (2026-09-27). The 241
    gitignored bundled GTFs were regenerated with it at 2026-09-28 22:46,
    during the a97cc00 final-panel work, after `viral_ref_final` was built.
    So cat42 and cat42b carry prefixed IDs for 124 RefSeq anellovirus genes,
    and `final` does not. Identity is unchanged, because t2g keeps the
    genome accession.
- [ ] `CAT-41` — catalogue display names that contradict their genome
  (found 2026-10-01 from the HPV16 "16,18" label, which is fixed). Still
  open:
  - all 7 Influenza D segments (NC_036615–21) display as "Influenza D virus
    segment 7";
  - HHV-6A (NC_001664.4) displays as "Human herpesvirus 6" (changing it
    changes sibling-note text);
  - XS2 (KC138720.1) displays as "Human papillomavirus 78";
  - smaller species/type mismatches: HPV61 for Alphapapillomavirus 3, HPV2
    for Alphapapillomavirus 4, HPV14 for 14D, HPV68 for 68a, HPV6 for 6b,
    "Human enterovirus 68, 70" for Enterovirus D.
  Re-run the scan in `CAT-41` after any catalogue merge.
- [x] `MECH-A2` — prefix-named rows lack catalogue decisions (residual of
  MECH-A step 4a, 2026-09-30). On a VIRTUS-style index (covid x213) HHV-6B
  resolves as `legacy_prefix`, key `name:Human herpesvirus 6b`, with no
  `sibling_group` or `risk_class`. The same virus is `taxid:32604` on the
  stored SRR20710641 index. So the HHV-6A/6B note can never fire there, and
  `eve_risk` takes the genus fallback. Candidate fix: adopt the catalogue
  row when the prefix name matches exactly one taxid's `common_name`.
  - Done 2026-10-01: a prefix name equal to a catalogue display name (these
    are unique per key) takes that key's row, keeping status
    `legacy_prefix`. The key's eve row is used when it has one, so
    `eve_risk` stays fail-closed. `virus_facts` uses the genus fallback
    only for `name:` keys, because `risk_class` is blank on 681 non-EVE
    rows. Covid x213 parity (job 25691885, against final4a): all 12
    multimap layers, X and `viral_summary.tsv` are identical. The 97
    `HUM_HERP6B_*` genes are now `taxid:32604` / `HHV-6`, and the
    genus-fallback warnings went from 1 to 0.
- [ ] `MECH-B` — virus-level Detection: group, then sum, then threshold.
  `accession_breadth` becomes coverage over reference genes (today it is always
  1.0), and `sensitivity.tsv` gets zero rows.
- [ ] `MECH-C` — a single Run Config writer; delete the Namespace → `k=v` →
  YAML round trip and `createconfig`. This is the root cause of `SW-13`.
- [~] `MECH-D` — (module landed 2026-10-03, see DEF-02) Chemistry module: one geometry for kb, STARsolo and the
  preflight, plus an R1-length / polyT check. The 10x v2 run as `-x 10xv3`
  gave 1.78M "cells" with no error.
- [ ] `MECH-E` — the reference pipeline: source adapters → one gate module →
  masking → GTF emitter → index with a recorded D-list. The prototype is
  `viral_panel_max_2026-09-28/01–04`.
- [ ] `MECH-F` — the corrected-BUS boundary, index-aware denominators, and one
  gene-role catalogue for every family.

### WP1E — Grill decisions for 3.0 (user-confirmed 2026-09-29)

The user answered a structured grill (41 questions, three rounds) and confirmed
the summary on 2026-09-29. These decisions are binding for the rows below and
for WP3–WP11. A `DEF-*` row is the implementation work a decision creates.

| # | Decision | Answer |
|---|---|---|
| Q1 | Release before measured defaults | Internal pre-release only; outputs stamped `defaults_status=provisional` |
| Q2 | Coordination | One orchestrator writes PLAN.md, `.living/`, inventory, claim registry, protocol |
| Q3 | MECH-A owner | This main session (user: "Continue MECH-A here") |
| Q4 | Default panel | Current default stays until MECH-A, MECH-B and the sibling experiment land |
| Q5 | 10x 5′ | In scope, as `-x 10xv2`/`10xv3` plus strand handling |
| Q6 | Chemistry declaration | Auto-detect from reads; fail on ambiguity or a mismatch with `-x` (override flag) |
| Q7 | Host strategies | Combined and STAR two-step both candidates; D2 decides |
| Q8 | Host genes in `--gtf` | Index build manifest of host/viral genes plus the column-5 structural guard; manifest wins; contradictions refused |
| Q8b | Anellovirus unit | Genus; unassigned = "Anelloviridae (genus unassigned)" (done, a566dad) |
| Q9 | Canonical environment | `env_full`, locked (REL-02/03); CI adds a snakemake 9 dry-run; retire py3.8 env for repo code |
| Q10a | Subagent commits | Own `v3/*` branches; only the orchestrator merges; every merge gated |
| Q10b | Concurrent implementers | 3; shared files (`menu.py`, `Snakefile`, `multimap.py`) serialised |
| Q10c | Downloads | Approved: chemistry subsamples, PBMC backgrounds, comparator envs; checksum-pinned manifests |
| Q10d | Tuning compute | 20k core-hours |
| Q11 | Covid anellovirus numbers (CAT-30) | Files frozen; retraction notice "poly-G no-signal reads ≈90 %, host-homologous ≈10 % (F-005, F-019); low-level divergent anellovirus not excluded"; claim retracted in registry |
| Q12 | Whole-genome features (CAT-38/39) | Virus-level detection only; tagged `whole_genome`; excluded from gene-level outputs and programmes; low-complexity filter required |
| R2.0 | Read-artefact filter | Default pre-count filter (homopolymer ≥15 nt, TSO/adapter, low complexity) with audit table, plus reference masking |
| R2.1 | Preregistration | `defaults_selection` in `protocol.yaml` (grid, objective, hard constraint, 1-SE conservative tie-break, sample bootstrap); reviewed, signed off, frozen before tuning |
| R2.2 | Already-seen data | Excluded from tuning and holdout; context only |
| R2.3 | Specificity constraint | Zero probable/strong calls on host-only and planted-homology negatives binds every default |
| R2.4 | Tuning reference | No tuning run until G4 passes; reference chosen after G4 |
| R2.5 | Sibling allocation (F-017) | Within-sibling-group variant (unique support or group EM) in D3; counting contract and conservation claim amended first |
| R2.6 | Strand | `--strand` option; otherwise inferred from a 1M-read pilot; recorded in manifest |
| R2.7 | Per-chemistry defaults | Allowed; live in the Chemistry module; preregistered per chemistry |
| R2.8 | Cell calling under two-step | `--called-cells` if given; else STARsolo host `GeneFull` counts into EmptyDrops (SCI-02) |
| R2.9 | Generator backgrounds | Real PBMC primary plus synthetic GRCh38-transcript control |
| R2.10 | Reruns now | All exploratory reruns now (see `EXPL-R2.10`) |
| R3.1 | Compute split | Shared defaults once on 10x v3 (~55 %); chemistry-specific per chemistry (~25 %); per-chemistry confirmation (~10 %); reserve (~10 %) |
| R3.2 | Decision rules | Per-default objective: recall (host strategy, STAR, read filter); weighted allocation L1 per Index Kind (multimap); cell Jaccard vs truth labels (cell calling, never viral outcomes); SCI-03 F1 (thresholds, tiers). No passing option → most conservative, reported "screening only" |
| R3.3 | Comparators | Full SCI-04 matrix, both index arms; ~10k core-hours, separate budget |
| R3.4 | Holdout blinding | Seeded partition; truth outside repo, owner-only, hash in protocol; only the H1 brief names the path; unblinding needs user's written sign-off after the defaults freeze |
| R3.5 | Releases | `3.0.0.devN` git tags only; `3.0.0rc1` after G2/G5/G6a; `3.0.0` after G6/G7. User sign-off at: protocol freeze (G3), defaults freeze, holdout unblinding, rc1 publication, 3.0.0 tag, any retraction |

Implementation rows:

- [ ] `DEF-00` — amend `analysis/v3_validation/protocol.yaml`:
  - add the `defaults_selection` section (R2.1, R3.1, R3.2);
  - add 5′ to the chemistry scope (Q5);
  - add the within-sibling D3 variant (R2.5);
  - add the per-chemistry defaults (R2.7);
  - add the both-background generator (R2.9);
  - add the blinding layout (R3.4);
  - decide the knee sensitivity arm: `sensitivity_only_callers` requires
    `knee`, whose estimator is defective (SW-23, 2026-10-01). Fix it or
    replace it in the amendment.

  Then an independent review, user sign-off, and the freeze (G3).
- [ ] `DEF-01` — read-artefact filter before `kb count`, with an audit table
  (R2.0, F-019). Reference homopolymer/low-complexity masking stays under
  CAT-17.
- [~] `DEF-02` — **priority after F-020**, because every 5′ run so far used
  about 7–9 % of its reads. The Chemistry module (`MECH-D`), covering:
  - auto-detection and fail-closed checks (Q6);
  - `--strand` and the pilot inference (R2.6, strand test job 25672276);
  - per-chemistry defaults, including SW-21 barcode correction (R2.7).
  - [x] `--strand {forward,reverse,unstranded}` plumbing: 2026-09-30, branch
    `v3/strand` (d8bf730). The option reaches `kb count` through
    `params.strand`; unset omits the flag, so kb keeps its per-technology
    default. The default stays unset until DEF-01, because unstranded
    re-admits the forward-only poly-G artefact (F-019).
    - `run_manifest.json` records `options.strand` only when it is set. So an
      old manifest still resumes when `--strand` is not given. An old manifest
      with `--strand X`, or a different explicit strand, refuses with a
      reason. Tiny e2e asserts that `kb_info.json` records `--strand
      unstranded`.
    - `reverse` fails on the 2-read fixture (kallisto exit 1), which is why
      the e2e uses `unstranded`.
  - [ ] Strand pilot inference and the per-chemistry default (the rest of
    DEF-02, MECH-D). Next: rerun the 5′ positives with `--strand reverse` and
    `unstranded`. These reruns are exploratory; they do not set the default.
    - **Reruns done (2026-10-03; full runs, emptyDrops, code d01043c;
      results in F-020's update).**
      - Covid x213/x216, reverse vs unstranded: 47.8/53.6 % and 46.2/54.4 %
        pseudoaligned, against 6.4/8.9 % with no strand.
      - HHV-6B: 2,405 molecules reverse vs 5,596 unstranded.
      - **Correction:** the HHV-6B baseline already ran unstranded in effect
        (60.2 %; kallisto 0.52.0, `-x 10xv2`, no strand flag). The covid
        baseline behaved as forward (kallisto 0.51.1, `-x 10xv3`, no flag).
        Whether version or chemistry sets the default is not separated. So
        the pilot must pass all three strands explicitly.
      - **New, unvalidated:** under reverse, covid x213 puts 55k + 25k
        molecules on two TTMDV placeholders (AB303552.1, AB303557.1), in
        64.9 % of called cells. These genomes have no homopolymer over 12 nt.
        Needs a read-level check before any claim. Logged, not started (no
        new tasks in this pass).
      - `--strand auto` check: reverse/unstranded = 0.85–0.89 on all three
        5′ libraries, forward/unstranded = 0.11–0.16. So τ = 0.8 picks
        reverse. For HHV-6B, reverse is the choice that loses 57 % of viral
        molecules.
    - **`--strand auto` landed (2026-10-03, 8ff6bec, opt-in; the default
      stays unset).**
      - `src/viralscan/strand.py`: the pilot runs on the first 1M pairs, with
        three **explicit** strands, and reads `p_pseudoaligned`.
      - Rule `infer_strand`, τ = 0.8: the larger qualifying ratio wins; a tie
        or unstranded ≤ 0 gives unstranded.
      - The manifest keeps `options.strand="auto"` and adds a per-sample
        `strand_inference` block (choice, rates, pilot_reads, tau) outside the
        fingerprint. Resume reuses the recorded choice; auto vs an explicit
        strand refuses.
      - Tests: `tests/test_strand_auto.py`. The three measured F-020 rate sets
        pick reverse; the forward and unstranded branches are synthetic only.
      - Still open: a measured 3′ ratio (done 2026-10-03, see (d) below) and
        a real-`kb` end-to-end pilot run.
  - **Chemistry module, 2026-10-03 (user decisions the same day).**
    - Decisions: SW-21 → bypass kb's barcode correction for a technology with
      no official on-list, so ViralScan's cell calling is the only caller.
      Drop-seq is inferred when no 10x on-list matches and R1 is 20 bp.
      Per-chemistry default *values* stay with DEF-00 (R2.7 preregistration);
      today's behaviour (strand unset) ships as provisional.
    - [x] (a) `src/viralscan/chemistry.py`: one table (`CHEMISTRIES`, kb name,
      CB/UMI length, ngs_tools on-list); `evidence.cb_umi_geometry` re-exports
      it (the unused `10xv3_5p` entry is gone). `detect()` samples 100k R1
      reads, streams each on-list once, and reads the UMI length from the
      TSO/poly-T start (26/28) or a trimmed R1 length. `resolve()` fails on
      ambiguity, on samples that disagree, or on a mismatch with `-x`;
      `--force-technology` keeps an explicit `-x`. Tests:
      `tests/test_chemistry.py`.
      Real libraries (19 s for all four): EBV SRR12682296 → 10xv2 (97.3 % v2
      list); HHV-6B SRR20710641 → 10xv2, 5′ (88.2 %, TSO at 26); covid x213
      → refused without `-w` (best bundled list v4 4.4 %), 10xv3 geometry
      with Cell Ranger's list (67.2 %); HSV-1 SRR8315713 → dropseq. EBV run
      as `-x 10xv3` is refused (the MECH-D failure). No bundled list covers
      GEM-X 5′. No local 10xv3 3′ library exists to check v3/v4 overlap;
      two lists above 0.5 fail closed. Bounded from the lists instead: v3 ∩ v4
      = 68,254 barcodes (1.0 % of v3), v2 ∩ v3 = 77,142 (10.5 % of v2), so a
      10xv3 library cannot clear 0.5 on two lists.
      Fix (same day): Drop-seq with its own `-w` list resolves to dropseq
      (a 20 bp R1 gives no UMI length, so it used to be refused).
    - [x] (b) run preflight, 2026-10-03: `-x` defaults to detection,
      `--force-technology` keeps an explicit `-x` (excluded from the
      fingerprint). Per-sample `chemistry_detection` blocks go into
      `run_manifest.json` outside the fingerprint (`record_manifest_block`,
      shared with `strand_inference`). The advisory `_whitelist_preflight`
      is gone; `check-whitelist` stays. Resume: the resolved `-x` is what is
      fingerprinted, so an old run made with the old `10xv3` default resumes
      when its reads are 10xv3, and refuses (with the reason) when they are
      not. Docs: `docs/faq.md`, `docs/cli_reference.md`. Tests:
      `TestResolveChemistry`, `test_force_technology_does_not_change_the_fingerprint`.
      Note: slurm scripts that run HHV-6B as `-x 10xv3` now stop; it is 10xv2.
    - [x] (c) SW-21 bypass, 2026-10-03: `chemistry.kb_whitelist_arg` gives
      `kb count -w None` for a chemistry with no on-list and no `-w`. kb 0.29.5
      then skips `bustools correct` and writes no `output.unfiltered.bus`, so
      multimap starts from the raw BUS (`correction="none"`, already logged).
      Tests: `TestKbWhitelistArg`, a Snakemake dry-run render. Not yet rerun
      on HSV-1 SRR8315713; expect ~1.8M barcodes and ~32k HSV-1 molecules.
    - [x] (d) measured 3′ ratio, 2026-10-03 (job 25695946,
      `scripts/slurm_strand_pilot_3p.sh`): EBV SRR12682296, 10xv2 3′, cat42b
      host+viral index, first 1M pairs, real `kb`. Pseudoaligned: forward
      72.6 %, reverse 2.8 %, unstranded 73.8 %. forward/unstranded = 0.98,
      so τ = 0.8 picks **forward**, as a 3′ library should. The forward
      branch is no longer synthetic-only.
    - [~] End-to-end check of (b)+(c): HSV-1 SRR8315713 run with **no `-x`**
      (job 25695959, `scripts/slurm_quant_sw21_hsv1.sh`, output
      `runs/combined_sw21/`). Pass: detection → dropseq, kb log shows
      `-w None`, and barcodes/HSV-1 molecules return to ~1.8M / ~32k (vs
      5,468 / 23,170 with the knee allowlist, job 25666447).
- [x] `DEF-03` — index build manifest of host and viral gene sets, overriding
  `--gtf` (Q8). Fix `reference_strategy.py:704` so it passes a viral-only
  GTF. Folds into MECH-A step 4.
  - 2026-09-30, branch `v3/def-03-build-manifest`. `build-ref` and
    `--reference` write `<index>.build_manifest.json` (sorted host and viral
    gene IDs, plus provenance sha256) after a successful `kb ref`. The
    `analysis` rule passes it to `build_identity_table(build_manifest=...)`.
  - Contradiction rule (user-confirmed 2026-09-30): compare de-versioned IDs,
    restricted to the index t2g.
    The viral set the no-manifest resolution would produce must equal the
    manifest's viral set, and every index gene must sit in one manifest set.
    Otherwise raise `BuildManifestContradiction` (a `ValueError`) in the
    `analysis` rule. `--gtf` is a union of panels, which is why the rule
    compares the resolved set and not the raw GTF set.
  - No manifest (every existing index: v2, final, max) → a logged warning,
    then the old path, byte-identical. `reference_strategy.py:704` now
    passes the all-virus GTF for both strategies.
- [ ] `DEF-04` — `whole_genome` tag in the catalogue and its exclusions (Q12).
- [ ] `DEF-05` — two-step cell calling: STARsolo `GeneFull` into EmptyDrops,
  `--called-cells` override, and fail closed when neither exists (R2.8,
  under MECH-F).
- [ ] `DEF-06` — counting-contract amendment and property tests for the
  within-sibling allocation (R2.5). Must land before the D3 grid is frozen.
- [ ] `DEF-07` — CAT-30 retraction notice beside the covid results, and the
  claim marked retracted in `claims/registry.json` (Q11). Needs the user's
  sign-off on the wording (R3.5).
- [ ] `DEF-08` — lock `env_full`, and add a snakemake 9 dry-run to CI (Q9).
- [ ] `CMP-06` (WP6B) — same-input comparison with evonk's original
  ViralScan 2.2.0 on the **same latest reference** as v3, so only the
  implementation differs (user request 2026-09-30).
- [ ] `DEF-09` — release tagging: `3.0.0.devN` after G1, with outputs stamped
  `defaults_status=provisional` in `run_manifest.json` (Q1, R3.5).

### WP1C — Simplification pass (ponytail audit, new 2026-09-27)

A read-only audit on 2026-09-27 ("ponytail") covered all 58 k tracked Python
lines. The user chose groups 1–3.

- **Kept out of scope, as provenance of registered numbers:** the
  `analysis/multimap_profiling` scripts (`fast_profile*.py`,
  `interpret_cprofile.py`).
- **Not taken:** `extras/build_anello_table.py`. It is the only regenerator of
  a frozen TSV.
- **Rejected:** replacing stdlib `urllib` with `requests`. That runs the wrong
  way on the stdlib-first rule.

- [x] `SIMP-01` — delete verified-dead code:
  - `multimap.normalize_barcodes`: no caller, and it references an undefined
    `output`.
  - `host_filter._kallisto_filter`: unreachable.
  - `sensitivity.capture_reference`: an alias with no caller.
  - `menu._config_value`: only its test calls it.
  - `KbCountOutputs.bus_txt`: the legacy path, only its test reads it.

  Kept: `evidence.blast_identity`. Its only caller is the BLAST integration
  test that `BASE-04` cites, and that test cannot run here.
- [x] `SIMP-02` — collapse duplicated helpers in `src/` to one copy each:
  - `sha256` of a file: `src/` keeps `run_safety.sha256_file`.
  - the packaged-TSV path resolver: 4 copies.
  - the gzip-aware text opener: 3 copies.
  - the TSV writer in `evidence_run.py`: 6 copies.
  - the `--verbose`/`--quiet` argparse flags in `menu.py`: 7 copies.
  - `detection._sum_axis0/1`.

  **The `scripts/` copies of `sha256`, the TSV readers/writers and the atomic
  JSON writers stay duplicated on purpose.** The legacy v2/v3 and fresh-control
  tools copy `benchmark_v3_multimap.py`, `compare_legacy_v2_v3.py` and
  `run_fresh_control.py` into hash-frozen packets and run those copies alone
  (`tests/test_benchmark_v3_multimap.py`, `tests/test_freeze_fresh_control_packet.py`).
  An import of `governance_utils` would break a frozen copy. The audit missed
  this.
- [x] `SIMP-03` — dependencies:
  - Drop `pyfiglet`. `menu.py` already falls back to plain text.
  - Drop `seaborn`. Its 2 plots move to matplotlib.
  - Swap `enrichment._bh_adjust` for `scipy.stats.false_discovery_control`,
    with `scipy>=1.11` pinned and an equivalence check against the old loop.
  - Apply the change in `pyproject.toml`, `environment.yml` (inventoried, so it
    needs a re-pin), `conda-recipe/meta.yaml` and CI.

## WP2 — Install, lock, and artifact parity

Objective: make the pip tier honest and make conda, OCI, and Apptainer execute
the same tested build. Estimated remaining effort: 5-8 engineering days plus
runner time.

### WP2A — Pip and locked environment

- [ ] `REL-01` — move full-workflow-only dependencies such as Snakemake out of
  mandatory pip runtime dependencies; make `doctor --profile pip` check only
  Python/API/reporting/validation capabilities.
- [ ] `REL-02` — define Python 3.11 on `linux-64` as the canonical full-workflow
  toolchain; advertise other platforms only after the same workflow passes.
- [ ] `REL-03` — generate and commit reproducible runtime/development lockfiles
  and an exact external-tool version manifest.
- [~] `REL-04` — build wheel and sdist once, run `twine check`, inspect packaged
  assets, install each in a clean environment, and write `SHA256SUMS`.
  - 2026-10-02: built once from a `git archive` export with the codex
    env (`build --no-isolation --skip-dependency-check`, no network).
    - Outputs: wheel 587,519 B and sdist 663,361 B. `twine check` passes
      both.
    - Both install offline into a `--system-site-packages` venv and report
      `viralscan 3.0.0.dev0`.
    - `SHA256SUMS`: whl `9061229c…`, sdist `d718456a…`; `sha256sum -c` is OK.
    - The first `check_ship_scope` run failed: 9 members were missing from
      `config/public_ship_scope.json` (`virus_identity.py`, `sensitivity.py`,
      `index_exclusions.tsv`, `run_complete.schema.json`). They were added in
      this merge; the re-check is recorded below.
    - Still open: clean installs from an index (needs the network; REL-05)
      and the `release.yml` steps (REL-13).
  - Re-check after the allowlist fix (2026-10-02, HEAD a968a63 export):
    the build succeeds, `twine check` passes on both artifacts,
    "ship-scope validation passed", and `sha256sum -c` is OK.
- [ ] `REL-05` — run `doctor --profile pip` and `validate-run` against a packaged
  v3 fixture from both clean installations.

### WP2B — Containers and parity

- [!] `REL-06` — Miniforge is digest-pinned; rebuild Docker from the committed
  lock and the exact tested wheel rather than from a mutable source install.
  - Blocked (2026-10-02): needs the committed lock (`REL-03`, `[ ]`) and the tested wheel (`REL-04`). No docker daemon here, and verification needs a CI build job.
- [ ] `REL-07` — publish a versioned OCI candidate, record its digest and tool
  manifest, and never assign `latest` to an RC.
- [ ] `REL-08` — build Apptainer/SIF from that OCI digest rather than performing
  a separate dependency solve; record the SIF SHA-256.
- [ ] `REL-09` — run the identical tiny workflow in locked conda, Docker, and
  Apptainer and compare validated molecule counts and deterministic output hashes.
- [ ] `REL-10` — save `analysis/release_parity/parity_report.json` with commands,
  versions, digests, hashes, and explained nondeterministic files.

### WP2C — Supply chain and release workflow

- [~] `REL-11` — security CI exists; audit the locked product environment rather
  than the scanner job, add OCI scanning, a full SBOM, dependency/data licence
  report, and reviewed vulnerability exceptions.
- [ ] `REL-12` — pin GitHub Actions by commit SHA and generate provenance/
  attestations for wheel, sdist, OCI, SIF, locks, and checksums.
- [ ] `REL-13` — make release publication depend on green CI for the exact tagged
  SHA; build distributions once and make all downstream jobs consume them.
- [ ] `REL-14` — add protected TestPyPI/PyPI environments, stable approval, GitHub
  release assets, and separate RC/stable container tag behavior.
- [!] `REL-15` — maintainers must provide final author, licence-holder,
  maintainer/contact, ORCID, affiliation, support-window, and trusted-publisher
  metadata.

- [ ] `REL-16` — **a standalone `kallisto` on `PATH` segfaults on a `kb ref`-built
  index, and `viralscan evidence` picks it up.** Found 2026-10-03 while running
  `viralscan evidence` in `conda/envs/viralscan_bench`: exit `-11` (SIGSEGV) on
  the cat42b `panel.idx`. The two binaries that both report **version 0.52.0**:

  | Binary | Size | md5 | On a `kb ref` index |
  |---|---|---|---|
  | `$ENV/bin/kallisto` (conda) | 2,795,304 B | `c9cb3b19f7dbf09e85b70b45e1d50c97` | SIGSEGV |
  | `site-packages/kb_python/bins/linux/kallisto/kallisto` | 8,518,336 B | `695d3d418a11263adbc451dd63ea1233` | works |

  The version string cannot be used to tell them apart, so a preflight
  `shutil.which('kallisto')` check passes and the failure surfaces as a crash
  mid-run. Related and already recorded: the mirror case where a *hand-built*
  index fed to `kb count` makes the bundled kallisto spin at 787 % CPU for
  minutes on one read pair (`.living/learnings.md`, 2026-10-03) — same
  version-lock, opposite direction, and neither fails loudly.
  - Fix: resolve kallisto through `kb_python`'s bundled path (or an explicit
    `--kallisto`), not `PATH`, wherever ViralScan shells out to it directly;
    record the resolved binary's md5 in the run manifest so a parity run can
    detect the swap.
  - Until then: run `viralscan evidence` with kb's `bins/linux/kallisto` first on
    `PATH`.

`G2` passes when a clean lock-created environment plus Docker and Apptainer run
the same tiny workflow with identical validated counts; installable schemas and
assets are present; security/licence reports contain no unresolved actionable
finding; and the exact-SHA release workflow is dry-run verified without
publishing.

## WP3 — Preregister scientific validation

Objective: freeze outcome-independent decisions before producing new v3
scientific results. Estimated effort: 1-2 days plus reviewer sign-off.

- [x] `SCI-01` — create schema-valid `analysis/v3_validation/protocol.yaml` and
  `README.md` with hypotheses, supported scope, datasets, factors, seeds,
  chemistry, input/reference hashes, and exclusions.
- [x] `SCI-02` — freeze the barcode universe/cell-calling rule, feature
  intersection, primary unique-only cross-tool layer, secondary ambiguity-aware
  layers, denominators, and failure-reporting policy.
  Frozen subsection SHA-256:
  `699854b71169222e74d26c2119f1c9d7742d5962ac4ed8c02dfa9f288e3339dc`;
  dependent hypotheses/endpoints/exclusions digest:
  `fd818661a8d06f28c5c23e78bfb1c5b48b151023f5f749453260d023adbff333`.
- [~] `SCI-03` — freeze training/holdout partitions, evidence-tier calibration
  rules, metrics, limit-of-detection method, and biological-sample bootstrap unit.
  All 14 round-1 specification findings against these sections are resolved:
  largest-remainder apportionment replaces ceiling rounding, the sort key is
  pinned, an empty eligible grid has a declared failure outcome, the tie-breaker
  names its standard error, the bootstrap refuses degenerate small samples, and
  two denominators no longer let failed rows flatter a result. Both sections
  remain **`pending`**: `F1` and `F2` need `VAL-01` and `REF-08` factor levels,
  and dormant validator rules now refuse a freeze without them.
  - Blocked (2026-10-02): F1/F2 factor levels need `VAL-01`/`REF-08`, and the freeze follows `DEF-00`. A review round before `DEF-00` would be wasted.

- [~] `SCI-04` — enumerate every ViralScan, STARsolo, traditional alignment,
  Venus, Viral-Track, and VIRTUS row with exact environment/reference
  requirements. Twelve workflows over 78 rows. Both fairness defects are fixed:
  the only exact-truth dataset now reaches all three dedicated comparators, and
  each runs a native-published and a matched-accession-index arm so a tool
  difference is separable from a reference difference. Breadth and tuning
  asymmetries are disclosed rather than removed. Sections remain **`pending`**
  until `SCI-05` round 2 passes; `tool_environments` still blocks on `REL-03`.
  - Blocked (2026-10-02): `tool_environments` needs `REL-03`. Freeze after a clean `SCI-05` round.

- [~] `SCI-05` — obtain an independent protocol review, resolve findings, then
  record the protocol SHA-256 and Git SHA before outcome-generating runs. Three
  rounds are complete and recorded, each verifying the previous round's fixes by
  tamper experiment rather than by reading resolution notes. Round 1: 24 findings.
  Round 2: found two round-1 fixes cosmetic and reproduced a hole in the amendment
  rail. Round 3: confirmed those fixes real, then found the ledger itself was not
  append-only in fact. Rounds 4-9 progressively hardened the digest scope, the
  hash chain, and the git anchoring, each round finding second-order defects in
  the previous round's fix. Rounds 10 and 11 returned no blocker, and two
  reviewers independently judged the specification buildable. All rounds are
  dispositioned and every reproduced tamper now fails closed. `verdict` stays
  `does-not-pass`: no round has yet passed clean with zero open Majors, which is
  the bar for flipping it. `R11-F2` — the last open item, the ledger checker's
  fail-open branches — is closed as of 2026-07-27.

  Standing limitation, established by the 2026-07-27 code review: every one of
  these rounds read `protocol.yaml`. None could have caught a frozen value that
  no call site consumes, because that is invisible from the protocol side and
  from the output side alike. Protocol review does not substitute for verifying
  that the code honours the contract.

`G3` passes when the protocol validates against its schema, has no unresolved
review finding, is hashed, and outcome-generating jobs have not preceded its
freeze commit.

## WP4 — Freeze production references

Objective: make reference contents reproducible and calibrate host-homology
safeguards without holdout leakage. Estimated effort: 3-5 engineering days plus
about 8 cluster hours per full GRCh38 build.

### WP4A — Profiles and provenance
  - Blocked (2026-10-02): run the next review round after the `DEF-00` amendment lands.

- [x] `REF-01` — profile names and expanded anellovirus behavior exist; the
  expanded panel is now the **default** rather than opt-in, and its gene IDs
  reach detection. Accession lists are frozen for `anellovirus-representative`
  and `anellovirus-expanded` via `src/viralscan/data/anellovirus_accessions.tsv`
  (2,042 accessions, 2,042 unique, 0 duplicates). `curated` and
  `broad-discovery` still carry **no** frozen accession list and remain
  label-only — see the new `REF-13`.
  - **2026-09-26 evidence for the default flip.** Measured 31-mer coverage of
    the bundled 20-genome RefSeq TTV panel against the 2,042 real human
    anellovirus genomes (`anellovirus.fa`): median **0.00 %**, and **85.8 %** of
    genomes share *zero* 31-mers with the panel. Positive control: EBV
    NC_007605.1 against itself = 100.0 % over 144,283 31-mers, so the method is
    sound. Per genus, zero-coverage genomes: Betatorquevirus 98.4 % (n=1,542),
    Gammatorquevirus/Samektorquevirus/Hetorquevirus/Gyrovirus/Memtorquevirus
    100 %. Leave-one-out capture for a *novel* strain: bundled panel 0.04 %
    (P(90 bp fragment captured)=0.024) vs expanded panel 20.5 % (P=1.0000).
    Cost of expanding: 5,994,773 bp total, 4,888,291 distinct 31-mers, only
    **1.21x** k-mer space inflation, 91.1 % of 31-mers unique to one genome,
    median genome 20.5 % redundant. CD-HIT-style downsizing does not help:
    20 genomes retain 1.1 % of the k-mer space, 204 retain 11.3 %, 2,042 retain
    100 %. Conclusion: keep all 2,042.
- [~] `REF-02` — fail-closed fetches and manifests exist; complete accession
  version, taxonomy, snapshot, retrieval date, SHA-256, length, licence, cluster,
  representative status, rationale, and missing-accession fields.
  - Blocked (2026-10-02): the licence fields need `REF-05`. Wiring the manifest into the production builder is `REF-03` (`[ ]`).
- [ ] `REF-03` — apply identical masking, duplicate-ID/sequence validation, and
  manifest generation to dedicated and combined build paths.
- [ ] `REF-04` — make frozen inputs rebuild byte-identical panel FASTA/GTF/t2g
  contents and save a reproducibility audit.
- [ ] `REF-05` — replace vague source-data licence text with reviewed terms for
  every redistributed or fetched reference source.
- [!] `REF-11` — publish the viral annotation panel archive and register its
  Zenodo DOI. `src/viralscan/data_fetch.py` pins
  `VIRAL_DATA_DOI = "10.5281/zenodo.20112332"`, but that identifier is **not
  registered**: `https://zenodo.org/api/records/20112332` returns
  `{"status": 404, "message": "The persistent identifier is not registered."}`
  and `https://doi.org/10.5281/zenodo.20112332` also returns 404, while an
  unrelated third-party DOI referenced elsewhere in the repo resolves normally.
  The 195 GTFs remain in `src/viralscan/data/` in the source tree but are absent
  from the installed package, so `viralscan data fetch` — and therefore every
  bundled-panel run from a clean install — fails for all users. Blocks `REL-05`,
  `DOC-05`, and `G2`; currently blocking the `GOV-06` attempt-3 cache pin.
  Owner action: publish the archive and register the DOI, or correct the pinned
  record identifier.

### WP4B — Full-genome competition and anellovirus

- [~] `REF-06` — the public `--genome-dlist` path and raw annotations exist;
  build the production curated human-plus-virus index with a checksum-pinned
  GRCh38 D-list.
  - **2026-10-03: grill approved, build and measure.**
    - Finding: the shipped cat42b has no genome D-list. kb's default D-list
      is the positional `combined.fa`, which is 465,769 host cDNA records
      plus the 2,343 viral genomes; that gives 632,261 D-list k-mers.
    - Literature: Sullivan et al., NAR 2025 (PMID 39657125) for the
      distinguishing flanking k-mer mechanism; Luebbert et al., Nat
      Biotechnol 2026 (PMID 40263451) use a host genome + transcriptome
      D-list, but for translated search.
    - Biomni review stalled with no output, so an internal review was used.
  - **Design (single variable):** `viral_ref_cat42d/` is cat42b's inputs
    and pinned code d01043c, with D-list = cat42b `combined.fa` + GRCh38-2024-A
    `genome.fa`, sha256 pinned in `dlist.sha256`.
    - Gates are one-sided: HPV16 may lose at most 1 %, EBV/HSV-1 at most
      5 %, HHV-2 must still be flagged as bleed, and per-gene losses are
      checked. Any increase is investigated.
    - Swap rule (default): cat42d replaces cat42b if every gate passes.
  - **Jobs:**
    - build 25695872;
    - regressions 25695873 (GSE189670) and 25695874 (EBV/HSV-1), both
      afterok;
    - covid x213 `--strand reverse` on cat42b (25695875) vs cat42d
      (25695876).
  - **Side check, tiny HPV16 index (kallisto 0.50.1):**
    - 90-nt reads crossing a CDS end into unmodelled sequence: 0/21
      aligned with the whole-viral-genome D-list (11 D-list k-mers), 1/21
      without it. Fully-inside controls: 18/21 both ways.
    - So these reads are lost by pseudoalignment itself, not by the D-list.
      Removing viral sequence from the D-list is not needed.
    - This supports F-022: gene models must reach the mRNA ends.
- [ ] `REF-07` — save the reference manifest, host-homology/low-complexity table,
  index/t2g/GTF, commands, versions, checksums, and build resource accounting.
- [ ] `REF-08` — calibrate homology/complexity exclusion thresholds using only
  preregistered training controls and freeze `thresholds.json`.
- [ ] `REF-09` — prove planted human-homology reads cannot reach probable/strong
  evidence on holdout; retain raw measurements and all excluded calls.
- [!] `REF-10` — an orthogonally confirmed anellovirus-positive sample is needed
  for real sensitivity claims; without it, ship screening support only and label
  every result accordingly.
- [~] `REF-13` — **new, 2026-09-26.** Detection-side reference visibility and
  name resolution. Two defects made most of the panel inert while the runs
  still looked clean:
  - `SENS-01` (**fixed**) `scripts/analysis.py` globbed only the packaged panel
    directory, but the expanded anellovirus GTFs are materialized into the
    *built index* by `build-reference`, so 2,022 of 2,042 genomes (91 %) were
    countable and never reportable. `anellovirus.candidate_gene_ids()` now
    derives the `{accession}_geneN` IDs the builder emits, and
    `--anellovirus-gene-ids` (default on) controls it. Over-inclusion is the
    safe direction: `detect_genes` only reports IDs that are real columns.
  - `SENS-02` (**fixed**) 49.2 % of gene IDs in the panel the covid runs
    actually used (`references/starsolo/.../viral_genome.gtf`, 4,650 genes)
    resolved to no virus name, so each became its own row in
    `viral_summary.tsv` — the covid run published 9 of 17 rows as bare gene IDs
    (`HHV1gp00p39`, `CeHV2gUL24`, `MPXV_gp132`), which silently broke
    `accession_breadth` (1.0 by construction), sibling cross-mapping and
    `eve_risk` for exactly the herpesvirus calls. Fixed by adding the
    underscore-delimited tokens to `VIRUS_NAME_MAP` and a new
    `VIRUS_GENE_ID_ALIASES` tier for concatenated schemes. 49.2 % -> 13.5 %,
    **0 regressions** across 7,443 real gene IDs. The boundary rule was *not*
    weakened (the `AICHIX`/`BORF1`/`BUNYAMW` guards are now regression tests).
  - Still open under this row: `curated`/`broad-discovery` accession lists;
    per-GTF SHA-256 + retrieval dates for the 195 packaged GTFs (no manifest
    row exists for any of them; *2026-10-03:* `scripts/write_gtf_manifest.py
    [-o OUT] [--check OUT]` now writes and checks a per-GTF sha256/size
    manifest. It hashes 324 GTFs on this machine. `source_date` is always
    `unknown`, because retrieval dates are not recoverable. No manifest is
    committed, since the GTFs are gitignored, and the check is not wired into
    CI); the 6 genes with neither CDS nor exon; the
    duplicated `NC_002076.2`; the malformed astrovirus feature column; and the
    2,520/2,692 gene IDs (93.6 %) with no `exon` record, which STARsolo
    comparators cannot count at all.

`G4` passes when panel contents reproduce byte-for-byte, manifests validate,
duplicates are absent, GRCh38 competition is operational, and holdout host-
homology negatives cannot become probable/strong calls.

### WP4B2 — Host-filtering design and reporting (new 2026-09-26)

- [x] `HOST-01` — pin every STAR parameter the host filter depends on
  (`STAR_FILTER_ARGS` in `scripts/host_filter.py`) and record each in
  `host_filter_audit.tsv`. The command previously set **no** alignment or filter
  options, inheriting whatever the installed STAR defaulted to. Two defaults
  were wrong for viral subtraction: `outFilterMismatchNmax 0` rejects any read
  with one host mismatch, so paralogues and allele variants escaped as
  "unmapped" and reached the viral index; `outFilterMultimapNmax 1` reports a
  multi-mapping read as unmapped, which is the dominant false-positive route
  for host repeats and EVEs. Now `4` / `20`, with `outFilterMatchNminOverLread`
  lowered `0.66` -> `0.9`. Closes the parameter half of G8 step 4; the
  per-fragment removal-reason half remains open (see `HOST-02`).
- [x] `HOST-02` — record `pct_retained` in the audit so a re-run that retains a
  different fraction is attributable to a parameter change. Per-fragment removal
  reasons are still **not** recoverable: `--outSAMtype None` discards the SAM and
  `fragment_lineage.tsv.gz` logs retained reads only.
- [x] `HOST-03` — add `pct_infected_comparable`, a strategy-independent
  denominator (absolute 200-molecule host-UMI floor intersected with the called
  set). `pct_infected_called` is **not** comparable across host-filtering
  strategies: on one covid PBMC sample, called cells fell 143,243 -> 28,921 while
  Alphatorquevirus UMI fell 1,167,103 -> 57,715, so `pct_infected_called` *rose*
  from 56.64 % to 62.89 % and inverted the comparison. `pct_infected_called`
  stays the within-run primary.
- [x] `HOST-04` — document the measured three-way host-control comparison in
  `docs/cli_reference.md` and **retract** the "~4x more sensitive" two-step claim
  in `BENCHMARK_COMPARISON.md`. That figure came from a 1M-read subsample of an
  implementation that no longer exists (kallisto CB-UMI-wide subtraction, since
  removed); at full depth on the same sample the arms differ by **3 %**
  (1,479,894 vs 1,434,619 UMI), and the original two-step arm was `blocked`, not
  completed. `--genome-dlist` is now documented as the *weakest* option: it cut
  `p_unique` 2.1 % -> 0.6 % (3.5x) to remove only 14 % of the anellovirus
  artifact, because a k-mer D-list cannot see diverged host sequence.
- [ ] `HOST-05` — **open.** The sensitivity cost of the v3 STARsolo path is still
  unmeasured. It blocks `MS-02`/`CMP-01`–`CMP-03` and the preregistered
  `D15`/`D16` endpoints; no number should be quoted for it until then.

## WP4C — Detection sensitivity and negative-result claims (new 2026-09-26)

Objective: make a negative result self-describing. Added after an audit that
found `limit_of_detection` present only as an unrun endpoint in
`schemas/v3/validation_protocol.schema.json` and in no code path, while
`min_counts`/`min_genes` gated only the UMAP and never detection. A run that
detected nothing carried no depth caveat, so "nothing there" and "did not look
hard enough" were indistinguishable from the output.

- [x] `SENS-01` — `src/viralscan/sensitivity.py`: Poisson detection
  probability, depth-only LOD95, and the exact-match k-mer capture curve
  `P = 1 - (1 - (1-d)^31)^(L-30)`. Depth is the sum of the count matrix, not
  raw reads, because only quantified molecules can be detected.
- [x] `SENS-02` — every run writes `results/sensitivity.tsv` (LOD95, band,
  `depth_sufficient`, `capture_measured`, `informative_negative`) and states the
  limit in `summary.txt` and `report.html`. `viral_summary.tsv` gains
  `pct_infected_comparable`.
- [x] `SENS-03` — `informative_negative` requires depth **and** a *measured*
  capture term, so it is `false` on almost every run by construction. This is the
  load-bearing design choice: depth was ample in every real run here (LOD95
  0.0003–0.0056 per 10k host UMI, all `informative`), while the covid samples
  called SARS-CoV-2 = 0 at 21.6 M quantified molecules. Reference capture — not
  depth — is the binding limit, and capture is unmeasurable without a control.
- [x] `SENS-04` — positive control: `--positive-control-gene` +
  `--positive-control-molecules` (required together), `--require-positive-control`
  to fail closed, `results/positive_control.json`, and bisection inversion of the
  capture curve to an implied divergence. `failed` (control invisible),
  `over-recovered` (not spike-in-specific) and `gene-not-in-reference` are all
  distinct, reported states.
- [x] `SENS-05` — LOD semantics calibrated against molecule-level downsampling
  of the bundled EBV LCL run (103,145,071 molecules, 1,636,934 EBV): P(detect)
  stayed 1.0000 down to 1,270 downsampled reads and first reached 0 at 127,
  i.e. the observed floor is the Poisson floor. `fragment_capture` cross-checked
  at 90 bp and 150 bp for 5–30 % divergence.
- [~] `SENS-06` — **open.** `SENS-01`–`SENS-05` make the limit *reportable*;
  they do not make a negative *certifiable* in practice, because no shipped
  workflow plants a control. Closing this needs either a spike-in recipe in
  `docs/vignettes/` or integration with `VAL-01`'s generator, after which
  `E8`/`D17` LOD95 can be estimated by the preregistered probit fit rather than
  reported as an analytic floor.

## WP4D — Gene-programme inference, layer 2 (new 2026-09-26)

Objective: for viruses layer 1 detected, distinguish latent from productive
expression per cell. Added after measuring that a per-gene comparison cannot do
this at all on this data.
  - Blocked (2026-10-02): needs a route decision (spike-in vignette vs `VAL-01`/`ANDET-08`). The probit fit needs VAL data.

- [x] `PROG-01` — catalogue generator `extras/build_gene_programs.py`. Hand-curated
  biology joined programmatically to bundled-panel attributes, with
  **overlap groups computed by exonic interval intersection** rather than
  hand-assigned. Two bugs were caught by generating rather than assuming: a
  monotonic sweep chained all 96 EBV genes into one group, and bounding boxes
  put LMP-2A (whose exons sit at both genome ends because LMP-2 is spliced
  across the origin and the genome carries terminal repeats) in a group with
  everything. EBV now resolves to 14 groups, largest 6.
- [x] `PROG-02` — ship `src/viralscan/data/gene_programs.tsv` (79 rows, 9
  viruses) and register it in `pyproject.toml` package-data,
  `config/public_ship_scope.json` (wheel + sdist) and `MANIFEST.in`.
  `include-package-data = false` means MANIFEST alone would not ship it.
- [x] `PROG-03` — `src/viralscan/gene_programs.py`. Evidence is the
  `counts_unique_viral` layer; breadth counts distinct non-overlapping overlap
  groups; the multimap-allocated breadth is reported alongside, never merged.
  `latent` and `mixed` are unreachable when `latency_observable_in_rna=false`.
- [x] `PROG-04` — integration: optional `gene_programs` Snakemake rule gated on
  `config["gene_programs"]`, depending on `log/detection.done` **and**
  `results/viral_summary.tsv` so it cannot run before layer 1; sentinel added to
  `rule all` only when enabled; `--gene-programs` / `--programme-min-breadth`
  CLI; `viralscan rerun-programs` operating in place (layer 2 changes no counts,
  so unlike `rerun-multimap` there is no reason to copy the run).
- [x] `PROG-05` — outputs `results/gene_program_summary.tsv` and
  `results/gene_program_cells.tsv`, plus a report section that surfaces
  `panel_completeness` and `latency_observable_in_rna` so a partial row is not
  over-read.
- [x] `PROG-06` — the EBV LCL regression test. Asserts the unique layer calls
  `latent` on a matrix built to reproduce the cross-mapping scenario, so the
  protection fails loudly if either the evidence layer or the overlap-group
  logic is removed.
- [x] `PROG-07` — **closed 2026-10-03 on the called-cell set (after `PROG-17`).**
  Same copy, `rerun-programs` with PROG-17: EBV scores **932** cells (exactly
  the called cells with marker evidence; cross-checked against
  `per_cell_viral.tsv`).
  - Unique layer: **526 latent**, 67 productive, 184 mixed, 155
    indeterminate.
  - Allocated layer: 339 latent, 12 productive.
  - Inversions (latent on unique → productive on allocated): **0**. Of the
    526 latent cells, 319 stay latent and 207 turn `mixed`.
  - The block below is the superseded first pass over 1,679 barcodes.
  Superseded first pass, 2026-10-03, on layer 2's own cell set:
  `rerun-programs` ran on a **copy** of the cat42b EBV run
  (`viral_ref_cat42b/runs_prog07/SRR12682296`, post-SW-13/SW-24, after
  PROG-11, code d3905e8). Layer 2 scored 1,679 cells. Only 932 of them are among
  the 2,763 emptyDrops-called cells, so these are **not** called-cell numbers:
  - Unique layer: **695 latent**, 78 productive, 185 mixed, 721
    indeterminate.
  - Allocated layer: 526 latent, 73 productive.
  - Inversions (latent on unique → productive on allocated): **0**. Of the
    695 latent cells, 478 stay latent on the allocated layer and 217 turn
    `mixed`.
  - Stays `[~]` until `PROG-17` lands. Then rerun on the called set; the
    command is cheap (minutes).
  - Every earlier number below came from `fresh12b`: another index, no
    barcode correction (SW-13) and no cell calling. They are not comparable
    and are kept as history only.
  - `docs/faq.md` and `docs/output_reference.md` carry a dated correction.
  - **Former status:** "under re-verification (2026-09-27); do not cite
    these numbers."
  An audit of the catalogue found lytic genes filed as latent markers for EBV
  (`BaRF1.1`, the ribonucleotide-reductase subunit, matched onto latent `BARF1` by a
  case-insensitive lookup; `BHRF1`, `BNLF2a/b`), so the latent counts below may be
  inflated by lytic reads. Re-measure with `viralscan rerun-programs` once `PROG-11`
  lands. Original entry, retained as history:
  - **Re-measured 2026-09-27 with `PROG-11`'s cited changes** (`rerun-programs`
    on `benchmark_runs/reference_strategy_2026-06-28_fresh12b/runs/ebv__viralscan__combined/SRR12682296`,
    code at `50253a6`). Unique layer: **895 latent** (was 2,240), 236 productive
    (was 102), 311 mixed (was 445), 3,094 indeterminate (was 2,968); cells with
    any marker evidence 4,536 (was 5,755). Allocated layer: **856 latent** (was
    1,277), 42 productive (was 5). So the unique layer's apparent latent gain
    shrinks from +75 % to +4.6 %: most of it was `BARF1.2`, which carried 13,668
    unique and 0 allocated molecules and is not a B-cell latency marker. About
    20 % of evidence-bearing cells reach a latent call; latent sensitivity, not
    direction, is now the open question, which is what `PROG-15`'s 3′/UTR
    measurement should explain. Still provisional: `BHRF1` and `BNLF2a/b` stay
    in the latent set until verified (`BNLF2a/b` carry 0 unique molecules, so
    they cannot inflate unique-layer calls; `BHRF1` could). The pre-change
    outputs are kept beside the new ones as `*.pre-PROG-11.tsv`.
  **measured on the real run, and the reason the design exists.**
  EBV LCL `SRR12682296`: aggregate LATENT 236,342 vs LYTIC 247,633 (ratio 1.15)
  in a cell line latently infected by construction, with `EBNA-1.1` at 920 UMI
  ~155x below `BHLF1` at 142,954 — so per-gene aggregate totals are uninformative.
  Per-marker, the uniquely-placing layer has **0** molecules on `BZLF1` and
  13,668 on `BARF1.2`, where the allocated layer has 9,308 and 0. Computed
  overlap groups: g5 = {EBNA-1, EBNA-2, EBNA-LP}, g50 = {BNLF2a, BNLF2b, LMP-1},
  g39 = {BTRF1, BcLF1}. End-to-end `rerun-programs` on that run: **2,240 cells
  latent, 102 productive, 445 mixed** on the unique layer versus 1,277 latent and
  5 productive on the allocated layer, with **0 cells inverted** between the two.
  The failure mode the unique layer fixes is lost sensitivity, not a wrong
  direction — recorded here so the claim is not overstated later.
- [x] `PROG-08` — **closed 2026-10-03** (user grill: "upgrade KSHV, close
  the rest").
  - **KSHV is now `complete`.** Its latent set:
    - ORF73/ORF72/ORF71 (LANA/v-cyclin/vFLIP): one mRNA family (PMID
      9733875), counted as one breadth unit via a new `CO_TRANSCRIBED`
      override in `extras/build_gene_programs.py`;
    - K12 kaposin, also lytic-induced (PMID 17913828);
    - vIRF-3/LANA2 (K10.5), B-cell latent (PMID 11119611).
  - K1 was removed: it is tied to lytic replication (PMID 27307571).
  - Only KSHV rows of `gene_programs.tsv` changed. Tests:
    `TestKshvLatency`, plus a calling test (cluster gives breadth 1, + K12
    gives 2).
  - The other 7 are documented as biology-limited in `docs/faq.md`. No real
    KSHV dataset exists in the repo yet; that validation is `PROG-18`.
  - Original text: **open.** Four of nine viruses have a genuine latency and
  reactivation split (EBV, CMV, HHV-6A, HHV-7). The other five are `partial`:
  HSV-1/2 latency is a single transcript (`LAT`), VZV's is inferred (ORF4), and
  HHV-6B's GTF carries no attributes at all. For those, `latent` is unreachable
  by construction — so the HSV-1 benchmark cannot demonstrate a latent call.
  Extending them needs either a fuller annotation source (the HHV-6B panel has
  no `product` text at all) or a decision to accept a weaker anchor set.
- [~] `PROG-09` — **open.** The catalogue is 9 viruses. Polyomaviruses (JC/BK/KI/WU,
  MCPyV), HPV, HBV, HDV, GBV-C, HIV-1 and HTLV-1 are deliberately excluded: they
  have no latency/lytic dichotomy representable from the panel's protein-coding
  genes, and for HIV/HPV/HBV the interesting state is DNA-level latency, which
  scRNA-seq cannot observe at all. Adding them would mean emitting
  `not_applicable` rows, which is what the code does for a detected virus with no
  model — a decision to make explicitly rather than by omission.

## WP4E — Named HPV ORFs for the oncogene-versus-capsid contrast (new 2026-09-26)

Objective: give HPV real, named genes so that transcriptional activity of the
E6/E7 oncoproteins can be told apart from passive L1/L2 capsid transcription in
oropharyngeal and tonsillar tissue. This is a **different axis from WP4D** and
does not change `PROG-09`: HPV still has no latency/lytic dichotomy representable
from protein-coding genes, so it remains out of the gene-programme catalogue and
`PROG-09`'s reasoning about DNA-level latency stands. What is new here is that
the ORFs are *named at all*, which WP4D did not require.

Trigger: the built index at
`covid_viralscan/viralscan_ref/` represents HPV with 4 accessions
(`NC_001526.4`, `NC_001356.1`, `NC_001352.1`, `NC_003461.1`) and 22 t2g rows whose
gene IDs are RefSeq `locus_tag` values — `HpV16gp1`…`HpV16gp8`, `HpV1agp1`…,
`HpV2agp1`…, `Hpv1gp01`…. Nothing in the index says which is E6 and which is L1.

- [x] `HPV-01` — **investigation first: the records already carry semantic
  names.** Every HPV complete-genome record examined annotates each CDS with
  `/gene="E6"` (RefSeq and most INSDC submissions) or
  `/product="transforming protein E6"` (the records that omit `/gene`).
  **No coordinate table is used or needed**, and the per-row `name_source` column
  records which qualifier each name came from. The bundled RefSeq GTFs already
  carried the answer in a `gene` attribute; the packager kept `locus_tag` as the
  ID and dropped `gene`, which is why the index cannot answer the question.
- [x] `HPV-02` — **a coordinate table was rejected on evidence, not taste.**
  Papillomavirus genomes are submitted as linearised circles cut at the
  submitter's chosen point, so the same E6 ORF sits at 7125-7601 in
  `NC_001526.4` (HPV16) and 105-581 in `NC_001357.1` (HPV18) — opposite ends of
  their records. `NC_001526.4` is cut inside E1 and therefore reports its ORFs as
  `E1, E2, E5, L2, L1, E6, E7` while `NC_001357.1` reports `E6, E7, E1, …`; both
  are the canonical order, differing only by where the circle was opened. The
  build's ORF-order check is rotation-tolerant for exactly this reason, and the
  rotation is asserted by a test.
- [x] `HPV-03` — **names independently confirmed against protein sequence.** The
  build aborts unless every E7 translation carries the LXCXE retinoblastoma-
  binding motif, every E6 translation carries its C-X2-C zinc fingers, E6/E7/L1/L2
  fall in their known length ranges (E6 ~150 aa, E7 ~100 aa, L1 504-569 aa,
  L2 474-525 aa), and L1 is at least 3x either oncogene. HPV16 and HPV18 E7 both
  read `…LXCYEQL…`; both E6 read `…IICVYCKQQL…`. Two apparent invariants were
  falsified by writing them first and watching all 16 records fail: **L1 is not
  the longest ORF** (E1 is, at ~650 aa — it is the replication helicase), and
  **L1 is not always longer than L2** (HPV-2 annotates L2 at 525 aa against L1 at
  511 aa; HPV-1 has them within one residue). Both are recorded in the code.
- [x] `HPV-04` — genotype coverage: **all 14 high-risk types** (16, 18, 31, 33,
  35, 39, 45, 51, 52, 56, 58, 59, 66, 69) plus HPV1 and HPV2, which are carried
  only so a rebuild does not *lose* the two types the current index already has.
  RefSeq has complete genomes for only 4 of the 14 (`NC_001526.4` REVIEWED,
  `NC_001357.1` VALIDATED, `NC_075191.1` and `NC_075233.1` PROVISIONAL); the
  other 10 are INSDC, chosen as the **oldest** complete genome carrying annotated
  CDS, since the earliest submission of a type is the prototype that genotyping
  assays and published amplicons target. Attempted and documented: 208 RefSeq
  papillomavirus complete genomes were enumerated (65 human, 57 types), and every
  one of the 10 remaining high-risk types was confirmed present in INSDC (37-572
  isolate records each) before one was selected.
- [x] `HPV-05` — `extras/build_hpv_reference.py` regenerates the TSV from NCBI,
  cache-first through `ncbi_fetch`'s own cache directory with its SHA-256 sidecar
  convention, so a record already fetched by any other ViralScan entry point is
  reused rather than re-downloaded. Refuses to write a partial catalogue: it exits
  non-zero on any unresolvable accession, unrecognised ORF symbol, failed
  invariant, or duplicate `(accession, gene)` pair.
- [x] `HPV-06` — `src/viralscan/data/hpv_genes.tsv` (124 ORFs, 16 genotypes) and
  `src/viralscan/hpv_genes.py`. Registered in `pyproject.toml` package-data and
  `MANIFEST.in`; `include-package-data = false` means MANIFEST alone would not
  ship it. Gene IDs are namespaced (`NC_001526.4_E6`), because the bare symbol
  `E6` occurs once in each of 16 genomes. Classes: 32 `oncogene` (E6+E7 per
  genotype), 1 `oncogene_locus` (HPV16 E6*), 32 `late_capsid`, 59 `early`.
- [x] `HPV-07` — **`E5` is classed `early`, not `oncogene`, and `E6*`/`E7*` are
  `oncogene_locus`, not `oncogene`.** E5 is a transforming protein several
  reviews call an oncogene, but it is not part of the E6/E7 axis and putting it
  in the oncogene class would make a positive call mean something it does not.
  E6* lacks the PDZ-binding motif and E7* lacks LXCXE, so neither is
  transforming, and their reads are indistinguishable from E6/E7 at the sequence
  level — folding them in would let E6* alone produce a confident positive
  oncoprotein call. Resolved by preferring an isoform symbol in `/product` over
  its own parent in `/gene`, which is what `NC_001526.4` requires.
- [x] `HPV-08` — `--emit-reference DIR` writes a merged FASTA (16 genomes,
  128 kB) + GTF (125 exon lines, 124 unique gene IDs) for `kb ref`, built from the
  parsed CDS rather than from `ncbi_fetch._genbank_to_gtf`. Sequences are
  generated on demand and **not committed**: reference size is a packaging
  decision governed by PR 8 and should be made once, deliberately. Structurally
  validated (FASTA/GTF seqnames agree, blocks within bounds); **not yet
  index-built**, as `kb`/`kallisto`/`bustools` are not on PATH in this
  environment.
- [x] `HPV-09` — **blocker, in `ncbi_fetch.py`.** `_genbank_to_gtf` set
  `gene_id` from `/gene=`, so in a merged multi-genome reference every
  genome's E6 collapses onto one row. Demonstrated: merging `NC_001526.4` and
  `NC_001357.1` through `fetch_reference()` yielded 16 ORFs on **9** gene IDs.
  A public `fetch_genbank()` accessor was also needed so this generator
  stopped importing `_efetch`, `_cache_valid` and `_write_cached`.
  **Both asks were implemented in WP4F (2026-09-26)** and **verified here
  (2026-09-26):** `_genbank_to_gtf` now emits genome-scoped `<accession>_<token>`
  gene IDs, and `fetch_genbank()` is public and caches the raw flatfile.
  Verification: the two cached `NC_001526.4`/`NC_001357.1` `.gtf` files did
  predate the fix (see `ANELLO-11`) — deleted and confirmed they regenerate
  from the retained `.gb` flatfile with no network call. Re-running the exact
  `fetch_reference(["NC_001526.4", "NC_001357.1"], ...)` merge that
  demonstrated the bug now yields **17 distinct gene IDs across 19 ORF lines**
  (the one repeat is a real intra-genome duplicate `locus_tag` on HPV16,
  correctly disambiguated `_dup2` by `_panel_gene_ids`, not a cross-genome
  collision). `extras/build_hpv_reference.py`'s own `fetch_genbank()` wrapper
  now delegates to the public `ncbi_fetch.fetch_genbank()` instead of
  reaching into its privates; `tests/test_hpv_reference.py` and
  `tests/test_ncbi_fetch.py` pass unchanged.
- [x] `HPV-10` — **the scientific limit, recorded so it is not over-read later.**
  What the catalogue supports: an **oncogene-versus-capsid contrast per
  genotype**. E6/E7 are early-region oncoproteins transcribed in
  carcinogen-driven HPV-positive oropharyngeal tumours, while L1/L2 are
  late-region and transcribed only in productive infection, so E6/E7 reports viral
  gene expression where L1 reports virion production. What it does **not**
  support: confident **per-genotype attribution of L1 signal**. L1 is the most
  conserved coding region in the genus — it is what pan-HPV PCR primers target —
  so L1 reads cross-map freely between all 16 genotypes and kallisto's
  multimapping will distribute them; a genotype label on an L1 count is not
  independent evidence of which type is present. E6/E7 are far more
  type-divergent and better behaved, but a cross-mapped count is still not a
  transcript count. Presence-of-HPV-transcripts and the oncogene-versus-capsid
  distinction are supportable; per-type L1 attribution is not. This is also a
  *transcriptomic* catalogue: a transcriptionally silent integrated genome — the
  common state in tonsillar crypt epithelium and the state that drives
  HPV-positive oropharyngeal carcinoma — produces no reads and is invisible
  here, not negative.
  - Closed 2026-10-02 (`v3/docs`): the limit is now in user docs
    (`reference_panel.md` and the FAQ). Presence and E6/E7 vs L1/L2 are
    supportable; per-genotype L1 attribution and low off-type calls are not.
    The measurement stays in `HPV-11` / `CAT-08`.
- [ ] `HPV-11` — next: rebuild the real index with the named HPV panel and
  measure the cross-mapping directly, rather than reasoning about it. Requires
  `HPV-09` or the `--emit-reference` path, plus `kb` on PATH. Until then the
  L1 claim in `HPV-10` is a literature-based expectation, not a measurement from
  this panel. Do not publish a per-type HPV number before this row is closed.
  - **2026-09-26, not evidence:** a local run tiled 1,837 error-free 90 bp reads
    from the panel's own E6/E7/L1/L2 CDS, pseudoaligned them against a
    panel-only `kb ref` index, and found **0** cross-genotype assignments (19
    within-genome only, E6/E8 and L1/L2 junctions). That result is expected by
    construction: reads drawn from indexed sequences match themselves, and HPV
    types are defined by ≥10 % L1 divergence. It says nothing about HPV-10's
    actual risk, which is reads from strains *not* in the panel. The valid test
    is `CAT-08`'s leave-one-out design: held-out isolates (HPV16 lineage variants,
    types outside the panel) scored for where their reads land.

## WP4F — Real gene structure for the Anelloviridae panel (new 2026-09-26)

Objective: replace the Anelloviridae panel's one-placeholder-gene-per-genome
representation with the CDS features NCBI actually annotates. This is a
**confirmed bug fix, not a research question**, and it is a different axis from
WP4D/WP4E: anelloviruses have no latency/lytic dichotomy and no oncogene/capsid
contrast, so they stay out of both catalogues. What they do need is for a hit to
be attributable to a *locus* rather than to a conservation rank.

### Root cause (confirmed in code, not inferred)

The panel never took the GenBank path, and the reason is a **discarded** GTF, not
a failed fetch:

| Site | What it does |
|---|---|
| `src/viralscan/scripts/ncbi_fetch.py` `_fetch_one` | fetches GenBank, runs `_genbank_to_gtf`, writes a real GTF, **discards the GenBank text** |
| `scripts/build_bundled_panel_ref.py:228` (was) | `anello_gtf_texts.append(_genome_as_transcript_gtf(text, acc))` — threw the real GTF away, rebuilt a placeholder from the FASTA |
| `src/viralscan/scripts/build_reference.py:1072` (was) | `_gtf_from_merged_fasta(final_fasta, final_gtf)` — same discard, in `build_anellovirus_reference` |
| `src/viralscan/scripts/build_reference.py:437` | `_genome_as_transcript_gtf` emits `gene_id = f"{accession}_gene{seq_idx}"` — the placeholder itself |

So the briefed hypothesis (that anelloviruses took
`_whole_genome_gtf_from_fasta` in `ncbi_fetch.py`) is **wrong for 97.7 % of the
panel**. Measured on the pre-existing NCBI cache: 1,995 of 2,042 cached `.gtf`
files already carried real CDS-derived exons and only 47 were placeholders. The
dominant emitter is `_genome_as_transcript_gtf` in `build_reference.py`, reached
from two callers. `ncbi_fetch._whole_genome_gtf_from_fasta` is real and does fire
for genuinely CDS-less records, but it is the minority path.

There was a **second, independent** defect that would have survived the first
fix: `_genbank_to_gtf` set `gene_id` from the bare `/gene=` or `/product=`, which
are not unique. Across the cached panel only 2,316 distinct gene IDs existed
across 1,995 annotated genomes, with `ORF1` shared by 150 genomes and `orf1` by
63. A merged index would have collapsed 1,995 genomes onto 2,316 columns with
cross-genome identity. `HPV-09` independently reports the same defect for HPV.

### Measured CDS coverage (full panel, 2,042 accessions)

Every accession was retrieved and audited; the numbers below are the generator's
own output, not a projection. **1,995 of 2,042 (97.7 %) carry real gene
structure — 2,515 genes — and 47 have no CDS feature in NCBI at all.**

| genus | panel | annotated | genes | coverage |
|---|---:|---:|---:|---:|
| Betatorquevirus | 1,542 | 1,517 | 1,699 | 98.4 % |
| Alphatorquevirus | 211 | 204 | 361 | 96.7 % |
| Anelloviridae (unclassified) | 185 | 175 | 264 | 94.6 % |
| Gammatorquevirus | 78 | 74 | 150 | 94.9 % |
| Hetorquevirus | 8 | 8 | 9 | 100 % |
| Samektorquevirus | 8 | 7 | 10 | 87.5 % |
| Gyrovirus | 6 | 6 | 17 | 100 % |
| Memtorquevirus | 4 | 4 | 5 | 100 % |
| **TOTAL** | **2,042** | **1,995** | **2,515** | **97.7 %** |

A 206-accession stratified pre-implementation survey across all eight genera gave
205/206 (99.5 %), consistent with the full run. One record was identified as
CDS-less in the survey (`KP343852.1`, a bare `source` feature); 47 are CDS-less in
full.

**Coverage is not the problem. The low CDS *count* is — and it is biology, not a
gap.** 1,740 of the 1,995 annotated genomes (87 %) carry exactly **one** CDS, and
1,713 of the 2,515 genes have the product `ORF1`. 75 % of the panel is
Betatorquevirus, whose ~2.8–3.0 kb TT-mini genomes genuinely have a single ORF
spanning the genome. There is no missing annotation to recover. The gain from
this fix is that the single gene is *named, product-labelled, and
genome-scoped* instead of anonymous, and that the 255+ multi-ORF genomes finally
get their real structure. **A further honest caveat on annotation quality:** 576
of 2,515 genes (23 %) carry the generic product `hypothetical protein`, and only
258 genes carry a `/gene` symbol at all — so "named gene" means
product-labelled for most of the panel, not functionally annotated.

- [x] `ANELLO-01` — the audit above, run **before** any code was written, so the
  decision to keep a placeholder fallback is evidence-based rather than assumed.
- [x] `ANELLO-02` — `ncbi_fetch.py` now caches the **raw GenBank flatfile** as
  `<accession>.gb` with the existing `.sha256` sidecar convention, and exposes it
  through a public cache-first `fetch_genbank()`. Re-deriving the annotation
  after a code change therefore costs zero NCBI requests. This also satisfies
  the public accessor `HPV-09` asked for; that row is left for its owner to
  close.
- [x] `ANELLO-03` — **gene naming: `<accession>_<token>`, token from
  `/locus_tag` → `/gene` → `/protein_id` → `cds<N>`.** `/locus_tag` is the
  submitter's stable locus name and wins when present — it is what makes the
  reference TTV record emit `NC_002076.2_TTVgp1/2/3`, i.e. the same
  `TTV_TTVgp1` names the bundled RefSeq GTF already uses. `/protein_id` is the
  fallback because it was observed to be reused by **zero** of the 206 sampled
  accessions, making it the only globally unique identifier NCBI offers for the
  238/275 sampled CDS features that carry no `/gene` and no `/locus_tag` at all.
  `orf2/5` is sanitised to `orf2_5`. The bare symbol is preserved in `gene_name`
  (GTF) and `gene_symbol` (TSV). **Missing and duplicated products:** a CDS with
  no identifier gets `cds<N>`; a token repeated within one genome gets `_dup2`,
  `_dup3` (not a bare ordinal, which is indistinguishable from a real `orf12`).
- [x] `ANELLO-04` — **the naming is genome-scoped on purpose.** The existing
  bundled convention is `{virusToken}_{locusTag}` (`TTV_TTVgp1`); for a
  2,042-genome panel the virus token must be the accession, because `TTV` would
  collapse the panel into one label. `{accession}_{token}` keeps the shape and
  works unchanged with `anello_name_map()`'s boundary-aware prefix rule.- [x] `ANELLO-05` — **circular topology, measured rather than assumed — and the
  panel does contain a wrap.** Anelloviridae are circular ssDNA, but NCBI
  annotates in a *linear* representation and only 906 of 2,515 genes (36 %) come
  from records that even declared `circular` on the LOCUS line, so the
  declaration is recorded (`topology`) and never used to interpret coordinates.
  Across the full panel there is **exactly one origin-spanning gene**:
  `KU243129.1` (2,824 bp, `ss-DNA`, `circular`) annotates
  `join(2677..2824,1..80)` — exon 1 at the end of the linear representation, exon
  2 back at the origin. NCBI writes the intervals in transcript order, so a
  coordinate sort would emit `1..80` first and silently transpose the gene's two
  exons into a scrambled transcript. The parser therefore **never re-sorts**; the
  only normalisation is reversing the interval list for minus-strand features,
  and `_origin_spans` flags a wrap structurally (in transcript order the first
  interval starts after the last interval ends) with `origin_spanning="true"`.
  A test asserts the real `KU243129.1` gene's exon *order*, not just its
  existence, so a regression here is a visible transposition rather than a
  silent one. No interval anywhere in the panel falls outside `[1, length]`.
- [x] `ANELLO-06` — **genogroup is not derivable and was not invented; the
  column ships empty except two verbatim NCBI values.** Corrected 2026-09-26,
  and re-verified directly against the cache rather than re-asserted: an
  earlier draft of this row claimed **zero** of 2,058 flatfiles carry a
  `/genotype` qualifier and cited two `/note`-only hits instead
  (`NC_002076.2`, `JN980171.1`). That "zero" claim was simply wrong — grepping
  the retained flatfiles for a literal `/genotype=` qualifier line
  (`grep -l '/genotype=' ~/.cache/viralscan/ncbi/*/*.gb`) finds exactly **two**
  records, `NC_014081.1` (`"6"`, `/organism` "Torque teno virus 3") and
  `NC_014094.1` (`"28"`, `/organism` "Torque teno virus 6") — both of which
  contradict their own organism species number, which is the concrete
  evidence that a genogroup must never be inferred from `/organism`. The
  `NC_002076.2`/`JN980171.1` pair is real too, and distinct: those two carry
  `genotype` only as free text inside a `/note` (confirmed separately with
  `grep -l '/note=.*genotype'`), never as a structured qualifier, so they stay
  out of `source_genotype` for the same reason the column is not back-filled
  from prose generally. Four different records, two different mechanisms —
  not a contradiction to reconcile, just two separate, now-verified facts.
  The shipped TSV column is named `source_genotype`, not `genogroup`,
  precisely because it is NCBI's own `/genotype` qualifier copied verbatim
  rather than a derived or inferred genogroup, and it is empty for all but
  the two `/genotype=` accessions. `test_source_genotype_is_never_invented`
  pins the exact pair so a future hand-fill or a broader NCBI regression
  cannot pass silently. Retained source fields, recomputed from the shipped
  TSV rather than copied from the generator's own (stale) docstring:
  `source/isolate` (2,338 genes across 1,899 accessions, laboratory sample
  codes such as `MDJHem2` or `SAfiA-468-6`) and `source/strain` (102 genes
  across 31 accessions), carried verbatim so a classifier can be fitted later
  without re-fetching.
- [x] `ANELLO-07` — `extras/build_anellovirus_genes.py` → `anellovirus_genes.tsv`
  → `viralscan.anellovirus.gtf_text_for()`, mirroring the `build_gene_programs.py`
  → `gene_programs.tsv` → `gene_programs.py` precedent. Cache-first and
  resumable (flatfiles cached, nothing written until every accession is
  attempted), `--accessions` / `--limit` / `--per-genus` for subset runs,
  `--min-coverage 0.95` so a silent NCBI regression cannot ship a
  mostly-placeholder catalogue, and a coordinate/uniqueness audit that refuses to
  write on any violation.
- [x] `ANELLO-08` — all three discard sites now consume the catalogue:
  `build_bundled_panel_ref.py` Step 4b, `build_anellovirus_reference` Step 4, and
  `build_combined_reference` Step 3 (**scoped to panel accessions only**, so the
  curated 195-genome bundled panel keeps its byte-identical whole-genome GTF and
  no existing index changes shape). Uncovered accessions still get a placeholder,
  because `kb ref` silently drops a sequence with no GTF row and the genome would
  then be neither quantified nor detectable.
- [x] `ANELLO-09` — `tests/test_anellovirus_reference.py`: 42 offline tests
  (panel integrity, no surviving `_gene1`, genome-scoped and unique gene IDs,
  spliced genes have >1 exon, coordinates within `[1, genome_length]`, strand
  `±`, ORF1/Rep present in all eight genera, genogroup stays empty, committed
  GenBank fixtures for the converter including the circular-wrap and
  minus-strand cases, generator run offline from a seeded cache) plus three
  `@pytest.mark.network` tests. The default selection needs no network and runs
  in ~4 s.
- [x] `ANELLO-10` — registered in `pyproject.toml` package-data, `MANIFEST.in`
  and `config/public_ship_scope.json` wheel **and** sdist allowlists, mirroring
  `gene_programs.tsv` / `hpv_genes.tsv`. Not added to `.dockerignore` /
  `docker_context`, matching how those two shipped.
- [x] `ANELLO-11` — **behaviour change to a private helper, recorded because it
  is observable.** `_genbank_to_gtf` gene IDs are now genome-scoped, so a
  cached `<acc>.gtf` written by an older build is stale. Delete the affected
  `~/.cache/viralscan/ncbi/<acc>/<acc>.gtf` (and its `.sha256`) to re-derive; the
  retained `.gb` means that costs no network. One assertion in
  `tests/test_ncbi_fetch.py` was updated to the new contract and the
  genome-scoping property given its own test.
- [~] `ANELLO-12` — **the honest limit, recorded so it is not over-read.** The
  bug is fixed and the 99.8 %-in-one-bucket artifact is now *explainable and
  testable*: with one anonymous gene per genome there was nothing else it could
  have been, and the decisive test is to rebuild the index and check that
  anellovirus UMI spread across many genome-scoped gene IDs instead of
  concentrating in the most conserved one. **That test is not run here** — it
  needs `kb ref` and a re-run of the COVID sample, neither available in this
  environment (`ANELLO-13`). **What remains true regardless, and is the
  important caveat: 1,740 of 1,995 annotated genomes (87 %) carry exactly one
  CDS spanning the whole genome.** For those genomes a "real gene" is still a
  whole-genome transcript, so the count is still a whole-genome count and still
  cross-maps against the rest of the panel in proportion to conservation. The fix
  makes that cross-mapping *attributable and visible* — a count can now be
  traced to a genome and a product, and a conservation-driven skew is
  distinguishable from a single-genome infection — but it does **not** make
  per-genotype anellovirus quantification defensible, and it cannot: a
  unique-sequence argument is needed, not a gene name. That is a different piece
  of work, closer to `SENS-06` (measured k-mer capture) than to this fix. The
  honest summary is: **genus-level anellovirus load becomes interpretable;
  genotype-level anellovirus load does not.**
  - Blocked (2026-10-02): the decisive test is `ANELLO-13` (`[ ]`, not started). Kept `[~]` as a recorded limit.
- [ ] `ANELLO-13` — next: rebuild the reference with the real-gene panel and
  **measure** the covid artifact rather than reasoning about it. Requires `kb` on
  PATH and a re-run of the COVID scRNA-seq sample; then compare the anellovirus
  UMI distribution across genome-scoped gene IDs against the 99.8 % single-bucket
  baseline. Until this row closes, the covid number in this section is a
  measurement of the *old* build and says nothing about the new one. Do not
  report a per-genotype anellovirus number before this row is closed.
  - **2026-09-27:** the bespoke covid build cannot run this test. It reads the
    static `references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.gtf`,
    which still carries 6,126 `_gene` placeholder IDs, so a rebuild through
    `slurm_build_ref_v2.sh` would reproduce the old artifact. Decided path: native
    `viralscan build-ref`, which first needs `CAT-01` (build-ref currently discards
    the real GTF). The retraction question is already answered by F-005
    (`.living/findings/`); this row only validates the WP4F code fix.
- [ ] `ANELLO-14` — the integration test
  `tests/integration/test_anellovirus_chain.py::TestAnellovirusLabelingChain::test_build_anellovirus_reference_produces_labelable_gtf`
  fails, and the failure predates WP1C (it fails at `af4d5b3`). The test still
  asserts the pre-WP4F placeholder `NC_002076.2_gene1`, but the builder now
  emits real gene IDs (`NC_002076.2_TTVgp1`…). Update the test's expected IDs
  to the WP4F gene structure. Do not revert the builder.

## WP4G — Anellovirus detection you can trust (new 2026-09-27)

Objective: make an F-005-type host-homology artifact impossible to publish by
default, and make anellovirus sensitivity a measured quantity. An audit
(2026-09-27) found no default step that would block the covid Alphatorquevirus
call today. "Reliable" here means a genus-level LOD95 from planted reads in 10x
geometry plus host-only negatives that never reach a reported call; without an
orthogonally confirmed positive sample (`REF-10`) every anellovirus result stays
`screening_only`, and no code change lifts that ceiling.

- [ ] `ANDET-01` — `accession_breadth` is always 1.0: it is computed over
  `found_genes`, which are already detected (`detection.py:166`, `:501-506`).
  Compute it over every index gene of the virus and add per-accession
  genome-coverage breadth, F-005's deciding gate (≤3.41 %).
- [ ] `ANDET-02` — read `host_homology_annotations.tsv` (written at
  `build_reference.py:768`, read by nothing) in detection; demote calls
  concentrated in host-homologous regions; surface `eve_risk` in the report.
- [ ] `ANDET-03` — `claim_scope` column (`screening_only` for Anelloviridae) in
  `viral_summary.tsv` and the report. `REF-10`'s label exists only in prose today.
- [~] `ANDET-04` — evidence replay reads the raw FASTQs (`evidence_run.py:154-155`,
  `:176`) instead of the host-filtered `kb_r1`/`kb_r2`; `--virus ttv` resolves to
  the 185 unclassified genomes only (`evidence.py:64`). Fix both; auto-run
  read-level host confirmation for detected anellovirus genera.
  Both bugs fixed 2026-09-27. `replay_fastqs()` returns the pair `kb count`
  quantified (and reconstructs `host_filtered/` for configs that predate
  `kb_r1`), both the replay and the extraction use it, and a missing input now
  fails closed instead of silently replaying different reads. `ttv` (or
  `Anelloviridae`) now selects every genus plus the unclassified group and the
  bundled "Torque teno virus" label; reproduced before the fix as
  `('Anelloviridae', ['AB303555.1_ORF1'])` on a three-genus fixture. Still open:
  the auto-run Snakemake rule.
  - Blocked (2026-10-02): the auto-run gate depends on `ANELLO-PRIOR.3` (`[ ]`). Default on vs opt-in is a user decision.
- [ ] `ANDET-05` — one genus name per genome: bundled `TTVgp1` IDs resolve to
  "Torque teno virus" while genome-scoped `NC_002076.2_TTVgp1` resolves to
  "Alphatorquevirus". Also fix the `UUKU` and `VARV` aliases.
- [ ] `ANDET-06` — correct three docs: `REF-01`'s "now the default" (the CLI
  default is off), the STAR-defaults claim at `host_filter.py:176-192`, and the
  genome-reference ranking in `docs/cli_reference.md:200`.
- [ ] `ANDET-07` — make the expanded panel the default quantification reference,
  **only after `ANDET-01`–`ANDET-04` land**. `protocol.yaml`'s
  `anellovirus_expanded` entry says "never the default" and sits under the frozen
  `frozen_inputs` digest, so this needs a `DEV-0xx` deviation record.
- [ ] `ANDET-08` — genus-level LOD95: plant held-out genomes from all 8 genera in
  10x v3 geometry into a checksum-pinned healthy-PBMC background (`VAL-04`),
  probit fit per the `SCI-03` method; the unplanted background must yield no
  reported call. A minimal `VAL-01` slice; closes `SENS-06` for this family.
- [ ] `ANDET-09` — STARsolo anellovirus alignment branch that does not depend on
  kallisto (grill 2026-10-03). The plan this replaces assumed bulk and rustar.
  The user chose scRNA inside v3, STAR/STARsolo, and default-on.
  - The evidence replay (`ANDET-04`) only sees reads in kallisto ECs, and
    31-mer capture is capped near 0.95^31 (F-013). Alignment is not.
  - Input is `host_filtered/R2`+`R1`. Without `--host-filter starsolo` the
    branch records `skipped_no_host_filter` (user).
  - [x] `ANDET-09a` index: `anello_star/` next to the kb index. Anelloviridae
    records of the panel's own `viral.fa`, one gene per contig, STAR
    SAindexNbases 10 / ChrBinNbits 11 (`scripts/build_bundled_panel_ref.py`).
    Built for cat42b 2026-10-03 (job 25696001): **2,040 contigs**, STAR
    2.7.11b, genome version 2.7.4a, `manifest.json` records the FASTA sha256.
  - [x] `ANDET-09b` Snakemake rule `anello_align`. STARsolo, Unstranded
    (F-020), GeneFull, `--soloMultiMappers Unique EM`, plan §11 filters,
    MultimapNmax = SAMmultNmax = 100.
  - [x] `ANDET-09c` evidence metrics merged into `viral_summary.tsv`, plus
    `detection_source`. Any genus with ≥1 unique alignment molecule gets a
    row (user). Labels, never filters.
  - [x] `ANDET-09d` `--anello-align/--no-anello-align`, **default False**.
    **Falls back to off if `ANDET-09e` fails** (user); the DEV record is
    written only on pass.
  - [x] `ANDET-09e` acceptance (`ANELLO-PRIOR.1` slice).
    - Plant: ~8 Modha BK genomes (absent from the panel), 10/100/1000
      molecules × 5′/3′/uniform windows, into 5 M x213 `host_filtered` pairs.
    - Negative: a synthetic set of host cDNA, GRCh38 and poly-G/poly-A/CAG/TSO
      artefact reads (user).
    - Pass: STARsolo planted−unplanted ≥ kallisto for every genome and >
      kallisto below 95 % identity; 0 alignment molecules in the negative; no
      planted read lost to the host filter.
    - **Result 2026-10-03** (jobs 25696086 plant, 25696105 arms, 25696114
      report). Criterion 2 **pass**: 0 alignment molecules in the synthetic
      negative, and 0 for kallisto too. Criterion 3 **pass**: 26,640/26,640
      planted reads survived the host filter. Criterion 1 **fails 1 of 8**, so
      the flag ships off (pre-registered user decision).
      Per genome, `uniform` window (molecules of 1,110 planted reads):

      | genome | genus / identity | kallisto | alignment |
      |---|---|---|---|
      | DRR140164_NODE_2 | alpha 0.969 | 184 | 1067 |
      | SRR7167047_NODE_1 | alpha 0.876 | 425 | 1071 |
      | SRR8862005_NODE_7 | beta 0.937 | 633 | 1059 |
      | SRR2037085_NODE_5080 | beta 0.452 | 499 | 861 |
      | SRR2037083_NODE_2038 | he 0.740 | 720 | 949 |
      | SRR6316308_NODE_9 | samek 0.182 | 50 | 121 |
      | SRR2037083_NODE_1170 | mem 0.108 | 29 | 107 |
      | **SRR2037085_NODE_7436** | **gamma, no alignment** | **3** | **0** |

      The one failure is the genome minimap2 could not place anywhere in the
      panel: the branch's ≥80 %-matched / ≤8 %-mismatch filters exclude it,
      while kallisto's EM still assigns a trickle. Whether those 3 molecules are
      real or EM mis-assignment is unresolved and decides whether the criterion
      was fair. **Next: relax the two filters on a sweep and re-measure
      (`ANDET-09f`), then revisit the default.**
    - Window effect, reported apart from criterion 1: the 5′ and 3′ windows lie
      outside the panel's CDS-only anellovirus models, so kallisto recovers ~0
      there regardless of identity (alpha 0.969: kallisto 0/2, alignment
      1110/1110). That is the UTR gap, not divergence, and it is why criterion 1
      uses the `uniform` window only.
- [ ] `ANDET-09f` — recalibrate the alignment filters against the held-out
  plant: sweep `--outFilterMatchNminOverLread` and `--outFilterMismatchNoverLmax`
  and report recovery vs the negative's false-positive count per setting. The
  acceptance set is built and reusable
  (`viral_ref_cat42b/runs_anello_plant/`, `--reuse-per-genome`).
  - 2026-10-03, code landed (09a builder step, 09b, 09c, 09d flag):
    - 09a still needs the cat42b index build (`--anello-star-only`).
    - 09d's default stays provisional until 09e.
    - Found by the real-STAR test: STARsolo drops homopolymer UMIs and then
      writes `CB:Z:-`/`UB:Z:-`, so those reads never count as molecules.
      Barcodes are corrected against kb's own on-list (`chemistry.onlist_path`)
      so molecules line up with the kallisto matrix. With no list, STARsolo
      keeps the raw CB.
  - **Held-out set (2026-10-03).** Modha et al. 2025, 829 genomes, TPA
    BK068993-BK069821. **0 of 829 are in the panel** (checked by accession),
    and whole-genome identity to the nearest panel genome puts **788 below
    85 %** — where an exact 31-mer cannot survive (0.85**31 ~ 0.007). Two are
    100 % identical to a panel genome under another accession.
    - Trap: the repository's `Modha_contigs.fas` is **ORF1 only** (all 829
      lengths equal the metadata's `ORF1_len`). Planting from it would sample
      just the hypervariable ORF1 and leave the 5'/3' windows undefined.
      Sequence and CDS coordinates both come from
      `Modha_genomes_annotated.gbk` instead.
    - Selection spans 6 genera x 5 identity bands, preferring full-length
      genomes: alpha >=95 (0.969) and 85-90, beta 90-95 and <85, gamma
      no-alignment, plus samek/mem/he <85 (the under-sampled genera).
  - Acceptance is scored by `scripts/anello_acceptance_report.py`: a
    per-genome assay (planted reads pulled back out by read ID, each
    (genome, window) set run through `kb count` and the branch's own STAR
    command, so recovery is a plain fraction) plus the three whole-arm runs
    (`viral_ref_cat42b/anello_acceptance.sbatch`).

## WP4H — Comprehensive human-virus catalogue (new 2026-09-27)

Objective: a generated, frozen human-host viral catalogue with real gene
structure and deterministic names. No comprehensive index exists: the widest
build (`_misc/viralscan_panel_ref_genomic`) holds ~2,216 genomes, ~88 % of them
anelloviruses, so ~99 other species.

- [x] `CAT-01` — `build-ref` discards the real GTF (`build_reference.py:609`) and
  turns every non-anellovirus accession into one `{acc}_gene1` gene (`:704`,
  `:930-947`), which would silently disable gene programmes on any natively built
  index. Use the `_genbank_to_gtf` output; emit `exon` rows for CDS-only
  features; placeholders only for CDS-less records.
  Done 2026-09-27. `index_gtf_by_seqname()` splits the merged GTF that
  `fetch_reference` already returns, and `viral_gtf_block()` picks the best
  annotation per accession: the packaged anellovirus catalogue, then the real
  NCBI CDS structure, and the whole-genome placeholder **only** for records
  carrying neither. The build logs the three counts and warns, naming
  accessions, whenever a placeholder is used, so a placeholder-heavy index is
  visible instead of silent. Measured on the real NCBI cache: HPV16
  `NC_001526.4` now yields **9** genes and HPV18 `NC_001357.1` **8**, where
  both previously collapsed to one `{acc}_gene1` bucket; `NC_002076.2` still
  takes the catalogue path with its 3 `TTVgp` genes. The existing
  `test_combines_mocked_host_and_viral_reference_without_kb_ref` asserted the
  placeholder `NC_045512.2_gene1`; it now asserts the real gene survives, which
  is the behaviour change.
- [x] `CAT-02` — **closed 2026-10-03.** The frozen `broad-discovery` list is
  now `src/viralscan/data/broad_discovery_accessions.tsv`: 4,397 rows
  (shipped 2,345, max 1,783, broad 268, legacy 1). It is a catalogue/list
  only. Building from it (`build-ref --profile broad-discovery`) is
  `REF-03`. Original row: `extras/build_virus_catalog.py` → `src/viralscan/data/virus_catalog.tsv`,
  cache-first via `ncbi_fetch.fetch_genbank()`: NCBI Virus RefSeq complete
  genomes with human host ∪ bundled panel ∪ 2,042 anelloviruses ∪ 16 HPV
  genotypes ∪ SARS-CoV-2. Becomes the frozen `broad-discovery` list (`REF-13`).
  Generator and seed catalogue done 2026-09-27; the Tier 1 union is `CAT-09`.
  Cache-first via `fetch_genbank()`, reusing `_locus_fields` and
  `_source_qualifiers`; taxonomy comes from the ORGANISM lineage. The four
  curation columns (`tier`, `persistence_class`, `risk_class`,
  `inclusion_rationale`) are carried forward on a re-run so a regeneration
  never drops a human decision, and `--check` reports drift without writing.
  **Seed catalogue: 2,215 accessions, 204 species, 30 families, 0 failures**
  (170 flatfiles newly fetched, the rest from cache). It detects 16 segmented
  species, including influenza A (8), influenza B (8), influenza C (7),
  Rotavirus C (11) and Rotavirus A (10) — the data `CAT-03` needs. Shipped in
  the wheel, sdist and Docker context (702 KB).
- [x] `CAT-03` — names from the catalogue (genome-scoped prefix → species/genus);
  segmented viruses grouped; target 0 % unnamed gene IDs (13.4 % today on the
  covid index).
  Done 2026-09-27. `src/viralscan/virus_catalog.py` loads the generated
  catalogue and exposes `catalog_name_map()` keyed by accession, versioned and
  bare, so the boundary-aware prefix rule already in `virus_grouping` resolves
  `{accession}_{gene}` with no new matching logic. `merged_name_map()` layers
  it under the anellovirus genus map and `VIRUS_NAME_MAP`, which keep
  precedence so no existing output is renamed; `detection.py`, `umap.py` and
  `evidence.py` now use it. A missing catalogue returns `[]` rather than
  raising, so an older checkout keeps working.
  **Measured on the 2,313-record reference:** unnamed gene IDs 173 → **0
  (100 % named, the target)**, and distinct groups 182 → 107 as segments and
  strains collapse into their species. Influenza A's 8 segments now group into
  one "Influenza A virus" instead of 8 separate viruses; before the change the
  same 8 gene IDs produced 8 rows.
- [ ] `CAT-04` — risk classes: exclude human endogenous retroviruses; flag
  integrated ciHHV-6, `EVE_RISK_GENERA`, and vector/reagent contaminants.
- [ ] `CAT-05` — duplicate guard (`validate_reference_records`) on every build
  path, including `scripts/build_bundled_panel_ref.py` (the `NC_002076.2`
  duplicate broke `kallisto index` once already).
  Second case found 2026-09-27 and still shipping: `NC_038359.1` and
  `AB303562.1` are the RefSeq and GenBank copies of one Gammatorquevirus
  genome. Both have canonical leave-one-out capture of exactly 1.000000, which
  is the detector — a genome whose every 31-mer is also contributed by some
  other panel member is by definition redundant. The guard should key on
  sequence identity, not accession, because these two differ by accession.
- [ ] `CAT-06` — the index manifest carries its own viral GTF and catalogue and
  `analysis.py` reads them, so `-gtf` no longer silently drops the panel and
  nothing depends on the unregistered Zenodo DOI (`REF-11`).
- [ ] `CAT-07` — diversity-aware representatives for high-diversity families,
  chosen by leave-one-out k-mer capture (the `REF-01` method), not one exemplar.
- [ ] `CAT-08` — build with native `viralscan build-ref` and measure leave-one-out
  confusability per family: where reads from non-indexed isolates land. Also the
  valid `HPV-11` test.

### WP4H expansion — breadth and measured diversity (new 2026-09-27)

Asked for: cover as many viruses as matter for infection and reactivation, and
give SARS-CoV-2, influenza and the torque teno viruses real strain diversity
instead of one exemplar.

**Census of the current reference**
(`references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.fa`,
measured 2026-09-27):

| | |
|---|---|
| Records | 2,313 |
| Anelloviridae + Gyrovirus | 2,041 (88 %) |
| Everything else | 272 records → 107 distinct species |
| …of which HHV-6B "pseudocontig" gene fragments | 97 records, not a genome |

Non-anellovirus content is therefore ~175 records over ~106 species, and nearly
every species has exactly one genome.

**Absent entirely** (checked by accession and by name):
- **SARS-CoV-2** — no `NC_045512`. Only SARS-CoV-1, MERS and HCoV-229E.
- **HIV-1, HIV-2, HTLV-1, HTLV-2** — the canonical latent retroviruses.
- **HCoV-OC43, NL63, HKU1** — 3 of the 4 seasonal coronaviruses.
- **hMPV, bocavirus, influenza D, TSPyV, HPyV6/7.**
- **HPV18** and every high-risk type except HPV16. The 16-genotype WP4E
  catalogue exists in `hpv_genes.tsv` but is in no index.

**Present but single-strain:** influenza A is 8 segments of one 1934 lab strain
(A/Puerto Rico/8 H1N1); influenza B one 1940 strain; influenza C one 1950
strain. No H3N2 and no circulating isolate.

**Decisions (2026-09-27, user):**
1. **Two-tier catalogue.** Tier 1 = every human-host RefSeq virus species, one
   genome each, the breadth floor. Tier 2 = a literature-curated
   persistence/reactivation set that gets deep strain diversity.
2. **Diversity is measured, not quota'd** — isolate counts set by leave-one-out
   31-mer capture (the `REF-01` method).
3. **Bundle sequences in the package** (self-contained, offline). This runs
   against PR 8's "move data off the package" direction; `CAT-13` is the gate.

- [x] `CAT-09` — **closed 2026-10-03** (user grill: eukaryotic, one per
  species, catalogue only). Virus-Host DB human-host RefSeqs (1,496 taxa)
  were filtered as follows:
  - phages out: 1,434 taxa left;
  - not yet catalogued: 1,026 taxa / 1,296 accessions, all fetched, 0
    errors;
  - host re-check: 1,191 kept and 105 sent to
    `analysis/cat09_sweep/review_host.tsv`;
  - collapse on NCBI Taxonomy species-rank taxids: 232 species, of which 68
    were already catalogued;
  - result: **164 species / 268 accessions added as `panel=broad`** (not
    indexed).

  Old rows are byte-identical. The script is
  `extras/cat09_human_host_sweep.py` (taxid→species cache in
  `analysis/cat09_sweep/species_taxids.tsv`); the merge is 357cb9e.
  - Defaults:
    - segments of a human-host taxon are all-or-none;
    - an organism or common name that is already catalogued is skipped, to
      avoid a display-name collision ("Norwalk virus");
    - phages are filtered on the GenBank lineage too.
  - Original row: Tier 1 breadth floor: NCBI Virus RefSeq complete genomes with
  human host, one representative per species. Closes SARS-CoV-2, HIV-1/2,
  HTLV-1/2, OC43/NL63/HKU1, hMPV, bocavirus, influenza D, TSPyV and HPyV6/7 in
  one step.
  **First batch done 2026-09-27**, closing every gap the census named. Each
  accession was resolved by a live NCBI `esearch`/`esummary` lookup, never from
  memory, per the warning in `docs/reference_panel_research_2026-09-27.md`:
  SARS-CoV-2 `NC_045512.2`; HIV-1 `NC_001802.1`, HIV-2 `NC_001722.1`; HTLV-1
  `NC_001436.1`, HTLV-2 `NC_001488.1`; HCoV-OC43 `NC_006213.1`, NL63
  `NC_005831.2`, HKU1 `NC_006577.2`; hMPV `NC_039199.1`; bocavirus
  `NC_007455.1`; TSPyV `NC_014361.1`, HPyV6 `NC_014406.1`, HPyV7 `NC_014407.1`;
  simian foamy `NC_001364.1`; influenza D `NC_036615.1`–`NC_036621.1` (7
  segments); and the 16 WP4E HPV genotypes including HPV18 `NC_001357.1`.
  Catalogue 2,215 → **2,249 accessions, 232 species, 31 families**; influenza D
  groups 7 segments into 1 virus. Remaining for this row: the systematic
  human-host RefSeq sweep beyond the ViralZone-derived list.
  **Cross-check against the Serratus screen**
  (`sources/viral_panel/Serratus_hits_all_viruses.tsv`, 129 viruses): all 17 of
  its `reactivation_candidate` viruses were already present before this batch —
  EBV leads at 6,597 high-confidence hits, then adenovirus 2,227, HSV-1 1,224,
  HCMV 1,119. That table marks HIV and HTLV as non-candidates and omits
  SARS-CoV-2 entirely, which is why the census, not the table, drove this list.
- [ ] `CAT-10` — Tier 2 persistence/reactivation set, curated from primary
  literature with a citation per entry and checked with `verify-references`
  before commit. Each row carries a `persistence_class` (`latent-episomal`,
  `latent-integrated`, `chronic-productive`, `persistent-commensal`,
  `recurrent-lytic`). Set: the 9 human herpesviruses; polyomaviruses (BK, JC,
  MC, WU, KI, TSPyV, HPyV6/7); HIV-1/2 and HTLV-1/2; HBV and HDV; HPV; the
  anelloviruses; parvovirus B19; adenovirus; measles (SSPE); enterovirus; HCV;
  pegivirus.
- [x] `CAT-11` — promote the 31-mer capture measurement to a real tool,
  `scripts/measure_kmer_capture.py`: per-genome coverage, zero-coverage
  fraction, leave-one-out capture, P(90 bp fragment captured), using
  `sensitivity.DEFAULT_K` and `fragment_capture`. The evidence behind `REF-01`
  currently lives only in the untracked `src/viralscan/scripts/kmer3.py`, which
  hard-codes institutional paths and must not be committed. Positive control
  already validated: EBV `NC_007605.1` against itself = 100.0 % over 144,283
  31-mers.
  Done 2026-09-27 as `scripts/measure_kmer_capture.py`, with 12 unit tests and a
  `--self-check`. It reproduces every `REF-01` coverage number exactly: EBV
  144,283 31-mers; bundled-panel median 0.00 % and 85.80 % zero-coverage;
  Betatorquevirus 98.4 %; the five minor genera 100 %; expanded-panel
  leave-one-out 20.55 %. **One number does not reproduce and is corrected:**
  `REF-01`'s "P(90 bp fragment captured)=1.0000" for the expanded panel is the
  analytic independence model; walking real windows measures **0.5097**,
  because shared k-mers cluster in conserved blocks instead of scattering.
  `.living/findings/` F-011 has the per-genus table. `REF-01`'s decision is
  unaffected — the expanded panel still takes zero-coverage from 85.8 % to
  0.15 % — but the panel misses about half the reads of a strain it does not
  contain, so 2,042 genomes is not the finish line, and `CAT-12` gates on the
  measured value.
  **Corrected 2026-09-27 (same day): the metric measured the forward strand
  only.** `kmers()` built single-strand k-mer sets, but kallisto indexes
  *canonical* k-mers (a k-mer folded with its reverse complement), so a panel
  genome deposited in the opposite orientation to a target scored as a miss
  when the real index would match it. `--strand canonical` is now the default
  and `--strand forward` reproduces the superseded numbers. Effect on the
  2,042-genome anellovirus panel: median leave-one-out 0.2055 → **0.2162**,
  median p_fragment 0.5097 → **0.5349**, panel 31-mers 4,888,291 → 4,843,779.
  The aggregate shift is small but the per-genome distribution is not — 350
  genomes improve, **0** degrade, 84 by more than 0.10 absolute p_fragment and
  55 by more than 0.25. The worst case, `MH649023.1`, went from 0.0329 to
  0.9514: the old metric called it essentially undetectable. The affected
  accessions cluster in the `MH648xxx`/`MH649xxx` submission block, which is
  deposited antisense to the rest of the panel. Every per-genus number below
  and in F-011 is therefore a *lower bound* until regenerated.
- [ ] `CAT-12` — set isolate counts greedily against a stated bar: **measured**
  leave-one-out P(90 bp fragment captured) ≥ 0.95 **and** zero-coverage genome
  fraction ≤ 5 %. Use the measured `p_fragment`, never the analytic
  `1-(1-c)^60` — see `CAT-11` and F-011. Measured baseline for the current
  anellovirus panel (2026-09-27): overall P = 0.5097, and per genus
  Alpha 0.788, Het 0.838, Beta 0.498, Gamma 0.332, Samek 0.255, Mem 0.247,
  Gyro 0.064 — so every genus is below the bar today and Gyrovirus is far below.
  - SARS-CoV-2: RefSeq plus one per WHO variant lineage.
  - Influenza A: one strain per relevant subtype, **all 8 segments each**, modern
    isolates beside the 1934 reference; influenza B both lineages; C and D one each.
  - Anellovirus: the 2,042 already clear P = 1.0, so the target here is the
    zero-coverage fraction and genus resolution. Add genomes for the
    under-represented genera — the 2026-09-27 held-out test put Gammatorquevirus
    at 8–40 % read recovery with no near neighbour, against 73–88 % for a
    Betatorquevirus with one at ~94 % identity.
  - HIV-1 per group/subtype; HPV, HBV, HCV, enterovirus, adenovirus by the same rule.
- [ ] `CAT-13` — **package size gate.** `src/viralscan/data/` is 8.1 MB today and
  195 of its 201 files are GTFs already gitignored (`.gitignore:99`) for size.
  Bundling sequences pushes the wheel toward PyPI's 60 MB project limit, and
  nothing in-repo would catch it: `check_ship_scope.py` is a path allowlist and
  `release.yml` has no size check. Bundle gzipped, add a wheel/sdist size
  assertion so CI fails loudly instead of the upload, and measure before `WP10`.
  If exceeded: keep Tier 2 bundled and move Tier 1 sequences to the fetch path,
  shipping the full TSV either way.
- [ ] `CAT-14` — host cross-talk gate. Every added genome is a fresh chance to
  call human reads viral. Measure per-accession host-homologous fraction for the
  whole catalogue into `host_homology_annotations.tsv`
  (`build_reference.py:768` writes it; nothing reads it — `ANDET-02`). A
  host-only negative must produce no reported call.
- [ ] `CAT-15` — gene programmes must survive the new index: the herpesvirus
  markers in `gene_programs.tsv` resolve to real index targets, and the EBV LCL
  `SRR12682296` re-run reproduces the `PROG-07` numbers (unique layer 895
  latent / 236 productive / 311 mixed / 3,094 indeterminate).
- [ ] `CAT-16` — scale check before the full build: gene count grows ~10–50×, so
  profile `analysis.py`/`detection.py` on the new GTF, and compute
  `accession_breadth` over *all* index genes of a virus rather than only
  detected ones (`ANDET-01`; `detection.py:166`, `:501-506`).

Blocking order: `CAT-01` and `CAT-03` first — without them a natively built
index has no gene structure outside the anelloviruses
(`build_reference.py:609` discards the fetched GTF; `:697`
`_genome_as_transcript_gtf` rebuilds one `{acc}_gene1` per accession), and
nothing groups segments (`virus_grouping.py:22-59` is a prefix match), so N
influenza strains would surface as 8 × N separate viruses.

### WP4H expansion correction and additions (new 2026-09-27)

Design document: [`docs/plans/2026-09-27-viral-reference-panel-expansion.md`](docs/plans/2026-09-27-viral-reference-panel-expansion.md).
It supplies the evidence for the rows above and corrects three of them.
Six failure modes it surfaces are not in any row above and are added here.

**Corrections to existing rows**

- `CAT-12` "SARS-CoV-2: RefSeq plus one per WHO variant lineage" is **not
  executable**. RefSeq holds **exactly one** SARS-CoV-2 genome (`NC_045512.2`,
  verified 3 ways via E-utilities) and **zero RefSeq records database-wide
  carry a `pango_lineage` attribute** — it is a BioSample attribute, not a
  nuccore one, and `NC_045512.2`'s flatfile contains no `pango`/`gisaid`/
  `lineage` string. The other 9.2 M sequences are unlabelled INSDC. Also
  Pango is depth-capped at 4 nodes (measured: 0 lineages deeper than 3 dots)
  while WHO counts >3,700 JN.1 descendants, so "one genome per lineage" is
  not well-posed at any depth. Superseded by `CAT-21`.
- `CAT-12` "Influenza A: one strain per relevant subtype" under-specifies. The
  axis carrying the diversity is the **clade within a subtype**, and the
  current 8 IAV RefSeq records are all 8 segments of *one 1934 lab strain*, so
  "one per subtype" still leaves H3N2 with no circulating isolate. Superseded
  by `CAT-22`. Clade names were renamed Feb 2023 (`3C.2a1b.2a.2a.3a.1` →
  `2a.3a.1`; current emergent `2a.3a.1 (J.2.4.1)` = "K") — **do not
  hard-code `3C` strings.**
- `CAT-13` is mis-scoped as a package-size problem. The full target panel is
  ≈30 Mbp against a 3.15 Gb host genome, <1 % by mass. The real build cost is
  the **host genome D-list** (~64 GB RAM, ~8 h,
  `scripts/build_genome_panel_ref.sh`), and the real risk is semantic, not
  byte count. See the document §7.3.

**New rows**

- [x] `CAT-17` — **low-complexity k-mer masking at index-build time.** A naive
  whole-genome anellovirus panel produced **1.44 % of R2 reads as false hits
  in the EBV LCL `SRR12682296`, 100 % attributable to homopolymer/tandem-repeat
  k-mers** — measured, 6,437/6,437 captured hit reads had *zero* genuine
  anellovirus k-mer after masking. 89.5 % of hit reads carried a homopolymer run
  ≥31 bp (median longest run 49, max 97). The k-mers a poly-A tail matched were
  literally `A`*31 and its near neighbours.
  **Done 2026-09-27** as a *k-mer-space* gate, which is deliberately not an
  N-masking property: our panel is 99.99 % unmasked (785 N in 9,925,822 bp)
  yet still carried **170 pure-homopolymer 31-mers across 9 records** and 7,236
  low-complexity k-mers total. The upstream hardmasked 2,023-rep panel carries
  **0** pure-homopolymer and 3 low-complexity k-mers — a 1,268× difference that
  the N-content check alone would have missed.
  Added `low_complexity_kmer_counts` (breaks k-mers down by
  `pure_homopolymer` / `long_run` / `few_bases` / `tandem`),
  `low_complexity_kmer_fraction`, `low_complexity_report`, and two gates on
  `validate_reference_records`: `max_pure_homopolymer_kmers` (absolute count,
  the shape the failure actually has — 170 k-mers over 0.06 % of the panel is
  >1 % of reads, because poly-A reads are not rare even though the k-mers are)
  and `max_low_complexity_fraction` (per-record). Both default to *off* so
  existing callers are unaffected; `build_anellovirus_reference` sets them to
  `0` when masking is requested and `2` / `0.05` under `--no-mask`.
  `write_reference_manifest` now emits `low_complexity_kmers` and
  `low_complexity_kmer_fraction` per sequence. `scripts/build_bundled_panel_ref.py`
  — the builder that produced the **shipped** covid index and which never called
  dustmasker — now writes a viral-only `viral.fa` and gates on
  `--max-pure-homopolymer-kmers` (default 0) *before* spending ~64 GB and ~8 h
  on `kb ref`. 15 new tests; 1,325 pass.
  The error text states that **a kallisto D-list cannot fix this**: a D-list
  filters host-homologous k-mers, not self-similarity inside a viral contig.
  **Open consequence:** the covid `viral_summary.tsv` files still publish
  `Alphatorquevirus` 1,167,103 / 1,605,631 UMI, and `main` @ `4fcd748` already
  closed F-005 as accession-level artifact. That retraction is a separate row
  (`CAT-30`) and is blocked pending explicit user approval.
  Same session's positive controls: EBV `NC_007605.1` = 29,207/2,000,000 R2
  reads (1.46 %), method validated by exact synthetic recovery at 1 %, 0.1 %
  and 0.02 % abundance.
- [ ] `CAT-18` — **two-index architecture.** Whole-genome pseudo-transcripts and
  real CDS transcripts for the same virus currently share one equivalence-class
  space, which is the mechanism behind the 99.8 % single-bucket failure
  (`MW455439.1_gene1`, 1,167,103/1,169,272 anellovirus UMI reported as
  "Alphatorquevirus", `ANELLO-12`). `CAT-01` fixed annotation, not the
  collision. Build `viral_gene` (real CDS only) and `viral_genome` (one
  pseudo-transcript per accession) and run two `kb count` passes on the same
  reads; write a `tier` column (`gene`/`genome`) into the v3 output schema so
  a reader can tell which space produced a number. Per-genome abundance and
  gene programmes must never read the `genome` tier.
- [ ] `CAT-19` — **`kb count` discards multimapping UMIs by default in the
  scRNA-seq path** (kallisto maintainers, pachterlab/kallisto#339), so growing a
  near-identical panel converts viral reads into *dropped* UMI rather than split
  counts. Panel expansion is therefore **not monotone-good**, consistent with the
  only benchmark that measured it (van Bemmelen et al. *BMC Genomics* 2026,
  doi:10.1186/s12864-026-12874-w: ρ = −0.73 F1 and −0.57 abundance vs set size
  for >99 %-identical genomes; the largest gain came from **geographic
  restriction**, +109 % abundance / +240 % F1, not from finer clustering).
  Concretely: use `kb count --multimapping` (not `--em`); record
  `max EC size` and `number of ECs discarded` from `kallisto inspect` on every
  build and fail on regression; and run **one experiment with kallisto
  `--distinguish`** (custom workflow, zero-indexed numeric target names) before
  committing to any per-virus genome cap. k is hard-capped at 31, so genomes
  differing by ~30–70 nt share k-mers in conserved regions.
- [ ] `CAT-20` — **anellovirus panel is probably structurally incomplete, and has
  no usable k-mer space.** `Anelloviridae[Organism]` = **42,755** nuccore
  records but genus-name queries reach only **21,283**; the ~21,472 genus-less
  legacy records (`Torque teno virus` 12,382, `Torque teno mini virus` 6,718,
  `Torque teno midi virus` 5,735, `SEN virus` 303) are structurally unreachable
  that way, and Entrez `[Organism]` matches the lineage, so a genus-derived
  collection missed them. **RefSeqViral for this family is 0**; only 178 RefSeq
  records exist and they are exactly the ICTV exemplars — everything else is
  INSDC bulk-submitted, so panel composition is currently dictated by whoever
  submitted last, not by ICTV or prevalence. Re-pull by **family taxid 687329**,
  filtering human host at the flatfile `/host` qualifier (`"Homo sapiens"
  [Organism]` returns 0 for anelloviruses and will mislead); seed with the **243
  ICTV MSL41 exemplars**; cross-check against SCANellome V2 (3,864
  representatives, Zenodo 10.5281/zenodo.7937276); add TTMV, TTMDV, **SEN virus
  (taxid 136966, outside every genus)**, and the 4 human circovirus genomes
  explicitly. Do **not** cluster (`REF-01`/`CAT-11` already showed dedup gives
  back k-mer space). Assign species by ORF1 identity to the exemplars — ICTV
  abolished genogroups, 0/2,042 records carry `/genogroup`, and the 2 that carry
  `/genotype` contradict their own organism. **Genus-level resolution is
  impossible from a nucleotide panel**: measured, **no 31-mer is shared by
  ≥1,000 of the 2,042 genomes** (max 54 k-mers reach 400–999); genera share at
  most ~44 % ORF1 *aa* identity. The only genus-sensitive route is translated
  search on ORF1 protein (kallisto `--aa`, PalmDB precedent; `--parity single`
  mandatory).
- [ ] `CAT-21` — **SARS-CoV-2 lineage policy.** RefSeq gives exactly one genome
  (§ `CAT-12` correction above). Decide: (a) ship `NC_045512.2` only, or
  (b) 4–6 genomes from INSDC/ENA with lineage assignment done locally
  (pangolin/UShER), which adds a GISAID DUA dependency and a CI
  reproducibility problem. Recommended set under (b): `NC_045512.2` (ancestral,
  carries UTR coverage later lineages lack) + one JN.1 descendant + one
  **BA.3.2** (only non-JN.1 branch in circulation; saltation with ~40 spike
  substitutions) + one XFG + one pre-Omicron VOC for re-analysis of 2021–2025
  data. Current WHO designations: VOI JN.1; VUMs XFG, NB.1.8.1, PQ.16.1.1,
  BA.3.2. **Blocked on a user decision (document §8 Q1).**
- [ ] `CAT-22` — **influenza clade-level budget**, replacing the "one per
  subtype" wording. IAV H1N1pdm09 6 (pre-2009 vs 2009; 3C.3a, 6B, 6B.1, 6B.2);
  H1N1 non-pdm09 2; H3N2 8 (3C.2a, 3C.2a1, 3C.2a1b.1a, 2a, 2a.1, 2a.3,
  2a.3a.1, 2a.3a.1 (J.2.4.1)); H5N1 6 (2.3.4, 2.3.4.4, **2.3.4.4b**, human-case
  genotypes B3.13 and D1.1); H7N9 3; H9N2 4; H1N2 2; H2N2 1; IB Victoria 4;
  **IB Yamagata 1 flagged `retired_from_surveillance`** (no confirmed detection
  since March 2020, dropped from NH 2026-27 CVVs, still in the Sept 2026
  Southern Hemisphere list and an LAIV component); IC and ID 1 each.
  **8 segments per strain, always.** ⚠ **Type by sequence, not organism string** —
  `OP212288` (A/Texas/61/2022, first dairy-cattle human H5N1, clade 2.3.4.4b
  B3.13) is deposited with `(H3N2)` in its organism field. Verified B3.13 =
  `PQ468757`–`PQ468764`; D1.1 = `PQ573551`–`PQ573557`. Clustering:
  **within-subtype only, 95 % nt / 85 % coverage** (`cd-hit-est -c 0.95 -aS
  0.85`) for RNA viruses, 98 % for adenovirus. Note CD-HIT was *excluded* from
  the van Bemmelen benchmark for runtime, so no benchmark-derived CD-HIT
  threshold exists; 95/85 is the closest documented viral-derep default
  (`votuderep` ANI ≥95 % + coverage ≥85 %, the CheckV method). The
  directionally relevant caution is PalmDB *Nat Biotechnol* 2025: 99 % aa
  clustering made **67.4 % of true taxa undetectable** vs 3.3 % unclustered, so
  they grouped by taxonomy instead — which is what `CAT-03` now does.
- [ ] `CAT-23` — **Serratus reactivation panel as an explicit acceptance
  criterion.** The "Serratus activation screen" is **Lareau CA, …, Satpathy AT.
  "Latent human herpesvirus 6 is reactivated in CAR T cells." *Nature*
  623(7987):608–615 (2023), doi:10.1038/s41586-023-06704-2, PMID 37938768,
  code `github.com/caleblareau/serratus-reactivation-screen` — *read*-based
  (Serratus petabase), not host-signature-based. Input panel 129 curated human
  viruses (ViralZone / Hulo 2011); reactivation criterion = DNA-genome-only,
  known latent cycles that reactivate *in vivo* (Traylen et al. *Future Virology*
  6:451, 2011), narrowed to 17 DNA viruses from Herpesviridae,
  Polyomaviridae, Adenoviridae, Parvoviridae; positive call = **≥100 reads AND
  ≥50 % mean mapped identity** per sample–virus. **Corrected 2026-09-27: the
  anchor is 17, not 18.** The in-repo source table settles it —
  `sources/viral_panel/Serratus_hits_all_viruses.tsv` has 129 rows and exactly
  **17** with `reactivation_candidate == TRUE`; **HHV-6B is not among them**
  (only `NC_001664`, HHV-6A). The upstream "17" was right and an earlier note in
  this row that said 18 was an enumeration slip on my part. **All 17 are
  already in the panel** — T2's latent-DNA
  core needs no new accessions, only real gene structure and strain diversity,
  so it becomes the *acceptance criterion* rather than new breadth. Make it a
  checkable row: `A1` in document §7.1. **T2 must be a superset, never a copy:**
  the Serratus panel contains **no SARS-CoV-2**, only 3 HPV genotypes, filed
  HHV-6 as HHV-6A only (needing a separate HHV-6B query), and TTV appears with
  206 incidental hits never flagged reactivation-relevant.
- [ ] `CAT-24` — **expected-negative arithmetic for anelloviruses and other
  low-prevalence latent agents.** TTV in an immunocompetent adult is ~10²–10³
  copies/mL plasma (*Viruses* 17(2):140, 2025, PMID 40143262) with >90 %
  prevalence (PMID 12721794), but scRNA-seq is 3′-biased and polyA-selected,
  TTV is ssDNA, and the panel has no conserved k-mer space. Measured on the EBV
  LCL: **EBV 1.46 % of R2 reads, zero genuine anellovirus reads.** So report
  anellovirus as presence-only, binary, `screening_only` (`REF-10`); attach the
  expected-value calculation to every negative instead of publishing a bare
  zero (the discipline already required for tonsil EBV, `SENS-06`); add a
  `detection_bound` column; and **never publish per-genotype anellovirus
  abundance** — 87 % of annotated genomes (1,740/1,995) still carry a single CDS
  spanning the whole genome and are competition buckets (`ANELLO-12`).

### WP4H external-asset adoption and adversarial corrections (new 2026-09-27)

**Vendored external references.** 41 files, 30.2 MiB, 0 failures, every
repo pinned by commit, under `extras/vendor_sources/` with
`VENDOR_MANIFEST.tsv` (per-file bytes + sha256) and `VENDOR_REPORT.md`.

| Source | Pin | What it gives us |
|---|---|---|
| `clareaulab/human_anellovirus_pangenome` | `3ed77e19` | **2,023 hardmasked reps** (5,927,006 bp, 284,762 N = 4.80 %, 0 pure-homopolymer k-mers), `anello_t2g.txt`, `anello_for_kallisto.gtf`, `simple_anello_metadata_V2.csv` (**3,545** rows), `vclust/{genus,species,clusters}.tsv`, ORF1 `.contree`/`.treefile`/`.faa`, and `example/expected_output/full_SRR32170409/` |
| `caleblareau/pan-viral-reactivation` | `74136de5` | `pan_virus_annotation_plain.tsv` (**724** rows, `EC/Gene/Nuccore/Virus`, 9 viruses) with standard nomenclature (`NC_006273.2`→`RL1`/`UL1`…/`US27`/`TRS1`; `NC_007605.1`→`BNRF1`/`EBNA-1`/`LMP-1`/`BARF0`/`BGLF1`…), `pan_virus_cds.fasta`, an HIV reference, and **6 real `*.kb.txt` outputs** as regression fixtures |
| `clareaulab/ad-hsv-mapping` | `53801369` | `HSV1-LATonly.fasta`, `HSV1-coding.fasta`, `VZV-coding.fasta`, `jg_NC_001806.2.gtf` |
| `yyoshiaki/VIRTUS3` | `b7873791` | `NC_007605.1_CDS_EBER12.fa` — 96 sequences = **94 CDS + EBER1 + EBER2** (named `rna-HHV4_EBER-*`) |
| `huangyh09/ViralScan` (Apache-2.0) | `d8279c37` | `Viral_GTF_maker.py` + `viral_reference/viruses_833.fasta` — **the provenance of our 195 GTFs** |

**Accession verification** against live NCBI/ENA/Datasets:
`extras/vendor_sources/ACCESSION_VERIFICATION.md` — 11 SARS-CoV-2 verified
(5 of 7 requested lineages + Wuhan-Hu-1), ~95 influenza segment accessions over
12 complete isolates, **9/9** endemic coronaviruses, **11/11** ssDNA viruses.

**Adversarial corrections to this work package (2026-09-27).** A review of my own
design found four things wrong, and two of them had already been found
independently upstream in this repo:

1. **Retracted — "genus queries reach only 21,283 of 42,755 records".** True of
   Entrez's `[Organism]` index but irrelevant here: the panel is already
   human-host filtered, and all 16 ICTV genera absent from the panel have
   **zero** complete genomes with `"Homo sapiens"[Host]` (Lambda 408/0, Eta
   256/0, Iota 139/0, Kappa 24/0 — swine, feline, canine, tupaia, pinniped).
   The ICTV 37-genus refresh is **deleted** from `CAT-20`; it buys no
   sensitivity and adds false-positive surface.
2. **Retracted — "our genus labels are unreliable".** Backwards. The 582 panel
   genomes with no GenBank genus are the ~584 the upstream ORF1 phylogeny
   **resolved**; our labels are better than GenBank's. **Do not reconcile them
   against GenBank lineages.**
3. **Retracted — my own collision measurement.** I grouped `t2g` by accession
   prefix and reported one accession carrying both real genes and a placeholder.
   Wrong key. F-013 / session `2026-09-27-001` found a second duplicate that
   also fixes *how* `CAT-05` must work: `AB303562.1` and `NC_038359.1` are the
   GenBank and RefSeq copies of one Gammatorquevirus genome, **both shipping**,
   an upstream dereplication failure. **The guard must key on sequence, not
   accession** — and it already does; see the note on `CAT-05` below.
4. **`OP212288` is not a mislabelled H5N1.** I reported it as H3N1 deposited
   `(H3N2)`. It is genuinely H3N2 (87.7 % identity to H3N2 vs 42.4 % to H5N1,
   567 aa CDS, monobasic HA0 cleavage site); the serotype lives in the
   DEFINITION isolate string. The real hazard is a **name collision**: both it
   and the A/Texas/61/2022 dairy-cattle H5N1 case are "A/Texas/61/2022". Type by
   sequence, not by organism *or* isolate string.

**Measured this session, `CAT-05` evidence.** The new gate rejected **both** the
upstream 2,023-rep panel and the shipped 2,312-record panel on the same
duplicate: `NC_038359.1` / `AB303562.1`, byte-identical sequences, present
upstream. This is an **upstream** dereplication failure we inherited, and it is
a live `kallisto index` hazard — `CAT-05` records a previous duplicate
(`NC_002076.2`) breaking the build.

- [x] `CAT-05` — **duplicate guard keyed on sequence, on every build path.**
  `validate_reference_records` already hashed sequences, and the new
  low-complexity gate made it fire in practice: run against
  `extras/vendor_sources/anello/ref/hardmasked_*.fa` it rejected the panel on
  `NC_038359.1` / `AB303562.1`. Partially landed with `CAT-17` — the guard is
  proven to catch a real duplicate. **Remaining:** decide the resolution
  (drop the RefSeq copy, or keep one and record why) and confirm
  `scripts/build_bundled_panel_ref.py` reaches the guard on every path.
- [ ] `CAT-25` — **pinned environment.** `environment.yml` declares
  kallisto 0.50.1 / bustools 0.43.2 / kb-python 0.28.2, but the running env is
  **0.51.1 / 0.45.1 / 0.29.5** — so every `CAT-11` number was measured off-pin.
  `blast=2.16.0` is declared (it provides `dustmasker`, **absent from PATH**) but
  `conda-recipe/meta.yaml` declares neither blast nor cd-hit. Decision: pin down
  to the declared versions and re-measure the capture bars and EC sizes.
- [ ] `CAT-26` — **freeze the current reference before the panel moves.**
  `docs/review-clear-execute-plan.md` reruns the EBV and HIV controls *"from
  identical FASTQs and the original combined reference"* and freezes reference
  hashes. A new panel must not perturb that arm; freeze a copy and treat the new
  panel as a third arm.
- [ ] `CAT-27` — **EBER1/EBER2 audit.** VIRTUS3 ships an EBV reference of 96
  sequences = 94 CDS + **EBER1 (167 nt) + EBER2 (173 nt)**, named
  `rna-HHV4_EBER-*`. Ours may omit the two non-coding RNAs, which are the
  highest-abundance latent EBV transcripts and arguably the best latent marker.
  Check, and add from RefSeq if absent. Cross-check already passes: VIRTUS3's
  94 CDS, `HSV1-coding.fasta` 77 and `VZV-coding.fasta` 73 exactly match the
  panviral HHV4/HHV1/HHV3 gene counts — three independent repos built against the
  same RefSeq release.
- [ ] `CAT-28` — **gene-nomenclature mapping.** Map
  `pan_viral_annotation_plain.tsv` (724 rows) into `gene_programs.tsv` and
  validate the markers. **Data hazards:** the file is EC-keyed so genes repeat
  (`BWRF1` ×6, `LMP-1` ×3, `US33A` ×3, `AAV2gp06` ×2, `K14` ×2), and **line 532
  is corrupted** — `530  Jvgp6  [NC_001699.1  Jcpolyomavirus` with a stray `[`.
  A join on accession drops JCV `Jvgp6` **silently**. Also records accession
  disagreements: HHV7 `U43400.1` (we use `NC_001716`), HHV8 `MK733606.1` (we
  use `NC_009333`), and HHV6B `AF157706.1` — the pseudocontig behind our 97
  "records" that are not genomes.
- [ ] `CAT-29` — **HSV-1 latency transcripts → Tier 2** (user decision).
  `HSV1-LATonly.fasta` is the thing `gene_programs.tsv` calls `partial` because
  "HSV-1's latent state is unreachable by construction". **But it is not
  spliced and not contiguous**: 3 separate records — exon 660, **intron 1,956**,
  exon 5,731 — all ICP0/LAT at 118805-127151. The intron is independently
  k-mer-countable, so despite the name this does *not* restrict to LAT mRNA and
  will add false-positive surface. Re-cut the exons before adopting, and
  re-evaluate HSV-1 `panel_completeness` afterwards.
- [ ] `CAT-30` — **retract the covid Alphatorquevirus claim.** `main` @ `4fcd748`
  closed F-005 — *"All 8 detected anellovirus accessions appear in Phase A
  (multi-chromosomal GRCh38 alignment) — accession-level artifact confirmation"*,
  *"no genuine viral infection in these COVID PBMC samples"* — but
  `covid_viralscan/results/*/results/viral_summary.tsv` **still publishes**
  `Alphatorquevirus` 1,167,103.0 UMI (x213-g, x216-g.jul2), 1,605,631.0
  (x216-g), 17.63–19.05 % of cells. The science is settled; the retraction never
  happened. **Blocked pending explicit user approval** — do not edit published
  results unilaterally. `CAT-17` makes this more urgent, not less: the shipped
  panel is the one carrying 170 pure-homopolymer k-mers.

**Sequencing note (superseded 2026-09-27):** `CAT-17` shipped. Next is
`CAT-30` (user approval), then the `CAT-05` duplicate decision, then
`CAT-20`'s adoption of the 2,023 reps — which stays **gated** on the F-013/F-011
reconciliation the previous session opened, per its own instruction *"do not act
on restore-the-1,522 until this is settled."* The full dependency-ordered
sequence is in document §9.

## WP4I — Complete latent/lytic state calling (new 2026-09-27)

Objective: make `gene_programs` biologically correct and measurable. Continues
`PROG-08`/`PROG-09`.

- [x] `PROG-10` — fix vacuous tests first: `test_measured_ratio` passes for any
  counts, and one EBV-regression assertion checks the opposite of its docstring.
  Done 2026-09-27. The directional guard was `not (state == productive and
  selected_state == latent)`, unreachable once `state == latent` is asserted two
  lines earlier; it now forbids the direction its comment names (allocated
  `productive`, unique `latent`), and
  `test_directional_guard_fires_when_allocation_is_productive_only` proves that
  combination is reachable, so the guard can actually fail. The ratio test only
  checked that the catalogue has at least as many productive groups as latent
  ones, whatever the masses; it is replaced by `test_calls_depend_on_breadth_not_mass`,
  which holds the support fixed and moves mass from 1.0/1.05 to 1000/0.01 and
  requires an identical call. Neither test reproduces the real-run numbers; that
  remains `PROG-07`'s job after `PROG-11`.
- [x] `PROG-11` — catalogue biology, each change verified against primary
  literature before editing: a kinetic-class column (`latent` /
  `immediate_early` / `early` / `late`); CMV UL122/123 and HHV-6A/7 U90/U86 are
  immediate-early, not latent; EBV `BaRF1.1`, `BHRF1`, `BNLF2a/b` are early lytic;
  KSHV's latent set lists ORF16 (vBcl-2, lytic) and `partial` is hand-set despite
  three independent latent groups; exact (not case-insensitive) name matching in
  `extras/build_gene_programs.py:888`, `:896`; `panel_completeness` derived by rule.
  Done 2026-09-27, each with a primary citation: EBV `BARF1` removed (latent only
  in epithelial cancers — PMID 32708965, 39329759), which also removes `BaRF1.1`
  (the lytic ribonucleotide reductase) that the case-insensitive match had pulled
  in; matching in `_resolve` is now case-sensitive, and regenerating dropped no
  other marker (79 → 77 rows, the diff is exactly the five intended changes).
  HCMV `UL122`/`UL123` are immediate-early, not latent, and HCMV is now
  `partial` with `latency_observable_in_rna=false`: single-cell HCMV latency has
  no restricted latency programme but mirrors a late-lytic one at much lower
  levels (Shnayder et al. 2018, PMID 29535194). KSHV `ORF16` is vBcl-2, lytic,
  not vGPCR (ORF74; PMID 20860481). The docs' BARF1.2 showcase numbers carry a
  dated correction. **Deviation:** "`panel_completeness` derived by rule" is
  dropped — HCMV has two independent latent anchors (`UL138`, `UL111A`), so the
  rule would call it complete while the biology says latency is unobservable;
  the facts stay hand-set, each with its reason. **Still unverified, so not
  edited** (literature search budget ran out): EBV `BHRF1`, `BNLF2a/b` as early
  lytic; the `BcLF1` and KSHV `ORF17` notes; KSHV additions `ORF72`, `ORF71`,
  `K12`, `K15`, `K8`, `K8.1`, `ORF57`, `ORF26`, PAN and `K1`'s class; HHV-6A/7
  `U90`/`U86`; HHV-6B `U95`; HHV-7's latency set.
  - Done 2026-10-02 (`v3/prog11`, c004f47; user confirmed both decisions).
    Sources were retrieved through Europe PMC.
    - **New columns:** `kinetic_class` (latent, immediate_early, early,
      leaky_late, late, unclassified) and `kinetic_pmid`. `unclassified`
      means no source was retrieved.
    - **EBV:** BHRF1, BNLF2a and BNLF2b are removed as markers; they are early
      lytic (29864140), and BNLF2a/b overlap LMP-1. This follows BARF1.
    - **HHV-6A, HHV-6B, HHV-7:** HHV-6A U90/U86, HHV-6B U95 and HHV-7 U90
      move to productive/IE. HHV-6A and HHV-7 are now `partial` (latency not
      observable in RNA).
    - **Note fixes:** KSHV ORF17 (protease, not MTA), VZV ORF4 (IE4, not
      IE62), EBV BcLF1 (major capsid protein).
    - **Result:** the catalogue goes from 77 to 74 rows. A new test forbids
      mixed-programme overlap groups.
    - **Flagged, out of scope:** EBV BBLF4/BBLF1/BGLF4/BALF5 notes contradict
      the CAGE table; KSHV ORF71/72 etc. were never in the catalogue.
- [ ] `PROG-12` — states `latent` / `reactivating` (immediate-early only) /
  `productive` / `mixed` / `indeterminate`; symmetric breadth thresholds (latent
  needs 1 group today, productive 2); a per-marker UMI floor. (Not "honour
  `non_overlapping`": `_breadth` ignores it deliberately — the overlap group is
  the unit.) HCMV needs its own handling here: because latency mirrors
  low-level late-lytic expression, a presence-based `productive` call is not
  specific either, so it needs a per-cell quantity threshold or
  `not_applicable`.
- [ ] `PROG-13` — merge 31-mer-identical repeat copies (HSV LAT/ICP0/ICP4 in
  TRL/IRL, VZV ORF62/ORF63 in TRS/IRS) to one gene_id in t2g, so they reach the
  unique layer.
- [ ] `PROG-14` — unresolved markers fail loudly instead of a log line
  (`scripts/gene_programs.py:180-187`); overlap groups from the active index's GTF.
- [ ] `PROG-15` — measure before adding antisense latency transcripts (VLT, LAT
  intron, LUNA): strandedness (`kallisto bus` runs with no strand flag), the share
  of EBV LCL reads outside annotated exons, and which markers are real index
  targets (the KSHV GTF has 26 exon rows for 96 genes).
- [x] `PROG-19` — **closed 2026-10-03.** `summarise_programs` skips the
  placeholder row when counting, and a no-model virus keeps only its "no
  programme model" caveat. Checked on the cat42b EBV LCL copy: no-model
  viruses report 0 and EBV still reports 932. Original text: (found 2026-10-03 with PROG-17) every `not_applicable`
  virus shows `n_cells_total=1` in `gene_program_summary.tsv`: the count is
  the placeholder row (`barcode=""`) that `run_one` emits. Next to
  `n_called_cells` it reads as one real cell. It should be 0. Not started.
- [ ] `PROG-18` — (found 2026-10-03 with PROG-08) validate the KSHV latent
  call on a real dataset: a public PEL or KS scRNA-seq set (BCBL-1/BC-3 PEL
  lines are latently infected). Expect mostly latent, with a lytic minority.
  Not started.
- [x] `PROG-17` — **closed 2026-10-03.** Detection writes
  `results/called_cells.tsv`; layer 2 scores only those barcodes and adds
  `n_called_cells` to the summary. Pre-PROG-17 run dirs fall back to the run's
  own `emptydrops_cells.tsv` when its method was emptyDrops, else re-call
  cells with the run config. A list that does not match the H5AD fails closed.
  The file is not a Snakemake `output:`, so resuming an old run does not
  re-run detection. Tests: `TestCalledCellsFile`.
  Original text: (found 2026-10-03 during PROG-07) layer 2
  (`rerun-programs` / `--gene-programs`) scores every cell in the multimap
  h5ad, not the called-cell set. On the cat42b EBV LCL run it scored 1,679
  cells, and only 932 of them are among the 2,763 emptyDrops-called cells
  (`per_cell_viral.tsv` `is_called_cell`). Restrict layer 2 to called cells,
  or report both sets, then rerun PROG-07. Not started (no new tasks in the
  2026-10-02 finish-pass).
- [ ] `PROG-16` — lytic acceptance test: KSHV `GSE190558` (`RUN-04`), induced vs
  uninduced; HSV-1 `SRR8315713` expected productive with no latent calls.

## WP4J — SFL tonsil TTV and HPV screen (new 2026-09-27)

Objective: answer whether torquetenovirus (TTV) and HPV are detectable in the
SFL tonsil CITE-seq pool (24 hashtagged donors, `s1`–`s24`), without publishing
an F-005-type host-homology call.

Checked before starting:
- **Only `x223` is gene expression.** `x225` is the antibody (ADT) library,
  according to cellranger `config.csv`; its R2 reads carry tag structure, not
  cDNA.
- **`x223` is 10x 5′ v3 R2-only, not 3′.** An earlier `.living` note said 3′.
  In this chemistry R2 is antisense to the transcript, and reads cluster near
  the transcription start, not the polyA site.
- The cellranger reference was GRCh38-2024-A.
- The existing index (`panel_ref_genomic`) carries only three HPV types: HPV16,
  HPV1 and HPV2. HPV18/31/33/45 are not in it. The opening note said four; the
  fourth accession, `NC_003461`, is not a papillomavirus.
- A native run is not safe yet:
  - kb's bundled 10xv3 whitelist is the 3′ list.
  - kb passes no strand flag unless given one, and `menu.py` exposes none.

- [x] `TONSIL-01` — host-subtracted screen from the existing cellranger BAM:
  1. Take the unmapped GEX reads (ADT reads excluded) and keep CB/UB.
  2. Align them with minimap2 to the panel's viral genomes plus the 16 WP4E HPV
     types.
  3. Re-check every candidate read against GRCh38.
  4. Remove UMI duplicates. Attribute cells to donors through the cellhashr
     singlets.
  5. Add a positive-control plant in 5′ geometry.

  Reads that aligned to the host never reach the unmapped set, so F-005 is
  excluded by construction.

  Decision rules, fixed before looking at results:
  - **TTV screening-positive:** ≥ 3 UMIs in ≥ 2 singlet cells of one donor. The
    call is family level unless the genus is unambiguous.
  - **HPV positive:** ≥ 3 UMIs on one type, with ≥ 1 read in URR/E6/E7.
  - **Anything else:** "not detected at this depth", never "absent".
  - The F-005 breadth gate (≤ 3.41 %) is not used to reject. 5′ capture
    concentrates true reads near the start site by design.
  - Every anellovirus result stays `screening_only` (`REF-10`).
  - Committed outputs name donors only as `s1`–`s24`.

  **Result (2026-09-27): no TTV and no HPV in any of the 24 donors — "not
  detected at this depth".** Full write-up: `.living/findings/` F-010.

  | Stage | Reads |
  |---|---|
  | Unmapped GEX reads with a valid barcode | 118,856,604 |
  | Pass the trim/complexity prefilter | 24,969,278 |
  | Any viral alignment | 34,208 |
  | Pass the viral filters | 1,511 |
  | Survive the host re-check | 23 |

  - **Anellovirus:** 1,852 raw hits, all low-complexity partial alignments
    (20–49 bp). 1,485 of them are 28-bp matches to a CAG trinucleotide repeat;
    the rest are poly(A) plus TSO-rc. None survive.
  - **HPV:** 2 raw hits, both shorter than 50 bp.
  - **The 23 survivors are not infections:**
    - HCV reads all fall on the 3′-UTR poly(U/UC) tract, at 9435–9505.
    - The macaque *Cercopithecine herpesvirus 2* reads fall on two GC-rich
      positions shared across donors.
    - HSV-1: 2 reads at one position, in one hashtag-negative barcode.
  - **HHV-6B:** 5 host-free reads at 5 genome positions and 5 barcodes. None of
    those barcodes is a singlet, so the reads cannot be attributed to a donor.
    This is below any call and is recorded, not claimed.
  - **Positive-control plant:** 1,000 reads each for HPV16 and TTV, into 1 M
    real reads. Recovery was 99.2 % for HPV16 and 99.4 % for TTV, all
    attributed to `s1`. No planted read was lost to the host re-check.
  - **Held-out anellovirus plant (sensitivity to strains not in the
    database):**
    - One genome per genus was removed from both references, and 1,000 reads
      were planted from each.
    - Recovery was 88 % (5′ window) and 73 % (uniform) for Betatorquevirus.
      Its nearest remaining genome is 94 % identical.
    - Recovery was 16 % and 35 % for Alphatorquevirus, nearest 84 %.
    - Recovery was 40 % and 8 % for Gammatorquevirus, which has no asm20
      alignment to any remaining genome.
    - **So the TTV zero rules out strains close to the 2,042 references. It
      does not rule out low-level divergent TTV.** Reads from those strains
      mostly fail to seed or fall below 85 % identity.
  - **Circular-record check:** 1 of the 1,853 dropped hits sits at a record end.
    It is G/C-run sequence, not a read spanning the origin (42 of 90 bases match
    across the junction).

  Non-obvious points:
  - **Prefilter.** Without it, 1 M unmapped reads gave 486,044 viral "hits"
    and 72 M alignments. These were TSO/poly(T) reads on the HCV poly(U) tract,
    poly(C) reads on EMCV, and poly(A) reads on A-rich anellovirus regions. The
    job ran out of memory.
  - **Missing reads.** The per-sample BAM holds only reads with a valid
    barcode. That is 1.35 B of the 1.61 B GEX reads, so the other ~16 % were
    never screened.
  - **HPV16 coordinates.** `NC_001526.4` is linearised at E1, so p97 is at
    position 7139, not 97.
  - **Limits of the result.**
    - The plant used genomes that are in the database, so this is not an LOD.
      Sensitivity to divergent anelloviruses is `ANDET-08`.
    - The EBV zero was predicted: E[EBV+ cells] ≈ 0.02–0.15.
    - The HPV prior in benign tonsil suspensions is low.
    - So none of these zeros is an informative negative (`SENS-06`).
    - BLAST spot-checks were skipped (no BLAST on the cluster). The competitive
      minimap2 re-check against GRCh38 plus the viral set stands in for them.
- [ ] `TONSIL-02` — native 5′ support so `viralscan` itself can run this library:
  - a `--strand` option passed to `kb count` (overlaps `PROG-15`);
  - the cellranger cell barcodes as `-w`;
  - `-x 0,0,16:0,16,28:1,0,0`.

  Measure on a 1 M-read subsample first: host mapping rate under forward,
  reverse and unstranded, then `viralscan check-whitelist`. The native counts
  must agree with `TONSIL-01`.

## WP4K — Catalogue↔index reconciliation and targeted reference adds (new 2026-09-27)

Objective: stop treating "catalogued" as "detectable". F-015 measured that the
built reference omits most of the catalogue, so the first task is not importing
genomes — it is making the build honest about what it contains.

- [x] `CAT-31` — **catalogue↔index reconciliation guard (highest value).**
  F-015: 13 of 16 catalogued HPV (incl. **HPV18, HPV31**) are absent from
  `viral_genome.dedup.fa`; **all 5 Retroviridae** (both HIV-1 and HIV-2) are
  absent; plus 7/30 influenza, 4/7 coronaviruses, 3/8 polyomaviruses. A
  catalogued-but-unindexed virus is *undetectable*, so this silently caps
  sensitivity. Add a build-time assertion in
  `scripts/build_bundled_panel_ref.py` that every catalogued accession intended
  for detection appears in the emitted FASTA, and emit an explicit
  `catalogued_not_indexed.tsv` rather than failing opaquely. Decide per family
  whether a miss is intentional (e.g. segment-only, partial CDS) and record the
  reason. **Do this before any import** — otherwise `CAT-32`/`CAT-33` reproduce
  the same invisibility at larger scale.
  - **2026-09-29 (`MECH-A`):** the guard's scope is now the catalogue's `panel`
    column. The catalogue grew to 4,128 rows so the Virus Identity table can
    name max-panel genomes. Only the 2,345 `panel=shipped` rows are detection
    targets of the shipped panel, so the 1,783 `max` rows do not turn into
    reported gaps.
  - **2026-09-28: the guard shipped; the 34 decisions did not, and that is the
    point.** `reconcile_reference_panel()` in
    `src/viralscan/scripts/build_reference.py` compares the assembled panel FASTA
    against `virus_catalog.tsv` on the version-stripped, underscore-normalised
    base accession and writes `catalogued_not_indexed.tsv`
    (`accession/family/species/status/reason`) as new Step 7/8, before `kb ref`
    so a low-sensitivity panel fails in seconds rather than after ~64 GB and
    ~8 h. Misses are `intentional` only if listed in the new
    `src/viralscan/data/index_exclusions.tsv` (`accession/reason/decided_by`),
    else `unexplained`. 45 tests in `tests/test_index_reconciliation.py`.
  - **`index_exclusions.tsv` is shipped deliberately empty.** All 34 misses are
    reported `unexplained` and the guard is red. Pre-allowlisting them would
    turn a red build green without adding a genome — the exact laundering this
    row exists to stop — and the 10 GenBank-only HPV in that set are the same
    high-risk types as the 3 RefSeq ones, so "we only index RefSeq" is not a
    real policy (2,027 non-RefSeq anelloviruses are indexed). The allowlist is
    where a reviewer records "we accept this loss"; the report is the work-list.
  - **Strict is opt-in** (`--strict-reconciliation`), off by default: the
    catalogue is 11× the bundled GTF set, so a default-on gate makes the build
    unusable rather than honest, and the report is written either way. CI turns
    it on. Unexpected records (indexed but uncatalogued) and stale allowlist
    entries are reported too, so the allowlist cannot rot into a blind spot.
  - **Gotcha for the follow-up:** `build_reference.py` is a governance-pinned
    artifact, so editing it makes `test_artifact_inventory` /
    `test_claim_registry` fail until the repo's usual 3-commit re-pin chain
    (`fd8cee1` → `be1c196` → `c31742d`) is replayed against the new commit.
  - **Left open:** all 34 rows of `catalogued_not_indexed.tsv` still need a
    genome or a named decision. `CAT-32`/`CAT-33` are the first consumers.
- [x] `CAT-32` — **add EBV type 2 (`NC_009334.1`).** F-015: the index carries
  only `NC_007605.1` (type 1 / B95-8). Measured 31-mer overlap with type 2
  (`NC_009334.1` / AG876) is **79.9% shared**, leaving 20.1% type-1-unique and
  20.4% type-2-unique — a type-2 infection is currently uncallable. One
  accession, enables EBV *typing* rather than only detection. Verify against
  `CAT-27`'s EBER audit so the two EBV entries do not double-count. *(Note: an
  earlier subagent report put the type-1-unique fraction at ~77%; that is the
  shared fraction. The 20% figure is the measured one — do not re-cite 77%.)*
- [x] `CAT-33` — **HPV expansion, sourced from RefSeq `NC_` not VIRTUS2 `lcl/`.**
  F-015: VIRTUS2 carries 92 HPV *types* vs our 16 → **76 genuinely new types**.
  The apparent "94 new accessions" is partly an artifact: VIRTUS2 stores HPV as
  `gi|…|lcl|HPV##REF.1`, which will never accession-match our `NC_` records even
  for the same type — so compare by **type**, not accession. `lcl|` records are
  RefSeq *local* submissions with weaker curation; source the same 76 types from
  curated `NC_` accessions so we get GTFs. Note this is largely redundant with
  `CAT-31`: fix the 13 already-catalogued HPV first, then decide how far to go.
- [x] `CAT-34` — **record the closed negative: no `chrEBV` in our host build.**
  VirDetect warns that some hg38 builds ship a `chrEBV` contig, which would make
  EBV silently invisible to host subtraction. F-015 checked both real D-list
  inputs — `refdata-gex-GRCh38-2024-A/fasta/genome.fa` (194 contigs) and
  `Homo_sapiens.GRCh38.116.gtf.gz` — and found **zero** viral contigs. Closed:
  the hazard does not apply. Re-run the check whenever the host reference is
  swapped rather than assuming it.
**OUTCOME 2026-09-28 — the reference is built. `CAT-31`…`CAT-35` shipped; see F-016.**

Final panel: **2,343 genomes** (323 GTF-backed + 2,020 anellovirus) from
`scripts/build_bundled_panel_ref.py` on the pinned toolchain. The reconciliation
guard ran inside that build and closed **32 of the 34** catalogued gaps; the only
two survivors are the CAT-05 duplicate pair, which is the correct outcome.

What the additions actually bought, and what it cost:

- **CAT-32** — EBV type 2 is indexed. EBER1/EBER2 are now emitted too: they were
  being dropped by a `key == "CDS"` filter in `ncbi_fetch.py`, so the highest-
  abundance latent EBV transcripts were unreachable. The regenerated EBV GTF has
  96 genes (94 CDS + 2 EBER), matching VIRTUS3's independent 96-record reference.
- **CAT-33** — HPV went from 16 catalogued types to **109 human-pathogen types**
  (all 15 IARC group-1 covered). The multimapping fear was **measured and
  disproven**: 98.68 % of k-mer space is type-discriminating, HPV16/HPV18 share
  *zero* 31-mers, and L1 is *more* type-specific than non-L1. But only 68/182
  types are curated RefSeq; **114 are INSDC-only, including 8 of the 15 group-1
  types**, so the catalogue cannot treat HPV rows as curated.
- **CAT-31** — the guard works and is wired before `kb ref`, so a lossy panel
  fails in seconds instead of after ~64 GB. Its allowlist ships *empty* on
  purpose: pre-allowlisting the 34 gaps would have turned a red build green
  without a genome being added.
- **CAT-34** — closed negative recorded; the D-list genome has 194 contigs and
  zero viral, so the `chrEBV` hazard does not apply to this build.

Three latent bugs found and fixed on the way, all in F-016: **20 duplicate
`transcript_id`s across the bundled GTFs** (a live `kb ref` crash, 301
occurrences namespaced), the EBER `misc_RNA` drop, and the CAT-17 gate itself —
**dustmasker masks only 0.01 % of this panel and cannot remove the
pure-homopolymer k-mers the gate fails on**, so the builder now masks runs of ≥31
identical bases directly (panel-wide pure-homopolymer 31-mers: 85 → 0 in the
worst record, 0 panel-wide) and the unachievable `0.0` fraction default became
`0.05`, with the absolute homopolymer gate as the real control.

**New rows, from what the build exposed:**

- [ ] `CAT-36` — **the catalogue is now behind the panel by 96 rows.** The build
  reports 96 genomes indexed but absent from `virus_catalog.tsv` (the INSDC HPV
  set, plus `AF157706.1` HHV-6B). They have no provenance, family or tier, and
  `CAT-31`'s reverse check will keep reporting them. Extend the catalogue, and
  propagate `source` / `refseq` / `oncogenic_class` for the HPV rows.
- [ ] `CAT-37` — **inherit a format-version stamp for the GTF cache.**
  `ncbi_fetch._cache_valid` only checks a file against its own sidecar, so the
  2,249 already-generated GTFs are silently reused and the EBER fix does not
  reach them until they are deleted. A cache-key change is the durable fix.
- [ ] `CAT-38` — **two whole-genome pseudo-transcripts entered the panel with
  the EBER fix.** HAV `NC_001489.1` (misc_RNA spans 100 % of the genome) and
  HTLV-2 `NC_001488.1` (94 % of 8,952 nt) are exactly the shape that collapsed
  99.8 % of anellovirus UMI into one bucket. Belongs to the `CAT-18` two-tier
  split; do not ship these two in a single-index panel without deciding.
- [ ] `CAT-39` — **`M12737` (HPV-8) and `NC_039089` (HPV-71) are thin.** Neither
  has usable CDS upstream, so they fall back to whole-genome pseudo-transcripts
  (`CAT-38`). Re-cut or exclude deliberately rather than by accident.
- [ ] `CAT-40` — **run `kallisto inspect` on the new index** (max EC size,
  discarded EC count) before trusting it. The k-mer-overlap metrics in F-016 are
  computed outside kallisto and cannot see the pseudo-inverse threshold effect,
  which is `CAT-19`'s actual failure mode.

- [x] `CAT-35` — **VIRTUS2/VirDetect/VIRTUS3 method extraction (research, done
  2026-09-27).** Conclusion recorded in F-015: take the **method**, not the
  genomes — host subtraction before viral quantification, per-strand counting,
  explicit multimap handling. Reject VIRTUS2's 37 influenza (legacy lab strains:
  PR8, H9N2 HK/97, H5N1 goose 1996, H3N2 NY/2004, B/Lee/1940) as no better than
  what we hold, and its 27 anelloviruses as already covered with better
  provenance by the upstream 2,023 reps. VIRTUS2 parsing hazard recorded: the
  list uses space-separated accessions (`NC 000883.2`) while its FASTA uses
  underscores, which silently drops 96 of 762 records under a naive regex.


## WP5 — Build and validate the truth panel

Objective: create deterministic read/molecule truth across supported chemistry
and ambiguity regimes. Estimated effort: 1-2 engineering weeks plus compute.

### WP5A — Generator

- [ ] `VAL-01` — implement a seeded generator spanning viral abundance, infected-
  cell fraction, host homology, sibling viruses, low complexity, PCR duplication,
  CB/UMI collisions, ambient/index hopping, and 10x v2/v3/Drop-seq geometry.
  - **2026-09-26:** `VAL-01` also unblocks `SENS-06` (WP4C). A generator that
    plants a target at a known abundance is exactly the positive control
    `SENS-04` consumes, so the two should be built together: it turns
    `informative_negative` from always-false into a measured quantity, and
    supplies the `E8`/`D17` probit input at the same time.
- [ ] `VAL-02` — emit paired FASTQs, `truth_manifest.tsv`, read/molecule truth
  tables, barcode/chemistry metadata, input hashes, and a run manifest.
- [ ] `VAL-03` — add synthetic host-only and adversarial GRCh38-homology
  conditions plus reagent/empty-droplet controls when available.
- [ ] `VAL-04` — add a checksum-pinned 10x healthy-donor PBMC v3 presumed-negative
  observational control; never label it absolute ground truth.
- [ ] `VAL-05` — add planted target and sibling-virus positives for EBV, HHV-6,
  HSV-1/2, KSHV, and a separately scored anellovirus panel.

### WP5B — Scorer and tiny gate

- [ ] `VAL-06` — implement molecule/cell precision, recall, F1, AUPRC, sibling
  confusion, host-homology false positives, burden concordance, calibration, and
  limit-of-detection scoring.
- [ ] `VAL-07` — validate the scorer against hand-computed fixtures and reject
  denominator, feature, barcode, or count-layer mismatch.
- [ ] `VAL-08` — run the golden tiny panel end to end and prove exact planted-
  molecule recovery, count conservation, determinism, row-order invariance, and
  chunk-size invariance.
- [ ] `VAL-09` — freeze generator/scorer version, seeds, manifests, and expected
  tiny outputs before cluster-scale execution.
- [ ] `VAL-10` — run the full training and untouched holdout panels, retaining
  failures rather than silently dropping conditions.

`G5a` passes when the golden tiny panel is exact and deterministic, the full
panel is manifest-complete, and every condition has a result or reproducible
failure record.

## WP6 — Run public positives and comparators

Objective: regenerate all eligible evidence with v3 and make cross-tool
comparisons denominator-, annotation-, barcode-, and layer-matched. Estimated
elapsed time: 2-4 weeks, dominated by downloads, queues, and comparator setup.

### WP6A — Public v3 reruns

- [ ] `RUN-01` — HHV-6B `SRR20710641` under measured 10x v3 geometry.
- [ ] `RUN-02` — EBV `SRR12682296` under measured 10x v2 geometry.
- [ ] `RUN-03` — HSV-1 `SRR8315713` under measured Drop-seq geometry.
- [ ] `RUN-04` — KSHV latent/lytic series `GSE190558`.
- [ ] `RUN-05` — KSHV plus EBV replicates `GSE154900`.

Each sample must produce a validated output tree plus input/reference hashes,
command manifest, tool versions, scheduler accounting, and failure log. Expected
virus recovery is contextual evidence, not ground truth.

### WP6B — Harmonized workflow matrix

- [ ] `CMP-00` — replace historical output-derived barcode unions/intersections
  with the frozen pre-outcome anchor/feature manifests and audited count-layer
  adapters; use `counts_unique`, retain structural zeros only from declared
  complete barcode domains, and record every missing/extra mapping.
- [ ] `CMP-01` — ViralScan combined and exact-fragment STAR two-step; keep
  kallisto two-step excluded unless exact fragment lineage becomes available.
- [ ] `CMP-02` — STARsolo combined and STAR host-filter/two-step.
- [ ] `CMP-03` — traditional host-genome subtraction followed by viral alignment,
  with and without CB/UMI retention.
- [ ] `CMP-04` — Venus in an isolated version/digest-pinned environment with no
  silent algorithm patch.
- [ ] `CMP-05` — Viral-Track in an isolated version/digest-pinned environment with
  no silent algorithm patch.
- [ ] `CMP-06` — (added 2026-09-30, user request) head-to-head with the
  **original ViralScan 2.2.0**, as released and run by evonk, on the **same
  reference** as v3.
  - **Scope.** The comparison measures *implementation* differences only. It
    runs 2.2.0 on the same inputs as v3. GOV-06, by contrast, compares
    *archived* 2.2.0 outputs.
  - **Reference (user decision 2026-09-30): one arm only.** Both versions use
    the **latest reference**: today the final 2,343-genome panel,
    `/exports/archive/hg-funcgenom-research/mdmanurung/viral_ref_final/build/{panel.idx,panel.t2g}`.
    *2026-10-03:* the latest reference is now
    `/exports/archive/hg-funcgenom-research/mdmanurung/viral_ref_cat42b/build/{panel.idx,panel.t2g}`
    (CAT-42: same 2,343 genomes, homopolymer and anellovirus low-complexity
    masks). Its anellovirus gene IDs carry the accession prefix.
    If the frozen G4 reference supersedes it, both versions switch to that.
    - Give both versions the same viral-only GTF, extracted from
      `combined.gtf`.
    - Never pass the combined host+viral GTF as `-gtf`: 2.2.0 would then count
      every host gene as viral (the Q8 issue).
    - There is no native-Serratus-index arm. Reference differences are out of
      scope.
  - **What 2.2.0 is.** The install is in env
    `/exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan`
    (`viralscan-2.2.0.dist-info`); evonk's run pattern is
    `/exports/archive/hg-funcgenom-research/evonk/viralscan/run_viralscan.sh`.
    Pin the version and the environment spec before any run. First check that
    2.2.0 accepts the 8-column kb-ref `panel.t2g` and its gene IDs; record
    whatever adapter is needed.
  - **Hold everything else equal:**
    - identical FASTQ hashes;
    - identical `-x` and whitelist;
    - the same strand. 2.2.0 has no `--strand`, so compare at kallisto's
      default and, once DEF-02 exists, report v3 at the inferred strand as a
      separate row;
    - the same cell-calling anchor where 2.2.0 allows one.
  - **Where it runs:**
    - **context** (already-seen data, R2.2, never tuning):
      - EBV SRR12682296 (10xv2);
      - HSV-1 SRR8315713 (Drop-seq);
      - HHV-6B SRR20710641 (5′, `-x 10xv2`);
      - evonk's own EBV SRR6825024;
      - covid x213/x216 and tonsil x223.
    - **confirmatory:** the truth-panel holdout, as an SCI-04 row.
  - **Report, per dataset:**
    - total and per-virus molecules;
    - called cells and the infected-cell denominator;
    - virus names and grouping;
    - false calls on the host-only and planted-homology negatives.
  - **Attribute every difference to a known implementation mechanism:**
    - barcode correction (SW-13/SW-20, F-018);
    - multimap allocation and siblings (F-017);
    - poly-G/low-complexity handling (F-019);
    - the strand default (F-020);
    - naming and grouping (MECH-A);
    - the summary headline (SW-16).

    Any unattributed difference is a finding in its own right.
  - **Dependencies:**
    - the context runs can start now as exploratory;
    - the confirmatory run needs G3/G4 and the frozen defaults (WP1E R2.1/R2.4);
    - the budget comes from the SCI-04 comparator allotment (WP1E R3.3).
Primary comparison rules: identical FASTQ hashes and viral sequences, same host
release, outcome-independent cell anchor, audited feature intersection, and
unique-only molecule parity. Recommended ambiguity-aware outputs are a separate
secondary analysis. Run at least EBV, HHV-6B, and HSV-1 for both dedicated
comparators.

`G5b` passes when every preregistered row is complete or transparently failed,
raw outputs and commands are retained, and an independent audit finds no
outcome-selected barcodes or mismatched denominator, annotation, feature, or
count layer.

## WP7 — Score, audit, and freeze results

Objective: turn completed runs into a single immutable scientific result bundle.
Estimated effort: 3-5 days after all runs finish.

- [ ] `RES-01` — compute all preregistered molecule/cell metrics, runtime, peak
  RAM, temporary storage, output size, calibration, and failure rates.
- [ ] `RES-02` — estimate uncertainty by bootstrapping biological samples, never
  by treating cells as independent experimental replicates.
- [ ] `RES-03` — calibrate thresholds on training data, freeze them, then evaluate
  the untouched holdout once; record every deviation.
- [ ] `RES-04` — independently audit inputs, references, barcodes, features,
  denominators, layers, manifests, scheduler records, and failed rows.
- [ ] `RES-05` — freeze `analysis/v3_results_bundle/` with a manifest containing
  every file digest, generation command, environment/tool version, Git SHA, and
  validation report.

`G5` passes only with 100% count-invariant compliance, complete ambiguous/
unresolved reporting, no probable/strong call in synthetic host-only or planted
host-homology negatives, positive recovery with sample-level uncertainty, and no
general superiority claim based on one virus or sample.

## WP8 — Regenerate documentation and claim evidence

Objective: make public instructions and claims match the actual v3 CLI, schemas,
and frozen results. Estimated effort: 4-7 days; quantitative pages wait for `G5`.

### WP8A — User documentation

- [~] `DOC-01` — reconcile README, installation, quickstart, CLI, outputs, API,
  FAQ, reference-panel, support, security, and migration docs with v3 contracts.
  - 2026-10-02 (`v3/docs`):
    - Added `docs/migration.md` (rebuild-only) to the toctree, and linked
      SUPPORT/SECURITY from the index.
    - Guard tests cover "2.5.0", faq, reference_panel and migration.
    - Removed the covid Alphatorquevirus counts from cli_reference, faq and
      output_reference, using the Q11 wording, and the "~4×" claim in
      `BENCHMARK_COMPARISON.md`.
    - Earlier, under DOC-02: the v2.5 container commands and the stale
      primary-call mode.
    - Still open:
      - Sphinx `-W`: nbsphinx isn't registered in the codex env and
        `sphinx_rtd_theme` is missing (installs not approved). A stripped
        build still has 37 pre-existing warnings: `csv` lexer, PLAN
        cross-refs in docs/plans, header jumps.
      - pip/full tier text (REL-01/02) and the data-fetch wording (REF-11).
- [x] `DOC-02` — generate CLI/default tables from the parser and test exact
  defaults; remove plain `em`, removed primary-call modes, silent knee fallback,
  v2.5 container commands, and raw-UMI language for fractional estimates.
  - Done 2026-10-02 (`v3/cli`, d2703eb).
    - `menu.build_parser()`; `scripts/gen_cli_reference.py [--check]` writes
      10 generated flag tables into `docs/cli_reference.md`.
    - `tests/test_cli_reference_parity.py` fails on any undocumented flag.
      The 11 that were missing are now covered.
    - The malformed 4-column rows are fixed.
    - Removed `--multimap-primary-call confidence` (showcase runbook) and the
      v2.5 container commands, and changed per-10k wording to molecule
      estimates.
    - A search found no "silent knee fallback" text, and plain `em` survives
      only in the notebooks, which belong to DOC-06.
- [x] `DOC-03` — clearly label every output as observation, model estimate,
  evidence tier, diagnostic flag, or biological interpretation.
  - Done 2026-10-02 (`v3/docs`, 33e78de).
    - Every column table in `output_reference.md` has a Kind column, with a
      legend of 5 labels: observation, model estimate, evidence tier,
      diagnostic flag, biological interpretation.
    - Four emitted viral_summary columns that were undocumented are added.
    - Tests check every label, plus viral_summary and
      `MULTIMAP_EVIDENCE_COLUMNS` coverage.
- [~] `DOC-04` — document combined versus two-step information loss, anellovirus
  screening limits, pip/full-workflow tiers, and legacy rebuild-only migration.
  - 2026-10-02 (`v3/docs`): FAQ sections on what combined vs two-step
    each lose, and on anellovirus screening limits (commensal prior, Kane
    et al.; F-019/F-021/F-022). Still open: the pip/full tier text, which
    needs REL-01.
- [ ] `DOC-05` — execute clean-install quickstart commands and validate the
  resulting run using only documented steps.

### WP8B — Vignettes and claims

- [x] `DOC-06` — eight notebooks exist but contain pre-v3 calls/values; rebuild
  them with negative and ambiguous examples and only v3 APIs/artifacts.
  - Done 2026-10-02 (`v3/vignettes`, c92b768). All 8 notebooks use v3
    APIs and columns:
    - `em-global`/`em-cell` and the `selected-method` primary call;
    - `viral_molecules_total_est`;
    - external cell calling, with knee as sensitivity-only.

    Legacy numbers are removed: 3.64×, 13–18 %/4,414, 28,922/30,849/19,920,
    and the hostresponse 0.87/0.967. Negative examples are in
    `specificity_true_negative`, `cell_calling_denominators` and
    `cell_type_enrichment`. Ambiguous examples are the HHV-6 sibling bleed and
    the rare-host EM case. The notebooks are committed output-stripped.
    - Verification: the 6 CI notebooks pass when their cells are run in a
      single exec namespace. nbclient could not run because no env has
      `ipykernel`, so CI `nbmake` is the kernel-path confirmation.
    - `results/hostresponse_ebv_matched/` is now unreferenced but still
      tracked; removing it is left to the user. DOC-07 is untouched.
- [ ] `DOC-07` — execute six lightweight notebooks in CI and the reference/full-
  workflow notebooks in the locked scheduled workflow; save logs and hashes.
  - **2026-09-27:** 3 of the 6 CI notebooks fail today, and they fail
    identically at `af4d5b3`. The breakage is notebook/API drift, not WP1C.
    Found while checking WP1C with nbclient, because `nbmake` is not installed
    locally.
    - `cell_calling_denominators`: `KeyError: "['viral_umi', 'total_umi'] not in index"`.
    - `qc_and_read_evidence`: `KeyError: 'viral_molecules_total_est'`.
    - `specificity_true_negative`: `KeyError: 'total_umi'`.
    - **Passing in both trees:** `multimapping_correction`,
      `cell_type_enrichment` (identical `padj` tables) and
      `host_response_depth_control`.
    - Fix the notebooks' column names against the current output schema before
      closing this row.
- [ ] `DOC-08` — add a balanced five-workflow pros/cons table generated from the
  harmonized benchmark rather than rhetorical claims.
- [ ] `DOC-09` — add a claim-registry schema, validator, stale-hash detection, and
  coverage check for README, docs, manuscript, tables, figures, and captions.
- [ ] `DOC-10` — require each quantitative, comparative, validated, performance,
  specificity, and installation claim to resolve to a `validated_v3` artifact.

`G6a` passes when Sphinx warnings are errors, links and shell examples pass,
all notebooks execute in their declared tier, parser/default parity is exact,
and claim coverage contains no missing or stale artifact.

Documentation gate:

```text
python3 -m sphinx -W -b html docs docs/_build/html
python3 -m sphinx -W -b linkcheck docs docs/_build/linkcheck
PYTHONPATH=src python3 -m pytest tests/test_docs_consistency.py -q
```

## WP9 — Regenerate and review the manuscript

Objective: create a submission package exclusively from the frozen v3 bundle.
Estimated effort: 1-2 writing weeks after `G5`, excluding author review.

- [ ] `MS-01` — replace the historical draft and hard-coded legacy figure script
  with generators for registered values, tables, figures, captions, and
  supplement sourced only from `analysis/v3_results_bundle/`.
- [~] `MS-02` — frame molecule-aware host-virus ambiguity, combined/two-step
  evaluation, read-level specificity/QC, reference provenance, and honest
  evidence tiers; remove every legacy count and unsupported superiority claim.
  - Blocked (2026-10-02): the prose waits on G5/MS-01. This pass's ship-doc removals of the covid counts are logged under `DOC-01`/`DOC-02`, not here.
- [ ] `MS-03` — rerun EBV host response using v3 labels with depth and
  mitochondrial controls; retain only if it replicates. Keep TTV solely as a
  host-homology case study absent orthogonal validation.
- [!] `MS-04` — authors must supply author order, affiliations, ORCIDs, CRediT,
  lead contact, funding, conflicts, acknowledgments, and ethics/consent/data-
  access text. Exclude private COVID libraries unless all requirements are met.
- [ ] `MS-05` — verify every citation, run independent methods/statistics and
  claim-to-artifact reviews, resolve all findings, and generate the final
  submission/availability package.

`G6` passes when every number/table/figure resolves to a frozen artifact digest,
there are no placeholders or ineligible private results, citations support their
sentences, and independent review has no unresolved correctness, denominator,
annotation, uncertainty, ethics, or unsupported-claim finding.

## WP10 — Publish and test `3.0.0rc1`

Objective: expose the exact candidate artifacts to real users without promoting
them as stable. Estimated release work: 1-2 days; testing window: 2-4 weeks.

- [ ] `RC-01` — freeze `3.0.0rc1` version, changelog, migration notes, CFF, locks,
  reference manifest, checksums, SBOM, licence report, and release notes on one
  green protected-main SHA.
- [!] `RC-02` — maintainers configure/approve TestPyPI or PyPI trusted publishing,
  GHCR/GitHub release permissions, and public artifact hosting.
- [ ] `RC-03` — tag `v3.0.0rc1` only after exact-SHA gates; publish prerelease
  wheel/sdist, versioned OCI, SIF, locks, checksums, attestations, reference
  archive, and GitHub prerelease. Do not update `latest`.
- [!] `RC-04` — recruit at least three external laboratories covering supported
  chemistries; give each a tester packet for clean install, tiny workflow,
  `validate-run`, and one supported real workflow.
- [ ] `RC-05` — track every correctness, data-loss, installation, documentation,
  and reproducibility failure to closure; rerun all gates after fixes.

`G7` passes when all three testers submit complete manifests, each supported
path works from public artifacts, and no release-blocking issue remains.

## WP11 — Publish stable `3.0.0` and submit

Objective: publish one exact reviewed build everywhere, archive it, then submit
the methods manuscript. Estimated release work: 1-3 days after approvals.

- [ ] `STB-01` — freeze final truth-panel/comparator/manuscript bundles and prove
  the release candidate's fixes did not change scientific outputs unexpectedly.
- [!] `STB-02` — reserve the Zenodo software DOI and benchmark/archive DOI, keep
  them distinct from the reference-data DOI, and provide them for metadata.
- [ ] `STB-03` — update CFF/metadata, pass exact-SHA gates, tag `v3.0.0`, and
  publish wheel, sdist, versioned plus `latest` OCI, SIF, locks, reference archive,
  release notes, migration guide, checksums, SBOM, licences, and attestations.
- [ ] `STB-04` — after the real PyPI sdist exists, insert its SHA-256 into the
  Bioconda recipe, lint/build/test it, submit the PR, and verify the installed
  package on the tiny workflow.
- [ ] `STB-05` — archive software and benchmark artifacts, verify every public
  DOI/link/download, update `CITATION.cff`, then submit the manuscript.

`G8` passes when the tagged SHA is the reviewed green SHA, public distribution
formats reproduce validated counts, all links and DOIs resolve, the software and
benchmark archives are durable, and the manuscript has been submitted. Journal
acceptance timing is external and is not a software completion condition.

## WP12 — Observe and maintain

- [~] `OPS-01` — publish supported-version, patch-release, deprecation, and
  security-response policies with named contact routes.
- [ ] `OPS-02` — establish reference-update cadence and manifest compatibility
  rules; never silently change a frozen reference under an existing version.
- [ ] `OPS-03` — triage false positives, false negatives, data loss, and reference
  drift as scientific incidents with reproducible packets.
- [ ] `OPS-04` — monitor install/usage failures during the first 90 days and ship
  patch releases from the same gate process.
- [ ] `OPS-05` — schedule independent truth-panel refreshes without changing v3
  thresholds retrospectively.

## Gate dashboard

| Gate | State | Required proof |
|---|:---:|---|
| `G0` governance | `[~]` | prescribed PyPA-frontend archive build and member check; governance validators otherwise pass |
| `G1` software | `[~]` | full unit/property/safety/tiny-workflow suite |
| `G2` distribution | `[~]` | clean installs, locks, Docker/Apptainer parity, supply-chain reports |
| `G3` preregistration | `[ ]` | reviewed, schema-valid, hashed protocol frozen before outcomes |
| `G4` references | `[~]` | byte-rebuild, GRCh38 D-list, calibrated holdout safeguard |
| `G5` science | `[ ]` | truth holdout, public positives, all comparators, uncertainty, frozen bundle |
| `G6a` docs | `[~]` | parser parity, executable docs, complete validated claim graph |
| `G6` manuscript | `[ ]` | artifact-generated manuscript plus independent audits |
| `G7` release candidate | `[ ]` | public RC artifacts and three-laboratory closure |
| `G8` stable/publication | `[ ]` | stable artifacts, archives/DOIs, Bioconda, manuscript submission |

## External inputs and authority

These items cannot be invented or completed by an implementation agent:

| Input | Needed by | Owner action |
|---|---|---|
| Author/ethics/funding/conflict metadata | `REL-15`, `MS-04` | authors approve final factual text |
| Publishing credentials and protected environments | `RC-02` | repository owner configures services |
| Three external laboratories | `RC-04` | maintainers recruit and coordinate testers |
| Zenodo DOI reservations | `STB-02` | archive owner reserves distinct DOIs |
| Orthogonal anellovirus-positive sample | `REF-10` | collaborator supplies lawful validated data, or claim stays screening-only |
| Published viral annotation panel archive and registered DOI | `REF-11` | archive owner publishes the panel and registers `10.5281/zenodo.20112332`, or the pinned record identifier is corrected |

## Stop rules

- Never tag or publish from a dirty, unreviewed, or non-green SHA.
- Never migrate or reinterpret a pre-v3 H5AD value as a v3 molecule count.
- Never select shared barcodes using ViralScan-positive outcomes.
- Never turn ambiguity/QC flags into biological conclusions automatically.
- Never claim real-anellovirus sensitivity, general superiority, or formal
  infection from a nonzero molecule without the required validation evidence.

## Evidence log

Append one line after each completed item:

```text
YYYY-MM-DD ITEM — command/result; artifact path(s); Git SHA; reviewer if required
```

- 2026-07-22 `BASE-03` — EBV v3 baseline conserved exact molecule mass;
  `analysis/v3_ebv_baseline/host_conservative.json`; SLURM 25316336.
- 2026-07-22 `BASE-04` — real-tool evidence integration passed;
  `tests/integration/test_exact_lineage.py` and
  `tests/integration/test_evidence_chain.py`; Git SHA `a146050`.
- 2026-07-22 `BASE-05` — full non-network suite: 720 passed, 21 deselected;
  Ruff check/format, data-governance check, and draft protocol validation passed.
- 2026-07-22 `SW-01` — wheel and sdist contain all six byte-matched v3 schemas;
  fresh wheel install loaded 6/6, missing-schema paths fail closed; 8 focused
  tests passed; Git SHA `6e1fe66`.
- 2026-07-22 `REL-04` partial — isolated wheel/sdist build and fresh-wheel schema
  smoke test passed; `twine check`, sdist clean-install test, and `SHA256SUMS`
  remain.
- 2026-07-22 `DOC-01`–`DOC-04` partial — public docs and tracked manuscript
  placeholder reject legacy quantitative/private claims and use v3 molecule and
  candidate-evidence terminology; 94 docs/CLI tests passed; Git SHA `3bb7d1b`.
- 2026-07-22 `SCI-01` — schema-valid, non-executable protocol draft with separate
  training/holdout gates; 16 focused validation tests and 4 docs-consistency
  tests passed; data-governance check passed; independent implementation review
  resolved four blockers and passed; `analysis/v3_validation/protocol.yaml`,
  `schemas/v3/validation_protocol.schema.json`, and
  `scripts/validate_v3_protocol.py`; Git SHA `d941261`.
- 2026-07-22 `SCI-02` — outcome-independent harmonization frozen with 24
  endpoint-mapped denominator contracts, canonical and dependent-field digests,
  49 focused adversarial tests, schema/CLI fail-closed checks, data-governance
  pass, and three independent reviewer passes; protocol remains non-executable;
  Git SHA `d941261`.
- 2026-07-25 `GOV-06` partial — all 44 retained BUS rows validated and
  aggregated; six downloaded EBV mates passed authoritative size/MD5 and local
  SHA-256 gates; pair/chemistry audits `25330878` and `25330881` completed all
  five controls; `analysis/legacy_v2_v3/control_inputs.tsv` is the sanitized
  shared manifest. Fresh attempt 1 (`25331024`, `25331026`) failed before
  execution because Slurm-spooled scripts could not resolve the packet root;
  attempt 2 uses a tested explicit packet root in arrays `25331035` and
  `25331037`. Full suite: 758 passed, 48 deselected; base Git SHA `a8f3287`.
- 2026-07-26 `GOV-06` partial — attempt-2 arrays `25331035` and `25331037` are
  terminal with all ten rows failed and every outcome retained. Three
  independent causes, each confirmed from primary evidence
  (`fresh_control_packet_attempt2/status/*.json`, `sacct`, on-disk output
  trees): (1) `run_fresh_control.py` `_v2_artifact_errors` checked a flat output
  layout while the legacy 2.2.0 CLI nests one level under the sample
  identifier, so `v2__SRR12682296`, `v2__SRR12682297`, `v2__SRR12682298`, and
  `v2__SRR6825025` recorded exit 65 with six phantom missing artifacts despite
  `workflow_exit_code` 0 and complete non-empty output trees; (2) all five v3
  rows raised `ViralScanDataError` during config creation because the Zenodo
  viral-annotation cache was never populated in the frozen packet environment;
  (3) `v2__SRR6825024` was genuinely out-of-memory killed in the legacy
  `multimap.py` at roughly 121.4 GiB peak resident set against the frozen
  128 GiB tier ceiling after 5 h 52 m. A fourth defect was found while
  verifying: the status payload satisfies none of `stage`, `attempt_id`, or
  `scientific_parameter_hash` from `protocol.yaml` `required_failure_fields`,
  so all ten records are non-compliant as failure records. Attempt 3 was
  prepared as a v3-only five-row packet and deliberately not submitted; the
  diagnostic is outcome-ineligible and yields priority to `SCI-03`.
  **Superseded 2026-07-28**: this entry originally read "Attempt 3 is frozen".
  It is not, and no attempt-3 packet exists on disk — only
  `benchmark_runs/legacy_v2_v3/fresh_control_packet_attempt2/`. The `GOV-06`
  work-package row and `analysis/legacy_v2_v3/TRACKER.md` both correctly record
  it as blocked on `REF-11`. Same drift class as the `SCI-03`/`SCI-04` entries
  corrected on 2026-07-27; found by the 2026-07-28 multi-agent review.
- 2026-07-27 `SCI-03` partial — `partitions` and `calibration` are frozen,
  schema-validated sections of `analysis/v3_validation/protocol.yaml`. Partitions
  allocate whole biological samples by deterministic stratified assignment over
  the declared factors at a 0.3 holdout fraction, with a single permitted holdout
  evaluation and explicit template/locus/molecule/cell-barcode leakage
  prohibitions. Calibration fixes a training-only threshold grid with a
  zero-false-positive constraint on host-only and planted-homology negatives, a
  conservative tie-breaker, nine endpoint-linked metrics, a probit LOD95 with
  extrapolation prohibited, and a 2000-replicate biological-sample bootstrap that
  forbids treating cells, molecules, reads, or technical repeats as independent
  replicates. Canonical digests: partitions
  `44173743fc7b0c3e9a48b165348faca196c535bd82eee64f4ea771724a3beb04`, calibration
  `252fa8ade5381266887527cde3ddc139a289e81d18fec9be1ce99da4fe4d6d7c`. Seeds
  `root`, `split`, `cell_calling`, `evidence_sampling`, and `bootstrap` are
  frozen; `generation` remains `VAL-01`. The `partitions_metrics` training
  blocker is closed; `data_hashes` stays open. Validator gains `_validate_sci03`
  with digest, seed, endpoint, and factor cross-checks. Draft gate valid;
  training gate still correctly blocked on `SCI-04`, `SCI-05`, `VAL-01`, and
  `REF-09`. Full suite 790 passed, 48 deselected.
- 2026-07-27 `SCI-04` partial — `workflow_matrix` and
  `failure_and_deviation_reporting` are frozen, schema-validated sections. The
  matrix enumerates eight workflows over 52 rows: ViralScan combined across all
  twelve datasets, ViralScan exact-fragment STAR two-step, STARsolo combined and
  host-filter two-step, traditional host subtraction with and without CB/UMI
  retention, and Venus and Viral-Track at published defaults in isolated
  digest-pinned environments. Both dedicated comparators cover EBV, HHV-6B, and
  HSV-1 plus the host-only and host-homology negatives. The kallisto two-step
  path is excluded with an explicit revisit condition rather than omitted.
  Seventeen per-row record fields are required, and five primary comparison rules
  fix identical FASTQ hashes, host release, the outcome-independent anchor, the
  audited feature intersection, and unique-only parity. Canonical digests:
  workflow matrix
  `fb4c6768e87a87f23fa554c1376cd6b3ce422d972d3f2e409d83985c154bde4c`, failure
  reporting `ce913b3ea10f23062ab0f86bbf81ce341ea9b02edad594d57cbab363c5e42975`.
  `_validate_sci04` enforces digest match, planned-section pairing, dataset and
  reference resolution, and row-count agreement, and at the training gate refuses
  any workflow lacking a pinned tool version and container digest.
  `environment_pinning` stays `pending` on `REL-03`, so the `workflow_rows`
  blocker is closed but `tool_environments` remains. Full suite 800 passed, 48
  deselected.
- 2026-07-26 `REF-11` opened — the Zenodo record pinned by
  `src/viralscan/data_fetch.py` is unregistered. `zenodo.org/api/records/20112332`
  returns `{"status": 404, "message": "The persistent identifier is not
  registered."}` and `doi.org/10.5281/zenodo.20112332` returns 404, verified from
  a network-capable host on which an unrelated third-party Zenodo DOI resolved
  200. The installed package under
  `benchmark_runs/legacy_v2_v3/env_full/.../viralscan/data/` contains only
  `anellovirus_accessions.tsv` and `__init__.py`, while 195 GTFs remain in
  `src/viralscan/data/` in the source tree. `viralscan data fetch` therefore
  cannot succeed for any user, which blocks the `GOV-06` attempt-3 cache pin and
  the clean-install paths behind `REL-05` and `DOC-05`.
- 2026-07-27 evidence-log correction — the `SCI-03` and `SCI-04` entries above,
  both dated 2026-07-27, state that `partitions`, `calibration`,
  `workflow_matrix`, and `failure_and_deviation_reporting` "are frozen". That was
  true when written and is no longer: `SCI-05` round 1 restored blockers and all
  four sections returned to `status: pending`, where they remain. The work-package
  rows carry the authoritative status. The entries are left in place because this
  log is append-only history, not current state.
- 2026-07-27 six-agent code review of `codex/viralscan-v3` vs `main` —
  `.living/outputs/reviews/2026-07-27-branch-codex-viralscan-v3.md`, with a
  behavioural tripwire audit alongside it. Two Majors, both fixed here. (1) The
  protocol froze `seeds.cell_calling: 20260727002` and named it the seed source
  for the shared cell anchor, but `call_cells` never passed `seed` to
  `emptydrops_cells`, so emptyDrops always ran at the signature default `100`;
  `emptydrops_seed` and `emptydrops_niters` are now declared config fields with
  CLI flags, and `emptydrops_cells` takes keyword-only required parameters so no
  future caller can omit one silently. (2) `check_git_sha_fields` accepted a
  record's before-digest unverified when its scoped section was absent at the base
  commit or when `git show` failed — the open `R11-F2`; both branches now error.
  Minors fixed: frozen FASTQ identity is re-derived from the bytes in
  `verify_frozen_fastq` rather than trusted from a weeks-old audit (kept out of
  `prepare_fresh_controls`, where it would re-read ~250 GB to close no additional
  window); `Rscript` is preflighted when the resolved cell caller is `emptydrops`,
  since `SW-11` made that path fail closed after `kb_count`, `analysis`, and
  `multimap`; `detect_cells` documents `viral_count_matrix`. One review finding
  was reversed on verification: `SIBLING_VIRUS_PAIRS` omits the protocol's
  EBV/KSHV pair *correctly* — the constant is a runtime EM-bleed heuristic
  requiring near-identity, the protocol list is an evaluation population spanning
  a relatedness gradient, and adding EBV/KSHV to the constant would annotate
  genuine co-infection as artifact. Both sides now say so.
- 2026-07-27 `DEV-019` — the protocol edit above changed
  `partitions.contract_sha256` to
  `5fd9b366001c0101dd174c79d2d446cbc957525b0a375e659937d4aaa4eb1a05`. An earlier
  draft of this entry claimed no ledger record was required because `partitions`
  is `pending`. That was wrong and unverified: `partitions` carries an
  eight-link chain that was reconciled at `0dae092a…` before the edit, so the
  edit orphaned it. `DEV-019` records the change, and declaring it re-digests the
  frozen `frozen_inputs` section, carried in the same record's
  `additional_digest_changes`. All five sections with chains now reconcile against
  the live protocol. The gap that let this happen is closed by
  `pending_section_ledger_drift`: `validate_amendment_ledger` binds only *frozen*
  sections, so a pending section that already has records could drift unnoticed
  until the moment it was frozen. That check now blocks the training phase, which
  must pass before any freeze.
- 2026-08-08 `GOV-03` — `analysis/v3_artifact_inventory.tsv` now has 30
  sanitized, stable-identity rows spanning claim-bearing inputs, references,
  intermediates, results, commands, environment, scheduler evidence, failure
  evidence, schemas, and documentation. The fail-closed schema/hash/cross-link
  validator passed; the provenance-incomplete retained EBV baseline remains
  private and ineligible for a public validated-v3 claim. Base Git SHA
  `26260e7cbc18cc0e7777379e1db12cc594778dd4`.
- 2026-08-08 `GOV-04` — `config/public_ship_scope.json` is the single positive
  wheel, sdist, Docker-context, public-documentation, claim-bearing, and
  governance-text allowlist. Packaging no longer recursively includes `docs/`,
  Docker no longer uses `COPY .`, the allowlisted-text institutional-path scan
  passed, and wheel/sdist/context members matched exactly when the distributions
  were built through the installed setuptools backend. Base Git SHA
  `26260e7cbc18cc0e7777379e1db12cc594778dd4`.
- 2026-08-08 `GOV-05` — the eight-record claim graph is schema-valid and covers
  every allowlisted public marker with checked artifact identities and hashes.
  Legacy counts are explicitly ineligible, reference homology is implemented but
  not calibrated, and the retained EBV baseline is provenance-incomplete rather
  than promoted. Focused governance tests: 18 passed. All three governance
  validators and the draft protocol validator passed. Base Git SHA
  `26260e7cbc18cc0e7777379e1db12cc594778dd4`.
- 2026-08-08 `G0` partial reconciliation — the default non-network suite passed
  958 tests with 73 deselected; Ruff check and changed-file format checks passed.
  The prescribed `python -m build --no-isolation` command is unavailable because
  the locked runtime has no PyPA `build` frontend; no network install was
  attempted. Repo-wide Ruff format remains red on ten protected/out-of-scope
  pre-existing files. A separate corrected-PATH integration run produced 19
  passed, 2 failed, and 17 errors: the shared real-workflow fixture requires a
  missing `conda` executable, and this Snakemake emits no DAG rule listing under
  the tests' `--quiet` invocation. `G0` remains `[~]` and no later gate is
  promoted.
- 2026-09-27 housekeeping — the WP4B2–WP4F work was squashed into `3379b7c`
  so every commit passes its own suite (the first half had registered
  `hpv_genes.*` and `anellovirus_genes.tsv` in the ship-scope allowlist before
  those files existed); the result tree is byte-identical to the pre-squash
  history, kept on `backup/wp4-pre-squash`. The schema 1.1.0 governance
  migration was finished in `f447ff3` (config at 1.1.0, redundant `text_files`
  dropped — all 86 of its entries are already covered by
  `ship_scope_source_paths()`) and its three self-referential inventory rows
  pinned in `7b879f5`. `471f889` makes the two bundled-GTF tests skip when the
  gitignored panel is absent. A clean-checkout worktree at `471f889` passed
  **1,260, skipped 6, failed 0**; `check_data_governance.py` and
  `validate_claim_registry.py --coverage` pass. Nothing pushed.
- 2026-09-27 `PROG-07` re-measurement — EBV LCL `SRR12682296` gene programmes
  after removing `BARF1.2`/`BaRF1.1` from the latent set: unique layer 895
  latent / 236 productive / 311 mixed / 3,094 indeterminate (was 2,240 / 102 /
  445 / 2,968); allocated layer 856 latent (was 1,277). Summary sha256
  `48d0a9d12fbc6ca0…` (pre-change `72218893ff71baa3…`, kept as
  `gene_program_summary.pre-PROG-11.tsv`). Git SHA `50253a6`.

- 2026-09-27 `TONSIL-01` — SFL tonsil x223 (5′ v3 GEX, 24 donors), host-subtracted
  screen from the cellranger BAM. Stages: 118,856,604 unmapped GEX reads →
  24,969,278 pass the prefilter → 34,208 viral hits → 1,511 pass the viral
  filters → 23 survive the host re-check. Anelloviridae 0 and HPV 0 in all 24
  donors. Plant recovery: HPV16 99.2 %, TTV 99.4 %. Scripts and outputs are in
  `benchmark_runs/sfl_tonsil_screen_2026-09-26/tonsil01/` (gitignored): SLURM
  25652114 (extract), 25652133 (screen), 25652132 (plant). `calls.tsv` sha256
  `6f853b9f36dbe1e3…`. All 4 FASTQ md5 checks pass (25652115). Opened at Git SHA
  `a86aa6e`.
- 2026-09-27 `TONSIL-01` sensitivity to held-out strains (SLURM 25652148):
  - Held out: MN770908.1 (Beta), MW679005.1 (Alpha), MW455373.1 (Gamma).
  - Recovery, 5′ window / uniform: Beta 87.7 % / 72.9 % (nearest 94.3 %),
    Alpha 15.5 % / 34.7 % (nearest 84.0 %), Gamma 39.5 % / 8.0 % (no asm20 hit).
  - `plant_ho/out/results/reads_final.tsv` sha256 `6dd488146ceceb78…`.
  - The scripts are gitignored, so this commit's SHA does not pin them. Their
    sha256 prefixes: `build_db.py` df022722c37c, `prefilter.py` 421964e9f693,
    `parse_hits.py` fca465b4fe5c, `plant.py` c1d2bd0b07cd, `heldout.py`
    cd98664c5004, `screen.sbatch` ebb86194b0fb, `extract.sbatch` 54a1b2c3d9a3,
    `plant.sbatch` d07115d49705, `heldout.sbatch` 5bac5aea2a8b.
  - `screen.sbatch` gained an optional reference-dir argument after the main
    run. Its default is unchanged.
- 2026-09-27 `SIMP-01` — dead code deleted. Unit suite 1,275 passed, 0 failed.
  A diff of collected test IDs against `af4d5b3` shows exactly one test
  removed, `test_config_value_serializes_none_as_empty_string` (its function
  was deleted). `tests/test_multimap.py` is inventoried, so its docstrings
  that name the removed `normalize_barcodes` were left unchanged; re-pinning
  it for wording alone is not worth it.
- 2026-09-27 `SIMP-02` — `src/` helpers collapsed, 186 lines removed and 84
  added. Unit suite 1,273 passed. The 2 failures are the expected stale
  `build_reference.py` pins, fixed by the follow-up pin commit. Collected test
  IDs are unchanged.

  Integration outcomes match `af4d5b3` exactly: 18 passed, 19 skipped, and 1
  failure that predates this work. The failure,
  `test_build_anellovirus_reference_produces_labelable_gtf`, still expects the
  pre-WP4F `_gene1` placeholder IDs.

  A before/after probe against an `af4d5b3` worktree gave identical results
  for:
  - the packaged-table paths and row counts;
  - sha256 of a `Path`, a `str`, and `hashlib`;
  - the NCBI cache validity check;
  - the gzip opener;
  - the matrix axis sums.

  All 8 `--help` screens exit 0. The `--verbose`/`--quiet` help text is now the
  same on every subcommand.
- 2026-09-27 `SIMP-03` — dependencies: `pyfiglet` and `seaborn` removed, and
  `scipy` pinned `>=1.11`.

  **BH swap (`enrichment._bh_adjust` → `scipy.stats.false_discovery_control`):**
  - Across 7 input sets, including 300 real `fisher_exact` p-values, old and new
    agree to 1.1e-16, which is one floating-point ulp.
  - NaN and p > 1 cannot reach `_bh_adjust`. `fisher_exact` clips to [0, 1],
    and a brute-force run over 1,296 tables with counts 0–5 found none outside
    that range.
  - The new `test_bh_adjust_matches_reference_step_up` pins the equivalence.

  **Plots moved to matplotlib:**
  - Bar heights, tick order and histogram bin counts are identical to
    seaborn's. Seaborn's error bars had zero length.
  - The histogram keeps seaborn's axis labels.

  **Banner:** the welcome banner is plain text unless `pyfiglet` happens to be
  installed. `menu.py` already fell back, which makes the `tests/conftest.py`
  stub redundant, so it was deleted.

  **Gates:**
  - Unit suite: 1,275 passed. The 1 failure is the expected `environment.yml`
    pin, fixed in the pin commit.
  - Collected test IDs: +1, the new BH test.
  - Integration outcomes: identical to `af4d5b3`.

  **Changed files:** `pyproject.toml`, `environment.yml`,
  `conda-recipe/meta.yaml`, `.github/workflows/ci.yml` and `release.yml`.
- 2026-09-27 WP1C regression follow-up. CI `mypy src/viralscan` found 3 new
  errors from SIMP-02/03, all missing or too-narrow annotations with no runtime
  change:
  - the `_bh_adjust` return type;
  - `_open_maybe_gzip` now accepts `str | Path`;
  - `_write_tsv` parameter types.

  After the fix, mypy matches `af4d5b3` exactly: 57 errors, 0 new. The
  untracked scratch file `kmer3.py` was excluded.

  Other checks:
  - **Names:** an AST check of all 606 `viralscan` names imported across 174
    tracked `.py` files found none missing.
  - **Smoke test:** the CLAUDE.md smoke test passes, including with `pyfiglet`
    blocked, where the fallback banner is used.
  - **Containers:** `Dockerfile` and `Singularity.def` build from
    `environment.yml`, so there is no parity drift.
  - **Unit suite at `e0a7b1f`:** 1,276 passed.

  Added `ANELLO-14` for the integration failure that predates this work.
- 2026-09-27 WP1C final gate at `ff86c3d`:
  - Unit suite: 1,276 passed.
  - Integration: identical to `af4d5b3` (18 passed; `ANELLO-14` fails in both).
  - mypy: 57 errors, the same as `af4d5b3`, 0 new.
  - CI vignettes run with nbclient: 3 pass in both trees; 3 fail in both trees
    with identical errors (logged under `DOC-07`).
  - `cell_type_enrichment` `padj` output is identical before and after the BH
    swap.
- 2026-09-27 `CAT-11` — `scripts/measure_kmer_capture.py` reproduces the
  `REF-01` evidence: EBV `NC_007605.1` 144,283 31-mers; bundled 20-genome panel
  median coverage 0.00 %, zero-coverage 85.80 %, Betatorquevirus 98.4 %
  (n=1,542), Gamma/Samek/Het/Gyro/Mem 100 %; expanded 2,042-genome panel
  leave-one-out 20.55 %, zero-coverage 0.15 %, panel k-mers 4,888,291.
  **Correction:** the expanded panel's fragment capture is a measured **0.5097**,
  not the claimed analytic 1.0000 (F-011). Population
  `viralscan_showcase/.../anellovirus.fa` (2,042 genomes); outputs
  `bundled20.{tsv,json}` and `expanded2042.{tsv,json}` under the session
  scratchpad; 12 tests in `tests/test_measure_kmer_capture.py`. Git SHA `db6c699`.
- 2026-09-27 `CAT-02`/`CAT-03` — seed catalogue and catalogue-driven naming.
  `extras/build_virus_catalog.py --from-fasta references/starsolo/all_virus_serratus_plus_anellovirus/viral_genome.fa`
  → `src/viralscan/data/virus_catalog.tsv`: 2,215 accessions, 204 species, 30
  families, 0 failures (170 flatfiles fetched, the rest cached). Naming measured
  on the same reference: unnamed gene IDs 173 → 0, distinct groups 182 → 107,
  and influenza A's 8 segments group into 1 virus instead of 8. 1,307 tests
  pass. Git SHA `f39e18d`.
- 2026-09-27 `CAT-09` first batch — 16 previously-absent viruses added to the
  catalogue, every accession verified by live NCBI lookup: SARS-CoV-2, HIV-1/2,
  HTLV-1/2, HCoV-OC43/NL63/HKU1, hMPV, bocavirus, TSPyV, HPyV6/7, simian foamy,
  influenza D (7 segments) and the 16 HPV genotypes. Catalogue 2,249
  accessions / 232 species / 31 families (was 2,215 / 204 / 30); 714 KB.
  1,307 tests pass. Git SHA `af97cf2`.
- 2026-09-27 `CAT-11` correction — **the capture metric measured the forward
  strand only.** kallisto indexes canonical 31-mers; `kmers()` did not fold
  reverse complements, so a panel genome deposited antisense to a target
  counted as a miss. `--strand canonical` is now the default;
  `--strand forward` regenerates the superseded numbers. Recomputed on the
  2,042-genome anellovirus panel:

  | | forward | canonical |
  |---|---|---|
  | median leave-one-out | 0.2055 | **0.2162** |
  | median p_fragment | 0.5097 | **0.5349** |
  | panel distinct 31-mers | 4,888,291 | 4,843,779 |

  Per genus (median p_fragment, forward → canonical): Alpha 0.788 → 0.808,
  Het 0.837 → 0.837, Beta 0.498 → 0.515, Gamma 0.332 → 0.336, Samek 0.255 →
  0.256, Mem 0.247 → 0.247, Gyro 0.063 → 0.063, unclassified 0.334 → 0.397.
  350 of 2,042 genomes improve and **none degrade**; 84 by > 0.10 absolute and
  55 by > 0.25. Worst case `MH649023.1` 0.0329 → 0.9514. Affected accessions
  cluster in the `MH648xxx`/`MH649xxx` block. 17 tests in
  `tests/test_measure_kmer_capture.py` (was 12); `--self-check` now asserts a
  reverse-complement-only panel still captures its target.
- 2026-09-27 `CAT-12` provenance — **the anellovirus panel's sensitivity
  ceiling is set by upstream dereplication, not by how much diversity exists.**
  The 2,022 non-RefSeq panel genomes come from
  `github.com/clareaulab/human_anellovirus_pangenome`, which took 3,545
  complete human-host anellovirus genomes (NCBI Virus taxon 687329, June 2025)
  and **CD-HIT clustered them at 95 % ANI / 85 % coverage down to 2,023
  representatives**, then resolved genus for ~584 NCBI-unclassified genomes by
  ORF1 protein phylogeny (MAFFT → trimAl → IQ-TREE, 40 clusters).

  Three consequences, each measured:

  1. A 31-mer survives only with zero mismatches, so at the 95 % ANI threshold
     the expected shared fraction is `0.95^31 = 0.2039`. Measured median
     canonical leave-one-out is **0.2162**. The agreement to ~0.01 says
     CD-HIT's threshold, not NCBI's holdings, set our capture. The 1,522
     genomes CD-HIT discarded sit in exactly the divergence band where 31-mer
     pseudoalignment fails, so **restoring them is the one lever that raises
     capture** (2,023 → 3,545, +75 % genomes, ~10.5 Mb). Caveat to test before
     acting: the binomial assumes substitutions are independent and uniform,
     while F-011 showed clustering matters; the two must be reconciled.
  2. **The panel's genus labels are phylogeny-derived and are better than
     GenBank's.** The 582/2,042 panel genomes I found with no GenBank genus are
     the ~584 the upstream ORF1 tree resolved. Upstream reports Betatorquevirus
     1,558 by phylogeny vs 1,360 by NCBI label; measured here independently,
     1,542 vs 1,368 — the same correction. Do **not** "fix" our labels against
     GenBank lineages; that would discard the better assignment.
  3. **The absent Anelloviridae genera are correctly absent.** Every genus
     missing from the panel has **zero** complete genomes with
     `"Homo sapiens"[Host]`, verified by live NCBI query 2026-09-27:
     Lambda 408 complete / 0 human, Eta 256/0, Iota 139/0, Kappa 24/0,
     Rho 19/0, Epsilon 18/0, Pi 14/0, Theta 10/0, Sigma 9/0, Upsilon 8/0,
     Mu 5/0, Xi 4/0, Zeta 2/0; Delta, Nu and Tau have none at all. They are
     swine, feline, canine, tupaia and pinniped viruses. Adding them to a human
     panel would add false-positive surface and no sensitivity. This closes the
     "add more genera" line of `CAT-12` and applies to Gyrovirus too.

  Masking checked and excluded as a confound: the shipped panel FASTA carries
  728 N bases in 5,998,625 (0.01 %), 16 genomes affected, worst 11 %.
