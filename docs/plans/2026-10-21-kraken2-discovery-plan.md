# ViralScan × Kraken2: cell-resolved taxonomic discovery — critical review and implementation plan

Status: **approved 2026-10-21** (user: "approve"; all nine §16 decisions as proposed). Implementation
proceeds ticket by ticket (§13.3), spikes DISC-04 and DISC-05 first.
Date: 2026-10-21. Branch: `codex/kraken2-discovery-plan` (from `origin/main` `b8e5146`).
PLAN rows: `DISC-00` (this document) and `DISC-01`…`DISC-14` (tickets, §13).

Evidence labels used throughout:

- **VERIFIED**: read in the code at the commit given.
- **INFERRED**: follows from code or tool semantics but was not executed.
- **PROPOSED**: a design choice made in this document.
- **UNRESOLVED**: needs a spike, a measurement, or a human decision.

**Classification is not infection. Read support is not molecule count. Host subtraction does not
prove non-host origin. A discovered candidate is not a validated viral call.**

---

## 0. Audited revisions

| Repository | Commit | Note |
|---|---|---|
| ViralScan | `b8e5146` (origin/main, 2026-10-20) | worktree `vs_work/kraken2-plan` |
| julianeweller/mg2sc (original) | `a8c01f0db2dd81ca13704b19c1882eff3e5519ad` | |
| sinanugur/mg2sc (Cellsnake fork) | `e744c2cc77098bcd1ff6a3a4876a5b6203c44ae6` | also fork HEAD |
| sinanugur/scrna-workflow | `e0b95334dcb4c00a105a2e5fe775056523a5be2c` (HEAD) | gitlink `workflow/mg2sc` → `e744c2c` (VERIFIED, `git ls-tree`) |
| sinanugur/cellsnake | `7adf99ddbcf7e69eb28dade218d062b265247434` | pins scrna-workflow `5d9279c`; microbiome files identical to HEAD |
| DerrickWood/kraken2 | HEAD `8e35df5` = `v2.18.0-RC1` (2026-10-19) | **pin stable `2.17.2` = `8c190b1` (2026-09-08)**, not the RC |

Paths below are relative to the repository named. For ViralScan, they are relative to `src/viralscan/` unless the
path says otherwise.

---

## 1. Executive recommendation

Build a **standalone post-hoc command**, `viralscan discover`. It reads a *completed* ViralScan run directory and writes
only into a separate output directory with its own manifest and fingerprint. The Snakemake integration is P1, considered
only after the M2 correctness gates pass. Even then, the primary run's `config.yaml`, options fingerprint and DAG stay
unchanged unless discovery is explicitly enabled.

The minimal pipeline:

```text
run_dir/ (complete: run_complete.json + run_manifest.json)
  │
  ├─ input resolution (fail closed)
  │    host_filtered/R1,R2  (primary already ran --host-filter starsolo)    → reuse
  │    else: explicit --host-subtract starsolo --host-index DIR             → discovery-only STAR pass into -o
  │    else: refuse (no implicit host subtraction, no unsubtracted mode in MVP)
  │
  ├─ stream the host-unmapped pair once:
  │    ordinal n, canonical read_id, raw_cb, raw_umi (R1, chemistry geometry)  → reads.tsv.gz
  │    R2 (cDNA mate only) → kraken_input.fa  with header ">n"
  │
  ├─ kraken2 --db DB --confidence C --minimum-hit-groups H --report ... --output kraken.out kraken_input.fa
  │
  ├─ join by integer ordinal; each ordinal exactly once (bool array)            → classifications.tsv.gz
  ├─ barcode correction: same bustools binary + the primary's resolved on-list → corrected_cb
  ├─ UMI groups: exact (sample, corrected_cb, raw_umi); taxon = LCA over classified reads (DB's own tree)
  ├─ cell universe: primary results/called_cells.tsv (barcode space checked, overlap reported)
  │
  └─ outputs: classifications, umi_groups, taxonomy, cell_taxon.h5ad (cells × assigned taxids),
              offcell_taxon.tsv.gz, candidates.tsv (observation state only in MVP), qc.json, provenance.json
```

Why this is the smallest credible design:

- It cannot alter primary inputs, counts, calls or fingerprints, because it never writes into the run directory and
  adds no option to the `run` parser.
- It reuses three verified ViralScan facilities instead of copying mg2sc's:
  - fragment-ID canonicalisation and paired-FASTQ validation (`scripts/host_filter.py:77-112`);
  - the STARsolo host pass (`scripts/host_filter.py:391-467`);
  - the exact bustools correction chain (`evidence.py:429-503`).
- The ordinal-header join removes every identity hazard found in mg2sc and Cellsnake: positional zip, `grep ^@`,
  `fgrep -w`, and last-record-wins awk.

What it does not do: it does not touch H5AD count layers, `detected`, evidence tiers, the v3 counting contract or the
reference panel. A Kraken2 result never creates an infection claim.

---

## 2. Repository audit

### 2.1 ViralScan (VERIFIED unless labelled)

**Host subtraction, `scripts/host_filter.py`**

- STARsolo is the only host-filter aligner. kallisto is rejected because BUS loses read identity (`:46-56`).
- The command (`_starsolo_filter`, `:391-467`):
  - `STAR --readFilesIn R2 R1` with `soloType CB_UMI_Simple` (`starsolo_barcode_args`, `:352-377`), so only the cDNA
    mate (R2) is aligned and R1 is the barcode read;
  - `--soloFeatures GeneFull --outSAMtype None --outReadsUnmapped Fastx`;
  - one of two pinned parameter sets (`STAR_FILTER_PARAM_SETS`, `:346-349`).
- `pinned` (`:295-319`): `outFilterMismatchNmax 4`, `outFilterMatchNminOverLread 0.9`, `outFilterMultimapNmax 20`,
  `outFilterMismatchNoverReadLmax 0.05`, `alignIntronMin 20`, `alignSJoverhangMin 8`, `alignSJDBoverhangMin 8`.
- `star-default` (`:322-343`) passes the STAR 2.7.11b defaults explicitly.
- `Unmapped.out.mate1` (cDNA) and `mate2` (barcode) become `host_filtered/R2.fastq.gz` and `R1.fastq.gz`. Pairing is
  checked on output: `validate_paired_fastq_ids`, and retained ≤ input (`:450-467`).
- Audit (`_write_filter_audit`, `:123-188`):
  - `host_filter_audit.tsv` holds input, retained, removed (= input − retained), `pct_retained`, the param set and every
    STAR argument;
  - `fragment_lineage.tsv.gz` holds `read_id, filter_decision, reason` for **every** input fragment, retained and
    removed;
  - there is **no CB/UMI** in the lineage, and the removal sub-reason is unrecoverable because the SAM is discarded.
- **Host filtering changes the kallisto input.** `runconfig.py:504-518` sets `kb_r1/kb_r2` to `read_filtered/`, then
  `host_filtered/`, then the raw sample, and `Snakefile:39-53` makes `kb_count` depend on them. Turning on the host
  filter in the primary run is therefore a quantification change.
- STAR's barcode on-list is the user's `-w` or the literal `None`, so no CB correction happens without `-w`
  (`starsolo_barcode_args`, `:376`). `kb count` instead uses `kb_whitelist_arg` (`chemistry.py`), which falls back to
  the bundled on-list.
- Under host filtering, EmptyDrops runs on `host_filtered/star_tmp/Solo.out/GeneFull/raw` (`scripts/cellcalling.py`
  `solo_raw_dir`, `:322`).

**Read-artefact filter, `scripts/read_filter.py`** (opt-in `--read-filter artefact|tso-trim`)

- It runs **after** the host filter and **before** kb (`Snakefile:57-83`).
- `artefact` mode drops these pairs:
  - `r1_tso`: TSO in the barcode span;
  - `r2_reagent`: TSO, polyT or TruSeq chimeras;
  - `r2_no_complex_body`: no complex body of 20 nt or more, which also removes low-complexity viral reads.
- Its lineage holds `read_id` and the decision only.

**Chemistry, `chemistry.py`**

- `CHEMISTRIES` (`:46-55`): 10xv1 14/10, 10xv2 16/10, 10xv3 16/12, 10xv4 16/12, dropseq 12/8.
- `cb_umi_geometry` (`:113-138`) accepts an explicit `bc:umi:seq` string only when CB is at R1 offset 0 and the UMI
  follows it directly. Anything else raises `ChemistryError` (P4-06).
- There is no reusable `(raw_cb, raw_umi)` extractor; callers slice R1 inline.
- 10xv1 keeps its UMI on a separate read (`:58-59` comment). Under `cb_umi_geometry`'s fixed R1 layout it would be
  mis-sliced, so discovery must reject it (PROPOSED).

**Barcode correction**

- ViralScan has no Python corrector. Correction is always `bustools correct`: inside `kb count`, and in
  `_replay_quant_bus` (`evidence.py:429-503`). The replay does `kallisto bus -n` (read number in the flag column), then
  `bustools sort → correct -w <plain on-list> → sort`.
- `whitelist_preflight.py` does exact set membership only (`whitelist_match_rate`).
- GEM-well handling is suffix stripping only (`scripts/cellcalling.py` `_strip_suffix`). The namespace is the
  per-sample run directory.

**Counting contract, `docs/v3_counting_contract.md`**

- A molecule is a corrected `(CB, UMI)` from the corrected BUS stream. BUS `count` is audit-only. There is no UMI
  error-correction step.
- Invariants: `unique + ambiguous + unresolved = input molecules`. Order and chunk invariance are enforced by an oracle
  property test (`tests/test_multimapping_properties.py`).
- The contract states that fragments discarded by STAR host filtering are unobservable in the BUS stream (SW-07 open).

**Read lineage**

- `parse_flagged_target_bus` (`evidence.py:303-335`) and `extract_exact_reads_by_number` (`:343-401`) map read numbers
  to *corrected* CB/UB, but only for reads captured from BUS, which means only reads that pseudoalign to the kallisto
  index.
- They write FASTA qnames `<CB>_<UMI>_<n>|<read_id>`, which `_cb_umi` (`:2003-2009`) parses.
- Reads absent from the index have no lineage.

**Cell universe**

- `results/called_cells.tsv` (`scripts/cellcalling.py:685-698`) has a `barcode` header and one suffix-stripped barcode
  per row.
- The method is `external`, `emptydrops` (default), `knee` or `none`, plus a provenance receipt
  (`_write_input_receipt`, `:518`).

**Fingerprint and run safety, `run_safety.py`**

- `_HASHED_KEYS` (`:111-120`) are: schema_version, viralscan_version, quantification_unit, allocation_method,
  reference_fingerprint, cell_calling, input_fingerprints and options.
- `options` are every scalar `vars(args)` entry, except entries in `omit_when_unset` whose value is `None` (`:200-239`;
  pattern comment "AUD-06/AUD-07/B-16: unset keeps the pre-existing fingerprint").
- `input_fingerprints` hash the sample FASTQs and references (`build_run_manifest`, `:148-177`).
- **The primary manifest does not hash `host_filtered/` outputs** (no `host_filtered` reference in `run_safety.py`),
  so a discovery run cannot verify them against the primary. It can only record its own hash at discovery time.
- `record_manifest_block` (`:388`) writes descriptive, unhashed manifest blocks. `MARKER_ARTIFACTS` (`:25`) limits what
  `run_complete.json` hashes.
- Precedents for a post-hoc command on a completed run: `evidence --run-dir -o`, `rerun-multimap --run-dir -o`,
  `rerun-programs --run-dir` (`cli/parsers.py:288-330, 454-541`).
- A new `RunConfig` field changes `config.yaml`. `config.yaml` is an input of `kb_count` and `host_filter`, so it can
  trigger primary reruns (SW-22 hazard). A standalone command adds no `RunConfig` field.

**Detection and tiers**

- `detected` = `infected_called >= detection_min_cells and called_umi >= detection_threshold·(1−tol)`
  (`scripts/detection.py:997-1000`).
- `evidence_tiers.assign_tier` has **no production caller** (grep).
- Ambient context (`soup_screen` `:738`, `_ambient_context` `:823`) is descriptive and never feeds `detected`.

**Evidence engine reuse**

- Read-set and FASTA agnostic (INFERRED from signatures):
  - `align_reads_to_viral` (`evidence.py:726`, minimap2 `sr` plus secondaries);
  - `write_competitive_fasta` (`:746`, `HOST|`/`VIRUS|` prefixes);
  - `coverage_table` (`:1284`), `alignment_qc_table` (`:1457`), `per_cell_alignment_qc` (`:1507`),
    `deduplicate_bam` (`:2253`), `competitive_blast_identity` (`:2372`);
  - `coverage_map_figure` (`:1694`; GTF optional).
- Coupled to kallisto or the panel: `resolve_viral_target`, `parse_flagged_target_bus`, `replay_*`,
  `viral_equivalence_classes`, `molecule_verdict`, `auto_evidence.candidate_genera` (Anelloviridae only).
- Adapter contract: reads FASTA named `<CB>_<UMI>_<n>` plus a reference FASTA with the `VIRUS|` prefix.

**Environment**

- Kraken2 is absent from `environment.yml`, `pixi.toml`, the Dockerfile, the Singularity definitions and the
  pyproject extras.
- `_check_required_tools` requires only `kb` and `snakemake`. Other tools are checked per flag.
- `validation.FULL_TOOLS` drives `doctor --profile full`. Kraken2 must **not** be added there.

### 2.2 mg2sc, original (`a8c01f0`, VERIFIED)

`k2sc.py` line numbers below use the **fork's** numbering (its added shebang shifts them by +1 versus the
original); the code is otherwise identical (§2.3).

- `src/scMeG-kraken.py:61` runs `samtools view -b -f 4`: read-level unmapped flag, no `-F 0x900`, no mate or MAPQ
  logic.
- `:68` runs `samtools fastq -n` (which drops 0x900 by default, so record counts can diverge; INFERRED). `:76` runs
  Kraken2 single-end with no `--paired`.
- Every command is a `Popen(..., shell=True).wait()` string with unquoted paths and **return codes ignored**
  (`:61-79`).
- `src/k2sc.py:71` joins `zip(pysam.AlignmentFile(bam), open(kraken))` **positionally**:
  - a read-ID mismatch is logged and skipped (`:77-80`), but both iterators have already advanced, so one
    missing or extra record misassigns everything after it;
  - `zip` truncates silently.
- `:83-89`: CB and UB come from BAM tags (corrected by Cell Ranger or STARsolo), and a bare `except` drops reads
  without them. INFERRED: Cell Ranger sets `UB` only on reads it can assign, so many unmapped reads may be dropped
  (untested).
- `:93-100`: a non-numeric taxid goes through a regex, and the fallback calls `sys.exit()` without `import sys`, so it
  raises `NameError`.
- `:119-121` `most_frequent = max(set(List), key=List.count)`: ties resolve by set iteration order. For `str`
  elements that order depends on per-process hash randomisation, so ties break **non-deterministically across runs**
  (INFERRED from CPython).
- Unclassified `'0'` takes part in the vote and becomes a matrix feature (`taxdict = {'0': 'unclassified'}`, `:180`).
- Votes are keyed on `(CB, UB)` only (`:103-112`), so two transcripts sharing a CB and UMI merge into one vote.
  Parent and child taxids are separate labels with no LCA.
- `dict2lists` (`:145-172`) builds the taxa × cells sparse matrix with an O(n²) `list.index`. Names come from a
  pre-generated `inspect.txt` in the DB dir (`:21, 174-195`), which Kraken2 does not create.
- `collapse_taxonomy.py`:
  - rank collapse is name-based;
  - UMIs assigned above the target rank are silently dropped;
  - rank sub-codes are compared as strings.

### 2.3 Cellsnake fork of mg2sc (`e744c2c`, VERIFIED by `diff` against the original)

- `k2sc.py` and `collapse_taxonomy.py` differ only by a shebang and one comment. Join, vote and aggregation are
  **unchanged**.
- `scMeG-kraken.py` adds:
  - `--confidence` 0.05 and `--minimum-hit-groups` 4;
  - `--complexity` 0, run as `fastp --stdin --stdout --complexity_threshold C -Q -L -A -y`;
  - `--bowtie`: a bowtie2 `--very-fast` host DB whose hits are excluded;
  - `--predb`: a second Kraken2 host DB at a hard-coded confidence of 0.1 and 12 hit groups, which silently overrides
    `--bowtie`;
  - `--report-minimizer-data` and `--classified-out`, both always on.
- The ID list `cat fq | grep ^@ | tr -d ^@` (`:99`) also matches FASTQ quality lines that start with `@`, and deletes
  `@` and `^` anywhere in the name.
- Re-selection uses `fgrep -w -f ids` against whole SAM and Kraken lines, then `awk 'NR==FNR{l[$1]=$0}…'`
  (`:127-133`):
  - substring matches are possible;
  - the **last record per read name wins**;
  - it needs `sponge` (moreutils), which is not declared.
- Return codes are still ignored. The fork's own comment on the zip-skip warning ("#due to fastp") acknowledges the
  desynchronisation that fastp causes.

### 2.4 scrna-workflow and Cellsnake (VERIFIED)

- Only `workflow/rules/microbiome.smk` is included (`workflow/Snakefile:429`). `microbiome2.smk` is dead code.
- Input is Cell Ranger `possorted_genome_bam.bam` (`microbiome.smk:2-8`), not STAR FASTQ output.
- `run_kraken` declares `threads: 10` and no resources, conda env or container. It does not use `--memory-mapping`,
  so every sample loads the database into RAM.
- Defaults conflict: `config.yaml:50-51` and `cellsnake/command_line.py:142-143` say 0.05 and 4; the
  `workflow/Snakefile:167-168` fallback says 0.01 and 3.
- The DB check `if not os.path.exists(db) and not os.path.isfile(db + "/inspect.txt")` (`command_line.py:187`) uses
  `and`, so a DB folder without `inspect.txt` passes. No DB version or checksum is recorded. fastp, bowtie2 and
  moreutils are not in `environment.yml`.
- `scrna-kraken2-collapse.py` collapses by rank name (default genus).
- `scrna-kraken2-data-parser.R`:
  - applies `min.cells=1, min.features=3`, so cells with fewer than 3 taxa are dropped;
  - removes `Homo*` **by name prefix**, so higher human ranks survive;
  - joins to the Seurat object by barcode string.

### 2.5 Kraken2 (VERIFIED at the pinned tag `2.17.2` = `8c190b1`; manual line numbers are identical at `8e35df5`)

- **`--confidence` is not a probability.** The manual (`docs/MANUAL.markdown:997-1003`) says "we have not yet
  developed a confidence score with a probabilistic interpretation". The score is C/Q, and the label climbs the tree
  until C/Q reaches the threshold.
  - The external-audit agent reports a code/doc discrepancy: the code's denominator includes ambiguous k-mers. This is
    UNRESOLVED and must be settled by the M1 synthetic test.
- `--minimum-hit-groups` (wrapper default 2) voids a call that lacks enough distinct minimizer hit groups.
  `--report-minimizer-data` adds per-taxon minimizer counts to the report, per run and not per cell.
  `--memory-mapping` avoids loading the DB into RAM.
- Output columns: `C/U`, sequence ID (first whitespace token; `/1`/`/2` stripped only with `--paired`), taxid, length,
  LCA hitlist.
- Output order equals input order regardless of threads: a priority queue reorders blocks by `block_id`
  (`src/classify.cc:499-686` at 2.17.2). The design does **not** rely on this.
- `k2 inspect --skip-counts` "only print[s] database summary statistics" (`MANUAL.markdown:307`), so a full taxonomy
  tree with ranks requires `k2 inspect` without `--skip-counts`, or parsing `taxo.k2d`.
- DB files: `hash.k2d`, `opts.k2d` (k, l, mask, db_version), `taxo.k2d`. There is no embedded build date or source,
  so provenance must be supplied by hashing.

---

## 3. Critical review: retain, modify, reject

| mg2sc / Cellsnake idea | Verdict | Replacement |
|---|---|---|
| Classify host-unmapped reads with Kraken2 and attribute them to cells | **Retain** (the core principle) | — |
| `samtools view -f 4` on a Cell Ranger BAM | Reject for MVP | STARsolo `--outReadsUnmapped Fastx` on FASTQ (fragment = cDNA mate) |
| Positional `zip(BAM, kraken)` join | **Reject** | integer-ordinal header and a bool-array exactly-once check |
| `grep ^@` / `fgrep -w` / awk last-wins / `sponge` | **Reject** | streamed Python; no shell pipelines |
| `shell=True`, return codes ignored | **Reject** | `subprocess.run([...], check=True)` via `evidence._run` |
| CB/UB from BAM tags | Reject | raw CB/UMI from R1 through `cb_umi_geometry`; correction with the same `bustools correct` and on-list |
| `max(set(labels), key=labels.count)` vote | **Reject** | LCA over classified reads using the DB's own tree (order-independent) |
| Unclassified (`0`) as a voter and as a feature | **Reject** | unclassified groups counted in QC; never a matrix variable |
| Name-based rank collapse; drop above-rank UMIs | **Reject** | taxid-based collapse; above-rank groups go to an explicit `above_rank` count |
| Remove `Homo*` by name | Reject | taxid lineage check against 9606 and its ancestors; reported as residual host |
| Dense matrices / `min.features=3` cell drop | Reject | sparse cells × taxa; called cells with zero counts kept |
| `--confidence 0.05 --minimum-hit-groups 4` defaults | **Do not adopt** | Kraken2's own defaults (0.0 / 2) as the baseline, tuned only by the M5 grid |
| fastp low-complexity filter | Modify → P1 | report `read_filter` decisions as an annotation; no extra tool in the MVP |
| `--predb` / `--bowtie` second host filter | Defer (P1/P2) | broad DB with human sequences gives a 9606 assignment, which is reported |
| `--report-minimizer-data` | **Retain** | per-run distinct-minimizer counts, used in candidate QC |
| Seurat/h5seurat conversion | **Remove** | AnnData only |

---

## 4. Architecture decision record (ADR-DISC-1)

**Context.** Kraken2 discovery needs host-subtracted reads. ViralScan's host filter is optional and *changes the primary
quantification input* when enabled. Primary fingerprints hash every `run` option.

**Alternatives.**

| | A. Parallel Snakemake branch now | B. Kraken2 on kallisto-unassigned reads | C. Standalone post-hoc command | D. C, then an optional Snakemake leaf |
|---|---|---|---|---|
| Sensitivity | same as C | **lower and unpredictable**: a microbial read that shares k-mers with a host or viral target is assigned and never reaches Kraken2; with a virus-only index almost everything is "unassigned" and nothing is subtracted | STAR host subtraction (measured loss, gate 12) | = C |
| Specificity | = C | worse: kallisto-unassigned ≠ host-absent (a k-mer D-list / pseudoalignment cannot see diverged host, HOST-04) | host reads removed by spliced alignment | = C |
| Barcode correctness | = C | needs BUS replay of unassigned reads, which BUS does not record | same bustools correct | = C |
| CPU/RAM/I/O | DB load inside the primary DAG; contends with kb | extra BUS replay | separate job; DB load once per invocation | = C |
| Backward compatibility | config.yaml and `run` options change → fingerprint and SW-22 rerun risk | same | **none touched** | none until P1, then via `omit_when_unset` |
| Testing | DAG + fingerprint + science together | hard | unit + CLI tests only | staged |
| Maintenance | highest | high | lowest | moderate, deferred |

**Decision (PROPOSED).** D, starting as C. M1–M2 ship `viralscan discover` only. M3 may add a Snakemake leaf rule
(template: `host_capture`) if a user needs it inside `viralscan run`, and it must pass gates G1–G3.

**Why STAR-unmapped and not kallisto-unassigned.** "Unassigned by kallisto" is defined against whatever index was
built: a combined, virus-only or host-only index each gives a different residual. It is not a statement about host
origin. STAR's spliced genome alignment is the only host-subtraction step in ViralScan with exact read identity and
pinned, audited parameters (HOST-01/02). Kraken2 discovery therefore consumes STAR-unmapped cDNA mates. The kallisto
result is used afterwards for comparison only (§9.4).

**Rejected.**

- A: premature coupling; a fingerprint and rerun risk with no scientific gain before M2.
- B: the denominator depends on the index, and subtraction is not host-based.

---

## 5. Data flow and boundaries

### 5.1 Host-unmapped input definition (P0)

Fragment = one `(R1, R2)` record pair with canonical ID `canonical_read_id(header)` (`scripts/host_filter.py:77-84`).
For every supported chemistry, R2 is the cDNA mate and R1 is barcode+UMI (`cb_umi_geometry`). **Only R2 is classified.**

A fragment is *host-unmapped* when STAR (2.7.11b, param set P ∈ {`pinned`, `star-default`}) emits it in
`Unmapped.out.mate1/2`. Under `soloType CB_UMI_Simple` alignment is single-end on the cDNA mate, so read-level and
fragment-level unmapped are the same thing (INFERRED from the two-file Solo convention; DISC-03 verifies it on a fixture).

| Case | Fate under `pinned` | Note |
|---|---|---|
| cDNA mate aligns ≥ 0.9 of length with ≤ 4 mismatches, ≤ 20 loci | removed | includes spliced and soft-clipped alignments that pass |
| host multimapper > 20 loci | **retained** ("too many loci") | repeats/EVEs reach Kraken2 → expect host-like calls |
| partial alignment < 0.9 of read | **retained** | host–virus chimeras and diverged host reach Kraken2 |
| host–virus homologous region / integrated virus (EVE) | removed if it aligns to host | **host subtraction loses true viral signal here**; measured, never assumed (gate 12) |
| barcode read | never classified | |
| barcode not on-list | still in `Unmapped.out` (alignment is CB-independent; INFERRED) | handled by the barcode stage, not dropped silently |

`--outSAMtype None` means the removal sub-reason is unknown. Discovery inherits this limit (HOST-02) and does not try
to fix it in the MVP.

### 5.2 Input resolution (fail closed, PROPOSED)

1. Require `run_dir/run_complete.json` and `run_manifest.json`. Refuse a running or incomplete run.
2. If `run_dir/host_filtered/{R1,R2}.fastq.gz`, `host_filter_audit.tsv` and `fragment_lineage.tsv.gz` exist:
   - check pair synchronisation;
   - check that the retained count equals the audit's `retained_host_unmapped` and the lineage's retained count;
   - record their sha256 as **discovery-time** hashes (`input_source: primary_host_filtered`,
     `primary_verified: false`, because the primary manifest does not hash them; §2.1);
   - use STAR param set P from the audit.
3. Otherwise, only with an explicit `--host-subtract starsolo --host-index DIR [--host-filter-star-params P]`:
   - run `_starsolo_filter` (reused, not copied) on the primary's sample FASTQs, with outputs under `-o/host_subtraction/`;
   - check that the sample FASTQ hashes equal the primary manifest's `input_fingerprints`, or refuse;
   - record `input_source: discovery_star_pass`.
4. Otherwise, refuse with a message naming both options. FASTQs deleted after a run (Track A policy) produce a clear
   "inputs absent" error, never a partial run.

**Pre- versus post-read-filter (PROPOSED).** Classify the **host-filtered pair (pre-read-filter)**. `read_filter`
`r2_no_complex_body` also removes low-complexity viral reads, and its purpose is kb artefact control, not discovery.
When the primary ran `--read-filter`, its `read_filtered/fragment_lineage.tsv.gz` decision is joined by `read_id` as
the annotation column `primary_read_filter`. A Kraken2 assignment carried only by fragments the read filter would have
dropped is a QC flag (`artefact_supported_only`), not a silent loss. The consequence: discovery sees a larger universe
than kb did, and the difference is reported.

### 5.3 Stages and their units

| Stage | Unit | Universe and audit |
|---|---|---|
| S0 input | fragment | `n_input` = host-unmapped pairs |
| S1 barcode parse | fragment | `n_input = n_parsed + n_bad_geometry` (short R1, N in CB/UMI) |
| S2 Kraken2 | fragment (R2) | `n_parsed = n_classified + n_unclassified`; every ordinal exactly once |
| S3 barcode correction | fragment | `n_parsed = n_cb_onlist + n_cb_corrected + n_cb_uncorrectable` (+ `n_cb_no_onlist` mode) |
| S4 UMI grouping | (sample, corrected CB, UMI) group | built from S3 on-list/corrected fragments only |
| S5 cell partition | group | `groups = called + offcell` |
| S6 matrix | resolved called-cell groups at their assigned taxid | one unit per group |

---

## 6. Identity and counting contracts

### 6.1 Record model (PROPOSED; simpler than the brief's list where the run directory already fixes a field)

`reads.tsv.gz` (one row per S0 fragment; primary key `ordinal`):

| field | type | note |
|---|---|---|
| `ordinal` | uint64 | 0-based record index in the resolved host-unmapped pair; Kraken2 sequence ID |
| `read_id` | str | `canonical_read_id`; unique (checked) |
| `raw_cb` | str | R1[0:cb_len]; empty if geometry failed |
| `raw_umi` | str | R1[cb_len:cb_len+umi_len] |
| `corrected_cb` | str | empty if uncorrectable |
| `cb_status` | enum | `onlist`, `corrected`, `uncorrectable`, `bad_geometry`, `no_onlist` |
| `kraken_status` | enum | `C`, `U` |
| `taxid` | uint32 | 0 if U |
| `primary_read_filter` | enum | `retained`, `removed:<reason>`, `not_run` |

`sample_id`, `chemistry`, `input_stage`, `database_id` and the parameters are constant per run and live in
`provenance.json`, not on every row. Full Kraken2 hitlists stay in `kraken.out.gz` (retained, keyed by ordinal).

### 6.2 Identity invariants (fail closed; each is a test)

- I1. Read IDs are unique across the resolved pair; mates match per record (`iter_paired_fastq_ids`).
- I2. Kraken2 output contains each ordinal in `[0, n_parsed)` **exactly once** (bool array; a duplicate, a missing
  ordinal or an out-of-range ordinal aborts).
- I3. A non-numeric Kraken2 sequence ID aborts (it means the wrong input or `--use-names`).
- I4. `retained` from the host-filter audit = `n_input` (when reusing primary outputs).
- I5. Truncated or gzip-corrupt input aborts. No partial outputs are published (§11.3).

### 6.3 Barcode spaces (P0; three exist)

| Space | Producer | Correction |
|---|---|---|
| raw | R1 bytes | none |
| kb/bustools | `kb count` with `kb_whitelist_arg` (user `-w`, else bundled on-list) | Hamming-1 to on-list |
| STARsolo | host filter with `-w` or literal `None` | none without `-w` |

**Policy (PROPOSED).**

- `corrected_cb` is produced in the **kb/bustools space**, using the primary's resolved on-list, read from the primary
  `config.yaml` and the `kb_whitelist_arg` resolution and recorded with its sha256, and the same `bustools` binary
  (version recorded).
- Mechanism (spike DISC-05):
  - write the S1 fragments as BUS text `(raw_cb, raw_umi, ec=0, count=1, flag=ordinal)`;
  - run `bustools fromtext → sort → correct -w → text -f` and read `corrected_cb` by flag.
  - If the flag column does not survive, fall back to a Python Hamming-1 corrector with a property test against
    `bustools correct`.
- Reconciliation gate (DISC-05). The primary `kb count` never runs with `-n`, so the run directory holds no read-number
  lineage. The gate runs `_replay_quant_bus` once on the kb input pair, which needs the primary kallisto index. It is a
  validation step, not part of every `discover` run:
  - on fragments that also appear in the primary BUS lineage, join by **canonical `read_id`**, not by read number,
    because read numbers shift when `read_filter` ran;
  - discovery `corrected_cb` must equal the primary corrected CB for ≥ 99.9 % of joined fragments;
  - every mismatch is listed.
- Called-cell space:
  - read the cell-calling receipt to determine whether `called_cells.tsv` came from kb space (EmptyDrops on the kb
    matrix, or external) or STARsolo space (two-step EmptyDrops on `GeneFull/raw`);
  - report the overlap rate of discovery `corrected_cb` with the called set;
  - if the spaces differ (STARsolo with no `-w`), report it and use the called set as given, with no re-correction.
    This is a documented limitation, not silent matching.
- **No `-w` and no on-list (Drop-seq):**
  - `cb_status = no_onlist` and `corrected_cb = raw_cb`, marked uncorrected in `qc.json`;
  - called-cell membership is by exact string.
- **Namespace:** one sample per run directory. Barcodes keep no `-N` suffix; suffixes are stripped exactly as
  `cellcalling._strip_suffix` does. Multi-sample aggregation is out of scope; a later tool must prefix with `sample_id`.
- **Unsupported chemistry:** anything `cb_umi_geometry` rejects, plus `10xv1`, aborts.

### 6.4 UMI policy (PROPOSED MVP)

- Group key: exact `(corrected_cb, raw_umi)`. There is no UMI error correction, which matches the v3 contract, where
  UMI identity is the exact string after bustools.
- Groups with an `N` in the UMI are `bad_geometry`.
- Group taxon is **taxonomy-aware**:
  - let `T` = the set of taxids of *classified* reads in the group;
  - `resolved_taxid = LCA(T)` on the DB's own tree (§6.6);
  - unclassified reads do not vote but are counted (`n_reads_unclassified`).
- Group classes:
  - `unclassified`: T empty.
  - `resolved`: `LCA(T) ≠ root(1)`, `LCA(T) ≠` cellular-organisms-level, and LCA ≠ 0.
  - `ambiguous`: the LCA is root, or `|T| > 1` and the taxids are on different lineages *and* the LCA is above the
    configured `--min-resolution-rank` (default: no floor; recorded).
  - Also recorded:
    - `lineage_consistent` (bool): all of T lie on one root-to-leaf path. For example, {561, 562} is consistent;
      `resolved_taxid` is then 561 (the LCA) and `deepest_consistent_taxid` = 562 is stored for audit only;
    - `n_distinct_taxids` and `conflict` (bool): `conflict=True` when T is not lineage-consistent. A conflicting
      group can still be `resolved` (for example, two strains under one species resolve to that species); the
      guarantee is only that it is never counted below its LCA.
- Two distinct transcripts sharing a CB/UMI, from different organisms: T spans lineages, the LCA climbs, and
  `conflict=True`. It counts as one unit at a high rank or as ambiguous, **never** as two species units. The MVP does
  not split groups.
- Determinism: LCA is a set function, so it is independent of read order, chunk order and thread count. No tie-break
  exists or is needed.

**Rejected for the MVP:** majority vote (non-deterministic or arbitrary ties, and it hides conflict); consensus with
fraction thresholds (an extra uncalibrated parameter). A `deepest_consistent` policy can be offered later as a declared
sensitivity variant (P1), never as the default.

### 6.5 Counting invariants

```text
n_input        = n_parsed + n_bad_geometry
n_parsed       = n_classified + n_unclassified
n_parsed       = n_onlist + n_corrected + n_uncorrectable  (+ n_no_onlist)
groups         = resolved + ambiguous + unclassified          (built from on-list/corrected/no_onlist fragments)
groups         = called_groups + offcell_groups
sum(cell_taxon.X) = resolved ∩ called groups                  (one unit per group, at assigned taxid)
collapse(rank R): sum over R-columns + above_rank_R = sum(cell_taxon.X)
```

No group contributes more than one unit at any resolution:

- the matrix holds groups at their `resolved_taxid`;
- a rank-R view maps each taxid to its rank-R ancestor, or to `above_rank_R` when the taxid sits above R;
- groups are never pushed down to a descendant;
- species and genus views are separate projections of the same groups, never summed together.

### 6.6 Taxonomy source

The tree, with ranks and names, comes **from the same DB**:

- `k2 inspect --db DB` (full, not `--skip-counts`), parsed by indentation;
- or a direct `taxo.k2d` parser (decided in DISC-04);
- cached by DB identity (a full `k2 inspect` scans the hash table), copied to `-o/taxonomy.tsv` and hashed into
  `provenance.json`.

`k2` exists at 2.17.2 (`scripts/k2`). For an older Kraken2, fall back to `kraken2-inspect`; the minimum
supported version is decided in DISC-01.

Never use a separate NCBI taxdump: merged or deleted taxids can disagree with the DB, and many DBs are built with
`--clean`, which removes `taxonomy/`. Ranks include `no rank` and intermediate codes such as `S1` or `G2`. Rank
collapse uses the canonical rank letters only.

---

## 7. CLI, config and API (PROPOSED)

```bash
# reuse the primary's STAR host-unmapped pair
viralscan discover --run-dir out/S1/ -o out/S1_discovery/ \
    --kraken2-db /shared/k2/k2_pluspf_20260XYZ/ \
    --confidence 0.0 --minimum-hit-groups 2 --memory-mapping --threads 16

# primary ran without host filtering → explicit discovery-only STAR pass
viralscan discover --run-dir out/S1/ -o out/S1_discovery/ \
    --kraken2-db DB --host-subtract starsolo --host-index /refs/star_grch38/ \
    --host-filter-star-params pinned --threads 16

# post-hoc rank view, no reclassification
viralscan discover-collapse out/S1_discovery/ --rank genus -o out/S1_discovery/genus.h5ad
```

**Rules**

- `-o` is required and must not be inside `--run-dir` (precedent: `rerun-multimap`).
- No option is added to `viralscan run` in M0–M2.
- Kraken2, and STAR only when `--host-subtract` is given, are checked with `shutil.which` at `discover` start.
- `kraken2` is called as `["kraken2", "--db", DB, "--threads", N, "--confidence", C, "--minimum-hit-groups", H,
  "--report", R, "--report-minimizer-data", "--output", O, (--memory-mapping), FASTA]` through `evidence._run`
  (list form, stderr surfaced).
- No `--use-names`, `--paired` or `--quick`.
- **DB validation:**
  - `hash.k2d`, `opts.k2d` and `taxo.k2d` must exist;
  - `k2 inspect --skip-counts` must succeed;
  - there is no download path at all, ever.
- **Environment:**
  - kraken2 goes in an optional pixi feature / `environment.discovery.yml`, not in `environment.yml`, `FULL_TOOLS` or
    the Dockerfile;
  - `doctor --profile discovery` lists it.

**Module layout:** `src/viralscan/discovery.py` (pure functions: parse, join, LCA, group, collapse, invariants)
plus `menu.run_discover`. It is not a Snakemake `script:` module, so the no-`__future__` rule does not apply. A third
file is added only if `discovery.py` passes about 600 lines.

---

## 8. Output schema (PROPOSED)

`-o/` layout. All files are written to `-o/.partial/`, then `os.replace`d. `discovery_complete.json` is written last.

| File | Key | Content |
|---|---|---|
| `reads.tsv.gz` | `ordinal` | §6.1 |
| `kraken.out.gz` | ordinal | raw Kraken2 output (audit) |
| `kraken_report.txt` | taxid | `--report --report-minimizer-data` |
| `umi_groups.tsv.gz` | (`corrected_cb`, `umi`) | `n_reads`, `n_reads_classified`, `n_reads_unclassified`, `taxids` (sorted, `;`), `resolved_taxid`, `group_class`, `lineage_consistent`, `deepest_consistent_taxid`, `conflict`, `called` |
| `taxonomy.tsv` | `taxid` | `parent_taxid`, `rank`, `name`, `depth` (DB-derived) |
| `cell_taxon.h5ad` | obs = called cells (all, zero rows kept), var = `resolved_taxid` | `X` = int32 sparse CSR of resolved groups; `var`: name, rank, parent, `is_host_lineage`; `uns`: provenance digest, invariant counts, `unit="resolved_umi_group"` |
| `offcell_taxon.tsv.gz` | (`corrected_cb`, taxid) | resolved groups on non-called barcodes (ambient and empty-droplet evidence) |
| `candidates.tsv` | taxid | §9; state ∈ {`observation`, `prioritized`, `sequence_supported`} |
| `qc.json` | — | every §6.5 count, barcode-space overlap, host-lineage fraction, `artefact_supported_only` |
| `provenance.json` | — | §11.2 |
| `discovery_complete.json` | — | output sha256 list + discovery fingerprint |

`X` is never added to, merged with, or written beside the primary H5AD layers. Rank views are produced by
`discover-collapse` from `cell_taxon.h5ad` plus `taxonomy.tsv`; Kraken2 is not run again.

---

## 9. False-positive control and candidate states

### 9.1 Known false-positive routes and MVP handling

| Route | MVP handling |
|---|---|
| Residual human reads (>20-locus multimappers, partials) | broad DB with human: assignments in the 9606 lineage → `is_host_lineage`, never candidates; **viral-only DB: candidates capped at `observation`** (residual host is invisible) |
| Human–microbe homology, EVEs, repeats, low complexity | LCA climbs; `deepest_consistent` audit; `primary_read_filter` annotation; minimizer-distinct ratio from the report |
| Adapter/TSO chimeras | `primary_read_filter` annotation; `artefact_supported_only` flag |
| Reference/reagent contamination, kitome | no built-in facility (VERIFIED: none in ViralScan); `--negative-control DIR` (P1) compares against another discovery output |
| Ambient RNA / empty droplets | off-cell versus called-cell group distribution per taxon (`offcell_taxon.tsv.gz`); a called-cell fraction not above off-cell is flagged |
| Barcode swapping / carry-over | per-taxon cell spread versus UMI concentration reported; no correction |
| Batch artefacts | out of scope for a single run; recorded as a limitation |

### 9.2 States (never equal to ViralScan tiers; `detected` is never touched)

1. **observation**: ≥ 1 resolved called-cell group at the taxon (or below it).
2. **prioritized**, which requires all of:
   - ≥ `k` resolved called-cell groups (`k` is reported, not calibrated);
   - not in the host lineage;
   - non-viral-only DB;
   - called-cell enrichment over off-cell evidence reported;
   - not `artefact_supported_only`;
   - negative-control comparison done, or explicitly `not_assessed`.
3. **sequence_supported**: a user-supplied reference FASTA for the taxon passed through the evidence adapter (§9.3),
   with recorded results. `not_assessed` ≠ `negative`; both are written explicitly.

### 9.3 Evidence adapter (P1, narrow)

For one prioritized taxid and a user-given `--candidate-fasta`:

- `extract_exact_reads_by_number` is called with a synthetic `lineage_by_read_number` built from `reads.tsv.gz` (the
  ordinals of the taxon's resolved groups' reads), producing a FASTA named `<CB>_<UMI>_<n>|<read_id>`;
- then `write_competitive_fasta(host_fasta, candidate_fasta)` → `align_reads_to_viral` → `coverage_table`,
  `alignment_qc_table`, `per_cell_alignment_qc`, `competitive_blast_identity`.
- The input pair must be the same resolved pair, so ordinals match.
- `molecule_verdict` and `resolve_viral_target` are **not** reused; they are kallisto- and panel-coupled.
- There is no automatic genome download. Validation is limited to `--max-candidates` (default 3).

### 9.4 Comparison with ViralScan calls (P1)

For taxa that map to a panel virus, cross-tabulate per called cell:

- discovery resolved groups versus primary `counts_unique`;
- `detected` versus the discovery state.

This is reported only and changes nothing.

---

## 10. Cell universe and interpretation

- Partitions: called cells (from the primary `called_cells.tsv`), off-cell on-list barcodes, uncorrectable or
  bad-geometry fragments. Each is counted in `qc.json`.
- Called cells with zero assignments are rows of `X`. Denominators are reported with every rate: called cells, called
  cells with ≥ 1 resolved group, and total groups.
- Not an absolute abundance. Detectability depends on:
  - poly(A) or 5′ capture: non-polyadenylated bacterial RNA is captured mostly by internal priming;
  - R2 read length;
  - sequencing depth;
  - DB content;
  - host-subtraction loss.
- A negative Kraken2 result is never evidence of absence.

---

## 11. Run isolation and HPC

### 11.1 Isolation guarantees

- Reads `run_dir` and writes only `-o`. `-o` inside `run_dir` is refused.
- No `RunConfig` field, no `config.yaml` write and no `run` option in M0–M2.
- A gate test (G2) compares a sha256 inventory of every file in `run_dir` before and after `discover`, and runs
  `viralscan validate-run run_dir` after.
- P1 Snakemake leaf (M3):
  - modelled on `host_capture`;
  - flags default `None` and are added to `omit_when_unset`;
  - a new test feeds `build_parser().parse_args([...])` into `build_run_manifest`, which closes the gap that existing
    fingerprint tests use hand-built Namespaces;
  - `_write_run_config` must not rewrite `config.yaml` when discovery keys are absent, so kb and STAR do not rerun
    (SW-22).

### 11.2 Discovery fingerprint (`provenance.json` → `discovery_fingerprint`)

The fingerprint is the sha256 of the canonical JSON of:

- `primary_run_fingerprint`, read from `run_manifest.json`; it is linked and never compared to or used as this
  fingerprint;
- `input_source` and the sha256 of the resolved R1/R2;
- `host_subtraction` (STAR version, param set and every argument, host index identity when run by discovery);
- `kraken2_version` (`kraken2 --version`);
- `db_identity`: sha256 of `opts.k2d` and `taxo.k2d`, plus `hash.k2d` size, mtime and **full sha256 by default**
  (the `--db-hash quick` shortcut is a human decision, §16);
- the taxonomy tree sha256;
- all parameters (confidence, min hit groups, group policy, rank floor);
- the on-list sha256 and `bustools --version`;
- the called-cells sha256 and the cell-calling receipt;
- the ViralScan version.

`discover` with an existing `-o`:

- same fingerprint and `discovery_complete.json` present → no-op;
- different fingerprint → refuse unless `--force`;
- `.partial/` present → clean `.partial/` only, then restart.

### 11.3 HPC behaviour

- Explicit opt-in; independent restartability (above); atomic publication (`.partial/` + `os.replace`).
- Streaming throughout:
  - one pass to write `reads.tsv.gz` and the FASTA;
  - Kraken2 streams;
  - the join uses a numpy bool/uint32 array of size `n_parsed` (8 bytes per fragment, i.e. 0.8 GB for 1e8 fragments;
    a `ponytail:` note sets 4e8 as the limit for sharding);
  - grouping uses `LC_ALL=C sort` on `(cb, umi)` with `-S` bounded, or CB-hash sharding. **Shard-count invariance is
    the chunking gate (G8).**
- Resources: `--threads` only for Kraken2. RAM is about the DB size without `--memory-mapping` (for example, a
  Standard-sized DB is tens of GB), or small with it on fast local storage.
- Copying the DB to `$TMPDIR` is a documented user choice (`--db` points at the copy), not something `discover` does.
- The shared DB is never written, deleted or cleaned.
- SLURM example: one `sbatch` per sample with `--cpus-per-task 16 --mem <DB size + 8G>`. Exit codes propagate.

---

## 12. Benchmark plan

### 12.1 Correctness (mandatory; synthetic; CI)

Each case is a fixture with an independent oracle (pure Python, re-derived from the raw records):

| Case | Expected |
|---|---|
| one read per molecule; multiple reads per UMI; multiple UMIs per cell | groups and units per §6.5 |
| conflicting reads (562 vs 1280) | LCA, `conflict=True`, never two units |
| parent/child (561, 562) | resolved 561, `deepest_consistent` 562 |
| cross-species UMI collision | one group, LCA above both |
| unclassified + classified in one group | resolved by the classified reads; unclassified counted |
| duplicate read ID; missing ordinal; extra ordinal; non-numeric ID | abort (I1–I3) |
| missing/short barcode; N in CB or UMI | `bad_geometry` |
| correctable / uncorrectable CB; Drop-seq (no on-list) | `corrected` / `uncorrectable` / `no_onlist` |
| unsupported chemistry (10xv1, offset geometry) | abort |
| host-lineage taxid (9606 subtree) | `is_host_lineage`, not a candidate |
| empty sample / all unclassified | valid empty outputs, invariants hold |
| truncated gzip / truncated FASTQ | abort, no published outputs |
| record order permuted; Kraken2 output permuted; shard count 1/2/7; threads 1/4 | byte-identical `umi_groups`/`X` |
| `-o` inside `run_dir`; incomplete run; missing DB file | refuse |

The Kraken2 step in the unit tests is a stub that writes output from a fixture table. One `@pytest.mark.integration`
test runs real `kraken2` 2.17.2 on a tiny custom DB built from three synthetic genomes plus one host decoy.

### 12.2 Scientific performance (frozen protocol; no claims before results)

These link to existing open rows rather than creating a parallel benchmark:

- **Host-subtraction loss (gate 12)** = `HOST-05` / `SW-07`. Spike known viral reads, including EVE-like and
  host-homologous segments, into a host background, run STAR `pinned` and `star-default`, and report the fraction of
  truth fragments removed. `lost_truth_counts` (`scripts/host_filter.py:191`) already computes D15/D16. Discovery reuses
  it.
- **DB comparison:**
  - viral-only (RefSeq viral) versus a broad DB with human (PlusPF-class), on the same inputs;
  - metrics: taxon sensitivity, false-positive taxa per sample, rank correctness against the truth lineage,
    host-lineage leakage;
  - RAM and time recorded per DB.
- **Threshold grid:** confidence {0.0, 0.05, 0.1, 0.2} × min hit groups {2, 3, 4}, which is 12 configurations. Only
  Kraken2 is rerun; the host subtraction is computed once. About 12 × 1 h CPU per sample for a mid-size library with a
  memory-mapped DB (UNRESOLVED; measured in DISC-12).
- **Positive controls:** existing ViralScan public positives (EBV, HHV-6B, HSV-1 Drop-seq, covid PBMC) give
  independent truth for panel viruses. Discovery must recover the panel virus as an observation at genus or species
  level in the cells where ViralScan has unique molecules.
- **Negative controls:** uninfected matched libraries (e.g. hsv1 37 °C uninfected; L-94) give the background taxa
  distribution.
- **Resources:** wall time, peak RSS and I/O for 1e6, 1e7 and 1e8 fragments.

Correctness gates are pass/fail. Scientific thresholds (`k` in §9.2, grid choice) are set only after a frozen
protocol run and recorded with their receipts.

---

## 13. Priority matrix, roadmap and tickets

### 13.1 Priority matrix

| Item | Priority | Depends | Effort |
|---|---|---|---|
| ADR, contracts, this document | P0 | — | done (DISC-00) |
| Input resolution + identity join + Kraken2 run + reads/kraken outputs | P0 | — | 2–3 d |
| DB-derived taxonomy + LCA + group classes | P0 | DB parse spike | 2 d |
| Barcode correction via bustools + reconciliation | P0 | spike | 1–2 d |
| Called/off-cell partition + `cell_taxon.h5ad` + qc/provenance | P0 | above | 1–2 d |
| Isolation gate tests (G1–G3, G8) | P0 | above | 1 d |
| `discover-collapse` | P1 | taxonomy | 0.5 d |
| Candidate states + negative-control comparison | P1 | P0 | 1–2 d |
| Evidence adapter | P1 | P0 | 1–2 d |
| Snakemake leaf + `omit_when_unset` + parser-fed fingerprint test | P1 | M2 gates | 1–2 d |
| Benchmarks (HOST-05 linkage, DB, grid, resources) | P1 | M2 + compute | 1–2 wk wall |
| Docs (`docs/discovery.md`, CLI reference) | P1 | M2 | 0.5 d |
| Bracken, multiple DBs, BAM/Cell Ranger input, assembly, pathways, interactive viz, tier promotion | **P2 / defer** | — | — |
| Seurat/h5seurat, name-based collapse, majority vote, fastp/bowtie2/predb pipelines, implicit DB download, dense matrices | **REMOVE** (never build) | — | — |

### 13.2 Milestones (reduced from six to five; M4 and M5 overlap)

**M0 — contracts (this document)**

- Done when §16 is approved and spikes DISC-04 and DISC-05 report.

**M1 — classification prototype**

- Files: `discovery.py` (input resolution, stream, ordinal FASTA, Kraken2 wrapper, join), `cli/parsers.py` (`discover`),
  `menu.run_discover`.
- Outputs: `reads.tsv.gz` with `cb_status=raw` pending M2, `kraken.out.gz`, `kraken_report.txt`, `qc.json`,
  `provenance.json`.
- Tests: identity cases, isolation (G2), refusal cases, stubbed Kraken2.
- Pass when: I1–I5 hold on all fixtures, and the primary-dir inventory is unchanged.

**M2 — single-cell aggregation**

- Files: `discovery.py` (correction, taxonomy, LCA, grouping, partition, h5ad), `discover-collapse`.
- Tests: all of §12.1 against the oracle; order, shard and thread invariance; real-Kraken2 integration on a tiny DB.
- Scientific check: the reconciliation gate of at least 99.9 % CB agreement on one real library (EBV SRR12682296).
- Pass when: every §6.5 invariant holds, the oracle matches, and output is byte-identical across permutations and
  shards.

**M3 — integration (optional)**

- A Snakemake leaf rule only if requested.
- Gates G1–G3 plus the parser-fed fingerprint test; `test_snakefile_dag` "absent by default" case.

**M4 — evidence and reporting**

- Candidate states, negative-control comparison, evidence adapter, `docs/discovery.md`.

**M5 — benchmark and release**

- §12.2 under a frozen protocol, with results in `results/discovery_benchmark/` and receipts.

### 13.3 Agent tickets

Each ticket is one branch `codex/disc-NN`, touches ≤ 3 files, and updates PLAN.md in the same commit.

| ID | Task | Files | Tests | Accept |
|---|---|---|---|---|
| DISC-01 | Optional discovery env: kraken2 `2.17.2` pin (bioconda availability check), `doctor --profile discovery` | `pixi.toml` feature or `environment.discovery.yml`; `validation.py` | test_validation profile | `doctor --profile full` unchanged |
| DISC-02 | `discover` parser + `run_discover` skeleton: run-dir completeness, `-o` refusal, tool preflight, `.partial/` publish | `cli/parsers.py`, `menu.py`, `discovery.py` | refuse incomplete / `-o` inside / missing tool | G2 inventory test passes |
| DISC-03 | Input resolution (primary `host_filtered` or explicit STAR pass via `_starsolo_filter`) + audit cross-checks | `discovery.py` | fixture run dirs: reuse, explicit pass, absent → refuse | I4; `input_source` recorded |
| DISC-04 | **Spike:** DB taxonomy (`k2 inspect` parse vs `taxo.k2d`) on a tiny built DB | scratch → `discovery.py` | parse equals the built DB's tree | decision recorded in PLAN |
| DISC-05 | **Spike:** `bustools fromtext→sort→correct→text -f` keeps flags; else Python Hamming-1 + oracle | scratch → `discovery.py` | correction property test | reconciliation ≥ 99.9 % on EBV |
| DISC-06 | Streaming pass: ordinal FASTA (R2 only) + `reads.tsv.gz` raw CB/UMI via `cb_umi_geometry`; reject 10xv1 | `discovery.py` | geometry cases | `n_input` invariant |
| DISC-07 | Kraken2 wrapper (`evidence._run`) + exactly-once join | `discovery.py` | stub Kraken2: dup/missing/extra/non-numeric abort | I2, I3 |
| DISC-08 | LCA + group classes + partition + `cell_taxon.h5ad` + `offcell_taxon.tsv.gz` | `discovery.py` | oracle property test (seeded, 40 fixtures) | §6.5 holds; G8 invariance |
| DISC-09 | `qc.json`, `provenance.json`, discovery fingerprint, resume/no-op/`--force` | `discovery.py` | fingerprint changes on each input; no-op on rerun | G9 |
| DISC-10 | Real-Kraken2 integration test on a tiny synthetic DB | `tests/test_discovery_integration.py` | `@integration` | output equals oracle |
| DISC-11 | `discover-collapse` | `discovery.py`, `cli/parsers.py` | rank sums + `above_rank` = total | gate 6 |
| DISC-12 | Benchmark protocol (frozen) linked to HOST-05/SW-07 | `docs/plans/…-discovery-benchmark.md` | — | approved before runs |
| DISC-13 | Candidate states + evidence adapter | `discovery.py`, `menu.py` | `not_assessed` ≠ `negative`; viral-only DB capped at observation | gate 10 |
| DISC-14 | Optional Snakemake leaf (only on request) | `Snakefile`, `runconfig.py`, `run_safety.py` | parser-fed fingerprint test; DAG absent-by-default | G1, G3 |

---

## 14. Acceptance gates

1. **G1**: discovery off → the `run` DAG, `config.yaml` bytes and the fingerprint equal `origin/main` for the existing
   fixtures. Test: parser-fed manifest equality.
2. **G2**: discovery on → the sha256 inventory of `run_dir` is unchanged, and `validate-run` passes.
3. **G3**: discovery-only STAR never changes `kb_r1/kb_r2` (its outputs live under `-o`).
4. **G4**: every row of `reads.tsv.gz` traces to an input record. Exclusions are counted, and missing or duplicate IDs
   abort.
5. **G5**: no barcode enters `X` unless it is on-list, bustools-corrected or `no_onlist` (declared). Reconciliation is
   reported.
6. **G6**: one unit per resolved group per resolution; rank projections conserve the total.
7. **G7**: no group is counted at a taxon more specific than the LCA of its classified reads. Root-level (ambiguous)
   groups never appear in `X`.
8. **G8**: output is byte-identical across record order, Kraken2 output order, shard count and threads.
9. **G9**: a missing or invalid DB, unsupported chemistry, an incomplete run, or corrupt or truncated input → a
   non-zero exit with a named cause and no published outputs.
10. **G10**: no code path writes `detected`, tiers or primary H5AD layers. A grep test asserts that `discovery.py`
    never imports `detection`, `evidence_tiers` or `multimap` writers.
11. **G11**: an existing `-o` with a different fingerprint is refused without `--force`; the primary run is never
    modified.
12. **G12**: host-subtraction loss is measured (HOST-05) before any sensitivity statement.
13. **G13**: performance claims cite frozen benchmark receipts.
14. **G14**: runs on Linux/SLURM with pinned kraken2 2.17.2, STAR 2.7.11b and bustools (recorded); no network on
    compute nodes.
15. **G15** (added): a viral-only DB never yields `prioritized` or above.
16. **G16** (added): `provenance.json` alone identifies the DB, inputs, tools and parameters well enough to rerun.

---

## 15. Risk register

| Risk | Sev | Mitigation | Status |
|---|---|---|---|
| Called cells in STARsolo uncorrected space (no `-w`) vs discovery bustools space | High | receipt read; overlap reported; documented limitation | UNRESOLVED until DISC-05 on a real run |
| bustools `fromtext` drops flags or rejects N | Med | Python Hamming-1 fallback + oracle | spike DISC-05 |
| `hash.k2d` sha256 cost (tens of GB) | Med | full by default; `quick` mode only by approval | decision §16 |
| Residual host via >20-locus multimappers and partials | High | broad DB with human; host-lineage flag; viral-only capped | designed |
| EVE / host-homologous viral loss in host subtraction | High | measured (G12); never claimed absent | open (HOST-05) |
| Kraken2 confidence denominator doc/code mismatch | Low | synthetic test pins the behaviour | M1 |
| `host_filtered/` not hashed by primary | Med | discovery-time hash; `primary_verified:false` stated | accepted limitation; optional P1 primary change |
| Track A FASTQ deletion removes inputs | Med | fail closed; explicit message | designed |
| DB on NFS without memory mapping → OOM / slow | Med | `--memory-mapping` documented; user stages to `$TMPDIR` | docs |
| Bioconda lacks kraken2 2.17.2 | Low | pin the newest available, record it; `kraken2-inspect` fallback if `k2` is absent | DISC-01 |
| Sample namespace collisions on later aggregation | Low | single-sample scope; `sample_id` in provenance | P2 |

---

## 16. Decisions requiring human approval before implementation

1. **Architecture:** a standalone `viralscan discover` first, with a Snakemake leaf only on request (ADR-DISC-1)?
2. **Input stage:** classify the host-filtered pair before the read filter, with `read_filter` decisions as an
   annotation (§5.2)?
3. **Barcode space:** correct into the kb/bustools space with the primary's resolved on-list, and only *report* a
   mismatch with STARsolo-space called cells (§6.3)?
4. **UMI policy:** exact `(CB, UMI)` groups; taxon = LCA of classified reads; no majority vote; `deepest_consistent`
   is audit only (§6.4)?
5. **Kraken2 defaults:** start at the Kraken2 defaults (confidence 0.0, min hit groups 2) and choose by the frozen grid,
   not Cellsnake's 0.05/4?
6. **DB identity:** full sha256 of `hash.k2d` by default (slow on large DBs), or size+mtime plus the small-file
   hashes?
7. **Viral-only DB policy:** cap at `observation` (G15)?
8. **Environment:** Kraken2 in an optional discovery environment only, never in `environment.yml` or the Dockerfile?
9. **Primary change (optional, P1):** start hashing `host_filtered/` outputs in the primary manifest via
   `record_manifest_block` (unhashed block; no fingerprint change)?
