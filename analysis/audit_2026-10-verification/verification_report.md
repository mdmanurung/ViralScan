# ViralScan audit — independent verification report

Verified against `$SRC` = `87c9e35ac89c73d668e8c5897b9a9cc920616766` (PR #67 merge, 2026-10-15; this is the
current `origin/main` tip after `git fetch`). Worktree: `<pinned-87c9e35>`.
Matrix: `verification_matrix.tsv` (51 rows). Repros: `repro/` (all assert `viralscan.__file__` is under `$SRC`).

How this was done: six domain agents re-derived each finding from source (own repros, plus audit repros patched to `$SRC`).
I then re-ran the repros for P4-02, P3-02, P3-03, P9-02 and P5-01 myself and traced P5-01's claim surface and the
P2-01 `host_homology` flag in code. Everything not listed in "verified by me" below is the agent's evidence, with its
repro on disk, but I did not re-run it.

## (a) Code-version situation

- Audit target `97a64f6`; `$SRC` is 4 commits later (`e793dae`, `4367adf`, `37cf4ca`, merge `87c9e35`).
- `git diff --stat 97a64f6 87c9e35` (whole tree, 9 files): `src/viralscan/scripts/cellcalling.py`, `scripts/emptydrops.R`,
  `tests/test_cellcalling.py`, `docs/faq.md`, `docs/output_reference.md`, `analysis/v3_artifact_inventory.tsv`,
  `PLAN.md`, `CHANGELOG.md`, `.living/decisions.md`.
- Effect on findings: `knee_cells` and `_in_barcode_set` (P4-03, P4-05) were **not** touched by the rewrite, so those are not STALE.
  The new CELLS-05 stop ("a kb matrix column the Virus Identity table does not list") covers only the *opposite* direction of
  P9-02, and only on the default emptyDrops path. The docs edits did not change the wording P3-03 relies on. **No finding is STALE.**
- Local/remote branch scan: `codex/umi-dedup-q137` (P3-01/P8-02) and `codex/multimap-3-options` (P9-01) are already merged to main;
  no unmerged `codex/*` branch touches the cited files (note-level only).
- Test baseline at `$SRC`: `4121 passed, 15 skipped, 138 deselected` (30m37s, `viralscan_test_full`). The skips are gitignored GTFs.
- **Branch fact:** `CLAUDE.md` and `PLAN.md` at `$SRC` still say all v3 work happens on `codex/viralscan-v3`, but that branch exists
  neither locally nor on `origin`. See `implementation_plan.md` (decision D1).
- **Audit-document inconsistency:** the report body says "47 rows", "High (5)" (lists 6), "Low (17)/Info (8)"; the TSV has 51 rows
  (1 Critical, 6 High, 20 Medium, 20 Low, 4 Info). The TSV is canonical here.
- Hygiene: one agent ran some repros with `python -I` before switching to `-B`, so two ignored `__pycache__` directories exist under
  `$SRC/src/viralscan/` and `.../scripts/`. They are gitignored; `git status` of `$SRC` is clean. The main checkout was not touched.

## (b) Verdict counts (51)

| Verdict | n | IDs |
|---|---|---|
| CONFIRMED | 11 | P3-02, P3-03, P3-07, P3-09, P4-02, P4-04, P4-06, P5-05, P8-01, P9-03, P9-05 |
| CONFIRMED-DOWNGRADED | 7 | P5-01 (Critical→High), P5-03 (High→Medium), P2-02, P3-06, P4-03, P1-02, P3-05 |
| CONFIRMED-UPGRADED | 0 | |
| PARTIAL | 21 | P2-01, P5-02, P7-01, P3-01, P9-01, P9-02, P2-03, P2-04, P2-05, P2-08, P4-01, P4-05, P5-04, P5-06, P5-07, P5-08, P7-02, P7-05, P1-01, P3-04, P8-02 |
| BY-DESIGN | 9 | P6-01, P6-02, P6-03, P6-04, P6-05, P7-03, P2-06, P2-07, P2-10 |
| REFUTED | 3 | P2-09, P3-08, P9-04 |
| STALE / UNVERIFIABLE | 0 / 0 | |

Severity movement: of the audit's 7 Critical/High, none stays at its label. P5-01 → High; P2-01, P5-02, P5-03 → Medium;
P3-01, P7-01 → Low; P6-01 → Info (documented state, not a defect). "Confirmed" in the audit's own confidence column held for only 11 of 51.

## (c) Per-finding notes (everything not plain CONFIRMED)

**P5-01 CONFIRMED-DOWNGRADED, Critical → High.** `scripts/umap.py:umap` does
`adata.var["highly_variable"] |= adata.var_names.isin(list(found_genes))` before PCA/neighbors/UMAP, then
`viral_neighbor_enrichment(coords, labels)` tests `labels = virus_positive` (derived from the same viral counts) in that embedding.
No flag excludes viral genes and the plot is not labelled viral-inclusive. *Verified by me:* on pure-noise host cells with random
virus-positive labels, `p<.05` in 15/15 runs (dominant viral counts) with viral genes in PCA, vs 27% (n=15) with them removed (agent: 10–20%, 95–100% / moderate).
So the *test*, not just the picture, is invalid. *Scope check (me):* nothing else consumes the embedding — `obsm`/`X_pca`/`X_umap` appear only in `scripts/umap.py`, which writes HTML plots and a `umap.done` sentinel; `enrichment.py` and `gene_programs.py` do not read it. Why not Critical: `--umap` is opt-in; the p-value is only a plot annotation (`conclusion_text`);
nothing in `claims/registry.json`, `viral_summary.tsv`, README or the vignettes consumes it; `docs/faq.md` describes it as clustering of "candidate-support labels" without claiming validation.
Real overclaim: `templates/report.html.j2:471` says significant clustering "suggests true infection rather than ambient contamination". The repo already
flags this sentence for review (`docs/plans/ambient-and-barcode-swapping-controls.md:65`), so it is a known, unfixed item. Unexplained: the viral-free arm still exceeds the nominal 5%
(7/60 for the agent, 4/15 for me) — a second, smaller null-calibration problem, not yet diagnosed.

**P2-01 PARTIAL, High → Medium.** Real: `validate_reference_records` rejects only exact-SHA-256 duplicates; 1-mismatch, 10%-diverged and reverse-complement host
homologues pass. `--genome-dlist` defaults to `None` on `build-ref`, `build_bundled_panel_ref.py` and every `--profile`; `reference_gates` host-crosstalk is a standalone CLI no build path calls.
Overstated: if the homologue is a host transcript in the same index, kallisto maps host reads to both (279/300, 0 viral-unique) and default `host-conservative` keeps that mass out of primary counts.
The audit's failure mode holds only when the homologue is *absent* from the host cDNA (295/300 viral-unique). Offline measurement in the repo: 51 of 2,343 genomes homologous (REF-07/F-023).
*Verified by me:* `evidence.py` has a `host_homology` flag (`flagged` at ≥25% host-preferred/competitive sampled reads), diagnostic-only — a read-level net exists but gates nothing.
Real exposure: user-supplied panels, and anything built without `--genome-dlist` (`host_homology_status=not_measured`, which docs say is not zero homology).

**P3-01 PARTIAL, High → Low.** The representative read is the best-viral read; BLAST/UMI FASTA, Q137 gates 1/2/4 and the ambiguity flag are built from it. This is a documented, user-approved design
(`.living/decisions.md` ~L1939; `docs/output_reference.md:943-944`). Gate 3 and `molecule_verdicts` use all reads and report `host_best` on a mixed molecule. alignment_qc's dedup layer uses a different rule, so the audit is wrong there.
No product gate consumes the representative read; only the DSR02 validation script's low-complexity verdict reads the BLAST-derived flag. A residual, honest gap: the flags are biased toward "virus", which is the conservative direction for a false-positive hunt but not for a read-level QC claim.

**P5-02 PARTIAL, High → Medium.** Reproduced (cell-CV AUC 0.787 vs group-CV 0.335; the audit quoted 0.413, seed-dependent). Default is `--cv cell`, no auto-switch. Group CV exists, fails closed, `cv_mode` is in every output row, the FAQ documents it, cell-mode results are labelled descriptive. A deliberately kept legacy default → needs user decision to change.

**P5-03 CONFIRMED-DOWNGRADED, High → Medium.** The audit's 0.774 reproduces (0.782) only at 400 samples × 200 features. At 1500 features the noise selection probability is 0.03/0.08/0.24/0.49 for 30/100/400/1000 positives (0%/0%/1.5%/25% of noise genes "stable"). Defect = fixed `C=1.0` not scaled with n. Docs call it descriptive and claim no error control.

**P9-01 PARTIAL, Medium → Low.** The stop rule does stop early (iteration 2), but the audit's 3.2 is θ; allocated mass is 2.56 and tracks pseudocount/tol (tol 1e-9 → 0.008). `em-global` is not default (default `host-conservative`; CLI alias `em` = em-cell), and end-to-end the screen passes but `detected=False` (max per-cell mass 5.9e-4) → "CANDIDATE ONLY", not a detection.

**P9-02 PARTIAL, Medium → Low.** *Verified by me:* a stale identity table listing `vB_g1` against a matrix holding `vB_g1.2` in 30/30 cells yields no Virus B row, summary line or sensitivity row, and no warning. Table genes merely absent from the matrix carry no counts, so nothing real is lost there. The CELLS-05 stop fires on the opposite direction only, and only for the default one-step emptyDrops path (not `--cell-calling none`, external, or two-step). `run_safety` fingerprinting blocks the normal `--resume` path. Cheap to fix and silent, so it keeps top priority despite the lower label.

**P7-01 PARTIAL, High → Low / P7-02 PARTIAL, Medium → Low.** Mechanism real: `kb_count` declares only `log/kb.done`. The audit's regex repro never evaluates `_kb_count_inputs` and so cannot show staleness. Agent's real snakemake 9.23.1 dry-run: in-place index edit and a purged `kb-python/` trigger no rerun; a changed index *path* rewrites `config.yaml` and reruns. Through `viralscan`, a non-empty output dir needs `--resume` (sha256 fingerprint of index, transcripts, FASTQs, mismatch refused) or `--overwrite` (wipes sentinels). Residual: direct `snakemake`, an unsupported entry point. **Could not verify:** whether `--resume` after a purged `kb-python/` fails loudly downstream (the agent did not run it, and I only grepped).

**P3-06 CONFIRMED-DOWNGRADED, Medium → Low.** IndexError on `SEQ=*` reproduced, but needs a viral-aligned primary record with no sequence; minimap2 primaries always carry SEQ, so only externally supplied BAMs trigger it.

**P4-01 PARTIAL, Medium → Low.** Silent only for two-step plus explicit `--strand forward` on a 5' library; default is unstranded and help text warns. Virus-only index without host filter is not silent (see false negatives).

**P4-03 CONFIRMED-DOWNGRADED, Medium → Low.** Not STALE (`knee_cells` untouched). At HEAD it calls 0/53000 on the audit curve (audit saw 1). Known, user-decided defect (PLAN SW-23, decisions 2026-10-01); knee is opt-in, default `auto`, warns every run, docs say not for reported numbers. The "opposite failure mode" docstring claim holds only for the audit's synthetic curve.

**P4-05 PARTIAL, Info.** No path yields suffixed or colliding obs_names.

**P2-02 CONFIRMED-DOWNGRADED, Medium → Low.** Wrap genes lose the first part (0/98 vs 98/98 reads pseudoaligned), but this is documented and test-pinned (`panel_integrity.md` L104, D-13). GenBank-fetched wrap CDS keep both exons; only junction order is wrong. Local cache: HBV (6 genomes, ~12 genes), EBV LMP-2, one TTV-like record. The shipped HBV GTF is absent from the pinned tree, so **not checked against shipped annotation**.

**P2-03/P2-04/P2-05 PARTIAL.** Placeholder gene is a logged fallback; shipped corpus has none (exposure: FASTA-only user references, 47/2,041 anellovirus accessions in an opt-in panel). The audit cited the legacy constant; the catalogue `sibling_group` has 4 groups (HSV, HHV-6, HHV-4/EBV types, KFDV/AHFV), no HPV row, and the guard is a note that never changes counts; HPV cross-mapping is predicted, not measured (PLAN HPV-11 open). `_parse_location` does drop nested `order()` intervals with a warning, but a scan of 11,557 features in 4,474 cached flatfiles found 0 dropped intervals — latent.

**P3-04, P3-05, P1-01, P1-02, P8-02, P5-04..P5-08, P7-05, P2-08.** Real but narrow or documented; see matrix notes. Examples: P3-05 needs `assign_tier`, which has no production caller; P1-02 needs a barcode file kb/bustools never emit; P5-04's confound is mitigated by `depth_alone_auc_mean` plus a `log.warning`; P7-05's location is wrong (`kb_outputs.py` does no I/O), "untested at scale" is true and recorded as deferred.

**BY-DESIGN.** P6-01 (protocol `draft`, `generated/` absent, claims registry already withholds accuracy claims; nothing under holdout/`generated/` opened), P6-04 (benchmark doc self-flags "NOT yet a fair comparison", `PUBLICATION_READINESS.md:122` excludes it; the 18×/3.7× figures still sit unretracted in the benchmark doc lines 23-25 and `.living/learnings.md`, not in README/docs/registry), P6-05 (refuses real backgrounds by design, blocker `real_background_anchor`), P7-03, P6-02, P6-03, P2-06, P2-07, P2-10.

**REFUTED.** P2-09 (the 2 taxids carry spelling variants of one species), P3-08 (`select_detection_matrix` returns `adata.X` unconditionally; `RunConfig` rejects any other `multimap_primary_call`), P9-04 (real bustools 0.45.1: output is exactly (barcode, umi) string-sorted; barcodes are 2-bit encoded, so numeric order equals lexicographic order).

**Confirmed, with changes in scope worth knowing.**
- **P3-02 (Medium):** worse than audited. Pinned minimap2 2.31 `-ax sr` emits no secondary alignments, so exact host/virus ties are *never* reported as ties. *Verified by me:* 200 perfect ties → 105 `host_best`, 95 `virus_best`, 0 `tie`; with `--secondary=yes` → 200 ties. Docs (`output_reference.md:984`, "score ratio 0.8, limit 5") and the verdict docstring describe the opposite. Only minimap2 2.31 tested. Consumers are diagnostic flags, which caps severity at Medium.
- **P3-03 (Medium):** *Verified by me:* `summary.txt` says "1/10 called cells (10.00%), 9/100 all barcodes (9.00%)", while `report.html` headlines "Infected cells (any virus) 9" and a 9.00% table row, never showing a called-cell rate. `output_reference.md:303-305` is unchanged and names `pct_infected_called` as the appropriate rate.
- **P4-02 (Medium):** *Verified by me:* same data gives `pct_infected` 0.04 / `viral_molecules_per_10k_est` 0.0156 combined, vs 100.0 / 10000.0 two-step, while `pct_infected_called` and `pct_infected_comparable` agree at 0.2. The bad values reach `viral_summary.tsv`, `summary.txt` and the HTML main table (badge above 1%, no called column). 10000.0 needs a single-virus toy matrix; the 10⁴× is a toy extreme, the 100%-style inflation is not.
- **P5-05 (Medium):** worse than audited. With exactly 10 positives and the default top-50%-depth filter each test fold is 1 positive × 1 negative (seed-sd 0.47); no low-support warning.
- **P9-03 (Low):** the audit's "never fires in float64" is wrong. Sums of k copies of 1/k fall below 1 for k=6 and 7, and `detect_genes` has no tolerance, so a virus ULPs under the threshold is silently dropped with no row and no CANDIDATE ONLY. Narrow trigger.

## (d) Audit false negatives (5 kept, evidence-backed)

1. **`--groups donor` is silently ignored under the default `--cv cell`.** `hostresponse.py`: `groups_column` is read only when `cv_mode == "group"`; no warning. Agent-verified, not re-run by me.
2. **minimap2 `-ax sr` never emits secondary alignments**, so "tie" is unreachable by default and docs say otherwise. Folded into P3-02. Reproduced by me (`repro/P3-02_minimap2_no_secondary.py`).
3. **Virus-only index without a host filter triggers an `error` from the strand/chemistry sanity gate on normal runs** (`analysis.py:_chemistry_sanity`, `host_in_index = not host_filter_aligner`; reproduced at p=0.3 by the agent). A false alarm, the opposite of P4-01's silence. Agent-verified only.
4. **`detect_genes` screen has no float tolerance** and its output alone feeds `compute_stats`, so P9-03's rounding drop also removes the virus from every downstream table. Extension of P9-03.
5. **Exact-hash dedup also passes reverse-complement duplicates** (folded into P2-01; reproduced by the reference agent). Also noted, Info: the DataFrame barcode path sorts before stripping `-1` (not kept in the five).

## (e) Positive claims tested

1. **"57/57 oracle matches; conservation holds" — survived.** An independent oracle (not importing the audit's) found 0 mismatches over 6 methods on ties, gene-less ECs, buffer sizes 1/2/3/7/64, unsorted-input rejection and 60 fuzz seeds; two deliberately broken mutants produced 948 and 344 mismatches, so the oracle can fail. Limits the audit did not state: conservation says nothing about EM *accuracy* (P9-01) nor about float rounding at thresholds (P9-03).
2. **"Phase 8: 17 prior majors, 12 fixed-verified, 0 still present" — not falsified.** Re-ran copies of the A-03 and B-03 repros against `$SRC` (both hold); A-04, C-02, C-05 hold by code/data inspection. A-01 end-to-end and the M-02..M-04 partial fixes were not rechecked.
3. **A "pass" that did not test what it claimed: the audit's P7-01 repro.** It regex-matches the literal text of `_kb_count_inputs` and never evaluates it. Its "passes" does not show staleness either way.
4. **Audit statement "never fires" for gate-leg slack (P9-03) — falsified** (float64 sums of 1/k for k=6,7), see above. **P9-04's premise — falsified** with the real bustools.
5. Not tested: the audit's "fail-closed chemistry autodetection" and "tier fail-closed" strengths; and anything needing a live `kb count`/STAR end-to-end run.

## Repro run table (every `repro/*.py`, run once from scratch on 2026-10-15, `viralscan_test_full`, each asserts `viralscan.__file__` under `$SRC`)

- 38 scripts, **38 exit 0** (36 in the batch run; `P6-05_production_refusal.py` and `P7-01_snakemake_staleness.py` need a scratch-dir argument, the batch run omitted it and they exited 1 with `IndexError`, then both exited 0 when re-run with one).
- Slowest: `P9-01_em_global_residue.py` 807 s, `P5-03_stability_null.py` 365 s, `P5-01_viral_in_embedding.py` 351 s, `P5-07_host_only_normalization.py` 155 s, `P7-01_snakemake_staleness.py` 127 s, `P5-02_audit_synthetic_confounding_demo.py` 120 s; all others under 90 s.
- Exit 0 means "ran to completion against `$SRC`", not "finding confirmed"; verdicts are in the matrix.
- The audit's own scripts that were previously unrun are now patched and run: `P1-multi_audit_run_oracle_vs_impl.py`, `P4-multi_audit_p4_demos.py`, `P9-multi_audit_p9_adversarial.py` (exit 0). I did not diff their outputs against the audit's saved `*_results.json`/`p9_output.txt`, so I make no claim here about whether they reproduce the audit's numbers.
- Raw per-script output: `<scratch>` (scratchpad, not preserved in the repo).

## Addendum (2026-10-15, after reading the frozen protocol sections)

Earlier drafts of `implementation_plan.md` over-stated the coupling to the frozen protocol. Facts: the frozen `harmonization` block defines the
cell universe (shared host-only called-cell anchor, `primary_cell_metric_population`, line 144), the feature universe and layers; it never
names `pct_infected*`, `viral_molecules_per_10k_est`, `detection_threshold` or `detect_genes`, and defines no detection tolerance.
`detected` is user decision AUD-06/07 (PLAN.md:5557, 5560). `utils.INFECTED_RELATIVE_TOLERANCE = 1e-6` already applies to `infected_masks`,
`tier_measures`, `evidence_tiers`; only `detect_genes` compares with a bare `>=` (P9-03). The all-barcode report headline (P3-03) contradicts the
frozen called-cell primary denominator. `em-global` is in the protocol's `multimap_allocation` grid (protocol.yaml:1143), so a G5 pick could inherit the P9-01
residue. `--cv cell` is not a recorded user decision. The 18×/3.7× figures appear only at `analysis/reference_strategy_benchmark/REFERENCE_STRATEGY_BENCHMARK.md:24-25`
(the earlier "also in `.living/learnings.md`" is unverified). New inconsistency: that doc says "COMPLETE — 12 of 12" (line 9) while protocol.yaml:2179-2180 records PFR-01 as
"4 of 12 rows completed, all EBV rows failed". The P5-01 glossary sentence ships in every HTML report, not only when `--umap` is on.

## Limits of this verification

No live FASTQ→report run. Local GTFs for shipped panel absent at `$SRC`. minimap2 behaviour tested at 2.31 only. P5-01's viral-free null excess is unexplained. Items marked "agent-verified only" were not re-run by me.

## Outcome (implemented 2026-10-15, branch `codex/aud3-integration`)

All decisions were settled with the maintainer before implementation. Frozen definitions (`analysis/v3_validation/protocol.yaml`
`harmonization` block) were not changed; `scripts/validate_v3_protocol.py` still passes.

| Findings | Unit | Commit | What changed |
|---|---|---|---|
| P9-02 | AUD3-01a | 145fc54 | A matrix column the Virus Identity table does not list now stops the run on every cell-calling path; table genes absent from the matrix log one warning |
| P9-03 (+ called-cell gate) | AUD3-01b | 8c235f7 | `detect_genes` and the called-cell gate use the existing relative tolerance 1e-6 |
| P3-03, P4-02 | AUD3-03 | e61b40c | report headline and table lead with the called-cell rate; two-step `pct_infected` and `viral_molecules_per_10k_est` are NA with a `denominator_note` |
| P5-01, P5-07, P8-01 | AUD3-02 | f402621 | neighbour-enrichment test runs on a viral-free host PCA; figure kept and labelled viral-inclusive; glossary overclaim removed; zero-HVG guard |
| P3-02 | AUD3-04 | e1410f8 | `minimap2 -ax sr --secondary=yes -p 0.8 -N 5` so exact host/virus ties are representable; version and options recorded in the evidence manifest |
| P2-01 (partial) | AUD3-05 | 7eba588 | warn-only near-duplicate screen on both strands (reuses the >= 95% shared 31-mer twin cut), visible "host homology not measured" warning, opt-in `--strict-reference`; a homologue diverged by about 8% is NOT screened (threshold belongs to REF-08) |
| P5-02, P5-03, P5-05, P5-06, `--groups` ignored | AUD3-06a, AUD3-06b | b4e7929, 1a2ca7a | warnings, "uncalibrated selection frequency" labelling, and `--cv auto` (group CV when a group column is named or a single donor-like column exists) |
| P9-01, P3-05, P1-01 | AUD3-07 | c701410 | documentation and pinning tests only (no numerical change) |
| P3-06, P4-04, P4-06, P4-03 (docstring) | AUD3-09a | 50e0bcd | `SEQ=*` guard, suffixed-whitelist diagnosis, explicit-geometry error on the read-level paths |
| P4-01, virus-only sanity false alarm, P3-01 (docs), P6-04 (banner) | AUD3-09b | 0f01d74 | sanity gate no longer errors on virus-only indexes; two-step `--strand forward` warning; do-not-cite banner on the 18x / 3.7x figures; representative-read bias documented |

Deferred or left open (tracked in `PLAN.md`): `AUD3-08` (P7-01/P7-02, `--resume` after a purged `kb-python/`; resume-time guard to follow as its
own unit), `AUD3-02-OPEN` (the viral-free neighbour test still gave p < 0.05 in 4 of 15 pure-noise runs), `AUD3-FLAKY` (one R cell-calling
test failed once under concurrent load), `AUD3-G3-EVIDENCE` (evidence made before AUD3-04 must be re-run before it feeds tier measure S),
`AUD3-REF08-DIVERGED` (P2-01 residual). Not fixed by decision: P6-01, P6-02, P6-03, P6-05, P7-03, P2-06, P2-07, P2-10, P2-09, P3-08, P9-04,
P2-05, P2-08, P1-01 (beyond docs), P1-02, P3-04, P7-05, P8-02, P5-08, P2-02, P2-03.

The `repro/` scripts that produced this evidence hard-code local scratch paths and are kept outside the repository.
