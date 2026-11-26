# Scope: `--quantifier simpleaf` backend (post-3.0)

Source: origin/main at 88034089 (2026-11-26), read only. Comparison numbers come from `round_2026-11-26_cmp_simpleaf.md` (repro, branch `codex/cmp-simpleaf`). I inspected no figure for them: `unsupported: no inspected figure`. Alevin-fry facts marked "run" were run on simpleaf 0.31.0 / alevin-fry 0.19.0 / piscem 0.25.0 by a mapping agent; "unverified" facts were not. Flags and defaults in sections 2 and 6 were rechecked against `--help` of those three binaries, and the resolution modes, barcode-correction defaults, chemistry registry keys, `run_safety.py:156,217` and the schema const line were rechecked again after review (2026-11-26). The `molecules.parquet` column list, the 2384/45477 figure and run timings were not re-run.

## 1. Recommendation

Do not build it yet. Run PR 0 (4-6 days, repro repo, no package change) and decide from its numbers. The default outcome after PR 0 is design C (side-by-side check, 4-6 days more). Promote to design A (simpleaf maps and corrects barcodes, ViralScan keeps its own molecule allocation and cell calling, kb stays the default) only if PR 0 passes the gate below.

What changed from the first draft: it said "build it" (30-40 days). Review showed the draft never stated what the backend buys once EM-CELL-AF-01 lands. That PR addresses the EBV-1/EBV-2 split, the main difference seen in the comparison, inside ViralScan; em-global already gives 0.4 % EBV-2. The speed gain is also not shown: the 5-14 min simpleaf figure is the mapper and quantifier alone, while the 20-63 min ViralScan figure is the whole pipeline. Design A still runs multimap, cell calling, detection and umap after simpleaf. The build decision therefore had to follow the evidence, not precede it.

Gate to build A. All of these must hold, and the thresholds are the user's call:

1. The molecule table from `alevin-fry quant --dump-molecules` agrees with ViralScan's BUS-intersection molecules (PR 0, with `--resolution` fixed, section 4).
2. Simpleaf co-solves with the current environment pins, or a separate environment is accepted knowing evidence, discover and translated cannot run in it.
3. A named gain that design C cannot deliver. Candidates: wall time saved per library of at least about 10 min (kb count alone against simpleaf map+quant at equal threads, plus ViralScan downstream time), or intron-aware host counts that change a downstream result (called-cell set, host-response outputs) in the S versus S+U+A measurement. If neither shows, stop at C.

Evidence replay, `discover` and `translated` stay on kb/kallisto in every outcome. If A is built, mark it experimental until the yardsticks in section 5 pass. Also agree the intermediate molecule format with the per-cell EM PR before PR 1, because both sit at the same seam (`multimapping.py:118-197`).

### What A buys after EM-CELL-AF-01

- EBV-1/EBV-2 split: nothing once EM-CELL-AF-01 lands, by the in-flight PR's own scope. PR 0 should check this against the planted truth rather than assume it.
- Viral specificity from USA expansion: nothing. Viral transcripts are all S, for an index built with `--spliced` (section 2A, limit 3).
- Intron-aware host counts: possible, but unmeasured. The host layer choice (S, S+U+A, S only) changes host UMI totals and the called-cell set, so it is an impactful choice, not a free gain. PR 0 measures called-cell overlap and host UMI per barcode, S versus S+U+A, on the three libraries.
- Speed: unmeasured for the pipeline. The comparison timed simpleaf alone; no kb-count-only or ViralScan downstream time exists. PR 0 measures both. If the saving is under about 10 min per library, drop speed as a reason.
- Cross-check against a second mapper: yes, but design C also delivers it.

## 2. Design options

### A. Mapper only; ViralScan allocates (build only if the gate passes)

Flow: `simpleaf index`, then `simpleaf quant --unfiltered-pl --dump-molecules --resolution <mode> --small-thresh 0` (it runs `piscem map-sc`, `generate-permit-list`, `collate` and `quant`; `--map-dir` skips mapping on reruns). `--resolution` is required and has no default. Simpleaf offers six modes: `cr-like`, `cr-like-em`, `parsimony`, `parsimony-em`, `parsimony-gene`, `parsimony-gene-em`. Alevin-fry adds `trivial`. The mode decides what `candidate_features` and `status` hold in the dump, so PR 0 must test it as a factor and the user must confirm it (section 6). The lower-level `piscem` and `alevin-fry` calls are an alternative. An adapter writes the files ViralScan already reads. `simpleaf quant` has no `--dump-eqclasses`; that flag exists only on `alevin-fry quant`.

- Data needed: per (cell, UMI) gene set, plus a barcode list, a gene list and a mapping rate.
- Does alevin-fry give it? Yes, with limits (run). `alevin/molecules.parquet` has `cell_barcode`, `rep_umi`, `n_umis`, `feature`, `candidate_features`, `reads`, `status` (`assigned` or `multi_gene`). `gene_eqclass.txt.gz` plus `geqc_counts.mtx` give per-cell class counts. A 60-cell dump took 9 s.
- Limit 1: the molecule table is after alevin-fry's parsimony collapse (`n_umis` > 1 on 2384 of 45477 rows). ViralScan intersects gene sets over all records of one CB+UMI (`multimapping.py:118-122`). The two can differ, and the difference depends on `--resolution`. Unverified how much.
- Limit 2: cells with fewer than `--small-thresh` reads (alevin-fry default 100) are resolved cr-like whatever `--resolution` says, so the molecule table mixes two resolution rules unless `--small-thresh 0` is passed. PR 0 must fix this setting.
- Limit 3: features are USA-expanded (S, `-U`, `-A` blocks, 255348 = 3 x 85116). The adapter must collapse them to genes. Viral transcripts are all S only for an index built with `simpleaf index --spliced viral_cdna.fa` next to a host splici build (`build_index.sh`, as in the comparison). A viral GTF passed through the splici builder would get U/A features for spliced viral genes.
- Keeps working: `multimapping.py` layers (`counts_unique`, `counts_multimap_*`, host-conservative default, em-cell, em-global), `cellcalling.py`, `detection.py`, `umap.py`, `gene_programs.py`, `hostresponse.py`, `virus_rollup.py`, `tier_measures.py`, `scautoqc.py`, `rerun programs`, `rerun multimap` layer swap, `validate-run`, `viralscan.compute`, the report template, `host_filter`, `read_filter`, `anello_align`, `star_confirm`.
- Breaks or needs a decision: evidence replay and `auto_evidence` (no `kallisto bus -n` analogue), `rerun multimap` full path (needs the retained molecule table), sensitivity (k-mer capture model), strand pilot, chemistry sanity thresholds, schema consts, golden digests, fingerprint of a directory index.

### B. Full backend on alevin-fry's gene matrix and resolution

Flow: take `quants_mat.mtx` and alevin-fry's own `--resolution` (cr-like, parsimony, parsimony-em), run ViralScan only on its output.

- Data needed: the final matrix only.
- Does alevin-fry give it? Yes (`quants_mat.mtx`, `quants_mat_rows.txt`, `quants_mat_cols.txt`; `--anndata-out` writes an h5ad).
- Keeps working: detection, umap and programs on a plain h5ad.
- Breaks: no `counts_unique`, no host-conservative default, no `adata_multimap.h5ad` v3 contract, no `viral_groups_*`, no molecule audit, no evidence tiers. `validate-run` skips the matrix contract. It conflicts with the project's host-conservative default. The comparison showed a different EBV-1/EBV-2 split (EBV-2 fraction 2.6 % vs 16.1 %), so results would not be comparable with kb runs.
- Verdict: cheap, but it drops the features that define ViralScan. Not recommended as a supported mode.

### C. Side-by-side secondary quantifier (default outcome after PR 0)

Flow: kb remains the primary run. A separate verb (for example `viralscan compare-quantifier`, or a repro-side script first) runs simpleaf on the same FASTQs and writes a comparison table (infected cells, per-gene rho, EBV-1/EBV-2 split, host UMI per barcode) beside the run.

- Cost: about 4-6 days after PR 0 (verb, tests, docs; or 2-3 days as a repro-side script). No change to the run contract, schemas, fingerprints, locks, containers or the touchpoint table below, which is where most of A's cost sits. It needs the three binaries on PATH, so the env solve check from PR 0 still applies, but a separate environment is enough.
- Keeps working: everything, because the primary run is untouched.
- What it delivers: the repeatable cross-check that the one-off comparison showed to be useful (infected cells within 6 %, per-gene rho 0.84-1.00), and the data for A's yardsticks if A is ever built.
- What it does not deliver: simpleaf results in the v3 layers, `validate-run` or the report. Promote to A only if a named benefit that C cannot give is shown (gate item 3 in section 1).

## 3. Touchpoint table (design A only)

Size: S under 0.5 day, M 0.5-2 days, L over 2 days. None of these rows applies to design C.

| Area | File:line | Change needed | Size |
|---|---|---|---|
| Quant rule | `src/viralscan/Snakefile:166-266` | Sibling rule chosen by `config["quantifier"]`; keep stage-and-swap, artifact check, `log/kb.done` sentinel (read at `:282`) | L |
| Output layout | `src/viralscan/kb_outputs.py:20-111` | Quantifier-aware layout, or write simpleaf products into `kb-python/counts_unfiltered/` names | M |
| Hard-coded kb paths | `scripts/analysis.py:159`, `scripts/cellcalling.py:719`, `detection.py:2269`, `run_safety.py:90-95,373`, `outputs.py:71,87,124-145`, `translated.py:232,241,570`, `run_context.py:4`, `evidence.py`, `scripts/evidence_run.py` | Route through the layout class first (own PR, no behaviour change); `translated.py` reads `counts_unfiltered/adata.h5ad` and stays kb-only | M |
| Molecule reader | `scripts/multimap.py:80-245,310-391` | New reader for `molecules.parquet`; replace `bustools correct/sort/text`; `read_ec` has no counterpart | L |
| Correction audit | `scripts/multimap.py:128-165` | Source from `generate_permit_list.json` `correction_stats` | M |
| Allocation | `multimapping.py:118-197,459-610` | Reuse; adapter feeds (cell, gene tuple) after S/U/A collapse | M |
| Cell calling | `scripts/cellcalling.py:275-312,525-568`, `emptydrops.R:49-75` | Unfiltered permit list matrix; USA-collapsed host layer; new receipt `count_layer` | M |
| Chemistry map | `chemistry.py:37-61,146-157,180-192` (on-lists come from kb-python `ngs_tools`), `kbparams.py`, `cli/parsers.py:1842-1890` | Table from kb `-x` names to simpleaf chemistry names and `--expected-ori fw/rc/both`; on-list source. The simpleaf registry has `10xv2`, `10xv3`, `10xv4-3p`, the `-5p` variants, ATAC, Flex and Visium entries, but no `10xv1` and no `dropseq` (both are in `CHEMISTRIES`, `chemistry.py:48-58`): custom geometry or refuse | M |
| Strand pilot | `strand.py:15-96` | New runner (`runner` is injectable at `:75`); recalibrate `tau=0.8` | M |
| Chemistry sanity | `chemistry_check.py:104-175`, `scripts/analysis.py:149-173` | Read `map_info.json`; re-measure 40 and 10 percent thresholds | M |
| Reference build | `scripts/build_reference.py:1654-1744`, `menu.py:1474-1560,1718-1812` | `simpleaf index` step; keep a kb-style `t2g.txt` with column 5 | L |
| Reference gates | `reference_gates.py:39-130` | No `kallisto inspect`; new index digest and planted-control gate | M |
| Virus identity | `virus_identity.py:430-435` | `manifest_path_for_index` puts `<index>.build_manifest.json` beside the index; define the path for a directory index | S |
| Config and CLI | `runconfig.py:147-186,639-669`, `defaults.py`, `cli/parsers.py:1799-1891` | `quantifier` field; reject or map kb-only flags | M |
| Required tools, env | `menu.py:72,1321-1333`, `validation.py:28-125,284-318` | Quantifier-dependent tool list; `simpleaf` profile; provenance for 3 binaries | M |
| Run contract | `run_safety.py:86-95,152-182,217,263,599-612` | Hash `quantifier`; fingerprint a directory index (`sha256_path` in `manifest_io.py` already hashes directories, and `:217` already admits a directory for `host_index`, so extend that rule to `index`); new resume check | M |
| Schemas | `src/viralscan/schemas/v3/run_manifest.schema.json:62` (const), `count_audit.schema.json:49,80`, `h5ad_contract.json:13` (`required_uns`), `table_columns.json:18,25`, `contract.py:7-16`; writers `scripts/multimap.py:149-164,590` and `run_safety.py:263` | Enum for `quantification_unit`, `bus_totals_source`, `barcode_correction`; version bump with kb runs unchanged | M |
| Sensitivity | `sensitivity.py:22-23,67-69,325,332` | `sha256_file(Path(config.index))` fails on a directory (use `sha256_path`); capture model not valid; mark "not assessed" | M |
| Detection labels | `scripts/detection.py:1840-1900`, `auto_evidence.py:97-121`, `docs/output_reference.md:374` | New `detection_source` value; stop printing "kallisto" for simpleaf | S |
| Discovery | `discovery.py:753,932-970,1205-1233,1499-1503` (`:753` requires bustools for barcode correction) | New barcode-space value; route (i) `not_assessed` without `ec_ids` | M |
| Evidence replay | `evidence.py:364-555`, `scripts/evidence_run.py:85-196,448-462` | Refuse for simpleaf with a clear message, or keep a kallisto index | M (refuse) / L (replay) |
| Rerun verbs | `menu.py:76-176,335-511` | Keep mtime guard for the new rule; refuse full path without the retained table | M |
| Inventory | `outputs.py:61-90,124-145`, `artifact_inventory.py:23-42` | Mark `af_map/`, `af_quant/` as intermediates (`map.rad` was 12.8 GB on EBV) | S |
| Demo and goldens | `demo.py:26-175` (`:113,141` need `kb` and `kallisto`), `tests/integration/` (`test_demo.py`, `test_golden_outputs.py`, `test_pipeline_counts.py`, `test_tiny_end_to_end.py`, `test_strand_replay_kallisto.py`), `tests/integration/golden/tiny_run.json` | Own simpleaf demo, tests and digests; kb digests unchanged; skip markers for missing simpleaf | M |
| Env and locks | `environment.yml:15-19`, `pixi.toml`, `conda-recipe/meta.yaml:28-30`, `environment.explicit.txt`, `environment.full.explicit.txt`, `environment.tools.lock.txt`, `environment.image.yml`, `environment.discovery.yml`, `environment.coverage.yml`, `uv.lock`, `tests/test_conda_lock.py` | Add simpleaf; separate env if it does not co-solve; new lock rows | L |
| Containers | `Dockerfile`, `Singularity.def`, `Singularity.pixi.def:41` (tool list text) | Add the three binaries to the images; image audit must see them | M |
| Licence | `conda-recipe/meta.yaml`, container notices | Check the licence of simpleaf, alevin-fry and piscem for redistribution in recipes and images; update licence fields and notices. Unchecked | S |
| CI | `.github/workflows/ci.yml:72,214,245,272`, `lock-audit.yml`, `image-audit.yml`, `release.yml`, `coverage-env.yml`, `notebooks-heavy.yml` | Separate simpleaf smoke job (ubuntu only); `ci.yml:72` runs pytest on `macos-latest`, so simpleaf-dependent tests need a skip marker there; lock and image audits and release checks must know the new binaries | M |
| Docs | `docs/v3_counting_contract.md:8-9`, `output_reference.md`, `faq.md`, `installation.md`, `limitations.md`, `CONTEXT.md:32-33,92-108` | Quantifier section; new unit definition | M |

## 4. Work breakdown

Effort is rough and assumes one developer. Order matters. Only PR 0 is committed; everything after it waits for the gate in section 1.

| # | PR | Days | Test or yardstick |
|---|---|---|---|
| 0 | Spike (repro repo, no package change): dump molecules on EBV, HHV-6B, KSHV and compare with BUS-intersection molecules, with `--resolution` as a tested factor (`cr-like`, `parsimony`, `parsimony-em` at least) and `--small-thresh 0` fixed; stage timing (kb count alone against simpleaf map+quant at equal threads, plus ViralScan downstream time); host S versus S+U+A (called-cell overlap, host UMI per barcode); EBV-1/EBV-2 split against planted truth after EM-CELL-AF-01; dry-run solve of `environment.yml` pins plus simpleaf 0.31.0, alevin-fry 0.19.0, piscem 0.25.0; test `simpleaf index --ref-seq` with `--decoy-paths` and a kb-style t2g passed at `simpleaf quant --t2g-map` (`simpleaf index` has no t2g option); test the toy demo index | 4-6 | Per-gene rho, infected-cell overlap, molecule counts per library, minutes saved per library. Stop conditions: molecule tables disagree materially, solve fails with no accepted fallback, or no named gain over C |
| C | Design C verb (or repro script) and its docs | 4-6 (2-3 as a script) | Comparison table reproduces the 2026-11-26 round |
| 0b | Fix REL-03 (`environment.explicit.txt:31` locks kallisto 0.52.0, `environment.full.explicit.txt:37` locks 0.51.1); needed before any new lock rows | 0.5-1 | `tests/test_conda_lock.py` |
| 1 | Agree the molecule format with the per-cell EM PR (one-page decision; depends on EM-CELL-AF-01 merge state) | 0.5 | Decision recorded in `.living/decisions.md` |
| 2 | Refactor: route all hard-coded `kb-python` paths through the layout class; no behaviour change. Start only after PR 0 passes, or ship it separately as a cleanup | 1-2 | Golden digests unchanged; full test suite in 3 chunks |
| 3 | `quantifier` option, config, fingerprint rule, schema enum bump, tool checks | 2-3 | `test_fingerprint_baseline.py` and `test_resume_fingerprint_unset_keys.py` unchanged for kb |
| 4 | Index builder: `simpleaf index`, kb-style `t2g.txt` beside it, directory fingerprint, reference gate | 3-4 | Planted control maps; index digest recorded |
| 5 | Chemistry and orientation tables (kb `-x` to simpleaf chemistry, `--expected-ori`), quant rule and converter: permit list, collate, quant, USA collapse, write the `counts_unfiltered` view | 5-7 | Matrix aligns by name with barcodes and genes; viral columns present |
| 6 | Molecule reader and adapter into `build_multimap_layers`; audit from `correction_stats` | 4-6 (may double if PR 0 shows the tables differ) | Molecule conservation (input = unique + ambiguous + unresolved); `validate-run` passes |
| 7 | Cell calling on the unfiltered simpleaf matrix; new receipt value | 2-3 | Called-cell overlap with the kb run on the three libraries |
| 8 | Strand and sanity threshold recalibration (needs approval) | 2-3 | Rates on the comparison libraries and a wrong-strand control |
| 9 | Refusals and labels: evidence, auto_evidence, discover route (i), sensitivity "not assessed", `detection_source` | 2 | Each verb exits with a clear message on a simpleaf run |
| 10a | Env, locks, containers, licence check | 3-4 | Lock and image audits pass; kb locks unchanged |
| 10b | Demo, golden digests, CI smoke job (ubuntu only), macOS skip markers | 3-4 | kb goldens unchanged; simpleaf goldens new |
| 11 | Docs and migration note | 2 | Docs drift test |

Totals. PR 0 plus C: 8-12 days (6-9 with a script). Design A after PR 0: 30-41 days summed from the rows above (PR 0b to 11). With a 15-20 percent contingency for the full-suite fix pass no PR budgets, A comes to about 40-57 days including PR 0. The ranges are judgment, not measured. They exclude review, merge and the waits for your approvals (impactful choices, thresholds). PR 0 can stop the project; it is also the only PR whose cost is committed now.

## 5. Validation plan

All pass criteria below are proposals. The thresholds are the user's call. The bounds in the second bullet come from the one comparison round on three herpesviruses that also motivated this scope, so they are not an independent yardstick. I could not open that note's saved figures, so they stay `unsupported: no inspected figure` until someone does.

- Planted truth panel (exact test in the repro): run kb and simpleaf on the same planted data. Proposed: no new false calls, and recovered molecules within a numeric tolerance the user approves (not "the kb run's tolerance", which is undefined).
- Planted low-abundance panel, including a non-herpes virus: explicit recall and false-call thresholds, set by the user. The three herpesvirus libraries do not cover low abundance.
- EBV, HHV-6B, KSHV comparison: rerun the round with the backend. The note reports the same viruses, infected cells within 6 %, per-gene rho 0.84-1.00. Proposed: keep those bounds, noting rho 0.84 is a permissive lower bound.
- Molecule conservation: the existing audit identity must hold for simpleaf runs. Proposed: exact.
- Allocation equivalence: counts from `molecules.parquet` against the BUS-intersection definition on one library, under each tested `--resolution`. Proposed: report the difference; decide the gate after PR 0.
- Golden digests for the kb path: byte-identical. Proposed: hard gate on every PR.
- EBV-1/EBV-2 split: gate it against the planted truth after EM-CELL-AF-01, and report it under ViralScan allocation and under alevin-fry `parsimony-em`. The note says the split depends on resolution.
- Strand and chemistry: recalibrated thresholds measured on the three libraries plus one wrong-orientation control.

## 6. Impactful choices for the user

Each needs its own approved plan row. "Recommendation" is the package default where one exists. Where none exists, no value is recommended. Rows tagged `inherited:` are choices an earlier stage already made; each needs a one-line confirm, not a new decision, unless the backend changes the data it acts on.

| Choice | Options | Recommendation |
|---|---|---|
| Build or stop after PR 0 | stop at C; build A | Stop at C unless the gate in section 1 passes. Needs your approval of the gate thresholds. |
| Index type | spliced+intronic via splici (simpleaf default); host splici plus viral FASTA via `--spliced` (used in the comparison, `build_index.sh`); spliced only; direct `--ref-seq` cDNA+viral | Package default for `simpleaf index`: `spliced+intronic` (simpleaf 0.31.0). The `--spliced` route is what built the comparison index, with viral all S. Direct `--ref-seq` matches the kb index but is not a simpleaf default. No recommendation among the non-default routes without a measured comparison (PR 0). |
| Host layer used for depth, emptyDrops, ambient screen, LOD | S, S+U+A, S only | None: no package default. This changes host UMI totals per barcode and the called-cell set. PR 0 measures the shift; count intron-aware counts as a benefit only if it moves a downstream result. |
| Dump `--resolution` mode | `cr-like`, `cr-like-em`, `parsimony`, `parsimony-em`, `parsimony-gene`, `parsimony-gene-em` (`trivial` is alevin-fry only) | None: `simpleaf quant --resolution` is required and has no default. It changes the molecule table that `build_multimap_layers` consumes, hence `counts_unique`. PR 0 tests it; you pick. |
| `--small-thresh` | alevin-fry default 100 (cr-like below it); 0 | A molecule dump for ViralScan needs a user-approved value. 0 is the value to test, with no yardstick yet beyond PR 0. |
| Allocation owner | ViralScan allocation on the molecule table; alevin-fry resolution as is | ViralScan allocation (design A), so the host-conservative layer exists. |
| Cell list | `--unfiltered-pl` with explicit `--min-reads`; `--knee` | `--unfiltered-pl` (emptyDrops needs all on-list barcodes). `--min-reads`: simpleaf 0.31.0 default is 10. Knee called 2137 cells against 2647 on EBV, so it is not a substitute. |
| Barcode correction policy | `--cell-bc-correction` `unique` or `frequency` | Package default `unique` (alevin-fry 0.19.0). `generate_permit_list.json` from the EBV comparison run shows `unique`. Audit meaning changes versus bustools. |
| Barcode neighbourhood | `hamming-1`, `substitution-or-shift-1`, `edit-1` | Package default: alevin-fry chooses per chemistry. The EBV comparison run resolved `substitution-or-shift-1`. |
| Barcode confidence | decimal or fraction | Package default 97.5 % for RNA (alevin-fry 0.19.0), used only under `frequency`. Inactive with the default `unique`. |
| Orientation map | kb `forward/reverse/unstranded/auto` to `fw/rc/both` (simpleaf defaults to `fw` for 10xv2/10xv3, `both` otherwise) | Table from chemistry geometry; `tau=0.8` recalibration needs approval. |
| Chemistry mapping and supported chemistries | kb `-x` names to simpleaf registry names or geometry strings | `10xv2` and `10xv3` map one-to-one; kb `10xv4` maps to `10xv4-3p`; 5' variants need a table. `10xv1` and `dropseq` are not in the registry: custom geometry or refuse. Unknown for tchem emit. You approve the supported list. |
| Decoys | pass `dlist.fa` to `--decoy-paths` or none | No recommended value until tested (PR 0); it changes viral-host homology calls. |
| Evidence under simpleaf | refuse; keep kallisto index; build piscem replay | Refuse (cheapest). Needs your approval. No yardstick for the alternatives. |
| Backend contract | simpleaf must feed multimap (v3 contract) or may run without | Must feed multimap. |
| Compat | keep `kb-python/` name and `log/kb.done`, or rename | Keep, with a doc note. A rename breaks digests and resume. |
| `quantifier` default entry | add to `DEFAULTS` (every kb fingerprint changes) or leave unset (no change) | Leave unset for kb. |
| Platforms | linux-64 only, or also macOS | linux-64 first, with simpleaf tests skipped on the `macos-latest` leg. |
| Version pin | simpleaf 0.31.0 or 0.32.0 | 0.31.0, pinned with binary digests; both force an index rebuild. |
| emptyDrops parameters and seeds on the simpleaf matrix, `inherited: src/viralscan/scripts/cellcalling.py`, `src/viralscan/scripts/emptydrops.R` | keep or change | Keep; the matrix changes (USA-collapsed, simpleaf barcodes), so confirm in one line. Seeds: no recommended value. |
| Viral-feature collapse (S only or S+U+A) | S only; S+U+A | `inherited: src/viralscan/multimapping.py`. Confirm. |
| Ambient screen and LOD thresholds, `inherited: src/viralscan/scripts/detection.py` | keep or change | Keep; the mapper changes from kallisto k=31 to piscem, so confirm in one line. |

## 7. Risks and open questions

- Molecule equivalence (high). Alevin-fry's table is post-parsimony. If it differs from the BUS intersection, `counts_unique` shifts. Unverified. PR 0 measures it.
- Value (high). The gain over a kb run once EM-CELL-AF-01 lands is not shown (section 1). PR 0 measures it before any build.
- Dump format stability (high). `--dump-molecules` and `--dump-eqclasses` have no documented contract. USA-mode parsimony prints "EXPERIMENTAL" in the 0.19.0 log. Simpleaf released about one minor per month in 2026 (0.25 to 0.32), and 0.27, 0.29 and 0.31 each raised minimum versions; 0.29 and 0.31 required rebuilt indexes. Pin all three binaries and record digests in the fingerprint.
- Evidence replay (high). Piscem `map-sc` writes RAD only. `alevin-fry view` has no read names or sequences. Read-level evidence and `counted_virus_best` (`molecule_verdict.py:16-17`) stay kallisto-bound.
- Sensitivity (high, scientific). The k=31 exact-capture model and the positive control certify negative calls. They do not transfer to piscem. Report "capture not measured" until measured.
- Disk (medium). `map.rad` was 12.8 GB and `map.collated.rad` 6.7 GB on EBV. `rerun multimap` needs the molecule table kept. Decide retention.
- Runtime and memory (medium, was low). Simpleaf took 12.7, 4.5 and 13.5 min per library, 5.2-5.6 GB peak (comparison round). ViralScan `run` times are for the whole pipeline, and design A still runs multimap, cell calling, detection and umap after simpleaf. The saving is unknown until PR 0 times kb count alone and the downstream stages.
- CPU features (medium). The earlier "Illegal instruction" crash on a non-avx2 node was not reproduced. Which of the three binaries needs avx2 is untested. Slurm jobs need `--constraint=avx2`.
- macOS (open). Only simpleaf's platform list (linux-64, linux-aarch64, osx-64, osx-arm64) was checked. Alevin-fry and piscem builds were not. `ci.yml:72` runs pytest on `macos-latest`; adding simpleaf to the environment files needs a skip path there.
- Environment and locks (medium). Whether simpleaf co-solves with `kb-python=0.30.2`, `snakemake=9.23.1`, `star=2.7.11b`, `r-base=4.5.3` is unchecked; PR 0 runs the dry-run solve. A separate `environment.simpleaf.yml` would mean evidence, discover and translated cannot run in it. REL-03 is still `[~]`; PR 0b fixes the kallisto lock mismatch first.
- Licence (open). Redistribution terms for simpleaf, alevin-fry and piscem in conda recipes and images are unchecked.
- Run contract (medium). `quantification_unit` is a schema const (`run_manifest.schema.json:62`) and is a fingerprinted field (`run_safety.py:156`). Replacing the constant changes every fingerprint and breaks `--resume`; adding an enum value does not. `count_audit.schema.json` is closed. Old kb runs must keep their fingerprints.
- Directory index (medium). `run_safety.py:217` fingerprints files and `host_index` directories only, `sensitivity.py:332` hashes a file, `virus_identity.manifest_path_for_index` expects `<index>.build_manifest.json`.
- Barcode space (medium). Alevin-fry's corrected space may not match bustools'. `discovery.py:959` would label simpleaf receipts "kb".
- CI cost (medium). A piscem index of the cat42d panel is about 5.3 GB. Use the demo fixture for CI. Whether the toy `transcripts.fasta` and `viral.gtf` build a piscem index is untested.
- Unverified: whether `simpleaf index --ref-seq` accepts the viral+host FASTA with a kb-style t2g at quant time; whether `--decoy-paths` matches kb's d-list masking (the comparison report said simpleaf has no d-list equivalent; the mapping agent found poison k-mer decoys as the nearest feature); whether tchem has a simpleaf emit mode.

## Review record

Review 1 (value versus cost):
- Changed "Build it" to "Do not build yet; run PR 0, then decide", and said so in section 1.
- Added the "What A buys after EM-CELL-AF-01" subsection and a written gate (molecule agreement, env solve, a named gain C cannot give).
- Costed design C (4-6 days, 2-3 as a script), made it the default outcome, and marked the touchpoint table as design A only.
- Moved speed from "low risk" to "medium"; added stage timing to PR 0 with a 10 min per library cutoff.
- Stated that viral specificity gains nothing from USA; added a PR 0 measurement of host S versus S+U+A.
- Added `--resolution` as a PR 0 factor and as its own row in section 6.

Review 2 (technical facts):
- Split the barcode-correction row into policy (default `unique`), neighbourhood and confidence (97.5 %, only under `frequency`). Checked against `alevin-fry generate-permit-list --help`.
- Corrected the resolution mode list (six in simpleaf, `trivial` alevin-fry only); checked against `simpleaf quant --help`.
- Added `10xv1`, `dropseq` (not in the registry) and `10xv4-3p` (rename) to the chemistry row; checked against `chemistries.json` and `chemistry.py:48-58`.
- Added the `--spliced` index route and qualified limit 3.
- Fixed three wording slips: constant versus enum fingerprint effect, `host_index` directory admission at `run_safety.py:217`, schema const at `:62`.
- Not re-run: the `molecules.parquet` column list, 2384/45477 and timings.

Review 3 (plan completeness):
- Moved the chemistry and orientation tables from PR 8 into PR 5, so PR 5 can run.
- Added the env dry-run solve to PR 0 as a stop condition, and PR 0b for REL-03.
- Made PR 2 start only after PR 0 passes; recorded the EM-CELL-AF-01 merge state as a PR 1 dependency.
- Split PR 10 into 10a (env, locks, containers, licence) and 10b (demo, goldens, CI, macOS skips); added a 15-20 percent contingency and a note that PR 6 may double. Total restated as about 40-57 days for A including PR 0.
- Section 5: said the bounds are circular, added a planted low-abundance panel, a numeric tolerance for the user to set, and an EBV split gate. The saved comparison figures are still uninspected.
- Section 6: added `inherited:` rows (emptyDrops, viral-feature collapse, ambient and LOD thresholds) and a "no recommended value" label where no default or yardstick exists.
- Added licence and macOS CI rows.
