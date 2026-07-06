# Analysis Manifest

<!-- One entry per registered analysis. See each analysis's UPPER_SNAKE_CASE.md doc. -->

### reference_strategy_benchmark
```yaml
name: reference_strategy_benchmark
question: Which reference strategy / aligner best detects each target virus (EBV, HHV-6B, HSV-1)?
inputs: [reference_strategy_refs, reference_strategy_fastqs]
outputs: results/reference_strategy_benchmark.tsv
doc: analysis/reference_strategy_benchmark/REFERENCE_STRATEGY_BENCHMARK.md
status: incomplete (4/12 rows complete; all EBV rows failed)
findings: "HHV-6B: ViralScan calls ~54% more positive cells than STARsolo (n=1 cond., count-layer + HHV-6A/6B caveats); HSV-1 undetected under both strategies"
blocked: Selectivity-Index (3b) and noisy-channel (7b) analyses need EBV rows + more complete conditions
tags: [reference-strategy, aligner, benchmark, viral-detection, incomplete]
```

Reference-strategy/aligner benchmark. Only 4/12 rows complete (EBV failed), so the
planned Selectivity-Index / channel analyses are blocked; documents the 2 valid
comparisons (HHV-6B aligner; HSV-1 strategy) and the completion needed. See doc.

### ideas/2026-07-01-cross-disciplinary-brainstorm
```yaml
name: 2026-07-01-cross-disciplinary-brainstorm
type: ideation
personas: [evolutionary-biologist, quantitative-geneticist, pharmacologist, stem-cell-biologist, causal-inference, ecologist, information-theorist]
ideas: 14
dir: analysis/ideas/2026-07-01-cross-disciplinary-brainstorm/
status: generated
tags: [ideas, brainstorm, ebv, host-response, reference-benchmark, anellovirus]
```

7-persona brainstorm on the EBV host-response finding + reference-strategy benchmark.
14 ideas grounded in real data; see `00_index.md` (grouped by feasibility). Convergent
themes: the specificity gap as a possible primed intermediate state; the 10-UMI
threshold sweep; three reframings of the reference benchmark; the underused
anellovirus panel.

### hostresponse_ebv_matched
```yaml
name: hostresponse_ebv_matched
question: Does the host transcriptome discriminate EBV+ from EBV− cells within one LCL sample?
script: scripts/hostresponse_ebv_matched.py
doc: analysis/hostresponse_ebv_matched/HOSTRESPONSE_EBV_MATCHED.md
inputs: [reference_strategy_fastqs (SRR12682296), GEO GSE158275 paper barcodes]
outputs: results/hostresponse_ebv_matched/
conventions: [robust-analysis, bioinformatics]
status: complete (depth-confounded headline; resolved effect AUC ~0.72)
headline: raw AUC 0.866 was depth-inflated (depth-alone 0.80–0.97); depth-MATCHED re-analysis gives the honest effect = AUC 0.72 (15 genes) / 0.68 (5 depth-robust genes) with depth-alone at chance (0.48). Real but moderate signal.
report: analysis/hostresponse_ebv_matched/reports/hostresponse_ebv_matched-report.pdf (comprehensive, 6pp, 2026-07-01; depth caveat + matched-design resolution)
followups: scripts/depth_confounder_check.py (confound diagnostics) + scripts/depth_matched_reanalysis.py (the fix) → results/hostresponse_ebv_matched/depth_{confounder,matched_reanalysis}.txt
tags: [ebv, host-response, scrna-seq, classifier, matched-cells, manuscript]
```

Matched-cell EBV host-response analysis (manuscript). Cross-validated host-gene
classifier separates EBV+ from EBV− cells at AUC 0.845. Registered values in
`analysis/hostresponse_ebv_matched/outputs/numbers.json`. scilintr clean.

### multimap_profiling
```yaml
name: multimap_profiling
question: Which of ViralScan's 4 multimapping methods is fastest, and what should be optimized?
scripts:
  - analysis/multimap_profiling/scripts/fast_profile.py        # main profiler (zip-loop + vectorised EM code)
  - analysis/multimap_profiling/scripts/profile_multimap.py    # RETIRED — old itertuples code
  - analysis/multimap_profiling/scripts/interpret_cprofile.py  # post-processing for cProfile text
  - analysis/multimap_profiling/scripts/register_values.py     # writes numbers.json
doc: analysis/ANALYSIS_MANIFEST.md (this entry)
inputs:
  - /exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_runs/reference_strategy_2026-06-28_fresh12b/runs/ebv__viralscan__combined/SRR12682296/kb-python/output.bus.txt (103M rows)
  - /exports/archive/hg-funcgenom-research/mdmanurung/viralscan_showcase/viralscan_showcase/fullrun/refs/merged/t2g_plus_anellovirus.txt
  - kb-python counts_unfiltered/adata.h5ad, matrix.ec, transcripts.txt, cells_x_genes.barcodes.txt
outputs: analysis/multimap_profiling/outputs/
conventions: [robust-analysis]
status: complete (2026-07-05)
headline: |
  _matrix_value scipy.__getitem__ dispatch = 86.4% of build_multimap_layers (cProfile on 1M rows,
  pre-fix code). Fixed in commit 3c53ad7 (direct CSR buffer access: indptr/indices/data +
  searchsorted). Three perf commits total: b7e9635 (itertuples→zip + EC-invariant hoisting),
  2c2e6f0 (vectorise em_gene_abundances sparse matvec), 3c53ad7 (kills __getitem__ bottleneck).
  Combined ~4× speedup. Old-code equal anchor: 19,964s full data. Methods are
  wall-time-equivalent (equal≈hc≈uw ±5%, all dominated by the same _matrix_value bottleneck).
  em wall-time NOT AVAILABLE (process killed; re-run would mix code versions).
numbers: analysis/multimap_profiling/outputs/numbers.json (27 values)
scilintr: 0 findings
code_version_note: |
  cprofile_equal.txt + cprofile_em.txt + fast_profile wall-time measurements (equal/hc/uw) were
  all taken on code with b7e9635 + 2c2e6f0 applied but WITHOUT 3c53ad7 (the major fix).
  The fast_profile_resume.py em timing (started but process killed) would have run on current
  code (all 3 commits applied) — mixing code versions in one table is misleading, so em is
  listed as NOT AVAILABLE. The cProfile outputs are the definitive before-fix provenance.
anchors:
  n_bus_records: 103145071
  n_cells: 848191
  n_genes: 43451
  n_ECs: 388677
  n_multi_gene_ECs: 317685
  multimapping_rate: 81.7%
  n_viral_genes: 4845
  rss_after_load_mb: 4984
  old_code_equal_wall_s: 19964.2    # OLD itertuples code on full data (profile_multimap.py)
  cprofile_matrix_value_fraction: 86.4%  # _matrix_value cumtime / build_multimap_layers cumtime
  walltime_equal_sub_s: 1139.75    # pre-3c53ad7 code, 5M rows (upper bound due to head-slice bias)
  walltime_hc_sub_s: 1078.88
  walltime_uw_sub_s: 1136.83
  fix_commit: 3c53ad7
  fix_speedup_approx: ~4× (all three commits combined)
tags: [multimapping, profiling, performance, em, benchmarking, complete]
```

Empirical performance profile of ViralScan's 4 multimapping methods on full-depth EBV LCL
(SRR12682296, 103M BUS records). Headline: `_matrix_value` scipy `__getitem__` dispatch = 86.4%
of `build_multimap_layers` runtime (cProfile, 1M rows, pre-fix code). The three methods
(equal/hc/uw) are wall-time-equivalent within ±5%, all dominated by the same bottleneck;
em adds ~9s overhead per 1M rows for the EM iteration loop. Bottleneck fixed in commit `3c53ad7`
(direct CSR buffer access: `indptr/indices/data` + `searchsorted`), preceded by `b7e9635`
(itertuples→zip + EC-invariant hoisting) and `2c2e6f0` (vectorise `em_gene_abundances`).
Combined ~4× speedup, byte-identical output, all 557 tests pass. Registered 27 values in
`analysis/multimap_profiling/outputs/numbers.json`. scilintr 0 findings.
