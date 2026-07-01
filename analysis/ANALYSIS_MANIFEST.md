# Analysis Manifest

<!-- One entry per registered analysis. See each analysis's UPPER_SNAKE_CASE.md doc. -->

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
