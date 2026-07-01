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
status: complete
headline: AUC 0.866±0.036, MCC 0.570±0.077 (n=1906; 1179 EBV+ / 727 EBV−; ≥10 UMI; leakage-corrected)
report: analysis/hostresponse_ebv_matched/reports/hostresponse_ebv_matched-report.pdf (comprehensive, 5pp, 2026-07-01)
tags: [ebv, host-response, scrna-seq, classifier, matched-cells, manuscript]
```

Matched-cell EBV host-response analysis (manuscript). Cross-validated host-gene
classifier separates EBV+ from EBV− cells at AUC 0.845. Registered values in
`analysis/hostresponse_ebv_matched/outputs/numbers.json`. scilintr clean.
