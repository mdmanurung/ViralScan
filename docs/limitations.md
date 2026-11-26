# What ViralScan cannot tell you

Read this page before you interpret a result. Each point below is stated elsewhere in the
documentation; this page keeps them in one place. No sensitivity, specificity or comparison with
other tools has been established for v3 yet (the truth-panel gates are open), so none is claimed.

## A nonzero count is not an infection

A viral molecule in a cell is evidence of viral sequence in that droplet. It is not proof of
infection. Host homology, endogenous viral elements, ambient RNA, barcode swapping and doublets can
all put viral sequence in a cell that is not infected. The `detected` column and the evidence flags
report what passed the gates; they do not turn an ambiguity flag into a biological conclusion. See
[Results](results.md) and the [FAQ](faq.md).

## A zero is a statement about depth and capture

"Not detected" means no molecule passed the gates at the depth sequenced. Whether that excludes a
virus depends on the detection limit at that depth and on how well the chemistry captures the
virus's transcripts (poly-A capture varies by virus and gene). `informative_negative` is false unless
a control or measured capture supports the negative. `indeterminate` is not a negative. See
"ViralScan found nothing. Is there really nothing there?" in the [FAQ](faq.md).

## Chemistry and sample type

- Single-cell RNA-seq with a cell-barcode and UMI geometry only. Bulk RNA-seq is not supported.
- 3' chemistries (10x 3', Drop-seq) are the supported default. For 10x 5' and unconfirmed
  libraries, `run` stops unless you state `--strand`; reads can be lost under the wrong strand, and no
  accuracy claim is made for 5'.
- Human host droplet data is the design target. Other hosts need `build ref --host` and have not been
  validated.

## Allocation and thresholds are model choices

Multimapping reads are allocated by a stated method (default `host-conservative`: host-virus
ambiguous reads get weight zero for viral genes). The infected-cell rule (`--infected-min-molecules`,
default 2.0) and the detection gate (`--detection-min-cells`, `--detection-threshold`) are defaults
that have not been validated against the truth panel; `defaults_status` in `run_manifest.json` says
`provisional` until they are. Report the values you used.

## Gene programmes (latent or lytic)

Programme states exist for a small set of viruses, the `reactivating` state is exploratory, and the
`mixed` state depends on viral load. They describe marker evidence in a cell, not a verified
biological state.

## References and annotations

Results depend on the reference index. The viral annotation panel is a curated working set, and its
licence and archive status are tracked in the [reference panel](reference_panel.md) page. A virus
absent from the index cannot be detected.
