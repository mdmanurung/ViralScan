# ViralScan 2.2.0's multimap correction bypasses UMI deduplication and is reported as "UMIs"

**Date**: 2026-07-27
**Branch**: `codex/viralscan-v3`
**Source**: `git show 8e3e3ed:src/viralscan/scripts/multimap.py` (the 2.2.0 commit)
**Status**: confirmed from source; mechanistic, not yet quantified against truth

## Finding

ViralScan 2.2.0's reported viral value is the sum of two quantities measured in
**different units**. The unique-mapping term is UMI-deduplicated; the
multimapping term is not.

### The unique-mapping term is correct

`counts_original` is `adata.h5ad` from `kb count`. `bustools count` collapses by
corrected cell barcode plus UMI, so these are genuine molecules.

### The multimapping term never deduplicates

`multimap.py` parses the BUS text with all four fields:

```python
bus_df = pd.read_csv(txt_file, sep="\t", header=None,
                     names=["barcode", "umi", "ec", "count"])   # line 322
```

and then never uses the UMI column again:

```python
bc, ec, count = row.barcode, row.ec, row.count                  # line 203
...
share = count / len(genes_in_ec)                                # line 225
```

`count` in a BUS record is the number of **reads** supporting that
(barcode, UMI, EC) triple. The correction therefore distributes read counts
fractionally across the genes of each equivalence class. No UMI collapse occurs
anywhere on this path.

### The two are then added

```python
adata.layers["counts_original"]  = adata_orig[:, adata.var_names].X.copy()
adata.layers["counts_corrected"] + adata.layers["counts_original"]
```

`counts_corrected` holds only the multimap shares — unique-mapping ECs are
deliberately skipped to avoid double counting, per the source comment. So the
reported total is *UMI-deduplicated molecules for unique mappers plus
read-weighted fractional shares for multimappers*.

### It is labelled UMIs

`final_results` computes `viral_counts` from `adata.X`, which is the corrected
matrix alone, and writes it to the summary as:

```
Total viral UMIs (corrected): <value>
```

That quantity contains no UMI collapse. The adjacent line, `Viral UMIs in
original (not corrected)`, is genuinely UMIs. The two labels use the same word
for different things.

## A second, smaller mixture

`multimap.py` reads `output.bus` — the raw kallisto output — rather than the
corrected and sorted unfiltered BUS that produced `counts_original`. So the
barcodes feeding the multimap term have not been whitelist-corrected, while the
term it is summed with came from the corrected path. Two barcode treatments in
one number.

`normalize_barcodes` only strips a `-1` suffix; it does not correct sequences.

## What this explains

- **Fractional reported values.** `share = count / len(genes_in_ec)` is why the
  archive contains values like `31,011.3333`.
- **The magnitude gap.** EBV `SRR12682296` reports 1,792,860 in the archive
  against 256,487 v3-unique molecules. The multimap term is inflated relative to
  the molecule term it is added to, by roughly the PCR duplication factor.
- **Why skin diverges most.** Legacy-to-v3 cell-set Jaccard is 0.975 for EBV
  controls but 0.118 for skin under host-conservative. Skin signal is low-mass
  and multimapper-dominated, which is exactly where the un-deduplicated term
  dominates the sum.

## What it does and does not establish

**Does:** the two versions' outputs are not on a common scale, and no rescaling
can put them on one, because the archive's value is a sum of two different
units whose mixing ratio varies per sample with duplication rate and multimapper
fraction. This is a firmer statement than "different counting rules".

**Does not:** show that v3 is more accurate. Demonstrating that still requires
planted ground truth — precision, recall, and F1 against known molecules, plus
false-call rate on negatives. `VAL-01` must build the truth panel first, and
2.2.0 is now a scored comparator there (`W10a`, `W10b`, added 2026-07-27).

## Provenance note

Established by reading the 2.2.0 source from git history. Emma's archived results
and environment were not touched, consistent with the read-only preservation rule
in `analysis/legacy_v2_v3/protocol.yaml`.
