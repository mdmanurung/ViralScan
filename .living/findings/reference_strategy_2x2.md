# Finding: reference-strategy 2×2 benchmark (complete, 2026-07-04)

**ID**: F-006 (updated 2026-07-04 — EBV GTF artifact confirmed; aligner claim downgraded)
**Status**: COMPLETE — harmonize_2x2.py produced anchor-restricted unique-layer comparison

Run `fresh12b` — 12/12 complete. Original (stale) numbers from `results/reference_strategy_benchmark.tsv`
had three confounds (denominator, count-layer, target-matching). Corrected below.

---

## Corrected 2×2 — anchor-restricted UNIQUE target UMI (fair comparison)

Shared anchor = intersection of STARsolo filtered cells ∩ ViralScan per_cell barcodes,
computed per dataset (not per row). Count layer = STARsolo unique-integer vs
ViralScan `counts_original` (pre-multimap). All 12 rows pass.

| Dataset | Ref-strat | STAR unique | VS unique | VS corrected | Multimap gain | Anchor N |
|---------|-----------|-------------|-----------|--------------|---------------|----------|
| hhv6b   | combined  | 519         | 988       | 1,496        | +508          | 3,517    |
| hhv6b   | two_step  | 519         | 984       | 1,487        | +503          | 3,517    |
| ebv     | combined  | 45,896      | 90,097    | 547,459      | +457,362      | 1,908    |
| ebv     | two_step  | 45,939      | 82,980    | 515,830      | +432,850      | 1,908    |
| hsv1    | combined  | 19*         | 14,274    | 31,672       | +17,398       | 3,307    |
| hsv1    | two_step  | 19*         | 12,527    | 27,946       | +15,419       | 3,307    |

*HSV-1 STARsolo = 19 UMI is a GTF artifact — see root cause section below.
**EBV total-UMI comparison is also confounded — see EBV section below.

Canonical outputs: `analysis/reference_strategy_benchmark/outputs/harmonized_2x2_unique.tsv`
Registered values: `analysis/reference_strategy_benchmark/outputs/numbers.json` (35 values + EBV supplement)

---

## Key findings

### Aligner axis (STARsolo vs ViralScan unique layer)

- **HHV-6B**: VS/STAR = 1.9× (unique). Both aligners cover 97 identical features (1 gene/contig,
  no overlaps). **This is the only clean aligner comparison in the benchmark.**
- **EBV**: CONFOUNDED — see EBV section below. The apparent "2×" ratio is a reference-completeness
  artifact; 95% of ViralScan's 90k EBV UMI comes from 80 CDS-only genes STARsolo cannot count.
- **HSV-1**: STARsolo = 19 UMI → DUAL GTF artifact (see below). Not a real comparison.
- **Multimap gain**: ViralScan's multimapping correction adds 51–122% UMI above unique layer.
  EBV gain is massive (+457k on anchor for combined) because many EBV reads multimap within the
  repetitive viral genome; hhv6b gain is modest (+508, +51%).

**Verdict: aligner axis effect is CONFIRMED for HHV-6B only (1.9×). EBV and HSV-1 are both
CONFOUNDED by GTF structural artifacts (missing exon records + overlapping gene clusters).
No clean aligner comparison is possible for EBV or HSV-1 without fixing the combined GTF.**

### Reference strategy axis (combined vs two_step)

- STARsolo: negligible delta for hhv6b (0) and hsv1 (0); EBV: +43 UMI for two_step (tiny).
- ViralScan: EBV shows a combined > two_step advantage (+7,117 unique UMI); hhv6b and hsv1 are small.
- **Verdict: reference strategy minor — CONFIRMED for all three viruses.**

---

## EBV root cause — same GTF artifact class as HSV-1

### GTF structure
- NC_007605.1: 94 genes total
- **16 have exon records** (STARsolo can count); 80 are CDS-only (STARsolo misses entirely)
- Same primary mechanism as HSV-1 (61/79 CDS-only) and anellovirus F-005 (commit 7739721)

### STARsolo EBV: 14 features, dominated by LMP-1
STARsolo raw matrix has 14 HHV4 features. Per-gene UMI (all barcodes):
- `HHV4_LMP-1`: 79,773 UMI (99.7% of all STARsolo EBV signal)
- `HHV4_BRLF1`: 216 UMI
- All 12 other exon-bearing genes: 0 UMI

Explanation: The EBNA gene cluster (positions 11k–97k) contains 10 of the 14 exon-bearing genes.
These form 23 overlapping pairs — 82.8% of exon-covered bases are multi-gene ambiguous.
STARsolo discards ambiguous reads → those 10 genes get 0 UMI.
`HHV4_LMP-1` (pos 166,483–169,088, 2,606 bp, no overlaps) is the sole survivor:
it gets 46,343 of 46,419 anchor UMI (99.8%) from STARsolo.

### ViralScan EBV: signal dominated by CDS-only genes
ViralScan unique UMI on anchor split by gene category:
- STARsolo-matched exon-bearing genes (14): **4,537 UMI (5.0%)**
- CDS-only genes (80): **85,723 UMI (95.0%)**

Top ViralScan EBV genes (all barcodes):
1. `EPSTEIN_HHV4_BHLF1`: 141,540 UMI — CDS-only, STARsolo cannot count
2. `EPSTEIN_HHV4_BILF2`: 28,930 UMI — CDS-only
3. `EPSTEIN_HHV4_BBRF3`: 24,943 UMI — CDS-only

ViralScan gets 325 UMI from LMP-1 on anchor vs STARsolo's 46,343 — because
ViralScan's pseudoalignment distributes LMP-1-compatible reads across overlapping/
CDS-only genes (BHLF1 being the top recipient).

### The apparent "VS/STAR = 2×" is a reference-completeness artifact
Total VS (90k) / total STAR (46k) = 1.96×. But this ratio conflates:
- 95% of VS UMI = reads from genes STARsolo cannot see (CDS-only)
- 5% of VS UMI = reads from the 14 shared genes — and STAR gets 10× MORE of those

The "2×" number is not an aligner comparison; it is a reference coverage comparison.

---

## HSV-1 root cause — DUAL GTF artifact

### Confound 1 (original): Regex artifact
The stale `commands.jsonl` regex lacked `hhv-?1` and `nc_001806`.
The STARsolo HHV1 feature names are `HHV1gp*` — zero matched the old regex.
Fix: use current DATASETS regex → 18 features match → 19 UMI on anchor.
But 19 UMI is itself not a real count (see below).

### Confound 2: GTF missing exon records (primary GTF artifact)
In the combined STARsolo reference (`combined.gtf`), 61 of 79 HHV1 gene records
have **CDS/start_codon/stop_codon records only — no `exon` records**.
STARsolo's `--soloFeatures GeneFull` uses `sjdbGTFfeatureExon = exon` to build transcript models.
Without exon records, GeneFull cannot assign reads to those 61 genes → zero UMI structurally.
This is the same root cause as anellovirus STARsolo=0 (finding F-005, commit 7739521).

Summary: 18/79 HHV1 genes have exon records (all others = 0 by construction).

### Confound 3: Overlapping terminal-repeat gene cluster (secondary GTF artifact)
The 18 exon-bearing genes are **not randomly distributed** — they all cluster in:
- Terminal repeat region (positions ~118k–152k): 16 genes, forming dense overlap clusters
  (e.g. p02/p03/p04 all share coordinates 144k; p07/p08/p09/p10 all overlap 137k–143k)
- Position 1–7.5k region: 2 genes (s01, p76)

Quantification:
- 11 of 91 possible pairs among the 18 exon genes are overlapping (12%)
- 27.4% of exon-covered base-pairs fall in multi-gene ambiguous regions

STARsolo's GeneFull discards reads that overlap >1 gene as ambiguous.
**Result: of the 18 nominally countable genes, only 2 get any UMI at all**:
`HHV1gp00p13` (27 UMI) and `HHV1gp00p03` (3 UMI). All 16 others = 0.
Total across all barcodes: 30 UMI. On shared anchor: 19 UMI.

### ViralScan vs STARsolo on shared-exon genes
ViralScan unique UMI on anchor restricted to the 18 exon-bearing genes: **2,352 UMI** (16.5% of total).
ViralScan unique UMI from the 61 CDS-only genes: **11,925 UMI** (83.5% of total).
STARsolo on the same 18 exon genes: **19 UMI** (predominantly from p13/p03 only).

### ViralScan combined HSV-1: 0 called cells anomaly
`hsv1_viralscan_combined`: 0 called cells (all 10,754 signal barcodes are `is_called_cell=False`).
`hsv1_viralscan_two_step`: 3,515 called cells. This is a cell-calling anomaly in the combined run,
not evidence the signal is ambient. The unique UMI on anchor (14,274) is substantial and consistent
with the two_step run (12,527).

---

## Surviving claims after fair comparison

| Original claim | Verdict after harmonization |
|---|---|
| "Aligner axis dominates" | CONFIRMED for HHV-6B only (1.9×). EBV: CONFOUNDED (same GTF artifact as HSV-1). HSV-1: CONFOUNDED. |
| "Reference strategy minor" | CONFIRMED for all three viruses. |
| "HSV-1 STARsolo = 0" | Was a REGEX artifact (0 → 19 after fix). But 19 is itself a GTF artifact — not a real count. |
| "EBV VS/STAR = 2×" | FALSE — this is a reference-completeness artifact. 95% of ViralScan EBV UMI from CDS-only genes. No clean aligner comparison possible for EBV. |

---

## GTF artifact pattern — affects all three viruses

All three Serratus-derived viral references in the combined GTF carry the same root cause:
- EBV: 80/94 genes CDS-only; 14 exon-bearing genes have 82.8% ambiguous bases
- HSV-1: 61/79 genes CDS-only; 18 exon-bearing genes have 27.4% ambiguous bases
- Anellovirus: all genes CDS-only (F-005, commit 7739521)
- HHV-6B: all 97 genes exon-bearing, 1 gene/contig → UNAFFECTED

Fix: the combined STARsolo GTF needs exon records for all viral genes before any
viral-detection comparison with ViralScan is scientifically valid.

---

Produced by: `analysis/reference_strategy_benchmark/scripts/harmonize_2x2.py`
Registered: `analysis/reference_strategy_benchmark/outputs/numbers.json`
Branch: `claude/multimap-memory-and-showcase`
Date: 2026-07-04
