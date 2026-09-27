# F-012 — The k-mer capture metric measured the forward strand only; 55 genomes were understated by >0.25

**Date**: 2026-09-27
**Status**: fixed the same day; `--strand canonical` is now the default
**Found by**: auditing `scripts/measure_kmer_capture.py` against what kallisto actually indexes

## The defect

`measure_kmer_capture.kmers()` built **forward-strand** k-mer sets:

```python
for i in range(len(seq) - k + 1):
    window = seq[i : i + k]
    if _ACGT.issuperset(window):
        out.add(window)          # no reverse-complement folding
```

kallisto and bustools index **canonical** k-mers — each k-mer folded with its
reverse complement, `min(kmer, revcomp(kmer))`. So a panel genome deposited
antisense to a target scored as a complete miss even though the real index would
match it. The metric was measuring something strictly harsher than the tool it
exists to model, and it was doing so silently.

This matters for Anelloviridae specifically because GenBank deposits of small
circular ssDNA genomes vary in both orientation and linearisation point; nothing
enforces a consistent strand.

## Measured effect

Recomputed across all 2,042 panel anellovirus genomes at k=31, 90 bp reads:

| | forward | canonical |
|---|---|---|
| median leave-one-out | 0.2055 | **0.2162** |
| median p_fragment | 0.5097 | **0.5349** |
| panel distinct 31-mers | 4,888,291 | 4,843,779 |

The aggregate shift is ~5 % relative and would be easy to dismiss. The
per-genome distribution is not:

- 1,692 of 2,042 unchanged
- **350 improve, 0 degrade** (canonical folding is a merge, so this is the
  expected direction and a sanity check that the fix is right)
- 84 improve by > 0.10 absolute p_fragment
- 55 improve by > 0.25

Worst cases — genomes the old metric declared near-undetectable:

| accession | panel label | p_fragment forward → canonical |
|---|---|---|
| MH649023.1 | Anelloviridae | 0.0329 → **0.9514** |
| MH649144.1 | Alphatorquevirus | 0.1539 → 0.9397 |
| MH648977.1 | Anelloviridae | 0.1028 → 0.8741 |
| MH648897.1 | Betatorquevirus | 0.0894 → 0.8171 |
| MH648946.1 | Alphatorquevirus | 0.1376 → 0.8455 |
| KP343847.1 | Betatorquevirus | 0.1988 → 0.8942 |

Per genus, median p_fragment (forward → canonical): Alpha 0.788 → 0.808,
Het 0.837 → 0.837, Beta 0.498 → 0.515, Gamma 0.332 → 0.336, Samek 0.255 →
0.256, Mem 0.247 → 0.247, Gyro 0.063 → 0.063, unclassified 0.334 → 0.397.

The affected accessions cluster hard in the **`MH648xxx`/`MH649xxx`** block — a
single submission batch deposited antisense to the rest of the panel. This is a
batch artifact of GenBank deposition, not biology, which is exactly the kind of
thing a strand-naive metric converts into a fake sensitivity gap.

## Why the median hid it

Nine genomes in ten were unaffected, so every summary statistic moved only
slightly and the bug survived the original 12-test suite and a `--self-check`
that used a self-vs-self positive control. A self-vs-self control cannot catch
this: the genome and the panel copy are the same strand by construction. The
test that catches it is a panel containing **only the reverse complement** of
the target, which must still yield coverage 1.0.

## The fix

`--strand canonical` is now the default; `--strand forward` regenerates the
superseded numbers so nothing already published becomes unreproducible. The
summary JSON records which mode produced it. Tests went 12 → 17, including the
reverse-complement-only panel case in both `--self-check` and the unit suite.

## What this invalidates

Every per-genus capture number recorded before 2026-09-27, including those in
`PLAN.md` `CAT-11`/`CAT-12` and in
[[kmer-capture-independence-model-overstates-detection]] (F-011), is a **lower
bound**. F-011's central correction — that the analytic independence model gives
1.0000 where the empirical walk gives ~0.51 — survives, because the empirical
value moved only from 0.5097 to 0.5349 and the analytic model is still wrong by
roughly a factor of two.

## Open question this raises

`kmers()` also does not wrap the genome origin. Anelloviruses are circular, so
the ~30 k-mers spanning the linearisation point are lost, along with every 90 bp
fragment straddling it. That is ~1 % of a 3 kb genome and probably immaterial,
but it is unmeasured, and it is worth knowing whether kallisto's index shares
the same blind spot — in which case it is a real detection gap to report rather
than a metric artifact to fix.

## Related

- [[anellovirus-capture-ceiling-is-set-by-cdhit-dereplication]] — F-013, which
  uses the corrected metric throughout.
- [[kmer-capture-independence-model-overstates-detection]] — F-011.
