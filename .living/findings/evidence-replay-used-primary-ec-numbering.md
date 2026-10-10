# F-029: the evidence replay selected reads with the primary run's EC numbering

**Date**: 2026-10-07
**Status**: cause confirmed on one run (GSM5725695 `combined_off`); fixed in code (`DSR-17`); affected evidence not yet regenerated
**Found by**: the `combined_off` KSHV evidence for GSM5725695 and GSM5725697 showed mostly host reads, and a rerun (job 25724752) gave different numbers

## Result

`viralscan evidence` replays the FASTQs with `kallisto bus -n`, then ran `bustools capture -e` with the primary
run's `kb-python/matrix.ec`. kallisto numbers equivalence classes in discovery order, which differs between
multithreaded runs, so the replay's records named different classes. For GSM5725695 `combined_off`:

| EC file passed to `bustools capture` | Records captured |
|---|---|
| primary `matrix.ec` (old behaviour) | 80,608 (host reads 77,212, KSHV 1,317) |
| the replay's own `lineage_bus/matrix.ec` | 66,671 |
| rerun, primary `matrix.ec` | 57,353 (host 55,756, KSHV 44) |

66,671 matches the 66,684 target reads in the `twostep_v2` evidence. `matrix.ec` of the primary run, the original
replay and the rerun have equal size but 93,708 to 224,345 differing lines.

## Implication

- Primary runs and `viral_summary.tsv` are unaffected (they use their own consistent EC file).
- Evidence verdicts for low-abundance targets made before 2026-10-07 are suspect: the DSR-02 `host_best` rows for
  HPV77, HPV29 and Cercopithecine herpesvirus, and the evidence behind F-028 (HPV77 in hpv16 SRR19537339).
  High-abundance KSHV in GSM5725696/GSM5725698 looked fine; why (early-discovered ECs have stable ids) is a
  hypothesis, not tested. GSM5725697 (abundant KSHV) was wrong, so the hypothesis is incomplete.
- Next: regenerate affected evidence with the fix (`DSR-17`), then re-judge F-028.
