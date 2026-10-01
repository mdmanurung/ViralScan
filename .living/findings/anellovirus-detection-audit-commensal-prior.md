# F-022 — Anellovirus detection audit: no true positive has ever been tested, and the gene models leave out the UTRs where 10x reads land

**Date:** 2026-10-01 · **Status:** preliminary (code audit + panel measurements; no new sequencing run)
**Tags:** anellovirus, commensal, sensitivity, specificity, UTR, gene-model, eve_risk, MECH-B, ANELLO-PRIOR
**Prompted by:** the user's commensal prior (ANELLO-PRIOR) and Kane et al., *Front. Microbiol.* (published 2026-01-15),
doi:10.3389/fmicb.2025.1716110. Volume not verified.

## What the paper does and does not support
- **The study.** Kane et al. used plasma **DNA** metagenomics: nuclease-treated virions, Phi29 amplification, NovaSeq sequencing.
- **Prevalence.** Anelloviridae was found in 79.1 % of the young cohort (0–16 y) and **100 %** of the older cohort (63–100 y).
- **Genera in the older cohort.** Alphatorquevirus 98.7 %, Betatorquevirus 86.7 %, Gammatorquevirus 83.5 %, Samektorquevirus 22.8 %.
- **Load.** Abundance rises with age.
- **Persistence.** The paper says anelloviruses persist "in leukocyte compartments (T cells and possibly granulocytes)".
- **What it does not show.** It measures viral DNA in plasma, not transcripts in cells. It supports "a true anellovirus call in human blood or tissue is expected and must not be filtered by family". It does not support "most cells, or most scRNA-seq libraries, should hold anellovirus mRNA".

## How ViralScan detects anelloviruses (default path)
1. **Reference.**
   - The panel holds 2,040 anellovirus genomes: 2,022 clareaulab CD-HIT 95 % representatives plus 20 RefSeq (F-013).
   - Genus labels come from `data/anellovirus_accessions.tsv`: Betatorquevirus 1,542, Alphatorquevirus 211, unassigned 185, Gammatorquevirus 77, others 26.
   - Gene models come from `anellovirus.gtf_text_for`. They are the NCBI **CDS** features in `data/anellovirus_genes.tsv` for 1,992 genomes, and a whole-genome placeholder for 48.
   - Homopolymer runs are masked at ≥20 nt. The build with this mask is in progress (CAT-42).
2. **Counting.**
   - kallisto pseudoaligns exact 31-mers.
   - There is **no read-level complexity filter** at count time.
3. **Multimapping (`host-conservative`, the default).**
   - Host/virus-ambiguous molecules give the virus 0.
   - Molecules shared between viral genes are split equally among them (`multimapping.py:540`).
4. **Detection.**
   - `detect_genes` (`detection.py:93`) keeps a gene when its total UMI across all cells is ≥ `detection_threshold`. The default is **1** (`defaults.py:23`).
   - Kept genes are then grouped by genus.
   - There is no anellovirus-specific gate.
5. **Flags.**
   - `eve_risk=True` on every anellovirus row: catalogue `risk_class=eve` on all 3,448 Anelloviridae rows, and `EVE_RISK_GENERA` in `constants.py:408`.
   - The rationale given is "EVEs in NALCN, LINC02742". That came from F-005, whose mechanism F-019 revised: ≥90 % of the signal was poly-G, ≤10 % host-best.
6. **Read-level evidence.** It exists only in the opt-in `viralscan evidence` subcommand: low-complexity, host-competitive and hotspot flags, all "diagnostic_only" (`evidence.py:879`).

## Measurements (2026-10-01)
| check | result |
|---|---|
| Human samples with an anellovirus call that has been read-validated | **0**. Covid x213/x216 = poly-G (F-019). Tonsil = low-complexity (F-010). GSE189670 normal raft = poly-A (F-021). |
| Where the raft artefact sits | In **both** rafts (SRR19537339 normal, SRR19537341 HPV16) every anellovirus molecule is on a `KP3438xx.1_gene1` whole-genome placeholder: KP343824, -47, -22, -42 (both rafts), plus -20 in the normal raft only. These are 5 of the 48 placeholders, the only models that include the non-coding genome ends, where the homopolymers are |
| How the raft artefact is reported | Gammatorquevirus in **27.9 %** of called cells (3,314/11,884). Anelloviridae unassigned 1.3 %, Betatorquevirus 2.1 %. Every row `evidence_tier=candidate_unique` |
| First modelled base, 1,995 CDS-model genomes | median nt **513** (q10 286, q90 576). The 5′ region upstream is in no model |
| TTV-1 NC_002076.2 | AATAAA at nt 3,073. ORF1/ORF2 models end at 2,901/2,875, only spliced gp1 reaches 3,077 |
| Gap from last modelled base to polyA site (signal + 20 nt), 1,922 genomes with an AWTAAA | median ~20 nt, **q90 ~220 nt**, max 598. 70 genomes have no AWTAAA near the model end |
| Per-gene threshold loss (MECH-B) on the raft runs | 0: the artefact molecules were all unique-EC. **Untested on a real, diffuse TTV signal** |
| Simulated sensitivity (F-011/F-012, leave-one-out, 90 bp) | median p_fragment Alpha 0.81, Beta 0.52, Gamma 0.34. **Measured on whole-genome FASTA** (`measure_kmer_capture.py --panel` reads FASTA), so it counts k-mers the current CDS-model index does not hold. It overstates current sensitivity |

## Interpretation
- **Specificity.** The default call has no defence against low-complexity reads. One UMI on one gene is enough to call a genus. F-019 and F-021 show that artefacts are reported as "infected" cell fractions, with nothing in `viral_summary.tsv` to tell them from a real call. Under the commensal prior this is the larger risk: a true call and an artefact look the same.
- **The tonsil null still stands, and the UTR gap does not explain it.** TONSIL-01 (F-010) screened with minimap2 against **whole** anellovirus genomes, UTRs included, accepting mismatches (≥50 bp aligned, ≥85 % identity), independently of the CDS models. It found 0 anellovirus reads in 25M prefiltered unmapped reads. So the best current reading is that **cellular anellovirus mRNA is rare in scRNA-seq libraries**. Plasma DNA prevalence (Kane et al.) does not imply mRNA in sequenced cells. This is one 5′ library from 24 donors, and the strand loss (F-020) does not apply to minimap2.
- **Sensitivity is unmeasured on real reads**, and the models leave out the regions 10x reads come from:
  - 3′ libraries sample the ~100–400 nt before the polyA site. By a first-AWTAAA heuristic, about 10 % of genomes have ≥200 nt of that unmodelled.
  - 5′ libraries sample from the cap. The upstream 5′ region (median 513 nt) is not modelled. TTV PCR assays target a conserved 5′ UTR region, so this probably drops cross-strain k-mers. That is literature lore, **not measured here** (see Action 2).
  - Caveat: circular genomes are linearised at different points, so "nt 1–513 = 5′ UTR" holds for RefSeq-style records, not necessarily for every deposit.
- **Genus balance vs. biology.** Alphatorquevirus is the most prevalent and abundant genus in blood, but it is 211/2,040 genomes. Gammatorquevirus (83.5 % prevalence) is 77 genomes and has the lowest simulated capture (0.34).
- **The `eve_risk` label is questionable for anelloviruses.** The flag fires on every anellovirus call. Its in-repo basis ("EVEs in NALCN, LINC02742") is the F-005 mechanism, which F-019 revised to ≥90 % poly-G and ≤10 % host-best. The EVE rule itself was user-confirmed (decisions 2026-09-30), so any relabel is a proposal for the user. The underlying risks are real: host homology, which the evidence subcommand measures, and low complexity. The label does not describe them.

## Action (proposed; not run)
1. **Positive control first.** Plant reads from held-out panel genomes into a real host library, as TONSIL-01 did. Do it as 3′ and 5′ reads drawn from the transcript ends, not uniformly. If possible, also find a public scRNA or bulk dataset with known high TTV load (transplant or immunosuppressed donors).
2. **Transcript models, measured first.**
   - Rerun `measure_kmer_capture.py` twice: on exon sequences cut from `viral_panel.gtf` (the current index), and on the same plus nt 1–513 and the last-exon→polyA gap. This gives current sensitivity and tests the 5′ UTR claim.
   - If the gain is real, extend models **to the mRNA extent only (cap to polyA site), never to the whole genome**. Whole-genome placeholders are exactly where every raft artefact landed, and the post-polyA GC-rich region (59 % GC on TTV-1) is non-transcribed.
   - Judge any extension on the negative controls as well (SRR19537339 normal raft, covid x213), not only on the plant.
3. **A default read-level gate (ANELLO-PRIOR).** Report the fraction of each virus's molecules with a ≥15-nt homopolymer or low entropy, and where reads land along the genome. Gate on that, never on family.
4. **Propose to the user:** rename or split `eve_risk` for Anelloviridae, for example into a "host-homology/low-complexity risk" note.
5. **Close MECH-B** (group → sum → threshold) before any anellovirus sensitivity claim.
