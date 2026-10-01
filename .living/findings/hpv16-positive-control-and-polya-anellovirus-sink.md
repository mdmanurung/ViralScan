# HPV16 positive control passes; a poly-T reference head sinks poly-A reads into an anellovirus (F-021)

Status: preliminary · Date: 2026-10-01 · PLAN: EXPL-HPV16, CAT-41

## Claim
On GSE189670 (isogenic NIKS keratinocyte rafts, 10x 3′ v3, NovaSeq), the
shipped 2,343-genome panel detects HPV16 in the infected raft (SRR19537341:
7,725 molecules, 2,442 called cells) and almost none in the parental raft
(SRR19537339: 8 molecules, 6 cells), about 970:1. In the same anellovirus-free
cell line it also calls Gammatorquevirus at 5,688 / 10,037 molecules. In the
normal raft all 10,037 sit on one whole-genome placeholder gene,
`KP343824.1_gene1`.

## Evidence
- Run: job 25684772, repo HEAD f40f72a, `-x 10xv3`, index
  `viral_ref_final/build/panel.idx`; outputs in
  `viralscan_work/hpv16_gse189670/runs/<SRR>/<SRR>/results/`.
- HPV16 per gene (unique + allocated), from `multimap_evidence.tsv`; names
  come from the bundled GTF `gene` attribute: E5 2,374; E1^E4 1,501; E7
  1,440; E1 1,100 (all allocated); E2 696; E6 451; L1 112; L2 51. Early
  genes dominate, as the authors report qualitatively (PMID 37031202). The
  genes overlap and 3′ reads pile up at polyA sites, so per-gene counts
  reflect polyA position.
- KP343824.1 ("UNVERIFIED: Torque teno virus isolate S57, complete genome",
  2,794 nt) starts with a 29-nt poly-T run. The builder masks only runs of
  ≥31 nt, so this one stays in the index.
- Read check: of the first 20M R2 reads of SRR19537339, **27,098 (0.14 %)**
  contain a 31-mer from the first 60 nt of KP343824.1 or its reverse
  complement. Examples are host mRNA 3′ ends running into poly-A tails, plus
  TSO (`AAGCAGTGGTATCAACGCAGAGTAC`) and Read-1 adapter reads. Only 2,094 of
  them are ≥80 % A/T over the whole read, so most are ordinary transcripts
  with a poly-A tail.

## Interpretation
This is the same class of failure as F-005/F-019: a homopolymer tract in a
reference genome becomes a single-gene sink. The feed differs. F-019 was
poly-G no-signal reads on the covid 5′ library. Here it is poly-A tails, which
every 3′ library contains. Any panel genome with a near-31-nt homopolymer can
do this. The HPV16 result is unaffected: its molecules are spread across 8
genes in the expected early-gene pattern, and it is ~0 in the isogenic
control.

## Caveats
- The knee cell call is not usable on these libraries (total ≥10 gives
  60,475 / 1.57M barcodes), so per-cell percentages are not reported.
- Only the first 20M of ~191M reads were scanned, and only for KP343824.1's
  head. The 5,688 Gammatorquevirus molecules of the HPV16 raft were not split
  by gene.

## Implications
- Lower the homopolymer mask threshold (for example ≥20 nt) or mask each
  genome's terminal homopolymer, and scan the whole panel for runs of 20–30
  nt. Owner: reference builder; not changed yet.
- The anellovirus calls on any 3′ dataset need this check before being
  reported.

Tags: hpv16, positive-control, anellovirus, homopolymer, poly-A, reference-artifact
