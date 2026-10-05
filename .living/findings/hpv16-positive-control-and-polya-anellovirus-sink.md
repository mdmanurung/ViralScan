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
poly-G no-signal reads on the covid 5′ library. Here it is poly-A tails. This
is library-dependent, not general to 3′ chemistry: KP343824.1 is also in the max
panel, and the 3′ EBV (10x v2) and HSV-1 (Drop-seq) runs put 0 molecules on it.
Any panel genome with a near-31-nt homopolymer can do this on a susceptible
library. The HPV16 result is unaffected: its molecules are spread across 8
genes in the expected early-gene pattern, and it is ~0 in the isogenic
control.

## Caveats
- The knee cell call is defective (PLAN `SW-23`): it lands at the 10-molecule
  floor. Cut instead at the barcode-rank steepest descent (total ≥952 / ≥739,
  computed outside viralscan), 1,869 / 8,023 HPV16-raft cells are HPV16+
  (23.3 %, holding 7,034 of 7,725 molecules). The normal raft has 0 / 11,045
  HPV16+ cells, and all 8 of its molecules are in empty droplets. These figures
  are indicative only.
- The read scan shows reads that contain the k-mer, not reads that kallisto
  assigned there.
- **Whole-genome test (2026-10-01), prompted by the commensal-virome prior
  (PLAN ANELLO-PRIOR).** In the same 20M reads, every 31-mer of KP343824.1
  (both strands, all 2,794 nt) was checked, with homopolymer-dominated k-mers
  (≥20 of one base) set aside. 27,098 reads hit only homopolymer k-mers; **0
  reads carry any other viral 31-mer**. A real anellovirus transcript would put
  reads on viral sequence, mostly near its 3′ end. So in this cell line the
  call is the poly-A sink, not commensal virus. This says nothing against
  anellovirus in human samples, where it is expected. Limits: 20M of ~191M
  reads, and exact 31-mers only (the same evidence kallisto uses to assign
  them).
- **Reportable cell-level numbers (emptyDrops, 2026-10-01):** the HPV16 raft
  has 2,012 / 16,079 called cells HPV16+ (12.5 %), or 1,927 / 9,086 (21.2 %) on
  the comparable (≥200 host UMI) denominator. The normal raft has 0 / 11,884.
  The anellovirus artefact reaches 3,314 / 11,884 (27.9 %) of the normal raft's
  called cells.
- Only the first 20M of ~191M reads were scanned, and only for KP343824.1's
  head. The 5,688 Gammatorquevirus molecules of the HPV16 raft were not split
  by gene.

## Implications
- Lower the homopolymer mask threshold (for example ≥20 nt) or mask each
  genome's terminal homopolymer, and scan the whole panel for runs of 20–30
  nt. Owner: reference builder; done (CAT-42, see the update below).
- Anellovirus calls on any dataset need a per-gene check (a single
  placeholder gene holding everything is the signature) before being
  reported.

## Update 2026-10-03 (CAT-42 closed)
- Two rebuilds:
  - `viral_ref_cat42`: homopolymer mask ≥20. It removed every
    Gamma/Betatorquevirus molecule (10,037 → 0).
  - `viral_ref_cat42b`: adds the anellovirus-only low-complexity 31-mer mask
    (`--lowcomplexity-kmer-mask`, 1,586 bases in 46 records). It removed the
    residual "genus unassigned" artefact: KP343822.1 220 → 0 and KP343842.1
    26 → 0 in the HPV16 raft; 273 → 0 in the normal raft.
- HPV16 stays at 7,725 molecules. EBV and HSV-1 are unchanged (±1 molecule).
  cat42b is now the current panel.
- The KP343822.1 sink also took 65 molecules in the EBV run and 19 in the
  HSV-1 run. So it was library-independent, unlike KP343824.1.
- Gene-ID rename (124 RefSeq anellovirus genes gain an accession prefix in
  cat42/cat42b). Cause: the `{accession}_{token}` rule (`ncbi_fetch.py:350`,
  3379b7c). The gitignored bundled GTFs were regenerated with it on
  2026-09-28 22:46, after `viral_ref_final` was built.

Tags: hpv16, positive-control, anellovirus, homopolymer, poly-A, reference-artifact
