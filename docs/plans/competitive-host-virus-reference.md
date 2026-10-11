# Competitive host + virus reference: design note (PLAN AUD-17, finding M-05)

Status: design only, written 2026-10-11. No code is proposed here that has been written, and nothing
below has been run. Claims about the current code cite the file they come from; claims about
performance are marked "to be measured".

Source of the finding: Biomni review, `analysis/audit_2026-10/findings.tsv` row `M-05` ("viral-only
matching cannot show uniqueness vs omitted host; host-first subtraction discards ambiguous viral
evidence; need explicit unresolved category"). Raised by Biomni, not independently audited. User
decision 2026-10-11 (PLAN WP4L `AUD-17`): record a design; M-02 (polyA capture columns) is
implemented separately.

## 1. Problem

A read set is quantified against some reference. Whether a viral call means "this read is viral"
depends on what the reference omits.

1. **Viral-only indices cannot show uniqueness.** A read that pseudoaligns to a viral gene in an
   index that holds no host sequence is "unique" only with respect to that index. If the read is
   also host-derived (EVE or repeat homology, low-complexity or TSO-led reads, ambient host RNA),
   the index cannot say so. `docs/v3_counting_contract.md` already states the boundary: `index_kind`
   is `host_only`, `combined` or `virus_only`, and a virus-only matrix "cannot provide host depth".
   The `F-028` HPV77 false positives (PLAN `DLIST-01`) are the concrete case: the match needs
   TSO-junction k-mers that a D-list does not mask.
2. **Host-first subtraction discards ambiguous viral evidence.** `--host-filter starsolo`
   (`scripts/host_filter.py`, `cli/parsers.py`) removes read pairs that align to the host genome
   before `kb count`. A removed pair is gone: the retained BUS stream cannot show that it was also
   viral-compatible. The counting contract says this explicitly ("SW-07's prefilter host-virus
   fragment boundary therefore remains open"). Genuine viral reads that share sequence with the host
   are therefore lost silently, and the loss is not reported anywhere.
3. **There is no explicit "unresolved between host and virus" outcome.** In a combined index,
   `multimapping.py` already counts host-virus ambiguous molecules
   (`host_virus_ambiguous_molecules`, `counts_host_viral_ambiguous`, the `host-conservative`
   allocation). That is a diagnostic layer. It does not reach the per-virus summary as a category a
   reader sees next to "detected" and "not detected". The word "unresolved" in the counting contract
   means something else (unresolved UMI collisions), so a new category needs a different name.

What already exists and what is missing, so the scope is not overstated:

| Capability | State in the repository |
|---|---|
| Combined host cDNA + virus `kb ref` index | Default of `viralscan build-ref` (`scripts/build_reference.py`, `build_combined_reference`). Host is **transcriptome only** (Ensembl cDNA, a cDNA-level GTF), no genome, no introns, no repeats. |
| Host genome as D-list | Opt-in `--genome-dlist` (passed to `kb ref --d-list`) plus `measure_host_homology` (minimap2) annotations written into the reference manifest. A D-list only stops extension at distinguishing flanks; it does not turn host genome sequence into countable targets (PLAN `DLIST-01`). |
| Host-virus ambiguity counting | `multimapping.py` layers `host_viral_ambiguous*`, audit count `host_virus_ambiguous_molecules`. Combined-index only. |
| Host subtraction | `--host-filter starsolo`, destructive, needs a STAR genome. |
| Competitive evidence alignment | `--host-fasta`, `--competitor-fasta` + `--competitor-manifest` for evidence alignment/BLAST (`evidence.py`), after calling, not during counting. |
| Truth data for host homology | `scripts/generate_truth_panel.py` dataset `synthetic_host_homology` with `host_homology_manifest.tsv` and `repeat_homology` planted reads. |

So the gap is narrower than "build a competitive index": (a) host sequence that is not in cDNA is
invisible to a combined index, (b) virus-only and host-filtered runs have no competition at all, and
(c) none of this surfaces as a reported category.

## 2. Design options

### Option A. Single competitive `kb ref` index (host cDNA + decoys + virus)

One index holds viral genes, host transcripts, and optionally extra host decoy targets. Every read is
pseudoaligned once; a read compatible with both host and viral targets lands in a mixed equivalence
class and is counted by the existing ambiguity machinery.

Variants, in increasing cost:

- **A1 (exists):** host cDNA + virus, `--genome-dlist` for the host genome. No new build code.
- **A2:** A1 plus host *non-transcript* decoys as countable targets: ncRNA, and the host
  genomic loci already identified by `measure_host_homology` (EVE/repeat hits), not the whole genome.
  The decoy FASTA entries are listed in the manifest as role `host_decoy`, so they never enter the
  viral identity table as viruses.
- **A3:** the whole host genome as extra targets. Rejected for now: index size and `kb ref` memory
  are the dominant risk (section 6) and a whole-genome target set turns every intronic or nascent
  host read into a "host" count, which changes the meaning of host depth.

Strength: one pass, one set of counts, no pair is discarded, ambiguity is observed per molecule.
Weakness: only as competitive as the decoy set; the k-mer index treats a decoy and a target
symmetrically, so a TSO-junction k-mer present in both host decoys and a viral CAG tract still
produces an ambiguous class rather than a clean call (the `DLIST-01` finding applies to D-lists; a
decoy *target* behaves differently, which is exactly what the validation plan must test).

### Option B. Two passes with an explicit `unresolved` outcome

Pass 1 quantifies against the viral index (or the combined index). Pass 2 re-examines the molecules
behind each viral call against a host reference and tags them instead of deleting them.

- Pass 1 is today's flow without `--host-filter`.
- Pass 2 either (i) re-pseudoaligns the extracted viral reads (`viralscan.reads.extract_virus_reads`
  already replays the exact molecules behind a call, PLAN `READS-01`) to a host-only index, or
  (ii) aligns them with STAR to the host genome, as `--host-filter` does, but writes the result as a
  per-molecule tag rather than a filter.
- Each viral molecule becomes `viral_unique`, `host_compatible` (also hits host), or
  `host_only` (host hit and no viral hit in pass 2, which signals a pass-1 index artefact).

Strength: reuses existing pieces, works with a viral-only prebuilt index, and gives the same
outcome vocabulary for users who need `--host-filter`. Weakness: a second alignment, STAR genome
needed for (ii), and the per-molecule tag must be threaded through the BUS-derived matrices, which
today carry no per-molecule host annotation. It also cannot replace the competition inside the
index: a molecule that pass 1 assigned to the wrong viral gene because host sequence was absent is
only caught if it still reaches pass 2 at all.

### Why not "just filter harder"

Tightening `--host-filter-star-params` trades host removal against viral-read loss and is already
parameterised (`pinned` vs `star-default`). It cannot make the loss visible, which is the actual
finding.

## 3. Required manifest and identity-table changes

Both options need the reference to say what kind of competition it provides. Proposed fields (names
are proposals, to be fixed when implemented):

`reference_manifest.json` (`schemas/v3/reference_manifest.schema.json`; a schema bump and a governance
re-pin, as for any schema change):

- per-sequence `role`: `virus` / `host_transcript` / `host_decoy`, replacing the implicit
  "everything not viral is host";
- `competition`: `none` / `host_cdna` / `host_cdna_plus_decoys` with the decoy FASTA sha256, source
  (genome release, `measure_host_homology` threshold) and count;
- the existing `genome_dlist` and `annotations` blocks stay; add `host_decoy_sha256` next to them.

`virus_identity.tsv` (`virus_identity.py`, one row per indexed gene_id):

- a `role` column that distinguishes `host` from `host_decoy`, so decoys are excluded from host depth
  and from every viral denominator by construction. Today rule 3 sends "everything else" to `host`;
  a decoy must not inherit host-gene semantics (for example in `host-conservative` allocation).
- `index_kind` gains a value such as `combined_decoy`, or is kept as is with `competition` carried
  in `uns` (the second is smaller; decide with the schema owner).

Run manifest (`run_manifest.schema.json`): record the reference `competition` value and, for option
B, the pass-2 host reference identity, so two runs with different competition cannot be compared
silently. `pct_infected_called` is already documented as not comparable across host-filtering
settings (PLAN, around the STAR host-filter rows); the same note applies here.

## 4. Outputs: an `unresolved` category

Per virus, `viral_summary.tsv` gains counts that partition the viral-compatible molecules:

- `n_molecules_viral_unique`: compatible with viral genes only;
- `n_molecules_host_viral_unresolved`: compatible with at least one viral and at least one host or
  decoy target (the existing `host_virus_ambiguous_molecules`, reported per virus instead of one
  total in the audit);
- under this proposed competitive-reference design, `detected` would use unique evidence only
  (the current pipeline supports other molecule-accounting modes). A virus whose only evidence is
  unresolved reads `detected = False` with a non-empty `unresolved` count, and the report says
  "unresolved", not "absent". This keeps `AUD-07` (infected rule) and `AUD-21` (programme rows only
  where `detected`) intact: nothing downstream starts treating unresolved mass as detection.

Name: `host_viral_unresolved`, deliberately not `unresolved`, to avoid the counting-contract term
for UMI collisions. Cell level: a per-cell count column in the same layer family as
`counts_host_viral_ambiguous`, not a new matrix. Documentation rows for `docs/output_reference.md`
(pinned, so a separate commit by the docs owner) follow at implementation.

## 5. Validation plan

Truth source: `scripts/generate_truth_panel.py`, dataset `synthetic_host_homology`
(`challenge_reads`, `repeat_homology`, truth in `host_homology_manifest.tsv` whose columns come from
`v3_score.required_columns`), plus `synthetic_host_only` as the negative control and the mixed
host-virus dataset for the positive control. The generator is deterministic and seeded; its output
is never committed (governance forbids tracked FASTQ).

Arms, all on the same planted FASTQ, one tiny panel:

| Arm | Reference | Host handling |
|---|---|---|
| V0 | virus-only | none |
| V1 | virus-only | `--host-filter starsolo` |
| C1 | combined host cDNA + virus | none |
| C2 | C1 + `--genome-dlist` | D-list |
| C3 | C1 + host decoys (option A2) | decoy targets |

Measured per arm, against planted truth: false-positive viral molecules from host-homology reads;
viral molecules lost to host handling (true viral reads present in the truth manifest but absent
from the counts); the size of `host_viral_unresolved`; recall on the mixed positive control.
Acceptance (to be fixed before running, not after): C3 or C2 must reduce host-homology false
positives versus V0 without losing more than a stated fraction of planted true viral molecules, and
every lost-or-demoted molecule must appear in the unresolved count rather than vanish (conservation:
unique + unresolved + discarded-by-design = planted viral molecules).

Known limit: `synthetic_host_homology` plants homology by construction. It shows the mechanism works
and bounds the cost; it does not show the false-positive rate on real data. That needs the real
control set and is outside this note (user-owned runs).

## 6. Risks

- **Index size, build time and `kb ref` memory.** All unmeasured here. `_run_kb_ref` already records
  resources, so the first experiment records index size, wall time and peak RSS for C1/C2/C3. A3
  (whole genome) is the case most likely to be infeasible on ordinary hardware, which is why it is
  not proposed.
- **Interaction with `--host-filter`.** Combining a competitive index with STAR subtraction removes
  exactly the pairs that would have become `host_viral_unresolved`. The two should be mutually
  documented: with a competitive index, `--host-filter` is not needed (its own help says so), and if
  both are used the unresolved count is a lower bound. A run manifest field must say which.
- **Interaction with DLIST-01/02.** `DLIST-01` found a D-list cannot mask the TSO junction, so
  decoy targets (option A2) must not be assumed to fix `F-028` either; `tso-trim` (`TSO-01`) remains
  the fix for that artefact and C3 is measured with and without it. `DLIST-02`
  (`bustools capture --complement` host-capture split, per-virus `reads_virus_only` /
  `reads_also_host`) is the same idea as option B(i) at BUS level; if option B is chosen, `DLIST-02`
  is its first step rather than a separate feature.
- **Host depth semantics.** Decoys must be excluded from host depth and from comparable-cell
  denominators (section 3), or `AUD`-era denominators silently change.
- **Equivalence-class growth.** More competing targets create more mixed classes; the `host-
  conservative` allocation then sends more mass to host. Whether that is acceptable is a measured
  question, not assumed.
- **Fingerprint and resume.** A reference `competition` field changes the reference manifest, so
  `--resume` behaviour for old runs needs the same caveat as the `C-13` fingerprint change.
- **Scope creep into hosts other than human/mouse.** `constants.ENSEMBL_SPECIES` limits hosts;
  decoy construction depends on `measure_host_homology` (minimap2) and on a host genome FASTA per
  species.

## 7. Recommendation

Choose **Option A (single competitive index) as the target design, delivered in the order A1, then
A2**, and make the **reporting category the first deliverable**, because it is valuable even with
today's indices and requires no new reference construction. Do not build A3. Keep Option B(i)
(`DLIST-02`) as the route for users who must stay on a virus-only prebuilt index, and do not extend
the STAR `--host-filter` beyond tagging.

Smallest first implementation step (not part of this task):

1. Report `host_viral_unresolved` per virus in `viral_summary.tsv` from the host-virus ambiguous
   molecules that `multimapping.py` already computes for combined indices, and for a virus-only
   index write the column as `not_applicable` (never 0, in line with `n_cells_reactivating`).
   No new build code, no schema change to the reference, one new documented output column.
2. Then run arms V0, V1, C1, C2 of the validation plan on one planted `synthetic_host_homology`
   panel to see whether the existing combined index plus D-list already meets the acceptance
   criteria, recording size/time/memory. Only if it does not, spend effort on A2 decoys and the
   manifest `role` change.

This ordering stops the work as early as the data allow: if C1/C2 already contain the false positives
and the unresolved count is informative, A2 is not needed.

## Documentation rows implied

For the docs owner (pinned files, not edited here): `docs/output_reference.md` `viral_summary.tsv`
row for `host_viral_unresolved` when implemented; `docs/cli_reference.md` only if a flag is added
(none proposed). The M-02 columns `polya_capture` and `capture_basis` are a separate row set.
