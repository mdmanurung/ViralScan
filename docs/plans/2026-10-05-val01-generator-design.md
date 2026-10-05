# VAL-01 truth-panel generator — design draft

**Status:** draft for review, 2026-10-05. Nothing is implemented.
**Covers:** `VAL-01` (generator), the generator side of `VAL-02`/`VAL-03`, and the factor levels that `SCI-03` needs.
**Contracts it must satisfy:**
- `analysis/v3_validation/protocol.yaml`, in particular:
  - `partitions` (`generated_sample_structure`, `sibling_virus_pairs`, the split algorithm, `leakage_prohibitions`);
  - `calibration.limit_of_detection`;
  - `datasets[synthetic_*]`;
  - `factors`;
  - `seeds`;
  - the four truth-manifest column lists in `audit_artifacts`.
- WP1E decisions R2.3, R2.9, R3.4, Q5 and Q10d (`PLAN.md:1682`).

Section 6 lists every place where the protocol, as written, contradicts itself or cannot be executed. Those items are inputs to the `DEF-00` amendment. This draft does not settle any of them quietly.

---

## 1. What the generator produces

The generator produces one directory per biological sample, as `file_layout` specifies:

- `{biological_sample_id}_R1.fastq.gz` and `{biological_sample_id}_R2.fastq.gz`;
- the four truth manifests (`truth_manifest`, `host_only_manifest`, `host_homology_manifest`, `mixed_host_virus_manifest`), each restricted to the datasets it belongs to;
- `barcode_pool.tsv`: the barcodes drawn for this sample, with a role of `cell` or `empty`. Drop-seq has no on-list, so this file is the only record of its barcode domain. It also feeds `cell_universe_manifest`;
- `run_manifest.json`, which records:
  - the generator version (git SHA plus `viralscan.__version__`);
  - the Python and numpy versions;
  - every input path with its sha256;
  - the stage seeds;
  - the factor levels;
  - the per-file sha256 of every output.

Datasets, all built by one generator with different switches:

| dataset id | viral plant | host source | extra population |
|---|---|---|---|
| `synthetic_factorial` | six targets × abundance grid | real healthy-donor PBMC, plus an unplanted twin (D1) | ambient, collisions, artefacts |
| `synthetic_host_only` | none | synthetic GRCh38 transcripts | artefacts, empty droplets |
| `synthetic_mixed_host_virus` | as factorial | as factorial | 10 % of viral molecules share a key with a host fragment |
| `synthetic_host_homology` | none | synthetic GRCh38 transcripts | reads from frozen GRCh38 loci homologous to the panel (REF-06/REF-08) |

## 2. Generation order (pinned)

1. **Enumerate conditions.** Take the cross-product of the stratification factors (`viral_abundance` × `host_virus_homology` × `chemistry`), with `minimum_samples_per_stratum` = 4 replicates per stratum. Recompute the sample floor, 4 × the stratum count, from the frozen level counts and refuse to run below it.
2. **Assign opaque sample IDs**: `s` + the first 12 hex characters of `sha256(seeds.split|dataset|stratum_key|replicate)`. This is where `seeds.split` is consumed. The largest-remainder split is otherwise deterministic, so a frozen split seed with no call site would be decoration. The split takes holdout members in ASCII order within each stratum, so an ID must carry no factor or replicate information. A readable ID such as `abund03_rep1` would put replicate 1 in holdout every time.
3. **Run the frozen split**: global largest-remainder apportionment with `holdout_fraction` 0.3, seeded by `seeds.split`.
   - No implementation exists yet; a grep of `scripts/` and `src/` for apportion/largest-remainder found nothing.
   - The split itself consumes no randomness: `seeds.split` enters only through the IDs in step 2.
   - Write it once, in a small `scripts/v3_partition.py` that the VAL-06 scorer imports too. If generator and scorer each had their own copy, two correct-looking copies could choose different holdout members.
4. **Assign covariates** (the six non-stratified factors) by a seeded permutation within each stratum, using stage seed `covariates` derived as in step 6 at the stratum level, and independent of sample-ID order. The levels cycle so that each stratum carries every level at least once, which needs ≥ 3 levels and ≥ 4 samples.
   - Drawing covariates after the split is deliberate: a covariate cannot then influence partition membership.
   - Drawing them independently of ID order matters just as much. Covariates tied to ID order would be confounded with the split.
5. **Partition the barcode pools.** The leakage rule ("no … UMI, or cell barcode is shared across partitions") is met by splitting each chemistry's on-list, or Drop-seq's random 12-mer pool, into disjoint training and holdout halves using a hash of the barcode. A sample draws only from its own partition's half.
6. **Derive stage seeds**: `int(sha256(protocol_version|dataset|sample|replicate|condition|"generator"|stage)[:16], 16)`, following `seeds.derivation`. The generator uses one numpy `Generator(PCG64)` per stage.
7. **Generate.** Write samples in sorted sample-ID order and molecules in a deterministic order.

## 3. Read and molecule model

**Geometry.** Barcode, UMI and read lengths come from `viralscan.chemistry.get()`, never from constants. The planter hard-codes the 10x v3 values (`CB_LEN, UMI_LEN, READ_LEN`).

- R1 is `CB + UMI`, plus a poly-T stub when the chemistry has one.
- R2 is 90 nt for every chemistry. One R2 length means the `chemistry` factor measures barcode geometry only. With per-chemistry R2 lengths (v2 98, v3 91, Drop-seq ~60), read length would be confounded with chemistry, and the LoD would change for k-mer reasons.
- Naming: the protocol says `drop-seq` and the package says `dropseq`. A single mapping table translates between them, and the manifests use the protocol's spelling.

**Molecule model, per cell:**
- host molecules: `reads_per_cell / mean_family_size`;
- transcripts drawn from a seeded log-normal expression vector over Ensembl `ENST*` records. This reuses `read_fasta(..., HOST_RECORD_PREFIX)` from the planter, the guard that keeps viral cDNA out of the host;
- there is no cell-type structure, because the cell anchor is the planted membership (`cell_universe.synthetic_rule`), not a clustering.

**Fragment model.** This is a 3′ chemistry. The read start is drawn uniformly in the 300 nt upstream of the cleavage site, and a read that reaches the site gets an untemplated poly-A tail. This is `polya_read()`, generalised over `read_len`.

- **Viral templates** are the frozen reference's viral gene models, so every planted molecule's `planted_gene_id` is in the common feature set (`eligibility_rule`).
- The cleavage site is the gene end, the same estimate as the planter's `cleavage_site()`. This keeps the planter's lesson: a clean genome slice would make any poly-A gate look free.

**Sequencing error.** One substitution model is applied to every read, host and viral alike, at 0.1 %/base, including the CB and UMI bases of R1. A 1-Hamming CB or UMI error therefore exists at every factor level, and the truth always records the source CB and UMI.

- This replaces the planter's 1 % viral rate. That rate stood in for strain divergence, which is not a factor here; siblings and homology carry divergence explicitly.
- Mutation is vectorised in numpy. The planter's per-base Python loop is too slow at the scale in section 5.

**FASTQ determinism.** Write through `gzip.GzipFile(fileobj=..., mode="wb", mtime=0, compresslevel=6)`.

- The planter's `gzip.open` stamps the current time into the header, so two identical runs give different sha256 values, which breaks VAL-02 hashing and VAL-08 determinism.
- The manifest also records the sha256 of the uncompressed stream, so a different zlib build cannot make an identical panel look different.

## 4. Factors — proposed levels

The table covers all nine factors. Two stratify (`viral_abundance` and `host_virus_homology`); `chemistry` is already frozen. The other six are covariates.

| factor | proposed levels | construction | truth columns | consumer |
|---|---|---|---|---|
| `viral_abundance` (stratifies) | **total planted molecules per target virus per sample**: 1, 3, 10, 30, 100, 300, 1000 (7 levels, half-log10 steps; see D3) | multinomial over that virus's designated infected cells | `viral_abundance_level`, `planted_molecule_id`, `is_planted_molecule` | E8 / D17 probit; LoD95 |
| `host_virus_homology` (stratifies) | **owned by REF-08, not VAL-01**; the tiny panel uses `none` | reads drawn from the frozen GRCh38 homologous-locus table | `host_virus_homology_level`, `host_source_locus`, `homology_identity`, `planted_as_host` | D9, D18, H2/H5, R2.3 constraint |
| `chemistry` (stratifies) | 10xv2, 10xv3, drop-seq (frozen); 5′ pending Q5/DEF-00 | `chemistry.get()` | `chemistry` | stratum, D17 per chemistry |
| `infected_cell_fraction` | per virus: 0.01, 0.04, 0.15 (20 / 80 / 300 of 2000 cells). Six disjoint sets at 0.15 use 1,800 cells; `co_infected` pairs share theirs | designated per virus; a cell is `is_infected_cell` only if it receives ≥ 1 planted molecule (see D3) | `infected_cell_fraction_level`, `is_infected_cell` | D12, D18, D19 |
| `sibling_virus_similarity` | `disjoint_cells`, `co_infected` (see D2) | pair members planted into non-overlapping or identical cell sets; per-pair ANI measured and recorded, not varied | `sibling_virus_similarity_level`, `planted_virus_id` | D13, D24, E5 |
| `low_complexity` | 0, 1 %, 5 % of library reads | equal parts of the planter's `artefact_read()` classes: polyG, polyA, CAG, TSO+polyA, with tract ≥ 40 nt; viral 3′ poly-A tails are always on and are not this factor | `low_complexity_level`; artefact reads get a manifest row with `origin = artefact:<class>`, an empty `planted_virus_id` and `is_planted_molecule = false`, so no scorer can count one as viral truth | R2.0 read filter, D9 |
| `pcr_duplication` | mean reads per molecule of 2, 5, 10 (1 + geometric) | each duplicate is an independent read with independent errors | `pcr_duplication_level`, one row per read | E1 conservation, UMI collapse |
| `cb_umi_collisions` | 0, 0.5 %, 2 % of viral molecules, half of them **same** and half **disjoint** | *same*: the exact CB+UMI of another molecule in the same cell, with a different gene; *disjoint*: a UMI 1 Hamming from another molecule's UMI in the same cell | `cb_umi_collision_level`; extra column `collision_group` | `collision_rule` (one corrected key → ≤ 1 truth molecule), E1 |
| `ambient_index_hopping` | 0, 1 %, 5 % of each virus's molecules | *ambient*: relocated to uninfected cells and empty droplets in the same sample; *hopped*: copied into another sample of the **same partition and chemistry**, under that sample's CB pool | `ambient_index_hopping_level`; extra column `origin` ∈ {cell, ambient, hopped, artefact:<class>}; `is_infected_cell = false` | D11, D18 |

Extra columns (`collision_group`, `origin`) are appended after the required ones. The validator checks only that the required columns are present. Check this during implementation.

**Empty droplets.** The protocol requires them (`synthetic_rule`, D11) but declares no count or depth. Proposal: 3,000 empty barcodes per sample, carrying ambient-only content at a log-normal median of 50 UMI. These values are generator parameters, so they must be frozen with the levels.

**Distribution over genes.** A virus's molecules are spread uniformly over its gene models. Many herpesvirus transcripts share 3′ ends, so a fragment can be compatible with more than one gene. For each planted fragment, the generator intersects the fragment interval with every gene model of that virus. It records the result as `truth_class`: `unique` for one compatible gene, `ambiguous` for more than one. The scorer then does not have to infer it, and D13 and the truth classes rest on a recorded fact. Uniform spreading is the simplest choice. Expression-weighted spreading (latent versus lytic) is a realism upgrade; add it only if review asks for it.

**Six targets in every sample.** E8 requires "every expected target at every abundance level within each chemistry". Each sample therefore plants all six targets at its abundance level, each into its own designated cell set. The sibling level decides whether pair members share cells.

## 5. Scale (arithmetic; assumptions labelled)

| quantity | value |
|---|---|
| pairs per sample | 2,000 cells × 25,000 reads = **50 M** (+ empty droplets) |
| `synthetic_factorial` samples | 4 × 7 abundance × *h* homology × 3 chemistry = **84 *h*** (252 at *h* = 3) |
| pairs, factorial alone | 12.6 B at *h* = 3; the other three synthetic datasets are roughly the same again → **~40 B** |
| FASTQ storage | ~75 B/pair gzipped (assumed: 28 + 90 nt records) → **~3.8 GB/sample, ~1 TB factorial, ~3 TB all four** |
| generation CPU | ~3–5 min per sample-core with vectorised numpy and gzip level 6 (estimate, unmeasured) → ~20 core-h per dataset; small |
| downstream CPU, the real cost | kb count ~2 core-h per 50 M reads; STAR two-step ~8 core-h (estimates). One ViralScan config × 252 samples ≈ 500 core-h; STAR ≈ 2,000 core-h |
| budgets | tuning 20 k core-h (Q10d), comparators 10 k core-h (R3.3) |
| truth manifest, as specified | one row per read → 12.6 B rows ≈ **2.5 TB TSV** for the factorial alone. Not feasible |

Generation is cheap; running the pipeline over the panel is not. At 50 M pairs per sample, a tuning grid of about 10 ViralScan configurations plus one STAR arm uses roughly 25–35 % of the 20 k budget on the factorial alone, before the other datasets. The options are in D5.

## 6. Decisions needed (inputs to DEF-00)

Each decision lists options with a recommendation first. None of them is settled here.

**D1. R2.9 background versus sample independence and exact truth.**
R2.9 locks a *real PBMC primary* background plus a synthetic GRCh38 control. A real background raises three problems:

- **Partition leakage.** Its barcodes and molecules are fixed, so the background cells themselves must be split between partitions: of the PBMC 10k, about 7 k cells go to training and 3 k to holdout.
- **Independence.** Every factorial sample would then reuse cells from that one pool. Samples are no longer independent biological samples, and the sample bootstrap (the `uncertainty` unit) would overstate precision.
- **No exact truth.** The background is presumed-negative, not exact truth, so the design needs planted-minus-unplanted pairs, as in `plant_anello_10x.py`. That doubles the runs.

Only the v3 PBMC exists today. v2 and Drop-seq backgrounds still need downloading (approved under Q10c).
- The only v3 PBMC on hand is also `pbmc_10x_healthy_v3`, the external-evaluation presumed-negative dataset (VAL-04). Planting into it for training would contaminate that evaluation.
**Decided (user, 2026-10-05): real PBMC is the primary background, so R2.9 stands.** Option (a) is rejected. The synthetic GRCh38 background stays as the control: `synthetic_host_only`, the host-homology negatives, and the host side of the mixed-molecule dataset.

Consequences the generator must implement:
- **Planted minus unplanted.** Every factorial sample has an unplanted twin: the same background reads with nothing added. The truth for real background reads is "unknown", never "host". Healthy PBMC can carry genuine EBV in B cells and genuine anellovirus, so a call in the unplanted twin is background signal and not a false positive. Recovery is scored on planted molecules only, specificity on the twin difference. The twins double the run count; this feeds into D5.
- **Background cells are the outcome-independent called cells of each library.** The anchor is the external called-cell list, else host-only EmptyDrops (`cell_universe.source_precedence`). Planted molecules go into those real CBs. The empty droplets are real, uncalled barcodes, which replaces the synthetic proposal in §4 for this dataset.
- **Partition by donor library, not by cell.** Each background library falls wholly in training or wholly in holdout. Canonical barcode keys are library-scoped (`library-id::CB`), so two libraries sharing an on-list sequence does not breach the leakage rule.
- **Depth comes from the background.** `reads_per_cell` is reached by subsampling each library, and a library shallower than the target is recorded as shallower, never topped up with synthetic reads.
- **`pbmc_10x_healthy_v3` is excluded as a background**, because it is the VAL-04 external-evaluation negative.
- **Downloads needed** (approved under Q10c): healthy-donor PBMC libraries for 10xv2, 10xv3 and Drop-seq, each checksum-pinned. They must not be already-seen data (R2.2).

**D1.1, the independence unit: decided (user, 2026-10-05) as option (b).** Samples drawn from one library are not independent biological samples.
- (b) **Chosen:** several donor libraries per chemistry, with the bootstrap resampling *donor libraries* (a cluster bootstrap) and generated samples nested within donors. Amend `partitions.unit` and `uncertainty.resampling_unit` to read "background donor library". The 0.3 holdout and `minimum_samples_for_interval` = 3 then need ≥ 10 donors per chemistry, which leaves 3 in holdout. That is about 30 public libraries in total.
- (c) Rejected: one library per chemistry, with the bootstrap unit declared to be that library. Every holdout interval is then not-estimable (one unit), so the LoD and recall claims carry point estimates only.

**D2. Sibling "absent" versus "every target at every level".**
E8 requires all six targets planted at every level in every chemistry. `sibling_pair_rationale` requires both pair members in the same condition. A sibling is therefore never absent, yet the factor text promises "absent/present sibling status".
**Decided (user, 2026-10-05): option (b).** The levels are `disjoint_cells`, `co_infected` and `one_member_only`.
- Construction constraints:
  - `one_member_only` falls on exactly one sample per stratum (the covariate cycling in step 4);
  - the dropped member of each pair is a seeded draw per pair, recorded as `sibling_virus_similarity_level` plus an absent `planted_virus_id`;
  - every virus therefore keeps ≥ 3 of the 4 limit-of-detection rows per stratum, and the E8 "every target at every level within each chemistry" requirement still holds.
- **DEF-00 must add the consumer.** Under the current `sibling_pair_rationale`, a condition contributes to D13/D24 only when *both* members are planted, so `one_member_only` samples would be generated and then scored by nothing. That is the frozen-value-with-no-call-site failure. The amendment adds an absent-sibling false-call metric: calls for the unplanted member in a `one_member_only` sample, denominated on that sample's anchor cells.
- (a) Rejected: the levels are `disjoint_cells` and `co_infected`. "Absent" is dropped from the factor text. Truth is per molecule, so sibling confusion is still measured exactly, and `disjoint_cells` gives cell-level absent-sibling truth: a cell planted with HHV-6A only, and an HHV-6B call in it, is a confusion. **Cost:** there is no sample-level absent-sibling test, meaning a sample in which the partner is truly absent. The host-only and homology negatives cannot stand in for one, because they contain no virus to bleed from.
- (b) Add `one_member_only`, where a seeded choice of one member per pair is planted. This loses LoD rows for the dropped member, so E8 needs ≥ 2 rows per virus × level × chemistry to stay identifiable.

**D3. Abundance units versus the detection knee.**
The protocol defines `viral_abundance` per infected cell. D17 fits *per-sample* detection, and MECH-B applies the threshold to the summed count. The knee therefore sits at roughly abundance × the number of infected cells, so a varying `infected_cell_fraction` would smear the probit x-axis.
**Decided (user, 2026-10-05): option (a).**
- (a) **Chosen:** redefine `viral_abundance` as total planted molecules per virus per sample (the levels above). The x-axis is then the quantity the threshold acts on. `infected_cell_fraction` then means "how spread out", not "how much".
- (b) Keep per-infected-cell units and hold the fraction fixed (1 %) inside the LoD grid. Fraction would vary only in a separate covariate arm.
- (c) Keep both as they are and report the smear.

Consequence of (a): at 1 molecule and a fraction of 0.25, at most 1 of 500 designated cells receives anything. The truth `is_infected_cell` must therefore mean "received ≥ 1 molecule", or cell recall would count undetectable cells as misses.

**D4. Leakage rule versus fixed viral templates.**
The leakage rule says "No template, locus, molecule, UMI, or cell barcode is shared across partitions". Both partitions plant the same six genomes and gene models.
**Decided (user, 2026-10-05): option (a).**
- (a) **Chosen:** amend the rule so that "template" and "locus" mean host and challenge (GRCh38 homologous) sources. Viral genomes are the estimand, not a leakage channel. Molecules, UMIs and barcodes stay disjoint, as in step 5.
- (b) Disjoint viral genes per partition. This breaks per-virus LoD comparability between partitions.

The rule also forces two generator behaviours:
- barcode pools are split by partition;
- index-hopped reads come only from a sample in the *same* partition.

Both are built into steps 3–5.

**D5. Scale and truth granularity.**
**Decided (user, 2026-10-05): option (b). Keep 25,000 reads per cell and change the truth granularity.** The run-cost risk against the 20 k core-hour tuning budget is accepted. Measure the real per-sample cost on the tiny and pilot panels before the grid is sized. Real libraries shallower than 25 k reads per cell are recorded at their own depth (D1).
- (a) Rejected: keep 2,000 cells; amend `reads_per_cell` to 5,000 (10 M pairs per sample, ×0.2 storage and CPU). Write truth **per read for every non-host read** (viral, ambient, hopped, artefact, challenge), because D15 needs exact paired-read truth. Write **per molecule for host molecules**, with the host molecule ID encoded in the read name (`h:<ENST>:<mol>`) **in training only**. Holdout uses opaque names (see D6). The truth stays exact, and the manifests shrink by about 10⁴.
- (b) **Chosen:** keep 25 k reads per cell with the same truth change. The run cost is about 5× that of (a).
- (c) Halve the cells to 1,000. This conflicts with the rationale "1 % infected → ~20 cells".

**D6. Holdout truth location (R3.4) versus the protocol locators.**
R3.4 requires holdout truth to live outside the repository, owner-only, with its hash recorded in the protocol. The asset locators point at `analysis/v3_validation/generated/*.tsv`.
- **Recommended:** the generator writes training truth to the repo locator and holdout truth to an owner-only path given on the command line, recording only the hash. Amend the locators to say this.
- The FASTQs must not carry truth either. For holdout samples the generator:
  - writes opaque read names, `r<n>`, mapped to truth only in the owner-only manifest;
  - writes reads in a seeded shuffled order, so file position cannot reveal host and viral blocks;
  - writes `run_manifest.json`, which holds the factor levels, planted abundance and stage seeds, to the owner-only path.

  Next to the holdout FASTQs it leaves only the sha256 values and the sample ID.
- **Seed secrecy.** A public `seeds.root` plus a deterministic generator lets anyone regenerate holdout truth. `seeds.generation` is still pending, so the choice is open now:
  - accept that blinding is procedural only, and say so in the protocol;
  - **chosen (user, 2026-10-05):** give holdout generation an owner-held secret salt, publishing only `sha256(salt)` in the protocol and revealing the salt at unblinding. Training stays reproducible from the public seed alone.

**D7. Dependencies that keep the production panel blocked.**
**Decided (user, 2026-10-05): do the REF work first.** VAL-01 coding pauses until the host-homology levels exist.
- Order:
  1. evaluate REF-06 (the cat42d D-list jobs 25695872–76 completed and are not yet assessed);
  2. REF-07's host-homology table (GRCh38 loci homologous to the panel: identity, aligned length, source);
  3. freeze the `host_virus_homology` levels from that characterization.
- REF-08's *threshold calibration* "uses only preregistered training controls". Those controls are the VAL-01 training panel's host-homology negatives, so threshold calibration has to follow the panel and cannot precede it. The factor levels need only the REF-07 characterization, which can come first.
- `host_virus_homology` levels come from REF-08, which needs REF-06 (GRCh38 D-list, not yet run). The stratum count, the sample floor, and therefore the panel cannot be enumerated until then.
- The generator accepts levels and a homologous-locus TSV as inputs. It refuses a production run while `factors.host_virus_homology.status != frozen`.
- The tiny panel runs with `none`.
- Viral gene models come from `curated_reference_manifest`, which is still pending. The tiny panel pins the current `viral_ref_cat42b` build by sha256.
- 5′ chemistry (Q5) needs its read model (reads at the TSS; R2 antisense under `-x 10xv2`) and the DEF-02 strand handling. It is out of the first build and added with DEF-00.

## 7. Tiny golden panel (VAL-08 target)

- One sample per chemistry, so three in total.
- 50 cells and 20 empty droplets per sample, at 200 reads per cell.
- All six targets at abundance 10, `co_infected`.
- One molecule of each collision kind, one ambient molecule, one hopped molecule, and one read of each artefact class.
- Homology level `none`.

Expected outputs and their sha256 values are committed. A pytest regenerates the panel twice and asserts byte identity, and checks that every manifest row maps to exactly one read.

## 8. Code layout (ponytail)

- `scripts/generate_truth_panel.py` holds the CLI and orchestration. It reads levels from `protocol.yaml`, so levels have one source of truth and the generator has no grid file of its own.
- The simulation primitives are lifted out of `scripts/plant_anello_10x.py` and parametrised over `read_len` and an `np.random.Generator`; the planter imports them back. The primitives are `read_fasta`, `circular_slice`, `polya_read`, `artefact_read`, `mutate` (vectorised), `umi` and `fastq`. The ANDET-09e planter keeps working, and `tests/test_plant_anello_tails.py` guards it.
- `scripts/v3_partition.py` holds the largest-remainder split, which the scorer shares.
- Skipped: a plugin or registry for factor constructions. Nine factors means nine functions.

## 9. What review should check

1. Are D1–D7 the full set of contradictions, and is each recommendation defensible?
2. Do the abundance levels bracket a knee? This rests on an unmeasured detection threshold; the grid assumes the knee lies at 3–30 summed molecules.
3. Is cycling covariates within strata after the split sound? Or should the covariates be fully crossed in a smaller sub-panel?
4. Does the `same` collision construction match `collision_rule` as the scorer will read it?
