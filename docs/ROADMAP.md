# ViralScan roadmap

Future work, grounded in (M) the manuscript's stated limitations
(`docs/manuscript_draft.md` Discussion), (R) code-review findings from this
session, and (T) existing `todo/TODOLIST.md` items. Each item notes rationale,
approach, effort, dependencies, and grounding tag. This is a planning document,
not a commitment — research items are marked exploratory.

Tiers: **P0** publication/release-blocking · **P1** high-value near-term ·
**P2** medium · **P3** long-term / research.

---

## A. Scientific & methods

### A1 · Read-start distribution + PCR-dup handling in `evidence` — **P2**
Per-position read-start histogram along the viral genome (3′ bias, subgenomic-RNA
junctions, EVE hotspots, coverage evenness). Adds `read_start_distribution()` +
`--read-start-profile`/`--dedup {umi,markdup,none}` to `viralscan evidence`.
Full plan: [`todo/read-start-distribution.md`](../todo/read-start-distribution.md).
*Grounding:* user request; extends the paper method (Chen et al. Nature 2024). *~1–2 d.*

### A2 · Per-cluster / per-cell EM allocation — **P3 (exploratory)**
The EM estimates **one transcriptome-wide** abundance vector (`em_gene_abundances`,
global pool) and allocates every cell's multimappers with it. In highly-infected
cultures, per-cluster (EM within Leiden clusters) or per-cell (alevin-fry-style)
allocation may be preferable. *Approach:* optional `--multimap-scope
{global,cluster,cell}`; cluster scope reuses the existing scanpy clustering, runs the
existing vectorised EM per group. Watch memory/runtime — global is the memory-cheap
default. *Grounding:* (M) "estimates one transcriptome-wide abundance vector … per-cell
allocation may be preferable." *Large; benchmark against the global default before
shipping.*

### A3 · Cell-level viral BAM with CB/UB tags — **P2**
`evidence` extracts reads to FASTA then realigns; produce a first-class per-cell BAM
with `CB`/`UB` tags so users can load viral reads in IGV **grouped by cell**. *Approach:*
carry the `(CB,UMI)` from the read name into BAM tags during/after
`align_reads_to_viral`; document an IGV "group by CB" recipe. *Grounding:* (M) "does not
produce cell-level BAM files … genome-browser inspection remains a separate alignment
step." *~1–2 d.*

### A4 · Combined-genome specificity mode / genome D-list — **P1**
The cDNA-only host reference cannot absorb reads from GRCh38 **non-coding** regions
that resemble viral references → spurious viral signal (the anellovirus / TTV
host-homology artifact: ~90% "viral" under cDNA-only vs 0/4.5M on a combined genome).
*Approach:* official `--genome-dlist genome.fa` in `build-ref` (mask host-homologous
k-mers via kallisto D-list; partially designed per PLAN 2026-07-06) **and/or** a
documented "genomic specificity control" workflow (orthogonal combined-genome STAR/kb
run) surfaced as a first-class QC. *Grounding:* (M) cDNA-only host-homology limitation;
"libraries with abundant non-coding reads … should include an orthogonal alignment."
*Medium; the D-list build is compute-heavy (~64 GB / 8 h).*

### A5 · Subgenomic-RNA junction detection (coronaviruses) — **P3 (exploratory)**
Once A1 lands, detect canonical leader–body junction pile-ups as a
genuine-transcription signature (vs genomic contamination). *Depends on A1.* *Research.*

---

## B. Validation & benchmarking

### B1 · Head-to-head vs dedicated viral scRNA tools — **P1**
The manuscript explicitly did **not** benchmark against Venus (Luebbert et al. 2024) or
ViralTrack (Bost et al. 2020) — the STARsolo comparison tests a general splice-aware
aligner only. *Approach:* extend the reference-strategy harness to run Venus/ViralTrack
on the same EBV/HHV-6B/HSV-1 samples; compare sensitivity, specificity, gene attribution,
runtime, memory. *Grounding:* (M) "did not complete a head-to-head comparison against
dedicated viral single-cell tools." *Large; strengthens the paper materially.*

### B2 · Gold-standard truth panel (planted-read simulation) — **P1**
No ground truth today — validation uses published rate ranges + matched comparisons.
*Approach:* spike known viral reads into a host library at known per-cell rates →
formal false-positive / false-negative rates, ROC, limit-of-detection. *Grounding:*
(M) "Formal FP/FN rates require negative controls and planted viral-read simulations …
left as a release-gated validation extension." *Large; the single biggest credibility
lever for the tool.*

### B3 · Complete the reference-strategy benchmark — **P1 (T)**
Only 4/12 rows done; all EBV rows failed, viralscan/two-step conditions blocked. Re-run
the failed SLURM conditions so `results/reference_strategy_benchmark.tsv` is analyzable.
*Grounding:* (T) TODOLIST open. *Unblocks B4/B5.*

### B4 · Selectivity Index + noisy-channel bit-cost — **P2 (T, blocked on B3)**
SI = on-target ÷ off-target viral UMI per (aligner × reference) → ranked reference
recommendation; per-pipeline I(true;measured) in bits (HHV-6A/6B BSC, HSV-1 denominator
artifact). *Grounding:* (T) TODOLIST, blocked on B3.

---

## C. Host-response module

### C1 · Gene-symbol annotation + pathway enrichment — **P2**
Host-response reports Ensembl feature IDs when the matrix lacks symbols, limiting
biological interpretation of the depth-robust gene set. *Approach:* map Ensembl→symbol
(gget/pybiomart), run enrichment (gseapy/gget enrichr — partially wired via
`viralscan[enrichment]`) on the **depth-controlled** labels, not the raw one. *Grounding:*
(M) "reports Ensembl feature IDs … limiting biological interpretation until annotation
and enrichment are added." *~1–2 d.*

### C2 · Additional depth-control designs — **P3**
The module already reports a depth-alone baseline, CPM-normalised and depth-matched
labels. Extend with propensity/covariate-adjusted designs and per-gene E-value reporting
surfaced in the HTML report. *Grounding:* (M) depth-confound discussion. *Exploratory.*

---

## D. Robustness / QC / correctness (review findings)

### D1 · Non-unique `var_names` handling — **P2 (R)**
`matrix_for_genes` uses `var_names.get_indexer`, which raises on duplicate accessions;
the codebase elsewhere anticipates non-unique names (`super_expressor` calls
`var_names_make_unique()`). *Approach:* validate/dedup var_names once at load, or make
`matrix_for_genes` first-occurrence-safe. Not a regression (old code also raised), but a
latent crash on duplicate reference accessions. *Grounding:* (R) primary-call review.

### D2 · EVE annotate robustness — **P3 (R)**
`_is_chromosome_subject` uses `re.match` anchored at start; old-style `gi|…|ref|NC_…|`
BLAST sseqids downgrade to `subject_not_chromosome` (safe failure, modern BLAST -outfmt6
unaffected). Add a `gi|`/`ref|` unwrap if legacy BLAST DBs are in scope. *Grounding:* (R)
EVE review. *Low.*

### D3 · `host_viral_ambig_fraction` >1.0 edge — **P3 (R)**
Under some data shapes the ambiguity numerator (separate layer) can exceed the adata.X
denominator → fraction >1.0. Pre-existing; document or clamp. *Grounding:* (R) detection
review. *Low.*

---

## E. Performance & scale

### E1 · numba the collapse loop — **P3**
The multimap loop is now ~3× faster (collapse) + memory-lean. For very deep samples, a
numba-jitted inner loop over the collapsed rows could give another 5–10×. *Approach:*
ragged CSR-style arrays for `ec_info`, `@njit` the per-record scatter. *Grounding:*
profiling; the loop is the runtime ceiling for 100M-record samples. *Medium; only if
deep-sample runtime becomes a user complaint.*

### E2 · Chunked/streaming BUS processing — **P3**
`build_multimap_layers` materialises the full `bus_df` (~5 GB RSS on deep samples).
Process the BUS file in chunks (aggregate `(cell,ec)` incrementally) to cap peak memory.
*Grounding:* profiling (peak-RSS on deep samples). *Medium.*

---

## F. Usability / packaging / release

### F1 · Open PR #7 → main — **P0 (T)**
The `claude/pub-readiness-hygiene` branch (25+ commits) needs its PR merged. Note the
git SSH commit-signing key mismatch (`user.signingkey`) — recent commits may be unsigned.
*Grounding:* (T).

### F2 · Tag v2.6.0 + publish (PyPI/conda) + Zenodo DOI — **P0**
Version is reconciled to 2.6.0; cut the tag, publish the wheel, mint a Zenodo DOI for the
GTF panel (already on Zenodo) and cite it in the manuscript data-availability statement.
*Needs owner (release action).* *Grounding:* release hygiene.

### F3 · Functional-script path hygiene — **P1 (T, R)**
~14 SLURM/infra scripts (`covid_viralscan/scripts/*.sh`, `scripts/*.sh`) and
`reference_manifest.json` embed `VS_CONDA_ENV=/exports/.../evonk/...` (third-party path)
as a default. Env-var-ize per script or exclude from the shipped artifact — decide ship
scope first. *Grounding:* (T, R) publication-hygiene branch goal. *Medium; per-script,
breaking-risk.*

### F4 · Stop tracking downloadable/generated data — **P2 (R)**
`ref/10x_version2_whitelist.txt` (3.8 MB, standard 10x file, referenced by 2 SLURM
scripts) and tracked `analysis/**/*.pdf` build artifacts. Fetch-instead-of-track the
whitelist; gitignore generated outputs. *Grounding:* (R) bloat analysis.

### F5 · CI gates: notebook execution + declared-deps — **P2**
Add an nbmake gate for the 6 CI-runnable vignettes (so `[skip-ci]` docs don't rot — the
enrichment vignette silently broke once) and a `deptry`/bare-venv import check to catch
undeclared deps despite CI's `--no-deps` install. *Grounding:* (R) vignette + CI-no-deps
learnings. *~1 d.*

### F6 · Docs: IGV recipe + read-start vignette — **P2**
Once A1/A3 land, add a "read-start profile & genome-browser inspection" vignette and an
IGV group-by-CB recipe to the 8-vignette suite. *Depends on A1, A3.*

---

## G. Reference data & annotation

### G1 · Reference annotation quality — **P2**
The EBV LMP-1/EBNA gene-attribution divergence vs STARsolo traces to **annotation**, not
the aligner (Serratus 96 EBV entries vs STARsolo's 16). Curate/version the viral GTF
panel; surface the reference annotation provenance in outputs (the manuscript argues for
"explicit reporting of viral reference annotations"). *Grounding:* (M) annotation-vs-model
discussion. *Medium.*

### G2 · Mouse / non-human host support — **P3**
Confirm/extend `build-ref --host mouse` (Ensembl) for non-human single-cell viral studies.
*Grounding:* generality. *Medium.*

---

## Suggested sequencing

1. **Ship the branch:** F1 (PR), F3 (path hygiene), F2 (tag/DOI) — publication-blocking.
2. **Strengthen the paper:** B2 (truth panel) and B1 (dedicated-tool benchmark) are the
   two highest-credibility additions; B3 unblocks the reference analyses.
3. **Feature depth:** A1 (read-start, requested) → A3 (cell BAM) → A4 (genome specificity).
4. **Polish:** C1 (gene symbols), F5 (CI gates), D1 (var_names robustness).
5. **Research:** A2 (per-cell EM), A5 (sgRNA junctions), E1/E2 (perf at scale).
