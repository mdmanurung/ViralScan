# ViralScan Showcase Runbook (kb-python edition)

**Purpose:** A self-contained, machine-executable spec for a coding agent to (1) resolve and
download a small curated set of public scRNA-seq datasets (plus one local sample), (2) point at a
ready combined host+virus reference, (3) run the **ViralScan** viral-load quantification workflow per
sample, and (4) apply calibration gates to produce a comparable cross-dataset report.

**Engine — read this first.** This repository's ViralScan is a **kb-python (kallisto | bustools) +
Snakemake** CLI (`src/viralscan/menu.py`), *not* the STARsolo tool of the same name. Its recommended
design is a **single pseudo-alignment pass against a combined host+virus reference**, where host and
viral targets compete in one quantification. There is **no host-pass / unmapped-read re-alignment
step**, **no manual whitelist or barcode-geometry handling** (kb derives geometry from `-x`), and
multimapping is corrected internally. Any STARsolo two-pass machinery from earlier drafts is gone.

**Intended consumer:** an LLM coding agent with shell access, the `viralscan` conda env on PATH,
≥16 cores, ≥64 GB RAM, ~1 TB scratch.

**Hard rule — raw FASTQ only:** `viralscan` pseudo-aligns **raw paired-end FASTQ** with `kb count`.
It does **not** ingest CellRanger BAMs or remap a processed matrix. A dataset is usable ONLY if raw
FASTQ is retrievable. Processed-matrix-only and BAM-only deposits are NOT usable (§5 gate).

**Read order matters:** for 10x, **`-s1` is the barcode+UMI read (R1)** and **`-s2` is the cDNA read
(R2)** — kb consumes them in that order. The per-sample output dir is named from the `-s1` filename
*before the first underscore* (`SRR12345_1.fastq.gz` → `output/SRR12345/`).

---

## 0. Global config
```bash
# Edit these, then `source` this block in every step.
export VS_ROOT=/exports/para-lipg-hpc/mdmanurung/viralscan_showcase   # writable workspace root
export VS_DATA=$VS_ROOT/data                  # downloaded FASTQ, per dataset id
export VS_OUT=$VS_ROOT/out                    # per-sample viralscan results
export VS_THREADS=16
export VS_SE_THRESHOLD=10                      # super-expressor call: >=10 viral UMI/cell (Lareau; this is the viralscan default)

# --- Ready combined host+virus reference (evonk; see §2). Used for EVERY sample. ---
export VS_REFDIR=/exports/archive/hg-funcgenom-research/evonk/old/intern/fasta_viruses/Serratus/transcriptome_human_hhv6_v2
export VS_INDEX=$VS_REFDIR/index_serratus.idx
export VS_T2G=$VS_REFDIR/t2g_serratus.txt
export VS_VIRAL_GTF=$VS_DATA/ref/viral_from_index.gtf   # generated in §2; passed via -gtf to every run

# --- Local test sample (evonk; read-only archive). ---
export VS_LOCAL=/exports/archive/hg-funcgenom-research/evonk/viralscan/data/healthy_samples

export NCBI_EMAIL=you@example.org              # only needed for the build-ref fallback (§2)
mkdir -p $VS_DATA/ref $VS_OUT
```

---

## 1. Environment setup (run once)
```bash
# ViralScan + its runtime engines. kb-python provides `kb`; snakemake is required.
# Download helpers (sra-tools/pysradb/ffq/seqkit) are for fetching public FASTQ only —
# they are NOT part of viralscan.
mamba create -y -n viralscan -c bioconda -c conda-forge \
  python=3.11 kb-python snakemake-minimal \
  sra-tools=3.1 pysradb ffq seqkit pigz pandas scipy numpy
mamba activate viralscan
python -m pip install -e /exports/para-lipg-hpc/mdmanurung/ViralScan   # this repo

# Sanity: both external binaries must be on PATH (viralscan preflights these).
command -v kb snakemake || { echo "kb / snakemake missing"; exit 1; }
viralscan --help >/dev/null && echo "viralscan OK"
# Require >= 2.3.0 — the --se-threshold / --multimap-method / --multimap-primary-call flags below
# do NOT exist in 2.2.0. (Verified failure mode: "unrecognized arguments: --se-threshold ...".)
python -c "import viralscan, packaging.version as v; assert v.parse(viralscan.__version__) >= v.parse('2.3.0')" \
  2>/dev/null || echo "WARN: viralscan < 2.3.0 — install this repo (pip install -e), do not rely on an older env copy"
```

> **Do not reuse a stale conda env's `viralscan`.** If you only have an older copy on PATH (e.g. an
> HPC env pinned at 2.2.0) but the repo checkout is newer, run the repo directly:
> `PYTHONPATH=<repo>/src python <repo>/src/viralscan/menu.py …` (this reuses the env's `kb`/`snakemake`
> while running the current CLI).

---

## 2. Reference — one ready combined index for everything
The reference is already built. `$VS_REFDIR` is a **single combined human + viral-panel kallisto
index** with host competition baked in:

- **226,005 human** transcripts (Ensembl ENST) + **~2,854 viral** transcripts across **187 viral
  accessions**.
- Contains every target this showcase needs (and more), so the **same index is used for all P0
  targets and for discovery**:

| Target | Accession in index | Rows |
|---|---|---|
| HHV-6B (Lareau strain) | `AF157706.1` | 97 |
| EBV (HHV-4) | `NC_007605.1` | 139 |
| HSV-1 (HHV-1) | `NC_001806.2` | 80 |
| HCMV (HHV-5) | `NC_006273.2` | 193 |
| KSHV (HHV-8) | `NC_009333.1` | 99 |
| + VZV, HHV-6A/7, HTLV, poxviruses, … | (≈180 more) | |

Because host (human cDNA) is in the same index, host and virus compete in one `kb count` pass — the
README's recommended design — so **no separate host-subtraction is needed**.

```bash
# Verify the index and its targets before running (read-only sanity check):
ls -lh "$VS_INDEX" "$VS_T2G"
for acc in AF157706 NC_007605 NC_001806 NC_006273 NC_009333; do
  printf "  %-12s rows=%s\n" "$acc" "$(grep -c "$acc" "$VS_T2G")"
done
```

### Required: derive the viral-gene GTF for this index (`-gtf`)
ViralScan's `analysis` step needs a GTF that tells it **which genes are viral**. With a *custom*
prebuilt index it will NOT use the bundled 195-virus panel (different gene naming → would match
nothing, and `viralscan data fetch` is irrelevant here). The authoritative list of viral genes for
this index is its own **t2g** — every non-host (non-`ENST`) row. Generate a minimal matching GTF
once and pass it via `-gtf` on every run:

```bash
# One-time: build a viral-gene GTF whose gene_id values exactly match the index's viral var_names.
grep -v '^ENST' "$VS_T2G" | awk -F'\t' '{print $2}' | sort -u | \
  awk 'BEGIN{OFS="\t"} {print $1,"ViralScan","exon","1","1",".","+",".","gene_id \""$1"\";"}' \
  > "$VS_VIRAL_GTF"
wc -l "$VS_VIRAL_GTF"   # ~2803 viral genes for transcriptome_human_hhv6_v2
# Spot-check the P0 targets are present:
for k in HUM_HERP6B EPSTEIN_HHV4 HUM_HERP1; do printf "  %-14s %s\n" "$k" "$(grep -c "gene_id \"$k" "$VS_VIRAL_GTF")"; done
```
> Why not `viral_genes.gtf.gz` / `all_gtf.gtf` shipped alongside the index? They don't match: the
> index gene names are prefixed (`ADENO_AAV2gp01`, `HUM_HERP6B_*`) while `all_gtf.gtf` is unprefixed
> (`AAV2gp01`, 0 overlap) and `viral_genes.gtf.gz` predates the HHV-6B addition (0 `HUM_HERP6B`).
> Deriving from the t2g is exact (2803/2803) and self-consistent with whatever index you point at.

### Portability fallback (only if evonk is unavailable)
Rebuild an equivalent combined index from public sources with the bundled `build-ref`:
```bash
viralscan build-ref --host human \
  --virus-accessions NC_007605.1 NC_001806.2 AF157706.1 NC_006273.2 NC_009333.1 \
  -o $VS_ROOT/ref_fallback --ncbi-email $NCBI_EMAIL
# -> $VS_ROOT/ref_fallback/{index.idx,t2g.txt}; then set VS_INDEX / VS_T2G to those paths.
# (Add more --virus-accessions to widen the discovery panel.)
```

---

## 3. Curated dataset registry (P0 subset)
Scope is intentionally small: one local smoke test + three public positive controls, **all run
against `$VS_INDEX`**. `access`: `LOCAL` (on-disk) · `PUBLIC_RAW` (raw FASTQ public) · `VERIFY`
(confirm raw reads in SRA first, §6a).

### 3a. Human-readable table
`-x` values below are **dry-run-verified** (1M-read subsample), not inferred — see §11.

| id | source | rep. SRR | access | host | `-x` (verified) | role | expected virus | dry-run |
|---|---|---|---|---|---|---|---|---|
| local_skn | evonk `$VS_LOCAL` (WS_SKN_KCL9369632) | — (S2 lane) | LOCAL | human | `10xv3` | **smoke test** (skin) | none | ✅ mechanics (empty, expected) |
| hhv6_carT_ref | GSE210063 | SRR20710641 | PUBLIC_RAW | human | `10xv3` | **positive control** (Lareau) | HHV-6B | ✅ HHV-6B @ 0.04% (Lareau band) |
| lcl_5lines | GSE158275 | SRR12682296 | PUBLIC_RAW | human | `10xv2` | positive control | EBV | ✅ EBV @ 3.34%, 285 SE |
| hsv1_fibroblast | GSE123782 (**Drop-seq**) | SRR8315713 | PUBLIC_RAW | human | `DROPSEQ` | positive control | HSV-1 | ✅ HSV-1 @ 0.31%, 4 SE |

> **GSE123782 mixes two assays** — *bulk* "stranded polyA RNA-seq" runs (2×76 bp, NOT usable) and
> *Drop-seq single-cell* runs (GSM3511317+, R1 = 20 bp = 12 CB + 8 UMI). Use only the **single-cell,
> infected** Drop-seq runs with `-x DROPSEQ`. The validated control is `SRR8315713` (single-cell,
> synchronous infection, **5 hpi** = peak lytic). Do **not** pick a bulk run (e.g. `SRR8315677`,
> which is also *uninfected*) and do **not** use a 10x `-x` here.
>
> The local sample's `-x 10xv3` is measured (R1 = 28 bp). EBV was *inferred* 10xv3 but is actually
> **10xv2** (R1 = 26 bp); HSV-1 is **Drop-seq**, not 10x — all caught by the dry run. **Always derive
> `-x` from the measured R1 length (§6b), do not trust the GEO tag.**

### 3b. Machine manifest (agent iterates over this — TSV; hsv1 row commented out)
```tsv
id	source	srr	access	host	technology	role	expected_virus	resolver
local_skn	WS_SKN_KCL9369632	-	LOCAL	human	10xv3	smoke_test	none	local_glob
hhv6_carT_ref	GSE210063	SRR20710641	PUBLIC_RAW	human	10xv3	positive_control	HHV-6B	pysradb_srp
lcl_5lines	GSE158275	SRR12682296	PUBLIC_RAW	human	10xv2	positive_control	EBV	pysradb_srp
hsv1_fibroblast	GSE123782	SRR8315713	PUBLIC_RAW	human	DROPSEQ	positive_control	HSV-1	pysradb_srp_dropseq_5hpi
```

---

## 5. Access-status gate (run per id)
```bash
gate_access () {  # $1=access
  case "$1" in
    LOCAL)      echo "PROCEED_LOCAL" ;;       # skip download; use on-disk FASTQ
    PUBLIC_RAW) echo "PROCEED" ;;
    VERIFY)     echo "VERIFY_FIRST" ;;        # §6a: if SRR has raw reads -> PROCEED, else SKIP_LOG
    *)          echo "SKIP_LOG: needs raw FASTQ; viralscan cannot scan processed/controlled data" ;;
  esac
}
```

---

## 6. Canonical per-sample workflow (reusable)

### 6a. Resolve + verify raw availability (public datasets only)
```bash
resolve_srr () {  # $1=GSE -> SRR list.  NB: pysradb has NO `gse-to-srr`; chain gse-to-srp | srp-to-srr.
  local srp
  srp=$(pysradb gse-to-srp "$1" 2>/dev/null | awk 'NR==2{print $2}')
  [ -n "$srp" ] && pysradb srp-to-srr "$srp" 2>/dev/null | awk 'NR>1{print $2}' | sort -u
}
# Before running, read the sample TITLE — confirm it is a 10x droplet run and the right condition
# (e.g. infected vs uninfected). The dry run found GSE123782 was bulk polyA RNA-seq, not 10x.
srr_title () { pysradb metadata "$1" 2>/dev/null | awk -F'\t' 'NR==2'; }
```

### 6b. Download FASTQ (ENA-first, SRA fallback) and identify R1/R2
```bash
fetch_fastq () {  # $1=id  $2=SRR  -> prints downloaded *.fastq.gz
  local d=$VS_DATA/$1/$2; mkdir -p "$d"; cd "$d"
  local urls=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$2&result=read_run&fields=fastq_ftp" \
               | awk 'NR==2{print $NF}' | tr ';' '\n')
  if [ -n "$urls" ]; then for u in $urls; do wget -q "ftp://$u"; done
  else  # SRA fallback; --include-technical keeps the 10x barcode read (R1)
    prefetch -O "$d" "$2" && fasterq-dump --split-files --include-technical -e $VS_THREADS -O "$d" "$d/$2/$2.sra"
    pigz -p $VS_THREADS "$d"/*.fastq
  fi
  ls -1 "$d"/*.fastq.gz
}
# viralscan wants: -s1 = barcode read (R1), -s2 = cDNA read (LONGER R2).
read_len () { seqkit stats -T "$1" | awk 'NR==2{print $7}'; }   # max read length

# CHOOSE -x FROM THE MEASURED R1 LENGTH — never trust the GEO chemistry tag (EBV was tagged v3 but
# is v2). R1 == 26 bp -> 10xv2 ;  R1 == 28 bp -> 10xv3.  If the two reads are EQUAL length (e.g.
# 76/76, 150/150) it is likely NOT droplet 10x at all — verify the title before proceeding.
pick_x () {  # $1=R1.fastq.gz  -> 10xv2 | 10xv3 | UNKNOWN
  case "$(read_len "$1")" in 26) echo 10xv2;; 28) echo 10xv3;; *) echo UNKNOWN;; esac
}
```

When streaming a subsample with `fastq-dump -X N --split-files --gzip`, assign reads by length:
**cDNA = longest read → `-s2`; barcode = the next-longest with len ≥ 20 → `-s1`** (an 8 bp `I1`
index read, if present, is ignored). Equal-length pairs ⇒ not 10x ⇒ skip (see §9).

### 6c. Per-sample ViralScan run — a single command runs the whole DAG
```bash
run_viralscan () {  # $1=id  $2=label(for out subdir via -s1 prefix)  $3=technology  $4=R1(barcode)  $5=R2(cDNA)
  viralscan \
    -i "$VS_INDEX" \
    -t "$VS_T2G" \
    -gtf "$VS_VIRAL_GTF" \
    -o "$VS_OUT/$1/" \
    -s1 "$4" -s2 "$5" \
    -x "$3" \
    -c $VS_THREADS \
    --se-threshold $VS_SE_THRESHOLD \
    --detection-threshold 1 \
    --multimap-method host-conservative \
    --multimap-primary-call confidence
  # Output lands in $VS_OUT/$1/<R1-prefix>/  (named from $4 basename before first '_').
}
# That single command runs: create_config -> kb_count -> analysis -> multimap -> detection -> umap,
# producing results/{viral_summary.tsv,per_cell_viral.tsv,multimap_evidence.tsv}, report.html,
# plots/, and kb-python/counts_unfiltered/adata_multimap.h5ad.
#
# Add --umap for UMAP HTML (slower). Add --cell-types labels.csv (barcode,cell_type) for the
# Fisher-exact per-cell-type viral enrichment table (results/cell_type_enrichment.tsv).
# host-conservative multimapping keeps host-virus ambiguous EC mass out of primary viral counts.
```

---

## 7. Per-dataset action plans (concrete invocations)

```bash
# ---- local_skn : smoke test (validate pipeline mechanics on real local data) ----
# Skin tissue; biologically we do NOT expect a strong viral signal. Success = the run completes
# and results/report.html are produced. A near-zero viral_summary is the EXPECTED outcome here.
for S in S2 S3 S4; do   # run script convention used S2/S3/S4 lanes
  R1=$VS_LOCAL/WS_SKN_KCL9369632_${S}_L001_R1_001.fastq.gz
  R2=$VS_LOCAL/WS_SKN_KCL9369632_${S}_L001_R2_001.fastq.gz
  run_viralscan local_skn WS_SKN_${S} 10xv3 "$R1" "$R2"
done

# ---- hhv6_carT_ref (GSE210063) : Lareau HHV-6B positive control -> sets the calibration ----
for SRR in $(resolve_srr GSE210063); do
  [ "$(srr_has_reads $SRR)" = OK ] || continue
  mapfile -t FQ < <(fetch_fastq hhv6_carT_ref $SRR)
  R1=$(printf '%s\n' "${FQ[@]}" | grep -E '_1\.|_R1' | head -1)   # barcode read
  R2=$(printf '%s\n' "${FQ[@]}" | grep -E '_2\.|_R2' | head -1)   # cDNA read
  run_viralscan hhv6_carT_ref $SRR 10xv3 "$R1" "$R2"
done

# ---- lcl_5lines (GSE158275) : EBV+ LCLs -> expect high EBV total_umi, large pct_infected ----
for SRR in $(resolve_srr GSE158275); do
  [ "$(srr_has_reads $SRR)" = OK ] || continue
  mapfile -t FQ < <(fetch_fastq lcl_5lines $SRR)
  R1=$(printf '%s\n' "${FQ[@]}" | grep -E '_1\.|_R1' | head -1)
  R2=$(printf '%s\n' "${FQ[@]}" | grep -E '_2\.|_R2' | head -1)
  run_viralscan lcl_5lines $SRR "$(pick_x "$R1")" "$R1" "$R2"   # dry-run verified: 10xv2
done

# ---- hsv1_fibroblast (GSE123782, Drop-seq) : 5hpi infected single-cell -> HSV-1 lytic ----
# Use only the single-cell INFECTED runs (GSM3511317+); SRR8315713 = 5hpi, validated. -x DROPSEQ.
for SRR in SRR8315713; do
  [ "$(srr_has_reads $SRR)" = OK ] || continue
  mapfile -t FQ < <(fetch_fastq hsv1_fibroblast $SRR)
  R1=$(printf '%s\n' "${FQ[@]}" | grep -E '_1\.|_R1' | head -1)   # 20 bp Drop-seq barcode
  R2=$(printf '%s\n' "${FQ[@]}" | grep -E '_2\.|_R2' | head -1)
  run_viralscan hsv1_fibroblast $SRR DROPSEQ "$R1" "$R2"
done
```

---

## 8. Validation & calibration gates
ViralScan already writes per-virus and per-cell tables; the gate reads them directly — no custom
matrix parsing.

```bash
score_sample () {  # $1=id  $2=run_subdir (the R1-prefix dir under $VS_OUT/$1)
  local r=$VS_OUT/$1/$2/results
  echo "== $1/$2 =="
  echo "-- viral_summary.tsv (virus_name total_umi infected_cells pct_infected) --"
  column -t -s$'\t' "$r/viral_summary.tsv" | sort -k2 -nr | head
  echo "-- super-expressors at >=$VS_SE_THRESHOLD UMI (from per_cell_viral.tsv) --"
  # per_cell_viral.tsv columns: barcode  virus_name  viral_umi  total_umi  viral_fraction
  awk -F'\t' -v T=$VS_SE_THRESHOLD 'NR>1 && $3>=T {c[$2]++} END{for(v in c) print v, c[v]}' \
    "$r/per_cell_viral.tsv" | sort -k2 -nr
}
```
**Gates**
- **Smoke-test gate (`local_skn`):** the run must complete and emit `results/viral_summary.tsv` +
  `report.html`. A near-zero viral signal is acceptable/expected — this validates mechanics, not
  biology. If the run *errors* (e.g. "no reads pseudoaligned"), fix `-x`/read order before anything
  else (§9).
- **Positive-control gate:** each of `hhv6_carT_ref` / `lcl_5lines` / `hsv1_fibroblast` MUST show its
  expected virus (HHV-6B / EBV / HSV-1) as a top entry in `viral_summary.tsv`. (Dry-run confirmed:
  HHV-6B @ 0.04%, EBV @ 3.34% / 285 SE, HSV-1 @ 0.31%.) If not → wrong `-x`, swapped R1/R2, or
  mislabeled sample. **Halt and
  report.**
- **Calibration (Lareau):** true reactivation cultures show super-expressors at ~**0.01–0.3 %** of
  cells at ≥10 viral UMI (`pct_infected`); `--se-threshold 10` is already the viralscan default.
  Use the HHV-6 control to confirm this rate is sane.
- **Discovery (free, same index):** because the index carries ~187 viruses, every run also reports
  off-target viruses. Flag a (sample, virus) hit only if super-expressors ≥ 2 AND it is
  biologically plausible; cross-check `multimap_evidence.tsv` `call_confidence` (prefer `strong`;
  treat `low_confidence`, i.e. host-virus ambiguous only, as suspect). Manually review hits.
- **Negative floor:** the `local_skn` smoke test doubles as a negative — any virus with
  super-expressors there marks a false-positive floor to subtract elsewhere.

---

## 9. Failure handling (autonomous rules)
- **`kb count` aborts with "no reads pseudoaligned"** (the Snakefile checks for this) → almost always
  wrong `-x` chemistry or swapped R1/R2. Retry the other 10x version (`10xv2`↔`10xv3`) and confirm
  R1=barcode (shorter) / R2=cDNA (longer) by length; if still empty, mark `chemistry_unresolved`
  and skip.
- **No SRR / `VERIFY` → `NO_READS`** → `SKIP_LOG`, continue.
- **ENA empty** → SRA `prefetch` fallback; if both fail twice → `SKIP_LOG` with the accession.
- **Two reads of EQUAL length** (e.g. 76/76, 150/150) → likely **not droplet 10x** (bulk/plate
  RNA-seq). Read the sample title; if not 10x, **exclude the dataset** (kb will fabricate barcodes
  and produce meaningless "cells" — as GSE123782 did).
- **Read order ambiguous** (R1/R2 similar length, e.g. some 5′ data) → cDNA = longer read = `-s2`;
  log the choice.
- **`pysradb gse-to-srr` "usage" error** → that subcommand does not exist; use
  `pysradb gse-to-srp <GSE> | …` then `pysradb srp-to-srr <SRP>` (fixed in §6a).
- **`unrecognized arguments: --se-threshold …`** → the `viralscan` on PATH is < 2.3.0; run the repo
  via `PYTHONPATH` (see §1).
- **`Viral reference annotations were not found … data fetch`** with a custom prebuilt index → you
  did not pass `-gtf`; supply the t2g-derived viral GTF from §2.
- **OOM in `kb count`** → lower `-c/--cores`. The combined index is host-dominated, so memory scales
  with the human transcriptome (~64 GB is comfortable for this index).
- **Output dir exists & non-empty** → `viralscan` prompts to overwrite; pre-clean `$VS_OUT/<id>/` in
  unattended runs.
- **build-ref fallback fails** → check `--ncbi-email` is set; downloads cache under
  `~/.cache/viralscan/`, so retries are cheap.

---

## 10. Output schema (one folder per sample)
`viralscan` writes a sample dir named from the `-s1` filename before the first underscore:
```
$VS_OUT/<id>/<R1-prefix>/
  config.yaml
  summary.txt
  report.html                                   # self-contained interactive report
  results/
    viral_summary.tsv                           # per-virus: total_umi, infected_cells, total_cells, pct_infected, umi_per_10k
    per_cell_viral.tsv                           # per barcode×virus: viral_umi, total_umi, viral_fraction
    multimap_evidence.tsv                        # unique / ambiguous / host-virus-ambiguous + call_confidence
    cell_type_enrichment.tsv                     # only if --cell-types given
  kb-python/counts_unfiltered/
    adata.h5ad
    adata_multimap.h5ad                          # layers: counts_original, counts_corrected, ...
  plots/                                          # *_histogram.png, SuperExpressor_*.png, umap_*.html (if --umap)
  log/                                            # create_config/kb/analysis/multimap/detection/umap .done
```
**Final deliverable:** a single `run_report.tsv` ranking every processed sample by super-expressor
count and top virus (built from §8 `score_sample` across all runs), plus `skipped.tsv`
(failed-verify / errored). Columns:
```
id  label  role  technology  n_cells  top_virus  top_virus_umi  pct_infected  n_superexpressors  call_confidence  status  notes
```

---

## 11. Dry-run validation (1M-read subsample, SLURM, 2026-06-21)
The whole chain was exercised on a 1M-read subsample per dataset against the evonk combined index
(repo viralscan 2.3.0, evonk env for `kb`/`snakemake`, t2g-derived `-gtf`). Runtime ~3–5 min/job,
peak RAM ~0.7 GB (index ≈ 0.4 GB) on 8 cores — so **full runs fit comfortably in `-c 8 --mem 32G`,
`-t 08:00:00`**, scaling roughly linearly with read count.

| id | sample | `-x` | result | notes |
|---|---|---|---|---|
| local_skn | WS_SKN S2 | 10xv3 | ✅ mechanics PASS | empty viral summary (healthy skin); 2 host-virus-ambiguous genes correctly excluded → no false positives |
| hhv6_carT_ref | SRR20710641 | 10xv3 | ✅ HHV-6B @ 0.041 % | within Lareau 0.01–0.3 % band |
| lcl_5lines | SRR12682296 | **10xv2** | ✅ EBV @ 3.34 %, 285 SE | strong; needed v2 not the inferred v3 |
| hsv1_fibroblast | SRR8315713 | **DROPSEQ** | ✅ HSV-1 @ 0.31 %, 4 SE | 5 hpi infected Drop-seq run; GSE123782's *bulk* runs are unusable |

All jobs ~3–5 min at 1M reads. Artifacts: workspace
`/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/dryrun/{sbatch,out,logs,dry_run_report.tsv}`.

**Gate verdict:** mechanics PASS on all four; **all three positive controls recover their expected
virus** (HHV-6B, EBV, HSV-1) at biologically sensible rates. The showcase is ready for full-depth
runs — note the per-dataset chemistry varies (`10xv3`, `10xv2`, `DROPSEQ`).

---

### Sourcing & honesty notes
- All pipeline commands map to the real `viralscan` CLI (default quantification mode + the
  `build-ref` fallback) and its documented outputs (`docs/output_reference.md`). The single-pass
  combined-reference design and built-in multimapping correction replace any STARsolo two-pass
  host-subtraction, which belongs to a different tool.
- The reference (`transcriptome_human_hhv6_v2`) is an **existing evonk asset** verified to contain
  226,005 human + ~2,854 viral transcripts over 187 accessions, including all P0 targets. No
  reference build is performed on the fast path.
- Calibration numbers (≥10 UMI; ~0.01–0.3 % of cells) are from Lareau et al. (HHV-6); `10` is the
  viralscan `--se-threshold` default.
- All three positive controls were **empirically validated** on a 1M-read subsample (§11), so the
  chemistries (`10xv3` HHV-6, `10xv2` EBV, `DROPSEQ` HSV-1) and SRRs are confirmed, not inferred.
  Still re-measure R1 length per run at full depth (§6b) — other runs in a series can differ.

---

## 12. KSHV (HHV-8) extension — co-infection showcase

**Rationale:** ViralScan's combined Serratus index carries KSHV (`HUM_HERP8_HHV8GK18*`, 86 genes)
and EBV (`EPSTEIN_HHV4_*`) in the same index used above. Two public KSHV scRNA-seq series are
processed here to demonstrate (a) **latent-vs-lytic viral-load contrast** and (b) **KSHV+EBV
co-infection detection** without rebuilding the reference.

### 12a. Dataset registry

| id | GSE | samples | SRR (lane-splits) | subsampling | expected |
|---|---|---|---|---|---|
| `kshv_islk219` | GSE190558 | GSM5725695 (latent iSLK.219) | SRR17180398–401 | none | HHV-8 (lytic ≫ latent) |
| `kshv_islk219` | GSE190558 | GSM5725696 (lytic iSLK.219) | SRR17180394–397 | none | HHV-8 |
| `kshv_islk219` | GSE190558 | GSM5725697 (lytic + casp-i) | SRR17180390–393 | none | HHV-8 |
| `kshv_islk219` | GSE190558 | GSM5725698 (lytic + casp-i + anti-IFN) | SRR17180386–389 | none | HHV-8 |
| `kshv_bc1_pel` | GSE154900 | GSM4682311 (BC-1 rep 1) | SRR12287751–752 | 75M spots/SRR | **HHV-8 + EBV** |
| `kshv_bc1_pel` | GSE154900 | GSM4682312 (BC-1 rep 2) | SRR12287753–754 | 75M spots/SRR | **HHV-8 + EBV** |

**GSE190558** (NextSeq 550, 10x): 4 SRRs per GSM condition are NextSeq lane-splits of the same
library — they **must be concatenated per-sample** before passing to viralscan to avoid UMI
collisions. Each lane alone covers <25 % of the barcodes.

**GSE154900** (NovaSeq 6000, 10x): BC-1 PEL (primary effusion lymphoma) is KSHV-positive and
EBV-positive — a double-infection, immunosuppression-relevant cell line. Raw depth is ~660M
paired-end reads per replicate, which is excessive for host-conservative multimapping. Subsample
to ~150M pairs by capping each SRR to 75M spots with `fasterq-dump -X 75000000`.

### 12b. Chemistry probe (always measure R1 — never trust the GEO tag)

```bash
# probe_kshv.sh: download 500k spots from one SRR per dataset, print read-length histogram.
# Submitted as SLURM job 25055294; results written to:
#   fullrun/logs/probe_kshv_25055294.out
# R1=26 bp → 10xv2  |  R1=28 bp → 10xv3
```

Confirmed chemistry (determined 2026-06-22 via ENA metadata — no probe needed):

| Dataset | SRR | Method | Evidence | TECH |
|---------|-----|--------|----------|------|
| iSLK.219 | SRR17180398 | bases/spot = 91bp = 28+8+55 (R1+I1+R3) | 10xv3 R1=28bp confirmed by arithmetic | `10xv3` |
| BC-1 PEL | SRR12287751 | ENA `_1.fastq.gz` byte-range fetch; reads = 28bp | 10xv3 barcode+UMI directly measured | `10xv3` |

Both are **10x Chromium v3** (16bp barcode + 12bp UMI = 28bp R1).

### 12c. Files (workspace: `fullrun/`)

| File | Purpose |
|------|---------|
| `manifest_kshv.tsv` | 6-row manifest: `id  sample  srr_csv  expect` |
| `sbatch/probe_kshv.sh` | chemistry probe (job 25055294) |
| `sbatch/full_public_kshv.sh` | array script: lane-concat, optional subsample, viralscan |
| `multi_virus_report.py` | multi-virus report; emits all viruses above noise floor per sample |

### 12d. Submitting the array

Chemistry is confirmed; `full_public_kshv.sh` already has `TECH=10xv3` hardcoded for both datasets. Submit directly:

```bash
cd /exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun
sbatch --array=1-6%2 sbatch/full_public_kshv.sh
```

`%2` keeps at most 2 concurrent tasks to avoid monopolising `highmem` while the primary array
(job 25053671) or other HPC jobs are running.

### 12e. Multi-virus co-infection report

After the KSHV array completes:

```bash
cd /exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun
python3 multi_virus_report.py   # KSHV samples only (default)
# Output: kshv_report.tsv  (id, sample, virus_name, total_umi, infected_cells, pct_infected, umi_per_10k)
```

Expected payoffs:
- `GSM5725695` (latent) KSHV UMI ≪ `GSM5725696` (lytic) — orders-of-magnitude contrast.
- `GSM4682311`/`GSM4682312` (BC-1 PEL): rows for **both HHV-8 and EBV** above floor.
  The script prints `*** CO-INFECTION ***` for samples with ≥2 viruses and a `PASS` line
  confirming KSHV+EBV co-detection.

### 12f. Gotchas

- **Lane concatenation:** skip per-lane partial runs; always cat all 4 lanes first. The array
  script does this automatically via the `srr_csv` column and a download-then-cat loop.
- **Subsampling:** `fasterq-dump -X 75000000` applied per SRR (not per sample); after concat
  BC-1 has ~150M paired reads — sufficient to detect both KSHV and EBV at >50 UMI with margin.
- **HHV-8 index prefix:** KSHV genes appear as `HUM_HERP8_HHV8GK18*` in the t2g-derived GTF
  and `viral_summary.tsv`. The `multi_virus_report.py` PASS check tests for `HHV8` or `HERP8`
  substring.
- **EBV prefix:** `EPSTEIN_HHV4_*`. The co-infection check tests for `HHV4`, `EPSTEIN`, or `EBV`.
- **No new reference needed:** the same `index_serratus.idx` and t2g-derived `viral_from_index.gtf`
  used for P0 are reused here.
