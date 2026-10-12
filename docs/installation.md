# Installation

## Requirements

Two tiers:

- **pip tier** (`pip install ViralScan`, Python ≥ 3.9): `build-ref`, `evidence`,
  `validate-run`, `doctor --profile pip`. No Snakemake or native tools needed.
- **Full tier** (quantification workflow): Python 3.11 on `linux-64` (the canonical,
  tested toolchain), [kb-python](https://www.kallistobus.tools/) (`kb` command) and
  [Snakemake](https://snakemake.readthedocs.io/en/stable/) ≥ 7.0, plus kallisto, bustools and STAR.

## Option 1 — Conda (recommended)

A ready-to-use conda environment file is provided:

```bash
conda env create -f environment.yml
conda activate viralscan
python -m pip install -e .
viralscan data fetch             # optional, see below
```

This installs Python, `kb-python`, `snakemake`, and all Python dependencies,
then installs the local ViralScan checkout. `viralscan data fetch` downloads
the bundled viral annotation panel, which a run needs only when its index has
no `<index>.build_manifest.json` (an index from `viralscan build-ref` has one
and carries its own viral gene list).

## Option 2 — pip (pip tier)

```bash
pip install ViralScan
viralscan doctor --profile pip
```

This does not install Snakemake. To add the Python workflow dependencies
(`snakemake`, `scanpy`, `scikit-learn`, `plotly`, `seaborn`, `pyfiglet`):

```bash
pip install "viralscan[full]"
```

The native tools (`kb`, `kallisto`, `bustools`, `STAR`) still come from conda:

```bash
conda install -c bioconda -c conda-forge kb-python kallisto bustools star
viralscan doctor --profile full
viralscan data fetch             # optional: only for indices without a build manifest
```

If `snakemake` or `kb` is missing, a workflow run stops with an error that names
the `full` extra. Prefer Option 1 for the full workflow.

### Optional extras

Core ViralScan does not require any optional extras for standard viral
quantification or host-response analysis. To enable **pathway enrichment**
in `viralscan hostresponse --enrichment` (via `gget.enrichr`), install the
`enrichment` extra:

```bash
pip install "viralscan[enrichment]"
```

This pulls in `gget>=0.27`. All other subcommands work without it.

## Option 3 — Container (Docker / Singularity)

Pre-built containers bundle every dependency including external tools.

Use a versioned image from the release notes when one is available. The v3
container workflow is prepared; a verified v3 image and Docker/Apptainer parity
receipt are still pending.

A local release build needs the tested wheel in `dist/` and the committed
`environment.explicit.txt` lock:

```bash
docker build --build-context wheel=./dist -t viralscan:release-test .
```

The image replays that lock and installs the supplied wheel. It does not build
the package from the checkout. For Apptainer, generate a definition from the
tested image's full digest:

```bash
python scripts/release_sif_definition.py \
  ghcr.io/mdmanurung/viralscan@sha256:<full-digest> > Singularity.release.def
apptainer build viralscan.sif Singularity.release.def
```

The root `Singularity.def` is a template and requires that digest.

## Verify the installation

```bash
viralscan --help
viralscan data fetch --help
viralscan build-ref --help
kb info
snakemake --version
```

You should see help text or version output for each command. If `kb` or
`snakemake` is missing, activate the conda environment or install those tools
before running a workflow.

## Download the bundled viral reference panel

ViralScan's default viral GTF annotation panel is distributed separately on
Zenodo to keep the Python package small. A run needs it only when its index has
no `<index>.build_manifest.json`: with a manifest (every index built by
`viralscan build-ref`) the viral gene list comes from the manifest and no fetch is
required. For an index without one, fetch the panel once:

```bash
viralscan data fetch
```

This downloads the archive for DOI `10.5281/zenodo.20112332`, verifies the
Zenodo checksum, and unpacks the GTF files under
`~/.cache/viralscan/data/`. Custom `-gtf` workflows do not require this cache.

Use a shared cache on HPC systems when many users or jobs should reuse the
same panel:

```bash
viralscan data fetch --cache-dir /shared/viralscan-cache
```

Then point quantification runs at the same cache:

```bash
viralscan --data-cache-dir /shared/viralscan-cache ...
```

Alternatively, set `VIRALSCAN_CACHE=/shared/viralscan-cache` before both the
fetch and quantification commands.
