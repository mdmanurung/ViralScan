# Installation

## Requirements

ViralScan installs in two tiers.

- **pip tier** (`pip install ViralScan`, Python ≥ 3.10): `build-ref`, `evidence`,
  `validate-run`, `doctor --profile pip`. No Snakemake or native tools needed.
- **Full tier** (quantification workflow): Python 3.11 on `linux-64` (the canonical,
  tested toolchain), [kb-python](https://www.kallistobus.tools/) (`kb` command) and
  [Snakemake](https://snakemake.readthedocs.io/en/stable/) ≥ 7.0, plus kallisto, bustools and STAR.

## Option 1 — Conda (recommended)

The repository provides a conda environment file:

```bash
conda env create -f environment.yml
conda activate viralscan
python -m pip install -e .
viralscan data fetch             # optional, see below
```

The environment contains Python, `kb-python`, `snakemake` and all Python dependencies.
The last command installs the local ViralScan checkout. `viralscan data fetch`
downloads the bundled viral annotation panel. A run needs the panel only when its
index has no `<index>.build_manifest.json`. An index from `viralscan build-ref`
has one and carries its own viral gene list.

## Option 2 — pip (pip tier)

```bash
pip install ViralScan
viralscan doctor --profile pip
```

This command does not install Snakemake. The `full` extra adds the Python workflow
dependencies (`snakemake`, `scanpy`, `scikit-learn`, `plotly`, `seaborn`, `pyfiglet`):

```bash
pip install "viralscan[full]"
```

The native tools (`kb`, `kallisto`, `bustools`, `STAR`) come from conda:

```bash
conda install -c bioconda -c conda-forge kb-python kallisto bustools star
viralscan doctor --profile full
viralscan data fetch             # optional: only for indices without a build manifest
```

A workflow run stops with an error that names the `full` extra when `snakemake` or
`kb` is missing. Option 1 is preferred for the full workflow.

### Optional extras

Standard viral quantification and host-response analysis need no optional extras.
**Pathway enrichment** in `viralscan hostresponse --enrichment` (via `gget.enrichr`)
needs the `enrichment` extra:

```bash
pip install "viralscan[enrichment]"
```

The extra installs `gget>=0.27`. All other subcommands work without it.

## Option 3 — Pixi

[Pixi](https://pixi.sh) installs the full tool suite (kb-python, kallisto, bustools,
STAR, samtools, Snakemake, minimap2, BLAST, CD-HIT, R/DropletUtils) from the same pins as
`environment.yml`, using the repository's `pixi.toml` (linux-64):

```bash
curl -fsSL https://pixi.sh/install.sh | bash   # once
git clone https://github.com/mdmanurung/ViralScan && cd ViralScan
pixi install                                    # default environment
pixi run doctor                                 # viralscan doctor --profile full
pixi run viralscan --help
```

Other environments: `pixi install -e dev` (pytest, ruff, mypy; `pixi run -e dev test`) and
`pixi install -e docs` (`pixi run -e docs docs` builds this site into `docs/_build/html`).
Commit `pixi.lock` once generated so every machine resolves identically.

## Option 4 — Singularity / Apptainer image

`Singularity.pixi.def` builds one image containing ViralScan and the whole tool suite from
`pixi.toml`:

```bash
apptainer build viralscan-pixi.sif Singularity.pixi.def   # from the repository root
apptainer run viralscan-pixi.sif --help
apptainer exec viralscan-pixi.sif kb info
apptainer run --bind /path/to/data:/data viralscan-pixi.sif \
  -t /data/t2g.txt -i /data/index.idx -o /data/out -s1 /data/a_1.fastq.gz -s2 /data/a_2.fastq.gz
```

This image re-solves its environment at build time. Use it for development and
convenience. For releases, build from a tested OCI digest (below), which replays the
committed lock.

## Option 5 — Container (Docker / Singularity release build)

Pre-built containers bundle every dependency including external tools.

A versioned image is listed in the release notes when one is available. The v3
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

Each command prints help text or version output. If `kb` or `snakemake` is missing,
activate the conda environment or install those tools before running a workflow.

## Download the bundled viral reference panel

The default viral GTF annotation panel is distributed separately on Zenodo to keep
the Python package small. A run needs it only when its index has no
`<index>.build_manifest.json`. With a manifest (every index built by
`viralscan build-ref`) the viral gene list comes from the manifest and no fetch is
required. For an index without one, fetch the panel once:

```bash
viralscan data fetch
```

This downloads the archive for DOI `10.5281/zenodo.20112332`, verifies the
Zenodo checksum, and unpacks the GTF files under
`~/.cache/viralscan/data/`. Custom `-gtf` workflows do not require this cache.

On HPC systems, a shared cache lets many users or jobs reuse the same panel:

```bash
viralscan data fetch --cache-dir /shared/viralscan-cache
```

Point quantification runs at the same cache:

```bash
viralscan --data-cache-dir /shared/viralscan-cache ...
```

Alternatively, set `VIRALSCAN_CACHE=/shared/viralscan-cache` before both the fetch
and quantification commands.
