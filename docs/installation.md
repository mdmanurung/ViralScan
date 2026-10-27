# Installation

ViralScan installs in two tiers.

- **pip tier** (`pip install ViralScan`, Python ≥ 3.10): `build-ref`, `evidence`,
  `validate-run`, `doctor --profile pip`. No Snakemake or native tools needed.
- **Full tier** (quantification workflow): Python 3.11 on `linux-64` (the canonical,
  tested toolchain), [kb-python](https://www.kallistobus.tools/) (`kb` command) and
  [Snakemake](https://snakemake.readthedocs.io/en/stable/) ≥ 7.0, plus kallisto, bustools and STAR.
  EmptyDrops cell calling additionally needs R with `DropletUtils`.

## Conda (recommended)

```bash
conda env create -f environment.yml
conda activate viralscan
python -m pip install -e .
viralscan doctor --profile full
```

`environment.yml` pins Python, `kb-python`, `snakemake`, kallisto, bustools, STAR,
samtools, minimap2, BLAST, CD-HIT and R/DropletUtils. The last two commands install
the checkout and check that every full-workflow tool is on `PATH`.

Optional extras, each installed with `python -m pip install -e '.[<extra>]'` in that
environment: `scautoqc` (SCTK and pybiomart for `--scautoqc`; the first opt-in run
queries Ensembl and caches the annotations in the sample's
`results/scautoqc_gene_annotations.tsv`) and `enrichment` (`gget>=0.27` for
`viralscan hostresponse --enrichment`). Standard quantification and host-response
analysis need neither.

## Pixi

[Pixi](https://pixi.sh) installs the same pins as `environment.yml` from the
repository's `pixi.toml` (linux-64):

```bash
curl -fsSL https://pixi.sh/install.sh | bash   # once
git clone https://github.com/mdmanurung/ViralScan && cd ViralScan
pixi install                                    # default environment
pixi run doctor                                 # viralscan doctor --profile full
pixi run viralscan --help
```

Other environments: `pixi install -e dev` (pytest, ruff, mypy; `pixi run -e dev test`),
`pixi install -e docs` (`pixi run -e docs docs` builds this site into `docs/_build/html`)
and `pixi install -e discovery` (adds Kraken2). Keep the checkout and the pixi cache on a
local disk: on a network filesystem the solve can stop with a uv cache lock timeout.
Commit `pixi.lock` once generated so every machine resolves identically.

## Containers

`Singularity.pixi.def` builds one image with ViralScan and the whole tool suite from
`pixi.toml`. It re-solves at build time, so it is a development and convenience image,
not a release artifact. Build and run it from the repository root; the build needs
`APPTAINER_TMPDIR` on a local disk and takes about 20 minutes:

```bash
export APPTAINER_TMPDIR=/tmp/$USER-apptainer && mkdir -p "$APPTAINER_TMPDIR"
apptainer build viralscan-pixi.sif Singularity.pixi.def
apptainer test viralscan-pixi.sif        # runs viralscan doctor --profile full
apptainer run --bind /path/to/data:/data viralscan-pixi.sif run \
  -t /data/t2g.txt -i /data/index.idx -o /data/out --r1 /data/a_1.fastq.gz --r2 /data/a_2.fastq.gz
```

Release images replay the committed `environment.explicit.txt` lock and install the
tested wheel instead of building from the checkout:

```bash
docker build --build-context wheel=./dist -t viralscan:release-test .
python scripts/release_sif_definition.py \
  ghcr.io/mdmanurung/viralscan@sha256:<full-digest> > Singularity.release.def
apptainer build viralscan.sif Singularity.release.def
```

The root `Singularity.def` is a template; the generator refuses a mutable tag and
needs the image's full digest. A versioned image is listed in the release notes when
one is available; a verified v3 image and Docker/Apptainer parity receipt are still
pending.

## pip only

```bash
pip install ViralScan                 # pip tier
pip install "viralscan[full]"         # adds snakemake, scanpy, scikit-learn, plotly, seaborn, pyfiglet
conda install -c bioconda -c conda-forge kb-python kallisto bustools star   # native tools
viralscan doctor --profile full
```

A workflow run stops with an error naming the `full` extra when `snakemake` or `kb`
is missing. `snakemake` has a transitive dependency (`connection_pool`) that fails
to build with recent setuptools; installing it from conda avoids that. Prefer the
conda route for the full workflow.

## The bundled viral reference panel

The default viral GTF annotation panel is distributed separately on Zenodo
(DOI `10.5281/zenodo.20112332`). A run needs it only when its index has no
`<index>.build_manifest.json`; every index from `viralscan build-ref` has one and
carries its own viral gene list. Custom `--gtf` and `--ncbi-accession` runs do not
need it either. When it is needed, fetch it once:

```bash
viralscan data fetch                                   # ~/.cache/viralscan/data/
viralscan data fetch --cache-dir /shared/viralscan-cache   # shared HPC cache
```

Point runs at a shared cache with `viralscan run --cache-dir /shared/viralscan-cache ...`
or set `VIRALSCAN_CACHE=/shared/viralscan-cache` for both commands. The command
verifies the Zenodo checksum before unpacking.
