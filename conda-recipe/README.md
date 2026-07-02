# Conda / bioconda recipe

`meta.yaml` is the [bioconda](https://bioconda.github.io/) recipe for `viralscan`.

## Submitting to bioconda

1. Publish the PyPI release first (the recipe sources the sdist from PyPI).
2. Fill the `source.sha256` in `meta.yaml` — the easiest way is to regenerate the
   recipe from the published sdist:

   ```bash
   pip install grayskull
   grayskull pypi ViralScan==2.4.0
   ```

   or copy the sha256 from the PyPI *Download files* page.
3. Fork [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes),
   copy this recipe to `recipes/viralscan/meta.yaml`, and open a PR. Bioconda CI
   builds and tests it (the `test:` block runs `viralscan --version` / `--help`).

## Building/testing locally

```bash
conda install -n base conda-build boa
conda build conda-recipe/ -c conda-forge -c bioconda
```

Notes:
- `noarch: python` — one build for all platforms.
- `snakemake-minimal` (not `snakemake`) keeps the closure small; `matplotlib-base`
  avoids the Qt GUI stack.
- The external binaries (`kb`, `kallisto`, `bustools`, `STAR`) come from the
  `kb-python` / `kallisto` / `bustools` / `star` run deps.
