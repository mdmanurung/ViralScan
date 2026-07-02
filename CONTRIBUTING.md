# Contributing to ViralScan

Thanks for your interest in improving ViralScan. This guide is for human
contributors. (`CLAUDE.md` in the repo root is guidance for AI coding agents and is
not required reading.)

## Development environment

ViralScan is a Snakemake/`kb-python` CLI. The runtime tools (`kb`, `snakemake`,
`kallisto`, `bustools`, `STAR`) come from conda:

```bash
conda env create -f environment.yml
conda activate viralscan
python -m pip install -e ".[dev]"      # editable install + dev tools
```

The `[dev]` extra provides `pytest`, `ruff`, `mypy`, and `pre-commit`. Install the
git hooks once:

```bash
pre-commit install
```

If `pip install -e .` fails while building `connection_pool` (a `snakemake`
transitive dependency) on older `setuptools`, make sure `snakemake` is provided by
conda first (the `environment.yml` flow does this), or run the test suite via
`PYTHONPATH=src` without an editable install (see below).

## Running the tests

```bash
# Fast unit suite (default: excludes network + integration markers)
PYTHONPATH=src pytest tests/ -q

# Network tests (hit the live NCBI API)
PYTHONPATH=src pytest -m network

# Integration tests (need kb/snakemake/minimap2/samtools/blastn on PATH)
PYTHONPATH=src pytest -m integration
```

`PYTHONPATH=src` lets the tests import the package without an editable install,
which is the most robust path across environments.

## Quality gates

Before opening a PR, make sure these pass (CI enforces them):

```bash
ruff check .
ruff format --check .
mypy src/viralscan
```

- Code style: `ruff` (line length 100) + `ruff format`.
- Type checking: `mypy` (strict on the core modules).
- Conventional Commits (`feat:`, `fix:`, `docs:`, `test:`, `chore:`) are preferred.
- Keep `CHANGELOG.md` updated under `## [Unreleased]`.

## Pull requests

1. Branch off `main`.
2. Add tests for new behavior; keep the suite green.
3. Update `CHANGELOG.md` and any affected docs (`README.md`, `docs/`).
4. Open the PR against `main`; CI runs lint + the test matrix.

## Reporting issues

Use the GitHub issue tracker:
<https://github.com/mdmanurung/ViralScan/issues>. For bug reports, include the
ViralScan version (`viralscan --version`), the command you ran, and the full error.

## License

By contributing, you agree that your contributions are licensed under the project's
MIT License (`LICENSE`).
