# Tester packet

This page is for laboratories testing a ViralScan candidate build (3.0.0b1 or an rc). It takes about
15 minutes plus one real run. Tags `b1`, `rc` and `dev` publish to TestPyPI and GHCR only, so the
install below uses test.pypi.org. Only a stable tag (`vMAJOR.MINOR.PATCH`) reaches pypi.org.

## 1. Install

Create the full-tier environment first. It supplies Python 3.11, `kb`, kallisto, bustools, STAR and
Snakemake 9, which pip cannot install ([Installation](installation.md)).

```bash
conda env create -f environment.yml
conda activate viralscan
```

`environment.yml` is in the source checkout and in the sdist. Then install the candidate from
test.pypi.org. The second index lets pip find the dependencies, which TestPyPI does not carry.

```bash
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ 'ViralScan==3.0.0b1'
```

For an rc, replace the version (`ViralScan==3.0.0rc1`). Install tchem as described in
[Installation](installation.md#tchem-chemistry-second-opinion); `viralscan check env` warns, and does
not fail, when it is missing.

## 2. Check the environment

```bash
viralscan check env --profile full --json check_env.json
```

The command exits 0 when every full-tier tool is on `PATH`. Fix what it lists as missing before the
demo.

## 3. Run the demo

```bash
viralscan demo
```

The demo copies a two-read-pair dataset (about 2 KB, bundled in the package) into `./viralscan_demo`,
builds a two-sequence kallisto index, runs `viralscan run` and `viralscan check run`, and compares 16
deterministic output files with pinned sha256 digests. It prints one line per file:

```text
  PASS      results/viral_summary.tsv
  PASS      check run
demo: PASS; report viralscan_demo/demo_report.json
```

It takes a few minutes, most of it `kb` start-up. Exit codes: 0 all files match and `check run`
passes, 1 a mismatch or a failed step, 2 a tool is missing or the output directory already holds files.
After a FAIL the run stays in `viralscan_demo/run` and the log in `viralscan_demo/demo.log`. After a
PASS only `demo_report.json` stays; add `--keep` to keep the run. `--output DIR`, `--cores N` and
`--force` are described in the [CLI reference](cli_reference.md).

The demo names `--strand unstranded` and `--tchem off` because its reads match no real kit. A real run
chooses its own strand ([Quickstart](quickstart.md)).

## 4. Run one real sample

Follow the [Quickstart](quickstart.md) on one sample of your own. Pick the chemistry you use most.
When it finishes:

```bash
viralscan check run --run-dir <run root> --json check_run.json
viralscan status --run-dir <run root>
```

## 5. What to send back

- `viralscan_demo/demo_report.json`, or the full terminal output of `viralscan demo` and, after a FAIL,
  `viralscan_demo/demo.log`.
- `check_env.json` from step 2.
- `check_run.json` from step 4.
- The output of `viralscan status --run-dir <run root>` for the real run.
- The chemistry, host, operating system, and how you installed (conda, pixi or container).

Do not send FASTQ files, count matrices or barcodes. If a step fails, send its full error text and the
command you ran.
