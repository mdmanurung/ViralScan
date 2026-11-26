# From install to the report

This page is one linear path for a routine host-aware run: install, build a
reference, check the library, quantify, validate, and open the results. Every
`viralscan` command in its code blocks is parsed by `tests/test_docs_sources.py`
against the current CLI, so a renamed or removed option fails the test suite. The
parse does not run any workflow. For the reasons behind each step and the full
option list, follow the links to [Installation](installation.md),
[Run ViralScan](quickstart.md) and the [CLI reference](cli_reference.md).

You need paired-end FASTQ files (R1 carries the cell barcode and UMI), Linux
x86-64, and a host species and viral accessions of your choosing.

## 1. Install

```bash
conda env create -f environment.yml
conda activate viralscan
python -m pip install -e .
python scripts/install_tchem.py --prefix "$CONDA_PREFIX"
viralscan check env --profile full
```

`environment.yml` and `scripts/install_tchem.py` are in the repository checkout.
tchem v0.1.0 is a separate download with a checked digest; its barcode lists and
the other install routes (pixi, containers) are in
[Installation](installation.md#tchem-chemistry-second-opinion). `check env` exits
nonzero when a required tool is missing.

## 2. Build the reference

```bash
viralscan build ref \
  --host human \
  --virus-accessions NC_045512.2 NC_002021.3 \
  --output ref_human_virus/ \
  --ncbi-email you@example.org
```

This writes `ref_human_virus/index.idx` and `ref_human_virus/t2g.txt`, which the
next steps take as `-i` and `-t`.

**No preset flag for the bundled panel.** `build ref` has no option that selects
the bundled viral panel (the GTF files kept as the bundled cache, see
[Reference panel](reference_panel.md)). It takes the viruses you name with
`--virus-accessions`, so the command line above is the documented form. Choose the
accessions for your question. ViralScan does not pick a panel for you, and this
page changes no default. Two opt-in builds exist: `--anellovirus` adds the
expanded Anelloviridae table to the host reference, and `--profile broad-discovery`
selects the frozen broad-discovery profile label. Neither is a routine choice; see
the [CLI reference](cli_reference.md). An index built here carries its own viral
gene list, so `viralscan data fetch` is not needed for this path.

## 3. Check the library

```bash
viralscan check chemistry \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  -i ref_human_virus/index.idx \
  -t ref_human_virus/t2g.txt
```

The command is read-only. It prints the `-x`, `-w` and `--strand` values it would
use and exits 1 when the chemistry is unresolved or the signals conflict. Fix that
before a long run.

## 4. Quantify

```bash
viralscan run \
  -t ref_human_virus/t2g.txt \
  -i ref_human_virus/index.idx \
  -o output/ \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --kb-params auto
```

`--kb-params auto` fills `-x` and `--strand` from the built-in detector, tchem and
the strand pilot, and stops when they disagree. Without tchem, pass `-x` and
`--strand` yourself. `run --help` lists the 16 essential options and
`run --help-all` the rest. The sample directory comes from the R1 file name before
the first underscore, so this run writes `output/sample/`.

## 5. Validate

```bash
viralscan check run --run-dir output/
```

A sound run prints `"ok": true` and exits 0. Validation checks the manifest, the
recorded input hashes, the completion marker and the count invariants. It does not
show that the reference covers your virus or that a call is infection.

## 6. What to open first

Open these three files, in this order. All paths are under `output/sample/`.

1. **`report.html`.** The sample README lists it as the place to start. It is
   self-contained and works offline.
2. **`results/viral_summary_core.tsv`.** One row per virus, 12 columns:
   molecule estimate, infected called cells, called-cell denominator and the
   `detected` flag. The full table is `results/viral_summary.tsv`.
3. **`results/sensitivity.tsv`.** The detection limit for each virus. Read it
   before you trust a negative: a virus that is absent or has `detected = false`
   was not detected, which does not show it is absent.

Then check what the run did and open it from Python:

```bash
viralscan status --run-dir output/
```

`status` lists the post-hoc steps run on the run so far. In Python, the same core table loads as a data frame:

```python
from viralscan import open_run

run = open_run("output")
run.summary(core=True)
```

A `detected` row is candidate evidence, not a diagnosis.
[Results](results.md) explains the units, the denominators and five checks before a
candidate becomes a claim.
