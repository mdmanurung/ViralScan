"""DAG-structure tests for the ViralScan Snakefile.

The pure-Python ``TestRuleOrdering`` class runs in the default suite without
needing the snakemake binary.  ``TestHostFilterDag`` runs ``snakemake -n``
and is marked ``@pytest.mark.integration``.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from viralscan.evidence import have_tools

SNAKEFILE = Path(__file__).parent.parent / "src" / "viralscan" / "Snakefile"

# Minimal config values needed to keep the Snakefile shell-block template
# from raising KeyError during a dry-run (kb_count's shell references these).
_BASE_CONFIG = [
    "output=/tmp/vs_dag_test/",
    "index=/fake/index.idx",
    "transcripts=/fake/t2g.txt",
    "technology=10xv3",
    "whitelist=",
    "host_index=",
    "host_filter_aligner=",
    "cores=2",
]


def _configfile(directory: Path, items: list[str]) -> str:
    """Write the typed ``config.yaml`` the Snakefile now reads (MECH-C) and return its path.

    ``items`` are ``key=value`` strings as the old ``--config`` wire took them; an empty value
    becomes ``None``, as snakemake's own parsing made it. ``output`` points at ``directory/out/``,
    where ``config.yaml`` must already exist because no rule produces it any more.
    """
    out = directory / "out"
    out.mkdir(exist_ok=True)
    cfg: dict[str, object] = {}
    for item in [*_BASE_CONFIG, *items]:
        key, value = item.split("=", 1)
        value = value.replace("{out}", f"{out}/")  # lets a test name files the DAG produces
        cfg[key] = yaml.safe_load(value) if value else None
    cfg["output"] = f"{out}/"
    path = out / "config.yaml"
    path.write_text(yaml.dump(cfg))
    return str(path)


class TestRuleOrdering:
    """rule all must be the first rule in the Snakefile.

    Snakemake uses the first rule it sees as the default target.  When
    ``host_filter`` was defined before ``rule all`` (conditional on
    ``config.get("host_index")``), running with ``--host-filter`` caused
    Snakemake to use ``host_filter`` as the default target — stopping the
    pipeline after host_filter.done with exit 0, no viral_summary.tsv.
    This test locks the fix in place.
    """

    def test_rule_all_is_first_rule(self) -> None:
        text = SNAKEFILE.read_text()
        matches = list(re.finditer(r"^rule\s+(\w+)\s*:", text, re.MULTILINE))
        assert matches, "No rules found in Snakefile"
        first = matches[0].group(1)
        assert first == "all", (
            f"First rule is '{first}', expected 'all'. "
            "When host_filter is defined first it becomes the default Snakemake "
            "target, halting the pipeline before kb_count (PLAN S2)."
        )

    def test_kb_count_shell_never_templates_raw_whitelist_config(self) -> None:
        """An unset whitelist must not reach kb as the literal ``-w None``.

        Snakemake parses ``--config whitelist=`` as None, and ``{config[whitelist]:q}``
        renders that as ``None``; kb reads ``-w None`` as "bypass barcode error
        correction". The whitelist has to be resolved in Python (a rule param).
        """
        assert "{config[whitelist]" not in SNAKEFILE.read_text()

    def test_kb_count_shell_never_templates_raw_strand_config(self) -> None:
        """An unset strand must not reach kb as the literal ``--strand None``."""
        assert "{config[strand]" not in SNAKEFILE.read_text()

    def test_kb_count_lists_filtered_fastqs_as_inputs_when_host_filter_set(self) -> None:
        """_kb_count_inputs must yield R1/R2 filter-FASTQ paths — not only the .done file."""
        import ast
        import textwrap

        # Exec the real function object lifted out of the Snakefile, with a
        # stub ``config`` global — a behavioural check, not a substring grep.
        # (The Snakefile as a whole is not parseable Python: ``rule all:`` etc.
        # are Snakemake DSL, so the function block is extracted line-wise.)
        lines = SNAKEFILE.read_text().splitlines()
        start = next(i for i, line in enumerate(lines) if line.startswith("def _kb_count_inputs("))
        block = [lines[start]]
        for line in lines[start + 1 :]:
            if line and line[0] not in " \t" and not line.startswith("#"):
                break
            block.append(line)
        func_src = textwrap.dedent("\n".join(block))
        func = ast.parse(func_src).body[0]
        assert isinstance(func, ast.FunctionDef)
        config: dict[str, str] = {}
        ns: dict[str, object] = {"config": config}
        exec(compile(ast.Module(body=[func], type_ignores=[]), str(SNAKEFILE), "exec"), ns)
        kb_count_inputs = ns["_kb_count_inputs"]

        config.update({"output": "/out/", "host_index": "/host.idx"})
        inputs = kb_count_inputs(None)
        assert "/out/host_filtered/R1.fastq.gz" in inputs
        assert "/out/host_filtered/R2.fastq.gz" in inputs
        assert "/out/log/host_filter.done" in inputs

        config["host_index"] = ""
        inputs_without = kb_count_inputs(None)
        assert not any("host_filtered" in path for path in inputs_without)

        # DEF-01: the read filter's outputs gate kb_count too; "off" adds nothing.
        config["read_filter"] = "artefact"
        assert "/out/read_filtered/R2.fastq.gz" in kb_count_inputs(None)
        assert "/out/log/read_filter.done" in kb_count_inputs(None)
        config["read_filter"] = "off"
        assert not any("read_filtered" in path for path in kb_count_inputs(None))


@pytest.mark.integration
@pytest.mark.parametrize("whitelist", [None, "gzip"])
def test_kb_count_rerun_preserves_metadata_and_previous_counts(tmp_path, whitelist):
    """Run the actual rule twice and prove failures cannot replace the last good count tree."""
    import gzip
    import textwrap

    if not shutil.which("snakemake"):
        pytest.skip("snakemake is not on PATH")
    work = tmp_path / "workflow with spaces"
    work.mkdir()
    extra = ["kb_r1=/fake/R1.fastq.gz", "kb_r2=/fake/R2.fastq.gz"]
    if whitelist:
        onlist = work / "on list.txt.gz"
        with gzip.open(onlist, "wt") as handle:
            handle.write("ACGT\n")
        extra.append(f"whitelist={onlist}")
    configfile = Path(_configfile(work, extra))
    out = configfile.parent
    (out / "notes.txt").write_text("user notes\n")
    (out / "whitelist_plain.txt").write_text("user-owned existing whitelist\n")
    (out / "run_manifest.json").write_text('{"user_metadata": true}\n')
    preserved = {p.name: p.read_bytes() for p in out.iterdir() if p.is_file()}
    bin_dir = work / "bin"
    bin_dir.mkdir()
    fake_kb = bin_dir / "kb"
    fake_kb.write_text(
        "#!/usr/bin/env python\n"
        + textwrap.dedent("""\
        import os
        import pathlib
        import sys
        out = pathlib.Path(sys.argv[sys.argv.index('-o') + 1])
        out.mkdir(exist_ok=True)
        (out / 'output.bus').write_text(os.environ['FAKE_KB_GENERATION'])
        if os.environ.get('FAKE_KB_MODE') != 'missing':
            counts = out / 'counts_unfiltered'
            counts.mkdir(exist_ok=True)
            if os.environ.get('FAKE_KB_MODE') != 'empty':
                (counts / 'adata.h5ad').write_text(os.environ['FAKE_KB_GENERATION'])
        sys.exit(1 if os.environ.get('FAKE_KB_MODE') == 'error' else 0)
        """)
    )
    fake_kb.chmod(0o755)
    real_mv = shutil.which("mv")
    assert real_mv
    fake_mv = bin_dir / "mv"
    fake_mv.write_text(
        "#!/usr/bin/env python\n"
        + textwrap.dedent(f"""\
        import os
        import pathlib
        import signal
        import subprocess
        import sys
        subprocess.run([{real_mv!r}, *sys.argv[1:]], check=True)
        if os.environ.get('FAKE_KB_MODE') == 'cancel-promotion' and pathlib.Path(sys.argv[1]).name == 'kb-python':
            os.kill(os.getppid(), signal.SIGTERM)
        """)
    )
    fake_mv.chmod(0o755)
    command = [
        "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        str(out / "log/kb.done"),
        "--cores",
        "1",
        "--forcerun",
        "kb_count",
        "--configfile",
        str(configfile),
    ]
    base_env = {
        **os.environ,
        "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
        "PYTHONPATH": str(SNAKEFILE.parent.parent),
        "XDG_CACHE_HOME": str(work / "cache"),
    }
    for generation in ("first", "second"):
        result = subprocess.run(
            command,
            cwd=work,
            env={**base_env, "FAKE_KB_GENERATION": generation},
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert (out / "kb-python/counts_unfiltered/adata.h5ad").read_text() == generation
        assert (out / "kb-python/output.bus").read_text() == generation
    for mode in ("error", "missing", "empty", "cancel-promotion"):
        result = subprocess.run(
            command,
            cwd=work,
            env={**base_env, "FAKE_KB_GENERATION": "failed", "FAKE_KB_MODE": mode},
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0, result.stdout + result.stderr
        assert (out / "kb-python/counts_unfiltered/adata.h5ad").read_text() == "second"
        assert (out / "kb-python/output.bus").read_text() == "second"
        assert (out / "log/kb_count.log").exists()
    assert {name: (out / name).read_bytes() for name in preserved} == preserved
    assert not list(out.glob(".kb-count.*"))
    assert not list(out.glob(".kb-previous.*"))
    assert not list(out.glob(".whitelist.*"))


@pytest.mark.integration
class TestHostFilterDag:
    """Snakemake dry-run DAG tests.

    These require the ``snakemake`` binary on PATH and are excluded from the
    default pytest run (use ``-m integration`` to include them).
    """

    @pytest.fixture(autouse=True)
    def _directory(self, tmp_path: Path) -> None:
        self.directory = tmp_path

    def _dryrun(self, extra_config: list[str]) -> str:
        """Return combined stdout+stderr of ``snakemake -n``.

        Skips when the binary is absent. Every other integration test here skips
        on a missing tool; these two failed instead, which made the required WP1
        gate red on any machine without snakemake and so made a green gate
        unavailable rather than merely inconvenient.
        """
        missing = have_tools(["snakemake"])
        if missing:
            pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")

        cmd = [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--dryrun",
            "--quiet",
            "all",
            "--configfile",
            _configfile(self.directory, extra_config),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return result.stdout + result.stderr

    def test_no_host_filter_plans_full_pipeline(self) -> None:
        """Without host_index all core rules must appear in the dry-run plan."""
        output = self._dryrun(
            [
                "sample1=/fake/R1.fastq.gz",
                "sample2=/fake/R2.fastq.gz",
                "kb_r1=/fake/R1.fastq.gz",
                "kb_r2=/fake/R2.fastq.gz",
            ]
        )
        for rule in ("kb_count", "analysis", "multimap", "detection", "umap"):
            assert rule in output, (
                f"Rule '{rule}' missing from dry-run plan (no host filter). Full output:\n{output}"
            )

    def _kb_shell(self, extra_config: list[str]) -> str:
        """Return the rendered kb_count shell command from ``snakemake -n -p``."""
        missing = have_tools(["snakemake"])
        if missing:
            pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")
        cmd = [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "all",
            "--dryrun",
            "--printshellcmds",
            "--forcerun",
            "kb_count",
            "--configfile",
            _configfile(
                self.directory,
                [
                    "sample1=/fake/R1.fastq.gz",
                    "sample2=/fake/R2.fastq.gz",
                    "kb_r1=/fake/R1.fastq.gz",
                    "kb_r2=/fake/R2.fastq.gz",
                    *extra_config,
                ],
            ),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return result.stdout + result.stderr

    def test_unset_whitelist_renders_empty_not_none(self) -> None:
        """No -w given → WL is empty, so kb uses its packaged on-list and corrects."""
        output = self._kb_shell([])
        # Snakemake's `:q` renders "" as nothing, i.e. a bare `WL=` assignment.
        assert re.search(r"^\s*WL=\s*$", output, re.MULTILINE), output
        assert "WL=None" not in output, output

    def test_dropseq_without_w_bypasses_correction(self) -> None:
        """SW-21: no official on-list, so kb gets -w None instead of a knee allowlist."""
        output = self._kb_shell(["technology=dropseq"])
        assert re.search(r"^\s*WL=None\s*$", output, re.MULTILINE), output

    def test_set_whitelist_is_passed_through(self) -> None:
        output = self._kb_shell(["whitelist=/fake/wl.txt"])
        assert re.search(r"^\s*WL=/fake/wl\.txt\s*$", output, re.MULTILINE), output

    def test_set_strand_is_passed_to_both_kb_branches(self) -> None:
        output = self._kb_shell(["strand=reverse", "whitelist=/fake/wl.txt"])
        assert re.search(r"^\s*ST=reverse\s*$", output, re.MULTILINE), output
        # Both the -w and the no -w kb invocation carry the flag, guarded on ST.
        flagged = [ln for ln in output.splitlines() if ln.lstrip().startswith("kb count")]
        assert len(flagged) == 2, output
        assert all('${ST:+--strand "$ST"}' in ln for ln in flagged), flagged

    def test_unset_strand_renders_empty_not_none(self) -> None:
        output = self._kb_shell([])
        assert re.search(r"^\s*ST=\s*$", output, re.MULTILINE), output
        assert "ST=None" not in output, output

    def test_host_filter_plans_full_pipeline(self) -> None:
        """With host_index ALL rules through umap must appear — PLAN S2 regression guard.

        Before the fix, only host_filter appeared and the run
        exited 0 without producing viral_summary.tsv.
        """
        output = self._dryrun(
            [
                "sample1=/fake/R1.fastq.gz",
                "sample2=/fake/R2.fastq.gz",
                "host_index=/fake/host.idx",
                "host_filter_aligner=kallisto",
                "kb_r1=/tmp/vs_dag_test/host_filtered/R1.fastq.gz",
                "kb_r2=/tmp/vs_dag_test/host_filtered/R2.fastq.gz",
            ]
        )
        for rule in ("host_filter", "kb_count", "analysis", "multimap", "detection", "umap"):
            assert rule in output, (
                f"Rule '{rule}' missing from host-filter dry-run plan (PLAN S2). "
                f"Full output:\n{output}"
            )


@pytest.mark.integration
class TestReadFilterDag:
    """DEF-01: kb_count and anello_align must wait for the read filter's output."""

    def _edges(self, tmp_path: Path, extra_config: list[str]) -> set[tuple[str, str]]:
        missing = have_tools(["snakemake"])
        if missing:
            pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")
        for mate in ("R1", "R2"):
            (tmp_path / f"{mate}.fastq.gz").touch()
        cmd = [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "all",
            "--dag",
            "--forceall",
            "--configfile",
            _configfile(
                tmp_path,
                [
                    f"sample1={tmp_path}/R1.fastq.gz",
                    f"sample2={tmp_path}/R2.fastq.gz",
                    "read_filter=artefact",
                    "kb_r1={out}read_filtered/R1.fastq.gz",
                    "kb_r2={out}read_filtered/R2.fastq.gz",
                    *extra_config,
                ],
            ),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        assert result.returncode == 0, result.stdout + result.stderr
        labels = dict(re.findall(r'^\s*(\d+)\[label = "(\w+)', result.stdout, re.MULTILINE))
        return {
            (labels[a], labels[b])
            for a, b in re.findall(r"^\s*(\d+)\s*->\s*(\d+)", result.stdout, re.MULTILINE)
        }

    def test_kb_count_waits_for_read_filter_without_host_filter(self, tmp_path) -> None:
        edges = self._edges(tmp_path, [])
        assert ("read_filter", "kb_count") in edges
        assert not any("host_filter" in edge for edge in edges)

    def test_read_filter_runs_after_host_filter_and_before_anello_align(self, tmp_path) -> None:
        edges = self._edges(
            tmp_path,
            [
                "host_index=/fake/host",
                "host_filter_aligner=starsolo",
                "anello_align=true",
                "anello_index=/fake/anello_star",
            ],
        )
        assert ("host_filter", "read_filter") in edges
        assert ("read_filter", "kb_count") in edges
        assert ("read_filter", "anello_align") in edges


@pytest.mark.integration
class TestAutoEvidenceDag:
    """ANDET-04: the opt-in auto_evidence rule is in the DAG only when enabled."""

    def _plan(self, tmp_path: Path, extra: list[str]) -> str:
        missing = have_tools(["snakemake"])
        if missing:
            pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")
        cmd = [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--dryrun",
            "--quiet",
            "all",
            "--configfile",
            _configfile(
                tmp_path,
                [
                    "sample1=/fake/R1.fastq.gz",
                    "sample2=/fake/R2.fastq.gz",
                    "kb_r1=/fake/R1.fastq.gz",
                    "kb_r2=/fake/R2.fastq.gz",
                    *extra,
                ],
            ),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        assert result.returncode == 0, result.stdout + result.stderr
        return result.stdout + result.stderr

    def test_rule_absent_by_default(self, tmp_path) -> None:
        assert "auto_evidence" not in self._plan(tmp_path, [])

    def test_rule_runs_after_detection_when_enabled(self, tmp_path) -> None:
        output = self._plan(tmp_path, ["auto_evidence=true"])
        assert re.search(r"^\s*auto_evidence\s+1\s", output, re.MULTILINE), output


@pytest.mark.integration
class TestGeneProgramsDag(TestAutoEvidenceDag):
    """AUD-08: ``gene_programs`` must plan on a fresh run (it once listed a phantom input)."""

    def test_rule_absent_by_default(self, tmp_path) -> None:
        assert "gene_programs" not in self._plan(tmp_path, [])

    def test_rule_runs_after_detection_when_enabled(self, tmp_path) -> None:
        # a fresh output dir: no results/viral_summary.tsv exists, and no rule produces it
        output = self._plan(tmp_path, ["gene_programs=true"])
        assert re.search(r"^\s*gene_programs\s+1\s", output, re.MULTILINE), output
        assert re.search(r"^\s*detection\s+1\s", output, re.MULTILINE), output


@pytest.mark.integration
class TestHostCaptureDag(TestAutoEvidenceDag):
    """DLIST-02: the opt-in host_capture rule is in the DAG only when enabled."""

    def test_rule_absent_by_default(self, tmp_path) -> None:
        assert "host_capture" not in self._plan(tmp_path, [])

    def test_rule_runs_after_the_count_when_enabled(self, tmp_path) -> None:
        output = self._plan(tmp_path, ["host_capture=true"])
        assert re.search(r"^\s*host_capture\s+1\s", output, re.MULTILINE), output
        assert re.search(r"^\s*kb_count\s+1\s", output, re.MULTILINE), output


@pytest.mark.integration
class TestCpuBudgetThreads(TestAutoEvidenceDag):
    """total_budget_v1: every rule reserves what its native tool uses, capped at --cores."""

    def _threads(
        self, tmp_path: Path, cores: int, override: str | None, scheduler: int | None = None
    ) -> dict[str, int]:
        cmd = [
            "snakemake", "--snakefile", str(SNAKEFILE), "--dryrun", "all",
            "--cores", str(scheduler or cores),
            "--configfile",
            _configfile(
                tmp_path,
                [
                    "sample1=/fake/R1.fastq.gz", "sample2=/fake/R2.fastq.gz",
                    "kb_r1=/fake/R1.fastq.gz", "kb_r2=/fake/R2.fastq.gz",
                    f"cores={cores}", "auto_evidence=true", "host_capture=true",
                ],
            ),
        ]  # fmt: skip
        if override:
            cmd += ["--set-threads", override]
        if have_tools(["snakemake"]):
            pytest.skip("snakemake not on PATH")
        done = subprocess.run(cmd, capture_output=True, text=True, check=False)
        assert done.returncode == 0, done.stdout + done.stderr
        # Snakemake omits the "threads:" line for the default of 1.
        return {
            name: int((re.search(r"\n\s+threads: (\d+)", block) or [0, 1])[1])
            for name, block in re.findall(r"rule (\w+):\n(.*?)(?=\n\[|\Z)", done.stdout, re.S)
        }

    @pytest.mark.parametrize("cores", [1, 6, 8, 32])
    def test_reservations_follow_budget(self, tmp_path, cores) -> None:
        got = self._threads(tmp_path, cores, None)
        assert got["kb_count"] == min(8, cores)
        for name in ("multimap", "auto_evidence", "host_capture"):
            assert got[name] == cores, got

    def test_reservations_come_from_the_rule_expressions_not_the_scheduler_clamp(
        self, tmp_path
    ) -> None:
        # Snakemake clamps threads to --cores, so only a scheduler larger than the
        # configured budget shows what the expressions themselves request.
        got = self._threads(tmp_path, 3, None, scheduler=64)
        assert got["kb_count"] == 3  # a literal `threads: 8` would give 8
        assert got["multimap"] == 3
        assert self._threads(tmp_path, 32, None, scheduler=64)["kb_count"] == 8  # the cap

    def test_external_set_threads_is_honoured(self, tmp_path) -> None:
        got = self._threads(tmp_path, 6, "auto_evidence=2")
        assert got["auto_evidence"] == 2 and got["kb_count"] == 6


@pytest.mark.parametrize("script", ["multimap", "auto_evidence", "host_capture"])
def test_scripts_use_scheduler_reserved_threads(script) -> None:
    source = (SNAKEFILE.parent / "scripts" / f"{script}.py").read_text()
    assert "with_cores(snakemake.threads)" in source


def test_multimap_pins_numerical_libraries_at_runtime() -> None:
    """Run the script's import-time preamble the way Snakemake would: BLAS vars end up 1."""
    source = (SNAKEFILE.parent / "scripts" / "multimap.py").read_text()
    preamble = source[: source.index("import anndata")]
    names = ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]
    code = f"snakemake = object()\n{preamble}\nprint(*[os.environ[n] for n in {names!r}])"
    env = {**os.environ, **dict.fromkeys(names, "6")}  # what Snakemake exports for 6 threads
    done = subprocess.run(
        [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
    )
    assert done.stdout.split() == ["1"] * 4


def test_snakefile_kb_count_uses_rule_threads_not_a_literal() -> None:
    text = SNAKEFILE.read_text()
    assert text.count("kb count -t {threads} ") == 2 and "kb count -t 8" not in text
    assert "threads: bounded_workers(KB_WORKERS" in text
