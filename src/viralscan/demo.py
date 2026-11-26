"""``viralscan demo``: the bundled tiny dataset and a one-command check (PLAN DEMO-01).

The demo takes the two-read-pair fixture the golden test uses through the real ``viralscan run``
path, runs ``check run`` on the result and compares the deterministic outputs with the pinned
sha256 digests in ``data/demo/expected_digests.json``. That file is a copy of
``tests/integration/golden/tiny_run.json``; ``tests/test_demo.py`` fails when the two differ.

Exit codes of :func:`run_demo`: 0 every file matches and ``check run`` passes, 1 a mismatch or a
failed step, 2 a required tool is missing.
"""

import hashlib
import json
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path
from typing import Any

from viralscan import __version__
from viralscan.run_safety import RunSafetyError, require_fresh_output
from viralscan.validation import missing_tools, tool_path

DEFAULT_OUTPUT = "viralscan_demo"
DEMO_INPUTS = ("R1.fastq", "R2.fastq", "transcripts.fasta", "t2g.txt", "viral.gtf")
REPORT = "demo_report.json"
EXIT_TOOLS = 2

#: The options of the integration tiny run (``tests/integration/test_tiny_end_to_end.py``). The
#: fixture matches no real kit, so tchem is off and the strand is named; its infected cells hold
#: one viral molecule each, so the pre-2.0 infected rule keeps a detected virus. No default changes.
RUN_OPTIONS = (
    "-x",
    "10xv3",
    "--no-visual",
    "--cell-calling",
    "none",
    "--yes",
    "--strand",
    "unstranded",
    "--tchem",
    "off",
    "--infected-min-molecules",
    "1",
)


def demo_data() -> Path:
    return Path(str(resources.files("viralscan.data").joinpath("demo")))


def expected_digests() -> dict[str, str]:
    digests: dict[str, str] = json.loads(
        (demo_data() / "expected_digests.json").read_text(encoding="utf-8")
    )
    return digests


def output_digests(run: Path) -> dict[str, str]:
    """Map ``results/<file>``, ``summary.txt`` and ``combined_summary.tsv`` of one run to a sha256.

    Only files two identical runs reproduce byte for byte are listed: the result tables and gene
    lists, ``combined_summary.tsv`` and ``summary.txt`` without its ``Runtime:`` line.
    """
    samples = sorted(p.parent for p in run.glob("*/config.yaml"))
    if len(samples) != 1:
        raise ValueError(f"expected one sample directory in {run}, found {len(samples)}")
    sample = samples[0]
    digests: dict[str, str] = {}
    combined = run / "combined_summary.tsv"
    if combined.is_file():
        digests["combined_summary.tsv"] = hashlib.sha256(combined.read_bytes()).hexdigest()
    for path in sorted((sample / "results").iterdir()):
        if path.is_file() and path.suffix in {".tsv", ".txt"}:
            digests[f"results/{path.name}"] = hashlib.sha256(path.read_bytes()).hexdigest()
    summary = (sample / "summary.txt").read_text(encoding="utf-8").splitlines(keepends=True)
    kept = "".join(line for line in summary if not line.startswith("Runtime:"))
    digests["summary.txt"] = hashlib.sha256(kept.encode("utf-8")).hexdigest()
    return digests


def compare_digests(got: dict[str, str], want: dict[str, str]) -> dict[str, str]:
    """``PASS``, ``FAIL`` (digest differs), ``MISSING`` (not produced) or ``UNEXPECTED`` per file."""
    verdict = {}
    for name in sorted(got.keys() | want.keys()):
        if name not in got:
            verdict[name] = "MISSING"
        elif name not in want:
            verdict[name] = "UNEXPECTED"
        else:
            verdict[name] = "PASS" if got[name] == want[name] else "FAIL"
    return verdict


def _step(label: str, cmd: list[str], cwd: Path, log: Path) -> bool:
    print(f"demo: {label} ...", flush=True)
    with log.open("ab") as handle:
        handle.write(f"$ {' '.join(cmd)}\n".encode())
        done = subprocess.run(cmd, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if done.returncode:
        tail = log.read_text(encoding="utf-8", errors="replace").splitlines()[-15:]
        print(
            f"demo: {label} failed (exit {done.returncode}); last lines of {log}:", file=sys.stderr
        )
        print("\n".join(tail), file=sys.stderr)
    return done.returncode == 0


def run_demo(output: Path, *, cores: int = 2, keep: bool = False, force: bool = False) -> int:
    from viralscan.validation import tool_provenance, validate_run

    missing = missing_tools(("kb", "snakemake")) + ([] if tool_path("kallisto") else ["kallisto"])
    if missing:
        print(
            f"viralscan demo: required tools are not on PATH: {', '.join(missing)}. "
            "Run `viralscan check env --profile full` to see what is missing; the demo needs "
            "the full tier (docs/installation.md).",
            file=sys.stderr,
        )
        return EXIT_TOOLS

    output = output.resolve()
    try:
        require_fresh_output(output, force=force)
    except RunSafetyError as exc:
        print(f"viralscan demo: {exc}", file=sys.stderr)
        return EXIT_TOOLS
    inputs, run = output / "inputs", output / "run"
    for stale in (inputs, run, output / REPORT):  # only reached with --force
        if stale.is_dir():
            shutil.rmtree(stale)
        else:
            stale.unlink(missing_ok=True)
    inputs.mkdir(parents=True)
    for name in DEMO_INPUTS:
        shutil.copyfile(demo_data() / name, inputs / name)
    log = output / "demo.log"

    index = inputs / "index.idx"
    kallisto = tool_path("kallisto") or "kallisto"
    steps_ok = _step(
        "build the kallisto index",
        [kallisto, "index", "-i", str(index), str(inputs / "transcripts.fasta")],
        output,
        log,
    ) and _step(
        "viralscan run",
        [
            sys.executable, "-m", "viralscan", "run",
            "-o", str(run),
            "-i", str(index),
            "-t", str(inputs / "t2g.txt"),
            "--r1", str(inputs / "R1.fastq"),
            "--r2", str(inputs / "R2.fastq"),
            "--gtf", str(inputs / "viral.gtf"),
            "-c", str(cores),
            *RUN_OPTIONS,
        ],
        output,
        log,
    )  # fmt: skip

    verdict: dict[str, str] = {}
    check: dict[str, Any] = {"ok": False}
    if steps_ok:
        print("demo: viralscan check run ...", flush=True)
        check = validate_run(run, verify_inputs=True)
        verdict = compare_digests(output_digests(run), expected_digests())
    for name, state in verdict.items():
        print(f"  {state:<10}{name}")
    print(f"  {'PASS' if check['ok'] else 'FAIL':<10}check run")

    ok = steps_ok and bool(check["ok"]) and all(state == "PASS" for state in verdict.values())
    report = {
        "ok": ok,
        "viralscan_version": __version__,
        "files": verdict,
        "check_run": check,
        "tools": tool_provenance(),
    }
    (output / REPORT).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if ok and not keep:
        shutil.rmtree(inputs)
        shutil.rmtree(run)
        log.unlink(missing_ok=True)
    print(f"demo: {'PASS' if ok else 'FAIL'}; report {output / REPORT}")
    if ok and keep:
        print(f"demo: kept the run in {run}")
    return 0 if ok else 1
