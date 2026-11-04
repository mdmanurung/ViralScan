#!/usr/bin/env python3
"""Run the tiny ViralScan workflow in one environment and compare environments (REL-07, 09, 10).

Three subcommands:

``versions``  Print the tool manifest (kb, kallisto, bustools, minimap2, samtools,
              snakemake, STAR, viralscan) read from inside the environment.
``run``       Run the identical tiny workflow (``tests/data/evidence_tiny``) in one
              environment and write ``result_<name>.json``: commands, versions, digest,
              validated molecule counts and a SHA-256 for every output file.
``compare``   Read two or more ``result_*.json`` files and write
              ``parity_report.json``. Exit 1 when molecule counts or any deterministic
              output hash differ.

An environment is a command prefix (``--runner``) that runs a program inside it:
empty for the current shell, ``micromamba run -n NAME``, ``apptainer exec --bind
$PWD x.sif`` or ``docker run --rm -v $PWD:$PWD -w $PWD --entrypoint env IMAGE
PATH=/opt/conda/envs/viralscan/bin:/usr/bin:/bin``. The work directory and the
repository must be visible at the same path inside the environment.

Run: python scripts/release_parity.py run --name conda --kind conda --out analysis/release_parity
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "data" / "evidence_tiny"

#: Tool -> command that prints its version. ``kb --version`` prints usage and ``kb info`` exits
#: non-zero, so read the kb-python distribution metadata instead.
TOOL_VERSION_COMMANDS = {
    "kb": [
        "python",
        "-c",
        "import importlib.metadata as m; print('kb_python', m.version('kb-python'))",
    ],
    "kallisto": ["kallisto", "version"],
    "bustools": ["bustools", "version"],
    "minimap2": ["minimap2", "--version"],
    "samtools": ["samtools", "--version"],
    "snakemake": ["snakemake", "--version"],
    "STAR": ["STAR", "--version"],
}

#: Output files whose bytes may differ between environments or between two runs in one
#: environment. Each entry is a path suffix from the sample directory (or the run root) and
#: the reason. They are listed in the report and compared by content elsewhere
#: (molecule counts), never by hash. Everything else must hash identically.
EXPLAINED_NONDETERMINISTIC = {
    "run_manifest.json": "start time, absolute paths and tool paths",
    "run_complete.json": "completion time and absolute paths",
    "README.txt": "run fingerprint, which includes the absolute output path",
    "chemistry_sanity.json": "SHA-256 of kb run_info.json, which holds absolute paths",
    "artifact_inventory.json": "SHA-256 of other explained files (h5ad, run_info.json)",
    "config.yaml": "absolute input and output paths",
    "summary.txt": "absolute paths and wall-clock timings",
    "report.html": "generation time and absolute paths",
    "kb_info.json": "kb call with absolute paths, run time and tool paths",
    "run_info.json": "kallisto/bustools call with absolute paths and run time",
    "inspect.json": "absolute paths",
    "reference_provenance.json": "absolute paths and timestamps",
    "adata_multimap.h5ad": "HDF5 container metadata; compared as molecule counts instead",
    ".h5ad": "HDF5 container metadata; compared as molecule counts instead",
    ".done": "completion-marker timestamps",
}
#: Whole directories that hold logs or scratch files.
EXPLAINED_NONDETERMINISTIC_DIRS = {
    "log": "timestamps and absolute paths",
    ".snakemake": "Snakemake scheduler state",
}


def explained_reason(relative: str) -> str | None:
    """Return why a relative path may differ between environments, or None."""
    parts = relative.split("/")
    for part in parts[:-1]:
        if part in EXPLAINED_NONDETERMINISTIC_DIRS:
            return EXPLAINED_NONDETERMINISTIC_DIRS[part]
    for suffix, reason in EXPLAINED_NONDETERMINISTIC.items():
        if relative == suffix or relative.endswith(
            suffix if suffix.startswith(".") else f"/{suffix}"
        ):
            return reason
    return None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_tree(root: Path) -> dict[str, str]:
    """SHA-256 of every file under ``root``, keyed by the path relative to it."""
    return {
        path.relative_to(root).as_posix(): sha256_file(path)
        for path in sorted(root.rglob("*"))
        if path.is_file() and not path.is_symlink()
    }


def first_line(text: str) -> str:
    return next((line.strip() for line in text.splitlines() if line.strip()), "")


def execute(runner: list[str], command: list[str], cwd: Path, log: list[dict]) -> str:
    """Run ``runner + command``; record it; return stdout."""
    full = [*runner, *command]
    result = subprocess.run(
        full, cwd=cwd, capture_output=True, text=True, check=False, timeout=3600
    )
    log.append({"command": shlex.join(full), "returncode": result.returncode})
    if result.returncode:
        sys.stderr.write(f"FAILED: {shlex.join(full)}\n{result.stdout}\n{result.stderr}\n")
        raise SystemExit(result.returncode)
    return result.stdout + result.stderr


def tool_versions(runner: list[str], viralscan: list[str], cwd: Path, log: list[dict]) -> dict:
    versions = {
        tool: first_line(execute(runner, command, cwd, log))
        for tool, command in TOOL_VERSION_COMMANDS.items()
    }
    versions["viralscan"] = first_line(execute(runner, [*viralscan, "--version"], cwd, log))
    return versions


def viral_molecules(sample: Path) -> dict:
    """Validated molecule totals from the run's own summary tables."""
    import csv

    with (sample / "results" / "viral_summary.tsv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        "viral_molecules_total_est": sum(float(r["viral_molecules_total_est"]) for r in rows),
        "viruses_detected": sorted(
            r["virus_name"] for r in rows if r["detected"].lower() == "true"
        ),
    }


def run_tiny(args: argparse.Namespace) -> dict:
    runner = shlex.split(args.runner)
    viralscan = shlex.split(args.viralscan)
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    log: list[dict] = []
    versions = tool_versions(runner, viralscan, work, log)

    index = work / "index.idx"
    out = work / "out"
    # The index must be built by the kallisto binary that ``kb count`` will read it with
    # (REL-16): kb-python bundles its own, and a different version reports "Corrupted index".
    kallisto = first_line(
        execute(
            runner,
            [
                "python", "-c",
                "from viralscan.validation import tool_path; print(tool_path('kallisto'))",
            ],
            work,
            log,
        )
    )  # fmt: skip
    versions["kallisto_used_by_kb_count"] = first_line(
        execute(runner, [kallisto, "version"], work, log)
    )
    execute(
        runner,
        [kallisto, "index", "-i", str(index), str(FIXTURE / "transcripts.fasta")],
        work,
        log,
    )
    execute(
        runner,
        [
            *viralscan, "run", "-o", str(out), "-i", str(index),
            "-t", str(FIXTURE / "t2g.txt"),
            "--r1", str(FIXTURE / "R1.fastq"), "--r2", str(FIXTURE / "R2.fastq"),
            "--gtf", str(FIXTURE / "viral.gtf"),
            "-x", "10xv3", "-c", "1", "--no-visual", "--tchem", "off", "--cell-calling", "none", "--yes",
        ],
        work,
        log,
    )  # fmt: skip
    report_path = work / "validate_run.json"
    execute(
        runner,
        [*viralscan, "validate-run", "--run-dir", str(out), "--json", str(report_path)],
        work,
        log,
    )
    validation = json.loads(report_path.read_text(encoding="utf-8"))
    samples = sorted(p.parent for p in out.glob("*/config.yaml"))
    if len(samples) != 1:
        raise SystemExit(f"Expected one sample directory under {out}, found {samples}")

    result = {
        "name": args.name,
        "kind": args.kind,
        "digest": args.digest,
        "versions": versions,
        "commands": log,
        "validate_run_ok": validation["ok"] is True and validation["issues"] == [],
        "molecules": viral_molecules(samples[0]),
        "files": hash_tree(out),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    target = args.out / f"result_{args.name}.json"
    target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {target}")
    return result


def compare(results: list[dict]) -> dict:
    """Compare per-environment results; the first one is the reference.

    A file is deterministic unless ``explained_reason`` names it. Deterministic files must be
    present in every environment with one hash. Molecule counts and the validation verdict must
    agree. Explained files are listed with their per-environment hashes.
    """
    if len(results) < 2:
        raise ValueError("Need results from at least two environments to compare")
    names = [r["name"] for r in results]
    if len(set(names)) != len(names):
        raise ValueError(f"Environment names must be unique: {names}")
    reference = results[0]
    problems: list[str] = []

    for result in results:
        if not result["validate_run_ok"]:
            problems.append(f"{result['name']}: validate-run reported issues")
        if result["molecules"] != reference["molecules"]:
            problems.append(
                f"{result['name']}: molecules {result['molecules']} differ from "
                f"{reference['name']} {reference['molecules']}"
            )

    every_path = sorted({p for r in results for p in r["files"]})
    deterministic: dict[str, str] = {}
    explained: dict[str, dict] = {}
    for path in every_path:
        hashes = {r["name"]: r["files"].get(path) for r in results}
        reason = explained_reason(path)
        if reason is not None:
            explained[path] = {"reason": reason, "sha256": hashes}
            continue
        if None in hashes.values():
            missing = [n for n, h in hashes.items() if h is None]
            problems.append(f"{path}: missing in {missing}")
        elif len(set(hashes.values())) > 1:
            problems.append(f"{path}: hashes differ {hashes}")
        else:
            deterministic[path] = reference["files"][path]

    return {
        "schema_version": 1,
        "ok": not problems,
        "problems": problems,
        "reference": reference["name"],
        "environments": [
            {
                "name": r["name"],
                "kind": r["kind"],
                "digest": r["digest"],
                "versions": r["versions"],
                "commands": r["commands"],
                "validate_run_ok": r["validate_run_ok"],
                "molecules": r["molecules"],
            }
            for r in results
        ],
        "deterministic_files": deterministic,
        "explained_nondeterministic_files": explained,
        "explained_rules": {
            **EXPLAINED_NONDETERMINISTIC,
            **{f"{d}/": why for d, why in EXPLAINED_NONDETERMINISTIC_DIRS.items()},
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    def add_env_args(p: argparse.ArgumentParser) -> None:
        p.add_argument("--runner", default="", help="command prefix that enters the environment")
        p.add_argument("--viralscan", default="viralscan", help="viralscan command inside it")

    versions = sub.add_parser("versions", help="print the tool manifest as JSON")
    add_env_args(versions)
    versions.add_argument("--digest", default="", help="image digest to record")
    versions.add_argument("--output", type=Path, help="write JSON here instead of stdout")

    run = sub.add_parser("run", help="run the tiny workflow in one environment")
    add_env_args(run)
    run.add_argument("--name", required=True, help="short environment name, e.g. conda")
    run.add_argument("--kind", required=True, choices=["conda", "docker", "apptainer"])
    run.add_argument("--digest", default="", help="image or SIF digest to record")
    run.add_argument("--work", type=Path, required=True, help="scratch directory for this run")
    run.add_argument("--out", type=Path, required=True, help="directory for result_<name>.json")

    cmp_ = sub.add_parser("compare", help="write parity_report.json from result files")
    cmp_.add_argument("results", nargs="+", type=Path)
    cmp_.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()
    if args.command == "versions":
        log: list[dict] = []
        manifest = {
            "digest": args.digest,
            "versions": tool_versions(
                shlex.split(args.runner), shlex.split(args.viralscan), Path.cwd(), log
            ),
        }
        text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(text, encoding="utf-8")
        else:
            print(text, end="")
        return 0
    if args.command == "run":
        run_tiny(args)
        return 0
    report = compare([json.loads(p.read_text(encoding="utf-8")) for p in args.results])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}: ok={report['ok']}")
    for problem in report["problems"]:
        print(f"  {problem}")
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
