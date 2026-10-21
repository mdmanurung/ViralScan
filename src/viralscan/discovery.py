"""Kraken2 cell-resolved taxonomic discovery on a completed run (PLAN DISC).

Design and contracts: ``docs/plans/2026-10-21-kraken2-discovery-plan.md``. Discovery is a
post-hoc layer: it reads a completed run and writes only into its own output directory, so it
cannot change primary inputs, count matrices, ``detected``, tiers or the run fingerprint. A
Kraken2 assignment is a taxonomic observation, never an infection call.

Implemented so far: run, output and database validation, atomic publication (DISC-02) and
host-unmapped input resolution (DISC-03). The classification stages land in DISC-06..09.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import shutil
from pathlib import Path
from typing import Any

from viralscan.run_safety import RUN_COMPLETE, RUN_MANIFEST, sha256_file

log = logging.getLogger(__name__)

#: The three files every Kraken2 database has; anything else (library/, taxonomy/) is optional.
KRAKEN2_DB_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")
#: ``kraken2`` classifies; ``k2 inspect`` gives the database's own taxonomy (DISC-04).
REQUIRED_TOOLS = ("kraken2", "k2")
#: Staging directory inside the output; published files are moved out of it with os.replace.
PARTIAL = ".partial"


class DiscoveryError(RuntimeError):
    """Discovery refused to run, or stopped, for a stated reason."""


def check_tools() -> None:
    missing = [tool for tool in REQUIRED_TOOLS if shutil.which(tool) is None]
    if missing:
        raise DiscoveryError(
            f"viralscan discover requires {', '.join(missing)} on PATH (Kraken2 >= 2.17). "
            "Kraken2 is optional and not part of the main ViralScan environment."
        )


def resolve_sample_dir(run_root: Path, sample: str | None) -> Path:
    """The sample directory of a *completed* run; refuses running or partial runs."""
    if not run_root.is_dir():
        raise DiscoveryError(f"Run directory does not exist: {run_root}")
    for marker in (RUN_MANIFEST, RUN_COMPLETE):
        if not (run_root / marker).is_file():
            raise DiscoveryError(
                f"{run_root} is not a completed ViralScan run ({marker} missing). "
                "Point --run-dir at the run root and let the run finish first."
            )
    samples = sorted(path.parent.name for path in run_root.glob("*/config.yaml"))
    if not samples:
        raise DiscoveryError(f"No sample directories (*/config.yaml) under {run_root}.")
    if sample is None:
        if len(samples) > 1:
            raise DiscoveryError(
                f"The run has several samples; choose one with --sample: {samples}"
            )
        sample = samples[0]
    elif sample not in samples:
        raise DiscoveryError(f"Sample {sample!r} not in run {run_root}; available: {samples}")
    return run_root / sample


def check_output_dir(output: Path, run_root: Path) -> None:
    """Refuse an output inside the run, or one that already holds results."""
    if output == run_root or run_root in output.parents:
        raise DiscoveryError(
            f"--output must be outside --run-dir: {output} is inside {run_root}. "
            "Discovery never writes into a primary run."
        )
    if output.exists() and any(p.name != PARTIAL for p in output.iterdir()):
        raise DiscoveryError(f"--output must be empty or absent: {output}")


def validate_kraken2_db(db: Path) -> None:
    """Fail closed on a missing or incomplete database; nothing is ever downloaded."""
    if not db.is_dir():
        raise DiscoveryError(f"Kraken2 database directory does not exist: {db}")
    bad = [
        name
        for name in KRAKEN2_DB_FILES
        if not (db / name).is_file() or (db / name).stat().st_size == 0
    ]
    if bad:
        raise DiscoveryError(f"Kraken2 database {db} is missing or has empty {', '.join(bad)}.")


def stage(output: Path) -> Path:
    """A clean staging directory; a leftover from an interrupted run is discarded."""
    partial = output / PARTIAL
    if partial.exists():
        shutil.rmtree(partial)
    partial.mkdir(parents=True)
    return partial


def publish(output: Path) -> None:
    """Move every staged file into *output* (same filesystem, so each move is atomic)."""
    partial = output / PARTIAL
    for path in sorted(partial.iterdir()):
        os.replace(path, output / path.name)
    partial.rmdir()


def _host_filter_audit(path: Path) -> dict[str, str]:
    """``category -> fragments`` column of ``host_filter_audit.tsv``."""
    with path.open(newline="") as handle:
        return {row["category"]: row["fragments"] for row in csv.DictReader(handle, delimiter="\t")}


def _lineage_retained(path: Path) -> int:
    with gzip.open(path, "rt", newline="") as handle:
        return sum(
            row["filter_decision"] == "retained" for row in csv.DictReader(handle, delimiter="\t")
        )


def read_host_filtered(hf: Path) -> dict[str, Any] | None:
    """A STAR host-unmapped pair in the host-filter layout, cross-checked against its own audit.

    Returns None when *hf* holds none of the files (the primary ran without ``--host-filter``).
    The primary manifest does not hash ``host_filtered/``, so the pair can only be hashed now
    (``primary_verified: False``).
    """
    files = {
        "r1": hf / "R1.fastq.gz",
        "r2": hf / "R2.fastq.gz",
        "audit": hf / "host_filter_audit.tsv",
        "lineage": hf / "fragment_lineage.tsv.gz",
    }
    present = {key for key, path in files.items() if path.is_file()}
    if not present:
        return None
    if present != set(files):
        missing = sorted(str(files[key]) for key in set(files) - present)
        raise DiscoveryError(f"Incomplete primary host filter output; missing {missing}.")

    from viralscan.scripts.host_filter import validate_paired_fastq_ids

    try:
        n_fragments = validate_paired_fastq_ids(str(files["r1"]), str(files["r2"]))
    except (ValueError, OSError, EOFError) as exc:
        raise DiscoveryError(f"Primary host-filtered FASTQs are unreadable: {exc}") from exc
    audit = _host_filter_audit(files["audit"])
    counts = {
        "fastq_pairs": n_fragments,
        "audit_retained": int(audit.get("retained_host_unmapped", -1)),
        "lineage_retained": _lineage_retained(files["lineage"]),
    }
    if len(set(counts.values())) != 1:
        raise DiscoveryError(
            f"Primary host filter outputs disagree on retained fragments: {counts}"
        )
    return {
        "source": "primary_host_filtered",
        "r1": str(files["r1"]),
        "r2": str(files["r2"]),
        "n_fragments": n_fragments,
        "star_param_set": audit.get("star_param_set"),
        "star_args": {
            k.split(":", 1)[1]: v for k, v in audit.items() if k.startswith("star_param:")
        },
        "sha256": {key: sha256_file(files[key]) for key in ("r1", "r2")},
        "primary_verified": False,
    }


def discovery_star_pass(
    run_root: Path,
    sample_dir: Path,
    out_dir: Path,
    host_index: Path,
    param_set: str,
    threads: int,
) -> dict[str, Any]:
    """Discovery-only STAR host subtraction of the primary's sample FASTQs into *out_dir*.

    Reuses the primary host filter (same STARsolo command and parameter sets); the primary's
    own inputs and outputs are only read. The sample FASTQs must match the primary manifest.
    """
    from viralscan.runconfig import RunConfig
    from viralscan.scripts.host_filter import _starsolo_filter, check_host_filter_tools

    config = RunConfig.from_yaml(str(sample_dir / "config.yaml"))
    manifest = json.loads((run_root / RUN_MANIFEST).read_text(encoding="utf-8"))
    recorded = {
        path: manifest.get("input_fingerprints", {}).get(key)
        for key, path in (manifest.get("input_paths") or {}).items()
        if key.startswith(("sample1:", "sample2:"))
    }
    for fastq in (config.sample1, config.sample2):
        path = Path(fastq).resolve()
        if not path.is_file():
            raise DiscoveryError(
                f"Sample FASTQ {path} is absent (deleted after the run?); a discovery-only "
                "host subtraction needs the original reads."
            )
        if str(path) not in recorded:
            raise DiscoveryError(f"{path} is not among the primary manifest's sample inputs.")
        if sha256_file(path) != recorded[str(path)]:
            raise DiscoveryError(f"{path} changed since the primary run (sha256 differs).")
    try:
        check_host_filter_tools("starsolo")
    except RuntimeError as exc:
        raise DiscoveryError(str(exc)) from exc
    if not host_index.is_dir():
        raise DiscoveryError(f"--host-index is not a directory: {host_index}")

    out_dir.mkdir(parents=True)
    r1, r2 = out_dir / "R1.fastq.gz", out_dir / "R2.fastq.gz"
    _starsolo_filter(
        config.sample1,
        config.sample2,
        str(host_index),
        config.technology,
        config.whitelist,
        out_dir,
        str(r1),
        str(r2),
        threads,
        param_set,
    )
    resolved = read_host_filtered(out_dir)
    if resolved is None:
        raise DiscoveryError(f"STAR host subtraction wrote no output to {out_dir}.")
    resolved.update(source="discovery_star_pass", host_index=str(host_index))
    return resolved


def resolve_input(
    run_root: Path, sample_dir: Path, partial: Path, args: argparse.Namespace
) -> dict[str, Any]:
    """Pick the host-unmapped pair (plan §5.2): the primary's, else an explicit STAR pass."""
    primary = read_host_filtered(sample_dir / "host_filtered")
    if primary is not None:
        if args.host_subtract:
            raise DiscoveryError(
                "The primary run already has host_filtered/ output; drop --host-subtract "
                "(discovery reuses it)."
            )
        return primary
    if not args.host_subtract:
        raise DiscoveryError(
            "The primary run has no host_filtered/ output. Kraken2 discovery needs host-"
            "subtracted reads: pass --host-subtract starsolo --host-index <STAR genomeDir> for "
            "a discovery-only STAR pass (the primary run is not changed)."
        )
    if not args.host_index:
        raise DiscoveryError("--host-subtract starsolo requires --host-index.")
    return discovery_star_pass(
        run_root,
        sample_dir,
        partial / "host_subtraction",
        Path(args.host_index).resolve(),
        args.host_filter_star_params,
        args.threads,
    )


def run_discover(args: argparse.Namespace) -> None:
    run_root = Path(args.run_dir).resolve()
    output = Path(args.output).resolve()
    db = Path(args.kraken2_db).resolve()

    check_tools()
    sample_dir = resolve_sample_dir(run_root, args.sample)
    check_output_dir(output, run_root)
    validate_kraken2_db(db)
    partial = stage(output)
    resolved = resolve_input(run_root, sample_dir, partial, args)
    (partial / "input.json").write_text(json.dumps(resolved, indent=2, sort_keys=True) + "\n")
    log.info(
        "Discovery input: %s, %d host-unmapped fragments.",
        resolved["source"],
        resolved["n_fragments"],
    )
    raise DiscoveryError(
        "Classification stages are not implemented yet (PLAN DISC-06..09); the resolved input "
        f"is staged in {partial} and nothing was published."
    )
