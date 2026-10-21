"""Kraken2 cell-resolved taxonomic discovery on a completed run (PLAN DISC).

Design and contracts: ``docs/plans/2026-10-21-kraken2-discovery-plan.md``. Discovery is a
post-hoc layer: it reads a completed run and writes only into its own output directory, so it
cannot change primary inputs, count matrices, ``detected``, tiers or the run fingerprint. A
Kraken2 assignment is a taxonomic observation, never an infection call.

This module currently holds the DISC-02 skeleton: run, output and database validation plus
atomic publication. The classification stages land in DISC-03..09.
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
from pathlib import Path

from viralscan.run_safety import RUN_COMPLETE, RUN_MANIFEST

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


def run_discover(args: argparse.Namespace) -> None:
    run_root = Path(args.run_dir).resolve()
    output = Path(args.output).resolve()
    db = Path(args.kraken2_db).resolve()

    check_tools()
    sample_dir = resolve_sample_dir(run_root, args.sample)
    check_output_dir(output, run_root)
    validate_kraken2_db(db)
    log.info("Discovery inputs valid: sample %s, database %s.", sample_dir.name, db)
    raise DiscoveryError(
        "Classification stages are not implemented yet (PLAN DISC-03..09); inputs were "
        "validated and nothing was written."
    )
