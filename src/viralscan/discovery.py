"""Kraken2 cell-resolved taxonomic discovery on a completed run (PLAN DISC).

Design and contracts: ``docs/plans/2026-10-21-kraken2-discovery-plan.md``. Discovery is a
post-hoc layer: it reads a completed run and writes only into its own output directory, so it
cannot change primary inputs, count matrices, ``detected``, tiers or the run fingerprint. A
Kraken2 assignment is a taxonomic observation, never an infection call.

Implemented so far: run, output and database validation, atomic publication (DISC-02),
host-unmapped input resolution (DISC-03) and the streaming identity pass (DISC-06). The
classification stages land in DISC-07..09.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import shutil
from array import array
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from viralscan.chemistry import cb_umi_geometry
from viralscan.evidence import _canonical_fastq_id, _open_maybe_gzip
from viralscan.run_safety import RUN_COMPLETE, RUN_MANIFEST, sha256_file

log = logging.getLogger(__name__)

#: The three files every Kraken2 database has; anything else (library/, taxonomy/) is optional.
KRAKEN2_DB_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")
#: ``kraken2`` classifies; ``k2 inspect`` gives the database's own taxonomy (DISC-04).
REQUIRED_TOOLS = ("kraken2", "k2")
#: Staging directory inside the output; published files are moved out of it with os.replace.
PARTIAL = ".partial"
#: 10xv1 keeps its UMI on a separate read, so R1[0:cb+umi] would mis-slice it.
UNSUPPORTED_CHEMISTRIES = frozenset({"10xv1"})
#: Kraken2 input: the cDNA mate (R2) of every fragment, named by its 0-based ordinal.
KRAKEN_INPUT = "kraken_input.fa.gz"
#: One row per fragment, keyed by ``ordinal`` (plan §6.1; later stages add their own tables).
READS_TSV = "reads.tsv.gz"
READS_COLUMNS = (
    "ordinal",
    "read_id",
    "raw_cb",
    "raw_umi",
    "geometry_ok",
    "cb_has_n",
    "umi_has_n",
    "primary_read_filter",
)


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


def barcode_geometry(technology: str) -> tuple[int, int]:
    """``(cb_len, umi_len)`` for chemistries with CB then UMI at R1 offset 0; others are refused."""
    if technology.strip().lower() in UNSUPPORTED_CHEMISTRIES:
        raise DiscoveryError(
            f"Chemistry {technology!r} is not supported by discovery (UMI not on R1)."
        )
    try:
        return cb_umi_geometry(technology)
    except ValueError as exc:  # ChemistryError is a ValueError
        raise DiscoveryError(f"Unsupported barcode geometry for discovery: {exc}") from exc


def _read_filter_rows(path: Path):
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            decision = row["filter_decision"]
            yield row["read_id"], decision if decision == "retained" else f"removed:{row['reason']}"


def stream_reads(
    r1: str,
    r2: str,
    cb_len: int,
    umi_len: int,
    out_dir: Path,
    read_filter_lineage: Path | None = None,
    read_filter_default: str = "not_run",
) -> dict[str, int]:
    """One pass over the pair: R2 to ``KRAKEN_INPUT`` as ``>ordinal``, identity to ``READS_TSV``.

    Every fragment is classified; barcode problems are handled later, by ordinal, so the
    Kraken2 output must cover ``[0, n_input)`` exactly once. ``N`` in CB/UMI is kept (kallisto
    and bustools encode it as ``G``, plan §6.3) and flagged. *read_filter_lineage* is the primary
    read filter's per-fragment lineage over this same pair; it is read in lockstep and must
    match every read ID. Duplicate read IDs, mate mismatches and truncation abort.
    """
    need = cb_len + umi_len
    counts: Counter[str] = Counter()
    # ponytail: 8-byte ID digests, 8 B/fragment (0.8 GB per 1e8); a false duplicate needs a
    # 64-bit collision (~3e-4 at 1e8 fragments) and fails closed. Shard by hash if RAM matters.
    digests = array("Q")
    lineage = _read_filter_rows(read_filter_lineage) if read_filter_lineage else None
    with (
        _open_maybe_gzip(r1) as fq1,
        _open_maybe_gzip(r2) as fq2,
        gzip.open(out_dir / KRAKEN_INPUT, "wt", compresslevel=1) as fasta,
        gzip.open(out_dir / READS_TSV, "wt", newline="", compresslevel=1) as reads,
    ):
        writer = csv.writer(reads, delimiter="\t")
        writer.writerow(READS_COLUMNS)
        ordinal = 0
        while True:
            rec1 = [fq1.readline() for _ in range(4)]
            rec2 = [fq2.readline() for _ in range(4)]
            if not rec1[0] and not rec2[0]:
                break
            if not all(rec1) or not all(rec2):
                raise DiscoveryError(f"Truncated or unpaired FASTQ record {ordinal} in {r1} / {r2}")
            read_id = _canonical_fastq_id(rec1[0])
            if read_id != _canonical_fastq_id(rec2[0]):
                raise DiscoveryError(f"FASTQ mate mismatch at record {ordinal}: {read_id!r}")
            digests.append(
                int.from_bytes(hashlib.blake2b(read_id.encode(), digest_size=8).digest(), "little")
            )
            seq1 = rec1[1].rstrip("\n")
            geometry_ok = len(seq1) >= need
            cb, umi = (seq1[:cb_len], seq1[cb_len:need]) if geometry_ok else ("", "")
            if lineage is None:
                annotation = read_filter_default
            else:
                step = next(lineage, None)
                if step is None or step[0] != read_id:
                    raise DiscoveryError(
                        f"Read-filter lineage is out of step with the input at record {ordinal}: "
                        f"{read_id!r} vs {step[0] if step else 'end of lineage'!r}"
                    )
                annotation = step[1]
            fasta.write(f">{ordinal}\n{rec2[1].rstrip()}\n")
            writer.writerow(
                [
                    ordinal,
                    read_id,
                    cb,
                    umi,
                    int(geometry_ok),
                    int("N" in cb),
                    int("N" in umi),
                    annotation,
                ]
            )
            counts["n_input"] += 1
            counts["n_bad_geometry"] += not geometry_ok
            counts["n_cb_has_n"] += "N" in cb
            counts["n_umi_has_n"] += "N" in umi
            counts[f"read_filter:{annotation.split(':')[0]}"] += 1
            ordinal += 1
    if lineage is not None and next(lineage, None) is not None:
        raise DiscoveryError("Read-filter lineage has more rows than the input pair.")
    ids = np.frombuffer(digests, dtype=np.uint64)
    ids.sort()
    if ids.size > 1 and bool((ids[1:] == ids[:-1]).any()):
        raise DiscoveryError(f"Duplicate read IDs in {r1}; fragment identity is ambiguous.")
    return {"n_input": 0, "n_bad_geometry": 0, "n_cb_has_n": 0, "n_umi_has_n": 0, **counts}


def read_filter_lineage_for(
    sample_dir: Path, resolved: dict[str, Any], read_filter: str
) -> tuple[Path | None, str]:
    """The lockstep lineage and fallback annotation for ``primary_read_filter`` (plan §5.2)."""
    if read_filter in (None, "", "off"):
        return None, "not_run"
    if resolved["source"] != "primary_host_filtered":
        # The primary filtered the raw pair; a discovery STAR pass is a different, reordered set.
        return None, "not_joined"
    path = sample_dir / "read_filtered" / "fragment_lineage.tsv.gz"
    if not path.is_file():
        raise DiscoveryError(f"The primary ran --read-filter {read_filter} but {path} is missing.")
    return path, "not_run"


def run_discover(args: argparse.Namespace) -> None:
    run_root = Path(args.run_dir).resolve()
    output = Path(args.output).resolve()
    db = Path(args.kraken2_db).resolve()

    from viralscan.runconfig import RunConfig

    check_tools()
    sample_dir = resolve_sample_dir(run_root, args.sample)
    check_output_dir(output, run_root)
    validate_kraken2_db(db)
    config = RunConfig.from_yaml(str(sample_dir / "config.yaml"))
    cb_len, umi_len = barcode_geometry(config.technology)
    partial = stage(output)
    resolved = resolve_input(run_root, sample_dir, partial, args)
    lineage, default = read_filter_lineage_for(sample_dir, resolved, config.read_filter)
    resolved["stream"] = stream_reads(
        resolved["r1"], resolved["r2"], cb_len, umi_len, partial, lineage, default
    )
    if resolved["stream"]["n_input"] != resolved["n_fragments"]:
        raise DiscoveryError(
            f"Streamed {resolved['stream']} but resolved {resolved['n_fragments']} fragments."
        )
    resolved.update(technology=config.technology, cb_len=cb_len, umi_len=umi_len)
    (partial / "input.json").write_text(json.dumps(resolved, indent=2, sort_keys=True) + "\n")
    log.info("Discovery input: %s, %s.", resolved["source"], resolved["stream"])
    raise DiscoveryError(
        "Classification stages are not implemented yet (PLAN DISC-07..09); the resolved input "
        f"and Kraken2 FASTA are staged in {partial} and nothing was published."
    )
