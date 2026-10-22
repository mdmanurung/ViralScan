"""Kraken2 cell-resolved taxonomic discovery on a completed run (PLAN DISC).

Design and contracts: ``docs/plans/2026-10-21-kraken2-discovery-plan.md``. Discovery is a
post-hoc layer: it reads a completed run and writes only into its own output directory, so it
cannot change primary inputs, count matrices, ``detected``, tiers or the run fingerprint. A
Kraken2 assignment is a taxonomic observation, never an infection call.

Implemented so far: run, output and database validation, atomic publication (DISC-02),
host-unmapped input resolution (DISC-03), the streaming identity pass (DISC-06) and Kraken2
classification with the exactly-once ordinal join (DISC-07), and UMI-group aggregation into a
cells x taxa matrix (DISC-08), provenance, fingerprint, publication and resume (DISC-09).
"""

from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import heapq
import io
import json
import logging
import os
import re
import shutil
import zlib
from array import array
from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np

from viralscan.chemistry import cb_umi_geometry, kb_whitelist_arg, onlist_path
from viralscan.evidence import _canonical_fastq_id, _open_maybe_gzip, _run
from viralscan.run_safety import RUN_COMPLETE, RUN_MANIFEST, sha256_file, sha256_path

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
#: Raw Kraken2 output (audit; hitlists kept) and report, and the ordinal-sorted calls.
KRAKEN_OUT = "kraken.out.gz"
KRAKEN_REPORT = "kraken_report.txt"
CLASSIFICATIONS_TSV = "classifications.tsv.gz"
TAXONOMY_TSV = "taxonomy.tsv"
#: Per fragment: corrected barcode, primary-encoded UMI and barcode status (plan §6.3).
BARCODES_TSV = "barcodes.tsv.gz"
UMI_GROUPS_TSV = "umi_groups.tsv.gz"
#: Cells x resolved taxa (one unit per resolved UMI group) and the off-cell remainder.
CELL_TAXON_H5AD = "cell_taxon.h5ad"
OFFCELL_TSV = "offcell_taxon.tsv.gz"
#: One row per resolved taxon: observation state only (plan R8); no calibrated tiers.
CANDIDATES_TSV = "candidates.tsv"
#: Candidate flag thresholds. Documented constants, not calibrated against a benchmark (plan §12.2).
#: Fewer distinct Kraken2 minimizers than this is a handful of reads or one repeated locus
#: (one 90-nt read carries roughly ten).
LOW_MINIMIZER_DISTINCT = 20
#: More than this fraction of a taxon's resolved groups outside the called cells.
OFFCELL_DOMINANT_FRACTION = 0.5
CANDIDATE_COLUMNS = (
    "taxid",
    "name",
    "rank",
    "groups_called",
    "groups_offcell",
    "n_cells",
    "reads",
    "distinct_minimizers",
    "minimizer_per_read",
    "host_relation",
    "db_has_9606",
    "state",
    "low_minimizer_diversity",
    "offcell_dominant",
    "host_lineage",
    "viral_only_db",
)
#: Provenance, invariant counts and the completion marker (written last, plan §11.2).
PROVENANCE_JSON = "provenance.json"
QC_JSON = "qc.json"
COMPLETE_JSON = "discovery_complete.json"
UMI_GROUP_COLUMNS = (
    "corrected_cb",
    "umi",
    "n_reads",
    "n_reads_classified",
    "n_reads_unclassified",
    "taxids",
    "resolved_taxid",
    "group_class",
    "lineage_consistent",
    "deepest_consistent_taxid",
    "conflict",
)
HUMAN_TAXID = 9606
#: R2 trimming before Kraken2 (--polya-trim): the first run of >= POLYA_MIN_RUN A (3' tail) is cut to
#: the end, and a leading run of >= POLYT_MIN_RUN T (5' of a reverse-strand tail) is cut from the start.
POLYA_MIN_RUN = 15
POLYT_MIN_RUN = 15
_POLYA = re.compile(f"A{{{POLYA_MIN_RUN},}}")
_POLYT = re.compile(f"^T{{{POLYT_MIN_RUN},}}")
READS_COLUMNS = (
    "ordinal",
    "read_id",
    "raw_cb",
    "raw_umi",
    "geometry_ok",
    "cb_has_n",
    "umi_has_n",
    "primary_read_filter",
    "r2_trim",
    "r2_len_classified",
)


class DiscoveryError(RuntimeError):
    """Discovery refused to run, or stopped, for a stated reason."""


def _gzip_out(path: Path) -> io.TextIOWrapper:
    """Text gzip writer with header mtime 0: identical content gives identical bytes (gate G8)."""
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", compresslevel=1, mtime=0), newline="")


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
    """Refuse an output inside the run, or one that already holds unfinished results.

    A published output (``COMPLETE_JSON``) is allowed through; ``run_discover`` then compares
    fingerprints to decide between a no-op and a refusal.
    """
    if output == run_root or run_root in output.parents:
        raise DiscoveryError(
            f"--output must be outside --run-dir: {output} is inside {run_root}. "
            "Discovery never writes into a primary run."
        )
    if (
        output.exists()
        and not (output / COMPLETE_JSON).is_file()
        and any(p.name != PARTIAL for p in output.iterdir())
    ):
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
    """Move every staged file into *output* (same filesystem, so each move is atomic).

    ``COMPLETE_JSON`` moves last: its presence means every other output is in place.
    """
    partial = output / PARTIAL
    for path in sorted(partial.iterdir(), key=lambda p: p.name == COMPLETE_JSON):
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


def _read_filter_rows(path: Path) -> Iterator[tuple[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            decision = row["filter_decision"]
            yield row["read_id"], decision if decision == "retained" else f"removed:{row['reason']}"


def trim_r2(seq: str) -> tuple[str, str]:
    """Cut a 3' polyA tail and a 5' polyT head; returns the sequence and none/polyA/polyT/both."""
    kind = []
    if (m := _POLYA.search(seq)) is not None:
        seq = seq[: m.start()]
        kind.append("polyA")
    if (m := _POLYT.match(seq)) is not None:
        seq = seq[m.end() :]
        kind.append("polyT")
    return seq, ("both" if len(kind) == 2 else kind[0] if kind else "none")


def stream_reads(
    r1: str,
    r2: str,
    cb_len: int,
    umi_len: int,
    out_dir: Path,
    read_filter_lineage: Path | None = None,
    read_filter_default: str = "not_run",
    polya_trim: bool = True,
) -> dict[str, int]:
    """One pass over the pair: R2 to ``KRAKEN_INPUT`` as ``>ordinal``, identity to ``READS_TSV``.

    Every fragment is classified; barcode problems are handled later, by ordinal, so the
    Kraken2 output must cover ``[0, n_input)`` exactly once. ``N`` in CB/UMI is kept (kallisto
    and bustools encode it as ``G``, plan §6.3) and flagged. *read_filter_lineage* is the primary
    read filter's per-fragment lineage over this same pair; it is read in lockstep and must
    match every read ID. Duplicate read IDs, mate mismatches and truncation abort. With
    *polya_trim* the R2 sent to Kraken2 loses its polyA/polyT ends (``trim_r2``); a read trimmed
    to nothing is sent as ``N`` (no k-mers, so it stays in ``[0, n_input)`` and comes back
    unclassified).
    """
    need = cb_len + umi_len
    counts: Counter[str] = Counter()
    # ponytail: 8-byte ID digests, 8 B/fragment (0.8 GB per 1e8); a false duplicate needs a
    # 64-bit collision (~3e-4 at 1e8 fragments) and fails closed. Shard by hash if RAM matters.
    digests = array("Q")
    rf_rows = _read_filter_rows(read_filter_lineage) if read_filter_lineage else None
    with (
        _open_maybe_gzip(r1) as fq1,
        _open_maybe_gzip(r2) as fq2,
        _gzip_out(out_dir / KRAKEN_INPUT) as fasta,
        _gzip_out(out_dir / READS_TSV) as reads,
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
            if rf_rows is None:
                annotation = read_filter_default
            else:
                step = next(rf_rows, None)
                if step is None or step[0] != read_id:
                    raise DiscoveryError(
                        f"Read-filter lineage is out of step with the input at record {ordinal}: "
                        f"{read_id!r} vs {step[0] if step else 'end of lineage'!r}"
                    )
                annotation = step[1]
            full2 = rec2[1].rstrip()
            seq2, trim = trim_r2(full2) if polya_trim else (full2, "none")
            fasta.write(f">{ordinal}\n{seq2 or 'N'}\n")
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
                    trim,
                    len(seq2),
                ]
            )
            counts["n_input"] += 1
            if trim != "none":
                counts["n_r2_trimmed"] += 1
                counts["n_r2_bases_trimmed"] += len(full2) - len(seq2)
            counts["n_bad_geometry"] += not geometry_ok
            counts["n_cb_has_n"] += "N" in cb
            counts["n_umi_has_n"] += "N" in umi
            counts[f"read_filter:{annotation.split(':')[0]}"] += 1
            ordinal += 1
    if rf_rows is not None and next(rf_rows, None) is not None:
        raise DiscoveryError("Read-filter lineage has more rows than the input pair.")
    ids = np.frombuffer(digests, dtype=np.uint64)
    ids.sort()
    if ids.size > 1 and bool((ids[1:] == ids[:-1]).any()):
        raise DiscoveryError(f"Duplicate read IDs in {r1}; fragment identity is ambiguous.")
    return {
        "n_input": 0,
        "n_bad_geometry": 0,
        "n_cb_has_n": 0,
        "n_umi_has_n": 0,
        "n_r2_trimmed": 0,
        "n_r2_bases_trimmed": 0,
        **counts,
    }


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


def kraken2_version() -> str:
    return _run(["kraken2", "--version"], capture=True).decode().splitlines()[0].strip()


def run_kraken2(
    db: Path,
    out_dir: Path,
    threads: int,
    confidence: float,
    minimum_hit_groups: int,
    memory_mapping: bool = False,
) -> Path:
    """Classify ``KRAKEN_INPUT``; returns the plain per-read output (``KRAKEN_REPORT`` beside it).

    No ``--use-names`` (IDs must stay integers), no ``--paired`` and no ``--quick``.
    ``--confidence`` is Kraken2's k-mer fraction threshold, not a probability (plan §2.5).
    """
    plain = out_dir / "kraken.out"
    cmd = [
        "kraken2",
        "--db",
        str(db),
        "--threads",
        str(threads),
        "--confidence",
        str(confidence),
        "--minimum-hit-groups",
        str(minimum_hit_groups),
        "--report",
        str(out_dir / KRAKEN_REPORT),
        "--report-minimizer-data",
        "--output",
        str(plain),
    ]
    if memory_mapping:
        cmd.append("--memory-mapping")
    cmd.append(str(out_dir / KRAKEN_INPUT))
    try:
        _run(cmd)
    except RuntimeError as exc:
        raise DiscoveryError(str(exc)) from exc
    return plain


def join_kraken2(plain: Path, n_input: int, out_dir: Path) -> dict[str, int]:
    """Check Kraken2's output covers every ordinal exactly once; write the ordinal-sorted calls.

    Kraken2 keeps input order, but the join does not rely on it: each sequence ID must be an
    integer in ``[0, n_input)`` seen once, ``C`` rows carry a taxid > 0 and ``U`` rows 0. Any
    violation aborts (plan §6.2, I2/I3). The raw output is kept gzipped as ``KRAKEN_OUT``.
    """
    seen = np.zeros(n_input, dtype=bool)
    taxid = np.zeros(n_input, dtype=np.uint32)
    with plain.open() as handle, _gzip_out(out_dir / KRAKEN_OUT) as raw:
        for line_number, line in enumerate(handle, 1):
            raw.write(line)
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise DiscoveryError(
                    f"Kraken2 output line {line_number} has {len(fields)} columns."
                )
            status, seq_id, tax = fields[0], fields[1], fields[2]
            if not seq_id.isdigit() or not tax.isdigit():
                raise DiscoveryError(f"Kraken2 output line {line_number}: non-integer ID or taxid.")
            ordinal, tax_int = int(seq_id), int(tax)
            if ordinal >= n_input:
                raise DiscoveryError(f"Kraken2 returned ordinal {ordinal} outside [0, {n_input}).")
            if seen[ordinal]:
                raise DiscoveryError(f"Kraken2 returned ordinal {ordinal} twice.")
            if (status, tax_int == 0) not in (("C", False), ("U", True)):
                raise DiscoveryError(
                    f"Kraken2 line {line_number}: status {status!r} with taxid {tax}."
                )
            seen[ordinal] = True
            taxid[ordinal] = tax_int
    if not seen.all():
        missing = np.flatnonzero(~seen)
        raise DiscoveryError(f"Kraken2 output lacks {missing.size} ordinals (first {missing[0]}).")
    plain.unlink()
    with _gzip_out(out_dir / CLASSIFICATIONS_TSV) as out:
        out.write("ordinal\tkraken_status\ttaxid\n")
        for ordinal, tax_int in enumerate(taxid.tolist()):
            out.write(f"{ordinal}\t{'C' if tax_int else 'U'}\t{tax_int}\n")
    n_classified = int(np.count_nonzero(taxid))
    return {"n_classified": n_classified, "n_unclassified": n_input - n_classified}


# ── Taxonomy (DISC-04: the database's own tree) ──────────────────────────────


def parse_inspect(text: str) -> dict[int, tuple[int, str, str]]:
    """``taxid -> (parent, rank_code, name)`` from ``k2 inspect --report-zero-counts`` output.

    Depth is the name's indentation (two spaces per level); the parent is the nearest
    preceding line one level up. The root's parent is itself. ``unclassified`` (0) is skipped.
    """
    tree: dict[int, tuple[int, str, str]] = {}
    stack: list[int] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 6 or not fields[4].strip().isdigit():
            raise DiscoveryError(f"Unexpected k2 inspect line: {line!r}")
        taxid, rank, raw_name = int(fields[4]), fields[3].strip(), fields[5]
        if taxid == 0:
            continue
        depth, rem = divmod(len(raw_name) - len(raw_name.lstrip(" ")), 2)
        if rem or depth > len(stack):
            raise DiscoveryError(f"k2 inspect indentation jumps at taxid {taxid}.")
        del stack[depth:]
        tree[taxid] = (stack[-1] if stack else taxid, rank, raw_name.strip())
        stack.append(taxid)
    if not tree:
        raise DiscoveryError("k2 inspect returned no taxonomy.")
    return tree


def load_taxonomy(db: Path, memory_mapping: bool = False) -> dict[int, tuple[int, str, str]]:
    # ponytail: a full inspect scans the hash table on every run; cache by DB identity
    # (DISC-09 provenance hash) if large databases make this slow.
    cmd = ["k2", "inspect", "--db", str(db), "--report-zero-counts"]
    if memory_mapping:
        cmd.append("--memory-mapping")
    try:
        return parse_inspect(_run(cmd, capture=True).decode())
    except RuntimeError as exc:
        raise DiscoveryError(str(exc)) from exc


# One memo per taxonomy tree (identity-checked; the tree is built once per run). Paths are
# root-first tuples filled lazily for the taxids actually queried, so a 2M-node tree costs nothing.
_PATHS: tuple[dict[Any, Any], dict[int, tuple[int, ...]]] = ({}, {})


def _path(tree: dict[int, tuple[int, str, str]], taxid: int) -> tuple[int, ...]:
    """Root-first ancestor path of ``taxid`` (root and ``taxid`` inclusive), memoised per tree."""
    global _PATHS
    if _PATHS[0] is not tree:
        _PATHS = (tree, {})
    memo = _PATHS[1]
    todo = []
    while taxid not in memo:
        parent = tree[taxid][0]
        if parent == taxid:
            memo[taxid] = (taxid,)
            break
        todo.append(taxid)
        taxid = parent
    path = memo[taxid]
    for node in reversed(todo):
        path = memo[node] = (*path, node)
    return path


def lineage(tree: dict[int, tuple[int, str, str]], taxid: int) -> list[int]:
    """``taxid`` up to the root, inclusive."""
    return list(reversed(_path(tree, taxid)))


def lca(tree: dict[int, tuple[int, str, str]], taxids: set[int]) -> int:
    paths = [_path(tree, taxid) for taxid in taxids]
    common = paths[0][0]  # the root, shared by every path
    for nodes in zip(*paths):
        if len(set(nodes)) > 1:
            break
        common = nodes[0]
    return common


def host_relation(tree: dict[int, tuple[int, str, str]], taxid: int) -> str:
    """``host`` (9606 or below), ``host_ancestor`` (a non-root ancestor of 9606) or ``other``."""
    if HUMAN_TAXID in _path(tree, taxid):
        return "host"
    if (
        HUMAN_TAXID in tree
        and taxid in _path(tree, HUMAN_TAXID)
        and not tree[taxid][1].startswith("R")
    ):
        return "host_ancestor"
    return "other"


def write_taxonomy(tree: dict[int, tuple[int, str, str]], out_dir: Path) -> None:
    with (out_dir / TAXONOMY_TSV).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxid", "parent_taxid", "rank", "name", "depth", "host_relation"])
        for taxid in sorted(tree):
            parent, rank, name = tree[taxid]
            depth = len(_path(tree, taxid)) - 1
            writer.writerow([taxid, parent, rank, name, depth, host_relation(tree, taxid)])


# ── Barcode correction (DISC-05: the primary's bustools and on-list) ─────────


def resolve_onlist(technology: str, whitelist: str | None) -> str | None:
    """The on-list ``kb count`` corrected against, or None when it skipped correction."""
    arg = kb_whitelist_arg(technology, whitelist)
    if arg == "None":
        return None
    if arg:
        return arg
    try:
        return onlist_path(technology)
    except ValueError:  # explicit bc:umi:seq geometry without -w: kb had no on-list either
        return None


def bustools_barcodes(out_dir: Path, onlist: str | None, threads: int) -> Path:
    """Run the staged BUS text through bustools; returns flag-sorted ``text -f`` output.

    Encoding (``N`` -> ``G``) and correction are bustools' own, so ``corrected_cb`` is in the
    primary ``kb count`` barcode space. Records ``correct`` cannot place are dropped by it.
    """
    from viralscan.validation import tool_path

    bustools = tool_path("bustools")
    if bustools is None:
        raise DiscoveryError("bustools (kb-python's) is required for barcode correction.")
    work = out_dir / "bus"
    raw, current = work / "raw.bus", work / "raw.bus"
    try:
        _run([bustools, "fromtext", "-o", str(raw), str(work / "fragments.txt")])
        if onlist is not None:
            plain = Path(onlist)
            if plain.suffix == ".gz":
                plain = work / "onlist.txt"
                with gzip.open(onlist, "rb") as src, plain.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
            _run([bustools, "sort", "-t", str(threads), "-o", str(work / "sorted.bus"), str(raw)])
            _run(
                [
                    bustools,
                    "correct",
                    "-w",
                    str(plain),
                    "-o",
                    str(work / "corrected.bus"),
                    str(work / "sorted.bus"),
                ]
            )
            current = work / "corrected.bus"
        by_flag = work / "by_flag.bus"
        _run([bustools, "sort", "--flags", "-t", str(threads), "-o", str(by_flag), str(current)])
        text = work / "by_flag.txt"
        _run([bustools, "text", "-f", "-o", str(text), str(by_flag)])
    except RuntimeError as exc:
        raise DiscoveryError(str(exc)) from exc
    return text


def _tsv_rows(path: Path) -> Iterator[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def correct_barcodes(out_dir: Path, onlist: str | None, threads: int) -> dict[str, int]:
    """Write ``BARCODES_TSV`` (ordinal order) with ``corrected_cb``, ``umi`` and ``cb_status``."""
    work = out_dir / "bus"
    work.mkdir()
    with (work / "fragments.txt").open("w") as handle:
        for row in _tsv_rows(out_dir / READS_TSV):
            if row["geometry_ok"] == "1":
                if int(row["ordinal"]) >= 2**32:
                    raise DiscoveryError(
                        "More than 2^32 fragments; BUS flags cannot hold the ordinal."
                    )
                handle.write(f"{row['raw_cb']}\t{row['raw_umi']}\t0\t1\t{row['ordinal']}\n")
    text = bustools_barcodes(out_dir, onlist, threads)
    counts: Counter[str] = Counter()
    with text.open() as bus, _gzip_out(out_dir / BARCODES_TSV) as out:
        out.write("ordinal\tcorrected_cb\tumi\tcb_status\n")
        pending = next(bus, None)
        for row in _tsv_rows(out_dir / READS_TSV):
            ordinal = int(row["ordinal"])
            cb = umi = ""
            if row["geometry_ok"] != "1":
                status = "bad_geometry"
            elif pending is not None and int(pending.split("\t")[4]) == ordinal:
                cb, umi = pending.split("\t")[:2]
                pending = next(bus, None)
                if onlist is None:
                    status = "no_onlist"
                else:
                    status = "onlist" if cb == row["raw_cb"].replace("N", "G") else "corrected"
            else:
                status = "uncorrectable"
            counts[status] += 1
            out.write(f"{ordinal}\t{cb}\t{umi}\t{status}\n")
        if pending is not None:
            raise DiscoveryError("bustools returned a flag that is not a staged ordinal.")
    if onlist is None and counts["uncorrectable"]:
        raise DiscoveryError("bustools dropped records without an on-list; encoding failed.")
    shutil.rmtree(work)
    return {
        f"n_cb_{key}": counts[key]
        for key in ("onlist", "corrected", "uncorrectable", "no_onlist", "bad_geometry")
    }


# ── UMI groups (plan §6.4) ───────────────────────────────────────────────────

ELIGIBLE = frozenset({"onlist", "corrected", "no_onlist"})


def _shard_fragments(out_dir: Path, n_shards: int) -> list[Path]:
    """Split eligible fragments by ``crc32(CB)`` so each shard holds whole UMI groups."""
    paths = [out_dir / f"shard{i}.tsv" for i in range(n_shards)]
    handles = [path.open("w") for path in paths]
    try:
        calls = _tsv_rows(out_dir / CLASSIFICATIONS_TSV)
        for barcode in _tsv_rows(out_dir / BARCODES_TSV):
            call = next(calls)
            if call["ordinal"] != barcode["ordinal"]:
                raise DiscoveryError(f"Ordinal mismatch at {barcode['ordinal']}.")
            if barcode["cb_status"] in ELIGIBLE:
                cb = barcode["corrected_cb"]
                handles[zlib.crc32(cb.encode()) % n_shards].write(
                    f"{cb}\t{barcode['umi']}\t{call['taxid']}\n"
                )
    finally:
        for handle in handles:
            handle.close()
    return paths


def _group_shard(path: Path, tree: dict[int, tuple[int, str, str]]) -> list[tuple[Any, ...]]:
    reads: dict[tuple[str, str], list[Any]] = {}
    with path.open() as handle:
        for line in handle:
            cb, umi, tax = line.rstrip("\n").split("\t")
            entry = reads.setdefault((cb, umi), [0, 0, set()])
            entry[0] += 1
            if tax == "0":
                entry[1] += 1
            else:
                entry[2].add(int(tax))
    path.unlink()
    rows = []
    for (cb, umi), (n, n_unclassified, taxids) in reads.items():
        resolved: int | str = ""
        consistent: bool | str = ""
        deepest: int | str = ""
        cls = "unclassified"
        if taxids:
            deep = max(taxids, key=lambda t: len(_path(tree, t)))
            consistent = all(t in _path(tree, deep) for t in taxids)
            deepest = deep if consistent else ""
            # Deepest when every classified taxon lies on one root-to-leaf path, else the LCA.
            top = deep if consistent else lca(tree, taxids)
            cls = "ambiguous" if tree[top][1].startswith("R") else "resolved"
            resolved = top if cls == "resolved" else ""
        rows.append(
            (
                cb,
                umi,
                n,
                n - n_unclassified,
                n_unclassified,
                ";".join(map(str, sorted(taxids))),
                resolved,
                cls,
                "" if not taxids else int(consistent),
                deepest,
                "" if not taxids else int(not consistent),
            )
        )
    rows.sort()
    return rows


def build_umi_groups(
    out_dir: Path, tree: dict[int, tuple[int, str, str]], n_shards: int = 16
) -> dict[str, int]:
    """Write ``UMI_GROUPS_TSV`` sorted by ``(corrected_cb, umi)``; independent of *n_shards*."""
    unknown = {int(r["taxid"]) for r in _tsv_rows(out_dir / CLASSIFICATIONS_TSV)} - set(tree) - {0}
    if unknown:
        raise DiscoveryError(
            f"Kraken2 taxids absent from the database taxonomy: {sorted(unknown)[:5]}"
        )
    shards = [_group_shard(path, tree) for path in _shard_fragments(out_dir, n_shards)]
    counts: Counter[str] = Counter()
    with _gzip_out(out_dir / UMI_GROUPS_TSV) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(UMI_GROUP_COLUMNS)
        for row in heapq.merge(*shards):
            writer.writerow(row)
            counts[f"n_groups_{row[7]}"] += 1
            counts["n_group_reads"] += row[2]
    return {
        "n_groups_resolved": 0,
        "n_groups_ambiguous": 0,
        "n_groups_unclassified": 0,
        "n_group_reads": 0,
        **counts,
    }


# ── Cell universe and matrix (plan §6.5, §8, §10) ───────────────────────────


def load_called_cells(sample_dir: Path) -> tuple[list[str], dict[str, Any]]:
    """The primary's called barcodes (suffix-stripped, sorted) and its cell-calling receipt."""
    from viralscan.scripts.cellcalling import CALLED_CELLS_TSV, _strip_suffix

    path = sample_dir / CALLED_CELLS_TSV
    if not path.is_file():
        raise DiscoveryError(f"{path} is missing; discovery needs the primary's called cells.")
    with path.open(encoding="utf-8") as handle:
        handle.readline()
        called = sorted({_strip_suffix(line.strip()) for line in handle if line.strip()})
    receipt_path = sample_dir / "results" / "cell_calling_input.json"
    receipt = json.loads(receipt_path.read_text()) if receipt_path.is_file() else {}
    # Two-step EmptyDrops on the host filter's STARsolo matrix calls cells in STAR's barcode
    # space, which is uncorrected without -w (plan §6.3); the overlap below shows the effect.
    if not receipt:
        space = "unknown"  # runs before the receipt existed; see the overlap counts
    else:
        space = "starsolo" if receipt.get("count_layer") == "GeneFull.raw" else "kb"
    info = {
        "called_cells_sha256": sha256_file(path),
        "n_called_cells": len(called),
        "cell_calling_method": receipt.get("method"),
        "called_cells_barcode_space": space,
    }
    return called, info


def build_cell_matrix(
    out_dir: Path, tree: dict[int, tuple[int, str, str]], called: list[str]
) -> dict[str, int]:
    """Write ``CELL_TAXON_H5AD`` (every called cell, resolved called groups) and ``OFFCELL_TSV``.

    ``X`` counts resolved UMI groups at their resolved taxid: one unit per group, never reads,
    never ambiguous or unclassified groups. Off-cell resolved groups go to ``OFFCELL_TSV``.
    """
    import anndata as ad
    import pandas as pd
    from scipy import sparse

    cell_index = {cb: i for i, cb in enumerate(called)}
    in_cells: Counter[tuple[int, int]] = Counter()
    off_cells: Counter[tuple[str, int]] = Counter()
    counts: Counter[str] = Counter()
    seen_called: set[str] = set()
    for row in _tsv_rows(out_dir / UMI_GROUPS_TSV):
        cb, cls = row["corrected_cb"], row["group_class"]
        where = "called" if cb in cell_index else "offcell"
        counts[f"n_groups_{where}"] += 1
        counts[f"n_groups_{where}_{cls}"] += 1
        if where == "called":
            seen_called.add(cb)
        if cls != "resolved":
            continue
        taxid = int(row["resolved_taxid"])
        if where == "called":
            in_cells[cell_index[cb], taxid] += 1
        else:
            off_cells[cb, taxid] += 1
    taxa = sorted({taxid for _cell, taxid in in_cells})
    col = {taxid: j for j, taxid in enumerate(taxa)}
    keys = sorted(in_cells)
    matrix = sparse.csr_matrix(
        (
            np.array([in_cells[k] for k in keys], dtype=np.int32),
            (
                np.array([k[0] for k in keys], dtype=np.int64),
                np.array([col[k[1]] for k in keys], dtype=np.int64),
            ),
        ),
        shape=(len(called), len(taxa)),
        dtype=np.int32,
    )
    if int(matrix.sum()) != counts["n_groups_called_resolved"]:
        raise DiscoveryError("Matrix total differs from the resolved called-cell groups (gate G6).")
    var = pd.DataFrame(
        {
            "name": [tree[t][2] for t in taxa],
            "rank": [tree[t][1] for t in taxa],
            "parent_taxid": [tree[t][0] for t in taxa],
            "host_relation": [host_relation(tree, t) for t in taxa],
        },
        index=pd.Index([str(t) for t in taxa], name="taxid"),
    )
    summary = {
        "n_called_cells": len(called),
        "n_called_cells_seen": len(seen_called),
        "n_called_cells_with_resolved": len({cell for cell, _t in in_cells}),
        "n_taxa": len(taxa),
        **{
            f"n_groups_{w}{c}": counts[f"n_groups_{w}{c}"]
            for w in ("called", "offcell")
            for c in ("", "_resolved", "_ambiguous", "_unclassified")
        },
    }
    adata = ad.AnnData(X=matrix, obs=pd.DataFrame(index=pd.Index(called, name="barcode")), var=var)
    adata.uns["discovery"] = {
        "unit": "resolved_umi_group",
        "note": "Kraken2 taxonomic observations; not viral molecule counts and not infection calls.",
        **summary,
    }
    adata.write_h5ad(out_dir / CELL_TAXON_H5AD)
    with _gzip_out(out_dir / OFFCELL_TSV) as handle:
        handle.write("corrected_cb\ttaxid\tn_groups\n")
        for cb, taxid in sorted(off_cells):
            handle.write(f"{cb}\t{taxid}\t{off_cells[cb, taxid]}\n")
    return summary


def parse_kraken_report(path: Path) -> dict[int, tuple[int, int]]:
    """``{taxid: (reads assigned directly, distinct minimizers)}`` from a minimizer-data report.

    Columns: pct, clade reads, own reads, minimizers, distinct minimizers, rank, taxid, name.
    """
    report = {}
    for line in path.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) != 8:
            raise DiscoveryError(f"{path}: expected 8 report columns, got {len(cols)}: {line!r}")
        if cols[6] != "0":
            report[int(cols[6])] = (int(cols[2]), int(cols[4]))
    return report


def build_candidates(
    out_dir: Path, tree: dict[int, tuple[int, str, str]], called: list[str]
) -> int:
    """Write ``CANDIDATES_TSV`` (one row per resolved taxon, ``state = observation``); returns its rows."""
    report = parse_kraken_report(out_dir / KRAKEN_REPORT)
    called_set = set(called)
    groups: dict[int, list[int]] = {}  # taxid -> [called groups, off-cell groups]
    cells: dict[int, set[str]] = {}
    for row in _tsv_rows(out_dir / UMI_GROUPS_TSV):
        if row["group_class"] != "resolved":
            continue
        taxid, cb = int(row["resolved_taxid"]), row["corrected_cb"]
        entry = groups.setdefault(taxid, [0, 0])
        if cb in called_set:
            entry[0] += 1
            cells.setdefault(taxid, set()).add(cb)
        else:
            entry[1] += 1
    has_human = HUMAN_TAXID in tree
    rows = []
    for taxid, (n_called, n_off) in groups.items():
        reads, distinct = report.get(taxid, (None, None))
        relation = host_relation(tree, taxid)
        rows.append(
            (
                taxid,
                tree[taxid][2],
                tree[taxid][1],
                n_called,
                n_off,
                len(cells.get(taxid, ())),
                "" if reads is None else reads,
                "" if distinct is None else distinct,
                "" if not reads or distinct is None else f"{distinct / reads:.4g}",
                relation,
                int(has_human),
                "observation",
                int(distinct is None or distinct < LOW_MINIMIZER_DISTINCT),
                int(n_off / (n_called + n_off) > OFFCELL_DOMINANT_FRACTION),
                int(relation != "other"),
                int(not has_human),
            )
        )
    rows.sort(key=lambda r: (-(r[3] + r[4]), r[0]))
    with (out_dir / CANDIDATES_TSV).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(CANDIDATE_COLUMNS)
        writer.writerows(rows)
    return len(rows)


# ── Candidate states (DISC-13) ───────────────────────────────────────────────

CANDIDATE_STATES_TSV = "candidate_states.tsv"
CANDIDATE_STATE_COLUMNS = (
    "taxid",
    "name",
    "rank",
    "state",
    "corroboration",
    "read_support",
    "reason",
)
EVIDENCE_COLUMNS = (
    "taxid",
    "route",
    "source",
    "n_adjudicated",
    "n_taxon_best",
    "n_shared_cells",
    "n_shared_umis",
)
#: Route (i) needs kb or STAR, route (ii) is kb-independent (any competitive aligner label).
CORROBORATION_SOURCES = frozenset({"kb", "star"})
ROUTES = ("corroboration", "read_support")
#: Documented constants, not calibrated (plan §12.2): at least this many molecules (reads for route ii)
#: scoring best for the taxon, at least this fraction of everything adjudicated (ties and other-reference
#: wins stay in the denominator), and at least this many of the taxon's called cells / molecules shared.
MIN_TAXON_BEST = 3
MIN_TAXON_BEST_FRACTION = 0.8
MIN_SHARED_CELLS = 1
MIN_SHARED_UMIS = 1


def read_evidence_table(path: Path) -> dict[tuple[int, str], dict[str, Any]]:
    """Read the documented evidence input table, one row per ``(taxid, route)``; fail closed.

    Columns ``EVIDENCE_COLUMNS``. ``route`` is ``corroboration`` (``source`` kb or star; counts are
    molecules from the competitive adjudication) or ``read_support`` (``source`` free label; counts are
    reads from a competitive alignment to the hit genome versus host + panel). ``n_adjudicated`` is all
    units tested, ``n_taxon_best`` those that scored best for the taxon alone; ``n_shared_cells`` and
    ``n_shared_umis`` are the cells and distinct (CB, UMI) molecules that are also resolved called-cell
    groups of the taxon in ``umi_groups.tsv.gz``. Ties and other-reference wins are not counted as taxon-best.
    """
    rows: dict[tuple[int, str], dict[str, Any]] = {}
    for raw in _tsv_rows(path):
        missing = [c for c in EVIDENCE_COLUMNS if c not in raw]
        if missing:
            raise DiscoveryError(f"{path}: missing columns {missing}")
        try:
            row: dict[str, Any] = {c: int(raw[c]) for c in EVIDENCE_COLUMNS[3:]}
            row["taxid"] = int(raw["taxid"])
        except ValueError as exc:
            raise DiscoveryError(f"{path}: non-integer value in {raw}") from exc
        row["route"], row["source"] = raw["route"], raw["source"]
        if row["route"] not in ROUTES or not row["source"]:
            raise DiscoveryError(f"{path}: bad route/source in {raw}")
        if row["route"] == "corroboration" and row["source"] not in CORROBORATION_SOURCES:
            raise DiscoveryError(f"{path}: corroboration source must be kb or star: {raw}")
        counts = [row[c] for c in EVIDENCE_COLUMNS[3:]]
        if min(counts) < 0 or row["n_taxon_best"] > row["n_adjudicated"]:
            raise DiscoveryError(f"{path}: inconsistent counts in {raw}")
        key = (row["taxid"], row["route"])
        if key in rows:
            raise DiscoveryError(
                f"{path}: duplicate evidence row for taxid {key[0]} route {key[1]}"
            )
        rows[key] = row
    return rows


def judge_evidence(row: dict[str, Any] | None, route: str) -> tuple[str, str]:
    """``(assessment, why)``: ``not_assessed`` (no row), ``negative`` or ``positive``. Never conflated."""
    if row is None:
        return "not_assessed", ""
    best, tested = row["n_taxon_best"], row["n_adjudicated"]
    if best < MIN_TAXON_BEST or best < MIN_TAXON_BEST_FRACTION * tested:
        return "negative", f"{route}: {best}/{tested} units best for the taxon"
    if row["n_shared_cells"] < MIN_SHARED_CELLS:
        return "negative", f"{route}: no shared called cell"
    if route == "read_support" and row["n_shared_umis"] < MIN_SHARED_UMIS:
        return "negative", "read_support: no shared (CB, UMI) molecule"
    return "positive", ""


def assign_state(
    candidate: dict[str, str], evidence: dict[tuple[int, str], dict[str, Any]]
) -> dict[str, str]:
    """State for one ``candidates.tsv`` row (Q12). Kraken2 alone never goes above ``observation``.

    ``corroborated``: route (i) positive. ``read_supported_observation``: route (ii) positive (and (i) not);
    an observation backed by reads, never an infection claim. Host-lineage taxa and a viral-only database
    (residual host invisible, approved G15) stay at ``observation`` whatever the evidence says.
    """
    taxid = int(candidate["taxid"])
    results = {r: judge_evidence(evidence.get((taxid, r)), r) for r in ROUTES}
    out = {c: candidate[c] for c in ("taxid", "name", "rank")}
    out["corroboration"], out["read_support"] = (
        results["corroboration"][0],
        results["read_support"][0],
    )
    if candidate["host_lineage"] == "1" or candidate["viral_only_db"] == "1":
        out["state"] = "observation"
        out["reason"] = "capped: " + (
            "host lineage" if candidate["host_lineage"] == "1" else "viral-only database"
        )
    elif out["corroboration"] == "positive":
        out["state"], out["reason"] = "corroborated", "same-taxon competitive corroboration"
    elif out["read_support"] == "positive":
        out["state"], out["reason"] = (
            "read_supported_observation",
            "kb-independent competitive read support",
        )
    else:
        why = [w for _, w in results.values() if w]
        out["state"] = "observation"
        out["reason"] = "; ".join(why) or "no orthogonal evidence supplied"
    return out


def build_candidate_states(out_dir: Path, evidence_path: Path | None = None) -> int:
    """Write ``candidate_states.tsv`` beside ``candidates.tsv``; returns its rows.

    ``evidence_path`` is the documented evidence table (``read_evidence_table``); without it every
    candidate is ``observation`` with both routes ``not_assessed``. Evidence for a taxid that is not a
    candidate is an error (a stale or mismatched table).
    """
    candidates = list(_tsv_rows(out_dir / CANDIDATES_TSV))
    evidence = read_evidence_table(evidence_path) if evidence_path else {}
    unknown = {t for t, _ in evidence} - {int(c["taxid"]) for c in candidates}
    if unknown:
        raise DiscoveryError(f"evidence rows for taxids not in {CANDIDATES_TSV}: {sorted(unknown)}")
    with (out_dir / CANDIDATE_STATES_TSV).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(CANDIDATE_STATE_COLUMNS)
        for candidate in candidates:
            state = assign_state(candidate, evidence)
            writer.writerow([state[c] for c in CANDIDATE_STATE_COLUMNS])
    return len(candidates)


# ── Rank views (DISC-11) ─────────────────────────────────────────────────────

COLLAPSE_RANKS = {"species": "S", "genus": "G", "family": "F"}


def collapse_rank(discovery_dir: Path, rank: str, output: Path | None = None) -> Path:
    """Project a published ``CELL_TAXON_H5AD`` onto ``rank`` ancestors (plan DISC-11).

    Each column moves to its nearest ancestor-or-self at ``rank`` (Kraken2 code S/G/F, so
    subspecies S1.. collapse to their species). Columns with no such ancestor (above the rank,
    or on a path that skips it) stay as themselves with ``above_rank = True``: nothing is
    reclassified, and every cell's row sum is conserved (gate G6).
    """
    import anndata as ad
    import pandas as pd
    from scipy import sparse

    if rank not in COLLAPSE_RANKS:
        raise DiscoveryError(f"--rank must be one of {sorted(COLLAPSE_RANKS)}.")
    code = COLLAPSE_RANKS[rank]
    tax_path, h5ad = discovery_dir / TAXONOMY_TSV, discovery_dir / CELL_TAXON_H5AD
    for path in (tax_path, h5ad):
        if not path.is_file():
            raise DiscoveryError(f"{path} not found; is {discovery_dir} a published discovery?")
    with tax_path.open(newline="") as handle:
        tree = {
            int(row["taxid"]): (int(row["parent_taxid"]), row["rank"], row["name"])
            for row in csv.DictReader(handle, delimiter="\t")
        }
    adata = ad.read_h5ad(h5ad)
    target: list[int] = []
    above: list[bool] = []
    for name in adata.var_names:
        taxid = int(name)
        if taxid not in tree:
            raise DiscoveryError(f"taxid {taxid} of {h5ad.name} is missing from {tax_path.name}.")
        hit = next((t for t in reversed(_path(tree, taxid)) if tree[t][1] == code), None)
        target.append(taxid if hit is None else hit)
        above.append(hit is None)
    new_taxa = sorted(set(target))
    col = {t: j for j, t in enumerate(new_taxa)}
    project = sparse.csr_matrix(
        (np.ones(len(target), dtype=np.int32), (np.arange(len(target)), [col[t] for t in target])),
        shape=(len(target), len(new_taxa)),
        dtype=np.int32,
    )
    matrix = (adata.X @ project).tocsr()
    if int(matrix.sum()) != int(adata.X.sum()):
        raise DiscoveryError("Collapsed total differs from the source matrix (gate G6).")
    above_by_taxon = dict(zip(target, above, strict=True))
    var = pd.DataFrame(
        {
            "name": [tree[t][2] for t in new_taxa],
            "rank": [tree[t][1] for t in new_taxa],
            "parent_taxid": [tree[t][0] for t in new_taxa],
            "host_relation": [host_relation(tree, t) for t in new_taxa],
            "above_rank": [above_by_taxon[t] for t in new_taxa],
            "n_source_taxa": np.bincount(project.indices, minlength=len(new_taxa)),
        },
        index=pd.Index([str(t) for t in new_taxa], name="taxid"),
    )
    out = ad.AnnData(X=matrix, obs=adata.obs.copy(), var=var)
    out.uns["discovery"] = {**adata.uns.get("discovery", {}), "collapsed_to_rank": rank}
    output = output or discovery_dir / f"cell_taxon_{rank}.h5ad"
    out.write_h5ad(output)
    return output


# ── Provenance, fingerprint and resume (DISC-09, plan §11.2) ────────────────


def _md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def db_identity(db: Path) -> dict[str, Any]:
    """sha256 of the three database files, plus the provider md5 check when a manifest ships.

    A manifest is any ``*.md5`` or ``md5sum.txt`` of ``<md5>  <file>`` lines; listed files that
    are present are verified and a mismatch is refused (a damaged download, not a setting).
    """
    identity: dict[str, Any] = {"sha256": {n: sha256_file(db / n) for n in KRAKEN2_DB_FILES}}
    manifests = sorted([*db.glob("*.md5"), *db.glob("md5sum.txt")])
    if manifests:
        checked: dict[str, str] = {}
        for manifest in manifests:
            for line in manifest.read_text(encoding="utf-8", errors="replace").splitlines():
                fields = line.split()
                if len(fields) == 2 and (db / fields[1].lstrip("*")).is_file():
                    name = fields[1].lstrip("*")
                    ok = _md5_file(db / name) == fields[0].lower()
                    checked[name] = "ok" if ok else "mismatch"
        if "mismatch" in checked.values():
            raise DiscoveryError(f"Kraken2 database {db} fails its md5 manifest: {checked}")
        identity["md5_manifest"] = {"files": manifests[0].name, "checked": checked}
    else:
        identity["md5_manifest"] = None
    return identity


def _input_identity(sample_dir: Path, run_root: Path, args: argparse.Namespace) -> dict[str, Any]:
    """What determines the classified reads, known before any expensive work."""
    if args.host_subtract:
        manifest = json.loads((run_root / RUN_MANIFEST).read_text(encoding="utf-8"))
        return {
            "source": "discovery_star_pass",
            "sample_input_fingerprints": {
                key: manifest.get("input_fingerprints", {}).get(key)
                for key in sorted(manifest.get("input_fingerprints", {}))
                if key.startswith(("sample1:", "sample2:"))
            },
            "star_param_set": args.host_filter_star_params,
            # ponytail: the index path, not its bytes (hashing a genome is minutes); hash it if
            # indices get rebuilt in place.
            "host_index": str(Path(args.host_index).resolve()) if args.host_index else None,
        }
    primary = read_host_filtered(sample_dir / "host_filtered")
    if primary is None:
        return {"source": None}  # resolve_input refuses with the actionable message
    return {
        "source": primary["source"],
        "sha256": primary["sha256"],
        "star_param_set": primary["star_param_set"],
        "star_args": primary["star_args"],
    }


def kb_input_check(config: Any, resolved: dict[str, Any]) -> str:
    """Whether the primary ``kb count`` read the pair discovery classifies.

    A run directory may have been copied, so ``kb_r1``/``kb_r2`` can name another path; then the
    contents must match. A different pair under ``host_filtered/`` is refused.
    """
    if resolved["source"] != "primary_host_filtered":
        return "not_applicable"
    if not (config.kb_r1 and config.kb_r2):
        return "not_recorded"
    kb = [Path(config.kb_r1), Path(config.kb_r2)]
    if any(p.parent.name != "host_filtered" for p in kb):
        return "kb_reads_other_pair"  # e.g. read_filtered/, a subset of this pair
    if [p.resolve() for p in kb] == [Path(resolved[k]).resolve() for k in ("r1", "r2")]:
        return "same_path"
    if not all(p.is_file() for p in kb):
        return "kb_input_absent"
    if [sha256_file(p) for p in kb] != [resolved["sha256"][k] for k in ("r1", "r2")]:
        raise DiscoveryError(
            f"The primary kb count read {kb[0]} / {kb[1]}, which differ from "
            f"{resolved['r1']} / {resolved['r2']}."
        )
    return "same_content"


def discovery_fingerprint(provenance: dict[str, Any]) -> str:
    """sha256 of the provenance subset that determines the outputs (not threads or paths)."""
    canonical = json.dumps(provenance["fingerprint_inputs"], sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def _forget_published(output: Path, complete: dict[str, Any]) -> None:
    """``--force``: remove exactly what a previous publication recorded, then the marker."""
    for name in complete.get("outputs", {}):
        if Path(name).name != name or name in ("", ".", ".."):
            raise DiscoveryError(f"{output / COMPLETE_JSON} lists a non-local output {name!r}.")
        path = output / name
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    (output / COMPLETE_JSON).unlink()


def build_qc(resolved: dict[str, Any]) -> dict[str, Any]:
    """Every invariant count of the run, plus barcode-space overlap and fractions."""
    kraken, cells = resolved["kraken2"], resolved["cells"]
    n_input = resolved["stream"]["n_input"]
    return {
        "n_fragments": resolved["n_fragments"],
        "stream": resolved["stream"],
        "kraken2": {k: v for k, v in kraken.items() if k.startswith("n_")},
        "barcodes": {k: v for k, v in resolved["barcodes"].items() if k.startswith("n_")},
        "groups": resolved["groups"],
        "cells": cells,
        "called_cells_overlap": cells["n_called_cells_seen"] / cells["n_called_cells"]
        if cells["n_called_cells"]
        else None,
        "fractions": {
            "kraken2_classified": kraken["n_classified"] / n_input if n_input else None,
            "kraken2_unclassified": kraken["n_unclassified"] / n_input if n_input else None,
        },
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


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

    from viralscan.run_safety import software_identity
    from viralscan.validation import tool_provenance

    called, cells_info = load_called_cells(sample_dir)
    onlist = resolve_onlist(config.technology, config.whitelist)
    bustools = tool_provenance(("bustools",))["bustools"]
    manifest = json.loads((run_root / RUN_MANIFEST).read_text(encoding="utf-8"))
    version = kraken2_version()
    identity = db_identity(db)
    software = software_identity()
    provenance: dict[str, Any] = {
        "schema_version": 1,
        "primary_run_fingerprint": manifest.get("run_fingerprint"),
        "viralscan": software,
        "kraken2": {
            "version": version,
            "confidence": args.confidence,
            "minimum_hit_groups": args.minimum_hit_groups,
            "memory_mapping": args.memory_mapping,
        },
        "polya_trim": (
            {"min_polya_run": POLYA_MIN_RUN, "min_polyt_run": POLYT_MIN_RUN}
            if args.polya_trim
            else None
        ),
        "db_identity": identity,
        "bustools": bustools,
        "onlist": onlist,
        "onlist_sha256": sha256_file(onlist) if onlist else None,
        "technology": config.technology,
        "input_source": _input_identity(sample_dir, run_root, args),
        "host_subtract": args.host_subtract,
        **{k: cells_info[k] for k in ("called_cells_sha256", "called_cells_barcode_space")},
    }
    # memory_mapping, threads, paths and the software build change speed or location, not output.
    provenance["fingerprint_inputs"] = {
        "primary_run_fingerprint": provenance["primary_run_fingerprint"],
        "viralscan_version": software["version"],
        "kraken2": {k: v for k, v in provenance["kraken2"].items() if k != "memory_mapping"},
        "polya_trim": provenance["polya_trim"],
        "db_sha256": identity["sha256"],
        "bustools_sha256": bustools["sha256"] if bustools else None,
        "onlist_sha256": provenance["onlist_sha256"],
        "technology": config.technology,
        "input_source": copy.deepcopy(provenance["input_source"]),
        "called_cells_sha256": cells_info["called_cells_sha256"],
        "called_cells_barcode_space": cells_info["called_cells_barcode_space"],
    }
    fingerprint = discovery_fingerprint(provenance)
    marker = output / COMPLETE_JSON
    if marker.is_file():
        try:
            previous = json.loads(marker.read_text(encoding="utf-8"))
        except ValueError as exc:
            raise DiscoveryError(f"{marker} is unreadable: {exc}") from exc
        if previous.get("discovery_fingerprint") == fingerprint:
            log.info("Discovery already published with the same fingerprint; nothing to do.")
            return
        if not args.force:
            raise DiscoveryError(
                f"{output} holds a discovery with a different fingerprint "
                f"({previous.get('discovery_fingerprint')} != {fingerprint}). "
                "Use a new --output, or --force to replace it."
            )
        _forget_published(output, previous)
    partial = stage(output)
    resolved = resolve_input(run_root, sample_dir, partial, args)
    resolved["kb_input_check"] = kb_input_check(config, resolved)
    rf_lineage, default = read_filter_lineage_for(sample_dir, resolved, config.read_filter)
    resolved["stream"] = stream_reads(
        resolved["r1"],
        resolved["r2"],
        cb_len,
        umi_len,
        partial,
        rf_lineage,
        default,
        args.polya_trim,
    )
    if resolved["stream"]["n_input"] != resolved["n_fragments"]:
        raise DiscoveryError(
            f"Streamed {resolved['stream']} but resolved {resolved['n_fragments']} fragments."
        )
    resolved.update(technology=config.technology, cb_len=cb_len, umi_len=umi_len)
    log.info("Discovery input: %s, %s.", resolved["source"], resolved["stream"])
    plain = run_kraken2(
        db, partial, args.threads, args.confidence, args.minimum_hit_groups, args.memory_mapping
    )
    resolved["kraken2"] = {
        "version": version,
        "confidence": args.confidence,
        "minimum_hit_groups": args.minimum_hit_groups,
        "memory_mapping": args.memory_mapping,
        **join_kraken2(plain, resolved["stream"]["n_input"], partial),
    }
    log.info("Kraken2: %s.", resolved["kraken2"])

    tree = load_taxonomy(db, args.memory_mapping)
    write_taxonomy(tree, partial)
    resolved["barcodes"] = {
        "onlist": onlist,
        "onlist_sha256": provenance["onlist_sha256"],
        "bustools": bustools,
        **correct_barcodes(partial, onlist, args.threads),
    }
    resolved["groups"] = build_umi_groups(partial, tree)
    eligible = sum(resolved["barcodes"][f"n_cb_{k}"] for k in ("onlist", "corrected", "no_onlist"))
    if resolved["groups"]["n_group_reads"] != eligible:
        raise DiscoveryError(
            f"UMI groups hold {resolved['groups']['n_group_reads']} reads, not {eligible}."
        )
    log.info("Barcodes: %s. Groups: %s.", resolved["barcodes"], resolved["groups"])
    resolved["cells"] = {**cells_info, **build_cell_matrix(partial, tree, called)}
    log.info("Cells: %s.", resolved["cells"])
    resolved["cells"]["n_candidates"] = build_candidates(partial, tree, called)
    provenance["taxonomy_sha256"] = sha256_file(partial / TAXONOMY_TSV)
    provenance["resolved_input"] = {
        k: v
        for k, v in resolved.items()
        if k in ("r1", "r2", "n_fragments", "primary_verified", "kb_input_check")
    }
    provenance["discovery_fingerprint"] = fingerprint
    _write_json(partial / PROVENANCE_JSON, provenance)
    _write_json(partial / QC_JSON, build_qc(resolved))
    _write_json(
        partial / COMPLETE_JSON,
        {
            "discovery_fingerprint": fingerprint,
            "outputs": {p.name: sha256_path(p) for p in sorted(partial.iterdir())},
        },
    )
    publish(output)
    log.info("Discovery published to %s (fingerprint %s).", output, fingerprint)
