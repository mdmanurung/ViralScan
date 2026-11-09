"""V3 run fingerprints and non-destructive output-directory handling."""

from __future__ import annotations

import atexit
import csv
import hashlib
import json
import logging
import os
import shutil
import signal
import socket
import subprocess
import sys
from collections.abc import Callable
from contextlib import suppress
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from viralscan import __version__
from viralscan.defaults import DEFAULTS
from viralscan.outputs import LAYOUT_VERSION

log = logging.getLogger(__name__)


def run_in_group(
    command: list[str], *, check: bool = True, capture_output: bool = False
) -> subprocess.CompletedProcess[bytes]:
    """Run *command* in its own session, like ``subprocess.run``.

    SIGINT/SIGTERM reaching viralscan is forwarded to that whole process group (the
    tool and every child it started, e.g. kallisto under kb), so none is orphaned.
    SIGTERM exits 143. Call from the main thread only (``signal.signal``).
    """
    pipe = subprocess.PIPE if capture_output else None
    proc = subprocess.Popen(command, start_new_session=True, stdout=pipe, stderr=pipe)

    def _on_sigterm(_signum: int, _frame: Any) -> None:
        sys.exit(143)

    previous = signal.signal(signal.SIGTERM, _on_sigterm)
    try:
        stdout, stderr = proc.communicate()
    except (KeyboardInterrupt, SystemExit) as exc:
        sent = signal.SIGINT if isinstance(exc, KeyboardInterrupt) else signal.SIGTERM
        try:
            os.killpg(proc.pid, sent)
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
        except ProcessLookupError:
            pass
        raise
    finally:
        signal.signal(signal.SIGTERM, previous)
    if check and proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, command, stdout, stderr)
    return subprocess.CompletedProcess(command, proc.returncode, stdout, stderr)


RUN_MANIFEST = "run_manifest.json"
RUN_COMPLETE = "run_complete.json"
LOCK_NAME = ".viralscan.lock"
COMBINED_SUMMARY = "combined_summary.tsv"
#: Per-sample artifacts whose sha256 the completion marker records (each only if
#: present). Deliberately narrow: the headline tables and the count matrix, not
#: logs, plots, or every intermediate file.
MARKER_ARTIFACTS = (
    "results/viral_summary.tsv",
    "results/virus_identity.tsv",
    "results/multimap_evidence.tsv",
    "kb-python/counts_unfiltered/adata_multimap.h5ad",
)


#: Per-sample count matrix every downstream rule reads; log/kb.done alone does not prove it.
KB_COUNTS = "kb-python/counts_unfiltered/adata.h5ad"


class RunSafetyError(RuntimeError):
    """Raised when output reuse would be unsafe or irreproducible."""


#: "provisional" until G5 validates the detection defaults; "final" from rc (DEF-09).
DEFAULTS_STATUS = "provisional"


def sha256_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_path(path: str | Path) -> str:
    """sha256 of a file, or of a directory (sorted relative paths + file hashes)."""
    path = Path(path)
    if not path.is_dir():
        return sha256_file(path)
    digest = hashlib.sha256()
    for child in sorted(p for p in path.rglob("*") if p.is_file()):
        digest.update(child.relative_to(path).as_posix().encode() + b"\0")
        digest.update(bytes.fromhex(sha256_file(child)))
    return digest.hexdigest()


def software_identity() -> dict[str, Any]:
    """Describe the executing package bytes and, in its own checkout, Git revision.

    Runtime files are hashed in relative-path order; ignored local GTFs and bytecode
    are excluded because they are reference inputs, not the installed software.
    Provenance failure is explicit and never prevents a run.
    """
    package = Path(__file__).resolve().parent
    identity: dict[str, Any] = {"version": __version__, "kind": "installed"}
    digest = hashlib.sha256()
    try:
        for path in sorted(package.rglob("*")):
            if path.is_file() and (
                path.suffix in {".py", ".tsv", ".json", ".j2", ".R"} or path.name == "Snakefile"
            ):
                digest.update(path.relative_to(package).as_posix().encode() + b"\0")
                digest.update(bytes.fromhex(sha256_file(path)))
        identity["build_sha256"] = digest.hexdigest()
    except OSError as exc:
        identity["build_status"] = f"unavailable: {type(exc).__name__}"
    root = package.parent.parent
    if (root / ".git").exists() and (root / "src" / "viralscan").resolve() == package:
        try:
            commit = subprocess.run(
                ["git", "-C", str(root), "rev-parse", "HEAD"],
                check=True,
                capture_output=True,
                text=True,
                timeout=5,
            ).stdout.strip()
            dirty = subprocess.run(
                ["git", "-C", str(root), "status", "--porcelain", "--", "src/viralscan"],
                check=True,
                capture_output=True,
                text=True,
                timeout=5,
            ).stdout.strip()
            identity.update(kind="checkout", git_commit=commit, source_dirty=bool(dirty))
        except (OSError, subprocess.SubprocessError):
            identity["git_status"] = "unavailable"
    return identity


# C-12: resource options are recorded in the manifest but never hashed, so a SLURM
# resubmission with another allocation can still --resume.
_INFORMATIONAL_OPTIONS = frozenset({"cores"})
_HASHED_KEYS = (
    "schema_version",
    "viralscan_version",
    "quantification_unit",
    "allocation_method",
    "reference_fingerprint",
    "cell_calling",
    "input_fingerprints",
    "options",
)


def _fingerprint(manifest: dict[str, Any], *, informational: bool = False) -> str:
    """Hash the fingerprinted fields; ``informational=True`` hashes the pre-C-12 way."""
    payload = {key: manifest[key] for key in _HASHED_KEYS if key in manifest}
    if not informational:
        payload["options"] = {
            k: v
            for k, v in (payload.get("options") or {}).items()
            if k not in _INFORMATIONAL_OPTIONS
        }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(canonical).hexdigest()


def _fingerprint_matches(previous: dict[str, Any], manifest: dict[str, Any]) -> bool:
    """True when ``previous`` is the same run; pre-C-12 manifests hashed ``cores``."""
    stored = previous.get("run_fingerprint")
    if stored == manifest.get("run_fingerprint"):
        return True
    # Old-style manifest: trust it only if its stored hash reproduces with cores
    # included, then compare it to this invocation with cores left out.
    return stored == _fingerprint(previous, informational=True) and _fingerprint(
        previous
    ) == manifest.get("run_fingerprint")


def build_run_manifest(args: Any) -> dict[str, Any]:
    """Build a canonical manifest from scientific options and input bytes."""
    path_fields = (
        "sample1",
        "sample2",
        "index",
        "transcripts",
        "gtf",
        "fasta",
        "f1",
        "whitelist",
        "called_cells_file",
        "host_index",
        "positive_control_receipt",
        # C-13: previously only path strings (host_index, a STAR directory, was skipped).
        "host_h5ad",
        "cell_types",
        "host_fasta",
        "viral_fasta",
    )
    input_fingerprints: dict[str, str] = {}
    input_paths: dict[str, str] = {}
    for field in path_fields:
        raw = getattr(args, field, None)
        if not raw:
            continue
        for i, value in enumerate(str(raw).split(",")):
            path = Path(value).expanduser().resolve()
            if path.is_file() or (field == "host_index" and path.is_dir()):
                input_fingerprints[f"{field}:{i}"] = sha256_path(path)
                input_paths[f"{field}:{i}"] = str(path)

    excluded = {
        "output",
        "yes",
        "resume",
        "overwrite",
        "verbose",
        "quiet",
        "_subcommand",
        # Overrides the chemistry check only; the technology itself is fingerprinted.
        "force_technology",
        "require_chemistry_sanity",
        # Credentials: never written to the manifest (AUD-09 / C-02).
        "ncbi_api_key",
        "ncbi_email",
    }
    # Options added after v3.0 manifests were first written. An unset (None) value
    # is omitted, so a manifest that predates the option still matches an
    # invocation that does not use it; any non-None value changes the fingerprint
    # and refuses --resume against such a manifest (old counts were made with
    # kb's default, not this value).
    omit_when_unset = {
        "strand",
        "host_filter_star_params",
        "read_filter",
        "auto_evidence",
        "host_fasta",
        "viral_fasta",
        "positive_control_scope",
        "positive_control_virus_key",
        "positive_control_receipt",
        "hostresponse_cv",
        "hostresponse_groups",
        "hostresponse_cv_folds",
        "hostresponse_cell_type_column",
        "hostresponse_cell_types",
        "hostresponse_panel_in_fold",
        "hostresponse_permutations",
        "hostresponse_permutation_unit",
        "hostresponse_permutation_block",
        "hostresponse_min_negative_cells",
        "hostresponse_min_groups",
        # AUD-06/AUD-07/B-16: unset keeps the pre-existing fingerprint.
        # `infected_min_molecules` is not here: the CLI always sets it (default 2.0,
        # INFECTED-MINMOL-02), so it is hashed and a manifest without it refuses --resume.
        "detection_min_cells",
        "hostresponse_min_umi",
        # PROG-MIXED-VETO-01
        "programme_mixed_veto",
        "programme_mixed_veto_alpha",
        # Label only (reference_manifest.json); default None. Omitted when unset so
        # manifests from before the flag still resume; a set value stays hashed.
        "reference_host_species",
    }
    # Options added later with a non-None default: the default value is omitted, so
    # a manifest that predates the option still matches (PROG-MINMOL-01, default 1).
    omit_when_default = {
        "programme_productive_min_molecules": DEFAULTS["programme_productive_min_molecules"],
    }
    options = {
        key: value
        for key, value in sorted(vars(args).items())
        if key not in excluded
        and not (key == "scautoqc" and value is False)
        and (
            isinstance(value, (str, int, float, bool, type(None)))
            or (
                key == "hostresponse_cell_types"
                and isinstance(value, list)
                and all(isinstance(v, str) for v in value)
            )
        )
        and not (key in omit_when_unset and value is None)
        and not (key in omit_when_default and value in (None, omit_when_default[key]))
    }
    # Pre-fix manifests hashed these two as null when unset; keep a null so a run made
    # without them still resumes. The values themselves are not recorded.
    options["ncbi_api_key"] = None
    options["ncbi_email"] = None
    reference_hashes = {
        key: value
        for key, value in input_fingerprints.items()
        if key.split(":", 1)[0] in {"index", "transcripts", "gtf", "fasta", "f1", "host_index"}
    }
    reference_fingerprint = hashlib.sha256(
        json.dumps(reference_hashes, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    payload = {
        "schema_version": "3.0.0",
        "viralscan_version": __version__,
        "quantification_unit": "bustools-resolved-cb-umi-molecule",
        "allocation_method": getattr(args, "multimap_method", None),
        "reference_fingerprint": reference_fingerprint,
        "cell_calling": {
            "method": getattr(args, "cell_calling", None),
            "external_file": getattr(args, "called_cells_file", None),
        },
        "input_fingerprints": input_fingerprints,
        "options": options,
    }
    payload["run_fingerprint"] = _fingerprint(payload)
    # Added after hashing so --resume still matches manifests written without them.
    payload["completion_marker"] = True
    payload["layout_version"] = LAYOUT_VERSION
    # Resolved at build time so validate-run finds the inputs from any cwd (C-04).
    payload["input_paths"] = input_paths
    payload["defaults_status"] = DEFAULTS_STATUS
    payload["software_identity"] = software_identity()
    from viralscan import virus_catalog
    from viralscan.virus_identity import load_build_manifest, manifest_path_for_index

    catalogue = virus_catalog.catalogue_path()
    observed = sha256_file(catalogue) if catalogue.is_file() else None
    expected = None
    index = getattr(args, "index", None)
    if index:
        # Analysis validates the build manifest; absence cannot establish a match.
        with suppress(OSError, ValueError):
            expected = load_build_manifest(manifest_path_for_index(index)).catalogue_sha256
    payload["reference_catalogue"] = {
        "expected_sha256": expected,
        "observed_sha256": observed,
        "matches": expected == observed if expected and observed else None,
    }
    return payload


def build_run_complete(run_root: Path) -> dict[str, Any]:
    """Describe a finished run: its fingerprint and the hashes of MARKER_ARTIFACTS."""
    run_root = Path(run_root)
    manifest = json.loads((run_root / RUN_MANIFEST).read_text(encoding="utf-8"))
    samples = sorted(p.parent.name for p in run_root.glob("*/config.yaml"))
    artifacts = {
        f"{sample}/{rel}": sha256_file(run_root / sample / rel)
        for sample in samples
        for rel in MARKER_ARTIFACTS
        if (run_root / sample / rel).is_file()
    }
    return {
        "schema_version": "3.0.0",
        "run_fingerprint": manifest.get("run_fingerprint"),
        "viralscan_version": __version__,
        "samples": samples,
        "artifacts": artifacts,
        "completed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def write_combined_summary(run_root: Path) -> bool:
    """Stack every sample's ``results/viral_summary.tsv`` into ``combined_summary.tsv``.

    A ``sample`` column leads; no value is recomputed. A valid no-virus run has
    summaries with a header and no rows; it still gets a header-only file. Returns
    False (and writes nothing) only when no sample has a viral summary file.
    """
    run_root = Path(run_root)
    rows: list[dict[str, str]] = []
    columns: list[str] = []
    seen = False
    for config in sorted(run_root.glob("*/config.yaml")):
        sample = config.parent.name
        path = config.parent / "results" / "viral_summary.tsv"
        if not path.is_file():
            continue
        seen = True
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns += [c for c in reader.fieldnames or [] if c not in columns and c != "sample"]
            rows += [{"sample": sample, **row} for row in reader]
    if not seen:
        return False
    staging = run_root / f".{COMBINED_SUMMARY}.tmp"
    with staging.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", *columns], delimiter="\t", restval="")
        writer.writeheader()
        writer.writerows(rows)
    staging.replace(run_root / COMBINED_SUMMARY)
    return True


def write_run_complete(run_root: Path) -> None:
    """Atomically (re)write ``run_complete.json`` at the run root."""
    from viralscan.artifact_inventory import write_layout_docs

    run_root = Path(run_root)
    write_layout_docs(run_root)  # README.txt files + results/artifact_inventory.json per sample
    staging = run_root / f".{RUN_COMPLETE}.tmp"
    staging.write_text(
        json.dumps(build_run_complete(run_root), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    staging.replace(run_root / RUN_COMPLETE)
    write_combined_summary(run_root)


def refuse_primary_tree(output: Path) -> None:
    """Refuse an ``output`` inside a run's or sample's tree unless it sits under ``posthoc/``."""
    from viralscan.outputs import POSTHOC_DIRNAME

    out = Path(output).resolve()
    for root in (out, *out.parents):
        # A bare config.yaml is common in any project dir; a sample also has kb-python/.
        sample = (root / "config.yaml").is_file() and (root / "kb-python").is_dir()
        if sample or (root / RUN_MANIFEST).is_file():
            if POSTHOC_DIRNAME not in out.relative_to(root).parts:
                raise RunSafetyError(
                    f"Output {output} is inside the ViralScan run {root}. Post-hoc results go "
                    f"under {root / POSTHOC_DIRNAME}/ or in a directory outside the run."
                )
            return


def require_fresh_output(output: Path, *, force: bool = False, forceable: bool = True) -> None:
    """The one overwrite guard of the post-hoc verbs: refuse a non-empty ``output``.

    ``force`` lets the verb write into it (same-named files are replaced). The verb
    itself decides what else it clears; ``hostresponse`` backs the old tree up and
    restores it on failure. ``forceable=False`` is for a verb that has no ``--force``.
    Any ``output`` inside a run's own tree (outside ``posthoc/``) is refused, forced or not.
    """
    refuse_primary_tree(output)
    if force or not output.exists() or (output.is_dir() and not any(output.iterdir())):
        return
    way = ", or pass --force to write into it (same-named files are replaced)" if forceable else ""
    raise RunSafetyError(f"Output {output} already holds files. Use a new --output{way}.")


def clear_run_complete(run_root: Path) -> None:
    (Path(run_root) / RUN_COMPLETE).unlink(missing_ok=True)


def restamp_run_complete(run_root: Path) -> bool:
    """Recompute the marker after an in-place mutation; no-op if the run had none."""
    if not (Path(run_root) / RUN_COMPLETE).is_file():
        return False
    write_run_complete(run_root)
    return True


def _write_manifest_atomic(output_dir: Path, manifest: dict[str, Any]) -> None:
    target = output_dir / RUN_MANIFEST
    staging = output_dir / f".{RUN_MANIFEST}.tmp"
    staging.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    staging.replace(target)


def recorded_strand_inference(output_dir: Path, sample: str) -> Optional[dict[str, Any]]:
    """Return the recorded ``--strand auto`` block for a sample, if any."""
    path = Path(output_dir) / RUN_MANIFEST
    if not path.is_file():
        return None
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if (manifest.get("options") or {}).get("strand") != "auto":
        return None
    return (manifest.get("strand_inference") or {}).get(sample)


def record_manifest_block(output_dir: Path, key: str, sample: str, block: dict[str, Any]) -> None:
    """Add a sample's ``key`` block to the manifest, outside the run fingerprint."""
    output_dir = Path(output_dir)
    manifest = json.loads((output_dir / RUN_MANIFEST).read_text(encoding="utf-8"))
    manifest.setdefault(key, {})[sample] = block
    _write_manifest_atomic(output_dir, manifest)


def record_resume_attempt(output_dir: Path, blocks: dict[str, dict[str, Any]]) -> None:
    """Log a ``--resume`` in ``resume_attempts`` and keep the first attempt's per-sample blocks.

    ``blocks`` maps a manifest key (``chemistry_detection``, ``tool_binaries``) to
    ``{sample: block}``. A block is only written where the manifest has none yet.
    """
    output_dir = Path(output_dir)
    manifest = json.loads((output_dir / RUN_MANIFEST).read_text(encoding="utf-8"))
    for key, per_sample in blocks.items():
        for sample, block in per_sample.items():
            manifest.setdefault(key, {}).setdefault(sample, block)
    manifest.setdefault("resume_attempts", []).append(
        {"at": datetime.now(timezone.utc).isoformat(timespec="seconds"), **blocks}
    )
    _write_manifest_atomic(output_dir, manifest)


def record_fetched_reference(
    output_dir: Path, accessions: list[str], fasta: Path, gtf: Path
) -> None:
    """Record the sha256 of an ``-acc`` run's fetched FASTA/GTF (C-13).

    These bytes do not exist when the fingerprint is computed. On ``--resume`` a
    hash that differs from the one recorded by the first attempt means the
    accession's content changed upstream, so existing counts are not reusable.
    """
    block = {
        "accessions": list(accessions),
        "fasta_sha256": sha256_file(fasta),
        "gtf_sha256": sha256_file(gtf),
    }
    manifest = json.loads((Path(output_dir) / RUN_MANIFEST).read_text(encoding="utf-8"))
    previous = (manifest.get("fetched_reference") or {}).get("run")
    if previous and (
        previous.get("fasta_sha256") != block["fasta_sha256"]
        or previous.get("gtf_sha256") != block["gtf_sha256"]
    ):
        raise RunSafetyError(
            "Cannot resume: the fetched NCBI reference differs from the one this run was "
            "started with (the accession changed upstream). Re-run with --overwrite."
        )
    record_manifest_block(output_dir, "fetched_reference", "run", block)


def record_strand_inference(output_dir: Path, sample: str, block: dict[str, Any]) -> None:
    """Add a sample's strand_inference block (outside the run fingerprint)."""
    record_manifest_block(output_dir, "strand_inference", sample, block)


def _is_litter(path: Path) -> bool:
    """Our own bookkeeping: the run lock and ``.<name>.tmp`` staging files left by a crash."""
    return path.name == LOCK_NAME or (path.name.startswith(".") and path.name.endswith(".tmp"))


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def acquire_run_lock(output_dir: Path) -> None:
    """Take ``<output_dir>/.viralscan.lock`` (pid + host); removed again at interpreter exit.

    A lock whose pid is dead on this host is replaced with a warning; any other existing
    lock is refused. A lock held by this very process is taken over (tests, re-entry).
    ponytail: a stale lock from another host needs a manual delete; no network liveness check.
    """
    output_dir = Path(output_dir)
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise RunSafetyError(f"Cannot create output directory {output_dir}: {exc}") from exc
    lock = output_dir / LOCK_NAME
    host = socket.gethostname()
    for _attempt in range(2):
        try:
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            holder = lock.read_text(encoding="utf-8", errors="replace").split()
            pid = int(holder[0]) if holder and holder[0].isdigit() else None
            same_host = len(holder) > 1 and holder[1] == host
            if pid is not None and same_host and (pid == os.getpid() or not _pid_alive(pid)):
                if pid != os.getpid():
                    log.warning("Replacing stale lock %s (pid %d is not running).", lock, pid)
                lock.unlink(missing_ok=True)
                continue
            raise RunSafetyError(
                f"Another viralscan run holds {lock} ({' on '.join(holder[:2]) or 'unreadable'}). "
                f"If no run is active, delete that file and retry."
            ) from None
        except OSError as exc:
            raise RunSafetyError(f"Cannot lock output directory {output_dir}: {exc}") from exc
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(f"{os.getpid()}\n{host}\n")
        atexit.register(_release_run_lock, lock, os.getpid())
        return
    raise RunSafetyError(f"Cannot take {lock}: it keeps reappearing.")


def _release_run_lock(lock: Path, pid: int) -> None:
    with suppress(OSError, IndexError):
        if lock.read_text(encoding="utf-8").split()[0] == str(pid):
            lock.unlink()


def prepare_output_directory(
    output_dir: Path,
    manifest: dict[str, Any],
    *,
    resume: bool,
    overwrite: bool,
    yes: bool,
    confirm: Callable[[str], str] = input,
) -> str:
    """Create, resume, or explicitly replace one resolved output directory."""
    output_dir = output_dir.resolve()
    if output_dir.exists() and not output_dir.is_dir():
        raise RunSafetyError(f"--output {output_dir} is a file, not a directory.")

    nonempty = output_dir.is_dir() and any(not _is_litter(p) for p in output_dir.iterdir())
    if not nonempty:
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            _write_manifest_atomic(output_dir, manifest)
        except OSError as exc:
            raise RunSafetyError(f"Cannot write to output directory {output_dir}: {exc}") from exc
        return "new"

    if resume:
        manifest_path = output_dir / RUN_MANIFEST
        if not manifest_path.is_file():
            raise RunSafetyError("Cannot resume: run_manifest.json is missing.")
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not _fingerprint_matches(previous, manifest):
            old_opts = previous.get("options") or {}
            new_opts = manifest.get("options") or {}
            old_strand = old_opts.get("strand")
            new_strand = new_opts.get("strand")
            reason = ""
            if old_strand != new_strand:
                reason = (
                    f" (--strand differs: previous run used {old_strand or 'the kb default'!r}, "
                    f"this invocation uses {new_strand or 'the kb default'!r}; "
                    "existing counts are not reusable)"
                )
            # C-12: name every option / input whose value changed, not just --strand.
            changed = sorted(
                k
                for k in old_opts.keys() | new_opts.keys()
                if old_opts.get(k) != new_opts.get(k) and k not in _INFORMATIONAL_OPTIONS
            )
            old_in = previous.get("input_fingerprints") or {}
            new_in = manifest.get("input_fingerprints") or {}
            changed_in = sorted(
                k for k in old_in.keys() | new_in.keys() if old_in.get(k) != new_in.get(k)
            )
            detail = ""
            if changed:
                detail += (
                    " Options that differ: "
                    + ", ".join("--" + k.replace("_", "-") for k in changed)
                    + "."
                )
            if changed_in:
                detail += " Input files that differ: " + ", ".join(changed_in) + "."
            raise RunSafetyError(
                f"Cannot resume: run fingerprint does not match this invocation{reason}.{detail}"
            )
        # AUD3-08: the Snakefile's kb_count declares only log/kb.done, so a purged
        # kb-python/ would be skipped silently and the run "completes" without counts.
        purged = sorted(
            p.parent.parent.name
            for p in output_dir.glob("*/log/kb.done")
            if not (p.parent.parent / KB_COUNTS).is_file()
        )
        if purged:
            raise RunSafetyError(
                f"Cannot resume: {KB_COUNTS} is missing for {', '.join(purged)} although "
                "log/kb.done exists (kb-python/ was deleted?). Delete that sample's "
                "log/kb.done to recount it, or use --overwrite."
            )
        # A resumed run is in progress again: the old marker no longer vouches
        # for it. (New and overwrite starts have no marker: empty dir / wiped.)
        clear_run_complete(output_dir)
        return "resume"

    if not overwrite:
        raise RunSafetyError(
            "Output directory is non-empty. Use --resume for an identical run or "
            "--overwrite to replace it."
        )
    if not yes:
        answer = (
            confirm(f"Overwrite all contents of {output_dir}? Type 'yes' to continue: ")
            .strip()
            .lower()
        )
        if answer != "yes":
            raise RunSafetyError("Overwrite cancelled; output directory was not changed.")

    for child in output_dir.iterdir():
        if child.name == LOCK_NAME:
            continue
        if child.is_dir() and not child.is_symlink():
            shutil.rmtree(child)
        else:
            child.unlink()
    _write_manifest_atomic(output_dir, manifest)
    return "overwrite"
