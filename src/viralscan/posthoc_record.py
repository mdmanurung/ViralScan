"""One JSON record per post-hoc invocation, and an index rebuilt from them (PLAN ``EVID-01``).

Every post-hoc verb (``evidence``, ``rerun-programs``, ``discover`` ...) runs inside
:func:`attempt`. It writes ``<attempt_id>.json`` with state ``started`` before the verb
runs and rewrites it with a terminal state (``succeeded``, ``failed``, ``interrupted``)
afterwards. A record left at ``started`` is reported ``interrupted`` once its process is
gone (same host) or ``unfinished`` otherwise. ``index.json`` is a view rebuilt from the
records after every attempt and by ``viralscan status``; it is never the source of truth.

Where records go:

* output inside the source run (the default ``<run>/posthoc/...``): ``<run root>/posthoc/records/``;
* output outside the run (``--output elsewhere``): ``<output>.records/`` beside it, so a
  run that the user keeps untouched (or cannot write) is never written to.

A verb pointed at a directory with no ``run_manifest.json`` above it leaves no record.
Recording never changes what a verb does: a record that cannot be written is logged and
the verb still runs. Records hold paths, settings, code identity and hashes; no reads.
"""

from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import json
import logging
import os
import secrets
import socket
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from viralscan.outputs import POSTHOC_DIRNAME

log = logging.getLogger(__name__)

SCHEMA = "viralscan.posthoc-record.v1"
RECORDS_DIRNAME = "records"
INDEX_NAME = "index.json"
EXTERNAL_SUFFIX = ".records"
TERMINAL_STATES = ("succeeded", "failed", "interrupted")
# outputs at most this large are hashed; larger files (BAMs, h5ad) record size and mtime only
MAX_HASH_BYTES = 64 * 1024 * 1024
MAX_OUTPUT_FILES = 10_000
# exit codes that mean the process was stopped, not that the verb failed (SIGINT, SIGTERM)
INTERRUPT_EXIT_CODES = (130, 143)

# verb -> argparse dest of the path it reads from (the run, a sample, or a post-hoc folder)
SOURCE_ARG: dict[str, str] = {
    "evidence": "run_dir",
    "evidence-plot": "evidence_dir",
    "rerun-programs": "run_dir",
    "host-loss": "run_dir",
    "rerun-multimap": "run_dir",
    "discover": "run_dir",
    "discover-states": "discovery_dir",
    "discover-collapse": "discovery_dir",
    "hostresponse": "run_dir",
    "translated": "run_dir",
}
VERBS = tuple(SOURCE_ARG)


def find_source_root(path: Path | None) -> Path | None:
    """The run root above *path*: the nearest directory (itself included) holding ``run_manifest.json``."""
    if path is None:
        return None
    path = Path(path).resolve()
    for candidate in (path, *path.parents):
        if (candidate / "run_manifest.json").is_file():
            return candidate
    return None


def _requested_output(args: argparse.Namespace) -> Path | None:
    value = getattr(args, "output", None) or getattr(args, "output_file", None)
    return Path(value).resolve() if value else None


def records_dir(source_root: Path | None, output: Path | None) -> Path | None:
    """Where an attempt's record goes; ``None`` (no record) when the source is not a run."""
    if source_root is None:
        return None
    if output is not None and not output.is_relative_to(source_root):
        return output.parent / f"{output.name}{EXTERNAL_SUFFIX}"
    return source_root / POSTHOC_DIRNAME / RECORDS_DIRNAME


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def new_attempt_id(verb: str) -> str:
    """Sortable, filename-safe and unique under concurrent starts: ``<UTC time>-<verb>-<8 hex>``."""
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}-{verb}-{secrets.token_hex(4)}"


def _sha256(path: Path) -> str:
    from viralscan.run_safety import sha256_file

    return sha256_file(path)


def _file_entry(path: Path, root: Path | None = None) -> dict[str, Any]:
    stat = path.stat()
    entry: dict[str, Any] = {
        "path": path.relative_to(root).as_posix() if root else str(path),
        "bytes": stat.st_size,
        "mtime": dt.datetime.fromtimestamp(stat.st_mtime, dt.timezone.utc).isoformat(
            timespec="seconds"
        ),
        "sha256": _sha256(path) if stat.st_size <= MAX_HASH_BYTES else None,
    }
    return entry


def _inputs(source_path: Path | None, source_root: Path | None) -> list[dict[str, Any]]:
    """The provenance files the verb starts from: run manifest and the sample's config."""
    found: list[Path] = []
    for directory in (source_root, source_path):
        if directory is None:
            continue
        for name in ("run_manifest.json", "config.yaml", "evidence_manifest.json"):
            candidate = Path(directory) / name
            if candidate.is_file() and candidate.resolve() not in found:
                found.append(candidate.resolve())
    return [_file_entry(path) for path in found]


def _outputs(output: Path | None) -> dict[str, Any]:
    if output is None or not output.exists():
        return {"files": [], "truncated": False}
    if output.is_file():
        return {"files": [_file_entry(output)], "truncated": False}
    files = sorted(p for p in output.rglob("*") if p.is_file())
    return {
        "files": [_file_entry(p, output) for p in files[:MAX_OUTPUT_FILES]],
        "truncated": len(files) > MAX_OUTPUT_FILES,
    }


def _settings(args: argparse.Namespace) -> dict[str, Any]:
    """The parsed options, JSON-safe; private (``_``) attributes are left out."""
    settings: dict[str, Any] = {}
    for key, value in sorted(vars(args).items()):
        if key.startswith("_") or callable(value):
            continue
        settings[key] = (
            value if isinstance(value, (str, int, float, bool, type(None))) else str(value)
        )
    return settings


def _write_json(path: Path, payload: Any) -> None:
    """Write *payload* to *path* atomically: a temp file in the same directory, then rename."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.parent / f".{path.name}.{os.getpid()}.{secrets.token_hex(4)}.tmp"
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    os.replace(tmp, path)


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def effective_state(record: dict[str, Any]) -> str:
    """The record's state, with an abandoned ``started`` read as ``interrupted`` or ``unfinished``."""
    state = str(record.get("state"))
    if state != "started":
        return state
    if record.get("host") == socket.gethostname() and not _pid_alive(int(record.get("pid", -1))):
        return "interrupted"
    return "unfinished"


def read_records(directory: Path) -> list[dict[str, Any]]:
    """Every readable record in *directory*, oldest first; unreadable files are skipped and logged."""
    records = []
    for path in sorted(directory.glob("*.json")) if directory.is_dir() else []:
        if path.name == INDEX_NAME:
            continue
        try:
            record = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            log.warning("Skipping unreadable post-hoc record %s: %s", path, exc)
            continue
        if isinstance(record, dict) and record.get("schema") == SCHEMA:
            records.append(record)
    return records


def rebuild_index(directory: Path) -> dict[str, Any]:
    """Rewrite ``index.json`` in *directory* from its records and return it."""
    entries = [
        {
            "attempt_id": r.get("attempt_id"),
            "verb": r.get("verb"),
            "state": effective_state(r),
            "exit_code": r.get("exit_code"),
            "sample": r.get("sample"),
            "target": r.get("target"),
            "output_path": r.get("output_path"),
            "started_at": r.get("started_at"),
            "ended_at": r.get("ended_at"),
        }
        for r in read_records(directory)
    ]
    index = {"schema": SCHEMA.replace("record", "index"), "built_at": _now(), "attempts": entries}
    _write_json(directory / INDEX_NAME, index)
    return index


class Attempt:
    """One post-hoc invocation's record; written at start and again at the end."""

    def __init__(self, verb: str, args: argparse.Namespace, argv: list[str]) -> None:
        source_value = getattr(args, SOURCE_ARG[verb], None)
        if isinstance(source_value, (list, tuple)):  # evidence-plot takes several --evidence-dir
            source_value = source_value[0] if source_value else None
        self.args = args
        self.source_path = Path(source_value).resolve() if source_value else None
        self.source_root = find_source_root(self.source_path)
        requested = _requested_output(args)
        self.directory = records_dir(self.source_root, requested)
        sample = getattr(args, "sample", None)
        if (
            not sample
            and self.source_path
            and self.source_root
            and self.source_path != self.source_root
        ):
            sample = self.source_path.relative_to(self.source_root).parts[0]
        self.record: dict[str, Any] = {
            "schema": SCHEMA,
            "attempt_id": new_attempt_id(verb),
            "verb": verb,
            "argv": argv,
            "state": "started",
            "exit_code": None,
            "error": None,
            "started_at": _now(),
            "ended_at": None,
            "host": socket.gethostname(),
            "pid": os.getpid(),
            "software": {},
            "source_path": str(self.source_path) if self.source_path else None,
            "source_root": str(self.source_root) if self.source_root else None,
            "sample": sample,
            "target": getattr(args, "target", None),
            "requested_output": str(requested) if requested else None,
            "output_path": None,
            "output_preexisting": bool(requested and requested.exists()),
            "settings": _settings(args),
            "inputs": [],
            "outputs": {"files": [], "truncated": False},
        }

    @property
    def path(self) -> Path | None:
        return self.directory / f"{self.record['attempt_id']}.json" if self.directory else None

    def _save(self) -> None:
        if self.path is None:
            return
        try:
            _write_json(self.path, self.record)
        except OSError as exc:
            log.warning("Could not write post-hoc record %s: %s", self.path, exc)

    def start(self) -> None:
        from viralscan.run_safety import software_identity

        try:
            self.record["software"] = software_identity()
            self.record["inputs"] = _inputs(self.source_path, self.source_root)
        except OSError as exc:
            self.record["error"] = f"provenance unavailable: {type(exc).__name__}"
        self._save()

    def finish(self, state: str, exit_code: int | None, error: str | None) -> None:
        output = getattr(self.args, "_posthoc_output", None) or getattr(self.args, "output", None)
        output_path = Path(output).resolve() if output else None
        self.record.update(
            state=state,
            exit_code=exit_code,
            error=error,
            ended_at=_now(),
            output_path=str(output_path) if output_path else None,
        )
        try:
            self.record["outputs"] = _outputs(output_path)
        except OSError as exc:
            self.record["outputs"] = {"files": [], "truncated": False, "error": str(exc)}
        self._save()
        if self.directory is not None:
            try:
                rebuild_index(self.directory)
            except OSError as exc:
                log.warning("Could not rebuild %s: %s", self.directory / INDEX_NAME, exc)


@contextlib.contextmanager
def attempt(
    verb: str, args: argparse.Namespace, argv: list[str] | None = None
) -> Iterator[Attempt]:
    """Run one post-hoc verb inside a record: ``started`` now, a terminal state when it ends."""
    current = Attempt(verb, args, list(sys.argv[1:] if argv is None else argv))
    current.start()
    try:
        yield current
    except KeyboardInterrupt:
        current.finish("interrupted", 130, "KeyboardInterrupt")
        raise
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
        if code == 0:
            current.finish("succeeded", 0, None)
        else:
            state = "interrupted" if code in INTERRUPT_EXIT_CODES else "failed"
            current.finish(state, code, None if isinstance(exc.code, int) else str(exc.code))
        raise
    except BaseException as exc:
        current.finish("failed", 1, f"{type(exc).__name__}: {exc}")
        raise
    else:
        current.finish("succeeded", 0, None)


def status_table(run_dir: Path) -> str:
    """Rebuild the run's index and render it as a plain table for ``viralscan status``."""
    root = find_source_root(run_dir) or Path(run_dir).resolve()
    directory = root / POSTHOC_DIRNAME / RECORDS_DIRNAME
    index = rebuild_index(directory) if directory.is_dir() else {"attempts": []}
    rows = [("attempt", "verb", "state", "sample", "target", "output")]
    rows += [
        (
            str(a["attempt_id"]),
            str(a["verb"]),
            str(a["state"]) + (f" ({a['exit_code']})" if a["exit_code"] not in (None, 0) else ""),
            str(a["sample"] or "-"),
            str(a["target"] or "-"),
            str(a["output_path"] or "-"),
        )
        for a in index["attempts"]
    ]
    if len(rows) == 1:
        return f"No post-hoc records under {directory}.\n"
    widths = [max(len(row[i]) for row in rows) for i in range(len(rows[0]) - 1)]
    lines = [
        "  ".join(cell.ljust(w) for cell, w in zip(row[:-1], widths, strict=True)) + "  " + row[-1]
        for row in rows
    ]
    return "\n".join(lines) + "\n"
