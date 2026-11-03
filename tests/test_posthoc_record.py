"""Post-hoc records (PLAN EVID-01): one JSON per invocation, an index rebuilt from them."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import socket
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from viralscan import menu, posthoc_record
from viralscan.outputs import posthoc_dir


def _run(tmp_path: Path) -> tuple[Path, Path]:
    root = tmp_path / "run"
    sample = root / "SAMPLE"
    (sample / "results").mkdir(parents=True)
    (root / "run_manifest.json").write_text('{"fingerprint": "abc"}\n')
    (sample / "config.yaml").write_text("technology: 10xv3\n")
    (sample / "results" / "viral_summary.tsv").write_text("virus\tdetected\nHSV-1\tTrue\n")
    return root, sample


def _tree_outside_posthoc(root: Path) -> dict[str, str]:
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file() and "posthoc" not in p.relative_to(root).parts
    }


def _args(sample: Path, output: Path | None = None) -> argparse.Namespace:
    return argparse.Namespace(
        _subcommand="rerun-programs", run_dir=str(sample), output=str(output) if output else None
    )


def _fake_verb(args: argparse.Namespace) -> None:
    """Stands in for a verb: resolve the default output, report it, write one file."""
    out = Path(args.output) if args.output else posthoc_dir(args.run_dir, "programs")
    args._posthoc_output = out
    out.mkdir(parents=True, exist_ok=True)
    (out / "gene_program_summary.tsv").write_text("virus\tn\nHSV-1\t3\n")


def _records(root: Path) -> list[dict]:
    return posthoc_record.read_records(root / "posthoc" / "records")


def test_success_writes_one_record_and_index_and_leaves_run_untouched(tmp_path: Path) -> None:
    root, sample = _run(tmp_path)
    before = _tree_outside_posthoc(root)
    args = _args(sample)
    with posthoc_record.attempt("rerun-programs", args, argv=["rerun-programs"]):
        _fake_verb(args)

    (record,) = _records(root)
    assert record["state"] == "succeeded" and record["exit_code"] == 0
    assert record["sample"] == "SAMPLE"
    assert record["source_root"] == str(root.resolve())
    assert record["output_path"] == str((sample / "posthoc" / "programs").resolve())
    assert [f["path"] for f in record["outputs"]["files"]] == ["gene_program_summary.tsv"]
    assert record["outputs"]["files"][0]["sha256"]
    assert {Path(i["path"]).name for i in record["inputs"]} == {"run_manifest.json", "config.yaml"}
    index = json.loads((root / "posthoc" / "records" / "index.json").read_text())
    assert [a["attempt_id"] for a in index["attempts"]] == [record["attempt_id"]]
    assert _tree_outside_posthoc(root) == before


def test_die_is_recorded_as_failed_and_reraised(tmp_path: Path) -> None:
    root, sample = _run(tmp_path)
    before = _tree_outside_posthoc(root)
    with pytest.raises(SystemExit):
        with posthoc_record.attempt("rerun-programs", _args(sample), argv=[]):
            menu._die("bad flag")
    (record,) = _records(root)
    assert (record["state"], record["exit_code"]) == ("failed", 1)
    assert _tree_outside_posthoc(root) == before


def test_exception_is_recorded_as_failed(tmp_path: Path) -> None:
    root, sample = _run(tmp_path)
    with pytest.raises(ValueError):
        with posthoc_record.attempt("rerun-programs", _args(sample), argv=[]):
            raise ValueError("boom")
    (record,) = _records(root)
    assert record["state"] == "failed" and "ValueError: boom" in record["error"]


@pytest.mark.parametrize(
    "stop, code",
    [(KeyboardInterrupt(), 130), (SystemExit(143), 143)],
    ids=["ctrl-c", "sigterm"],
)
def test_stop_signals_are_recorded_as_interrupted(
    tmp_path: Path, stop: BaseException, code: int
) -> None:
    root, sample = _run(tmp_path)
    with pytest.raises(type(stop)):
        with posthoc_record.attempt("rerun-programs", _args(sample), argv=[]):
            raise stop
    (record,) = _records(root)
    assert (record["state"], record["exit_code"]) == ("interrupted", code)


def test_abandoned_started_record_reads_as_interrupted_or_unfinished(tmp_path: Path) -> None:
    root, sample = _run(tmp_path)
    attempt = posthoc_record.Attempt("rerun-programs", _args(sample), [])
    attempt.start()  # never finished: the process "died"
    dead = subprocess.run(
        [sys.executable, "-c", "import os; print(os.getpid())"],
        capture_output=True,
        text=True,
        check=True,
    )
    record = json.loads(attempt.path.read_text())
    record["pid"] = int(dead.stdout)
    assert posthoc_record.effective_state(record) == "interrupted"
    record["host"] = "another-host"
    assert posthoc_record.effective_state(record) == "unfinished"
    record.update(host=socket.gethostname(), pid=os.getpid())
    assert posthoc_record.effective_state(record) == "unfinished"


def test_external_output_keeps_the_run_untouched(tmp_path: Path) -> None:
    root, sample = _run(tmp_path)
    out = tmp_path / "elsewhere" / "programs"
    args = _args(sample, out)
    with posthoc_record.attempt("rerun-programs", args, argv=[]):
        _fake_verb(args)
    assert not (root / "posthoc").exists() and not (sample / "posthoc").exists()
    (record,) = posthoc_record.read_records(tmp_path / "elsewhere" / "programs.records")
    assert record["state"] == "succeeded" and record["source_root"] == str(root.resolve())


def test_large_outputs_record_size_but_no_hash(tmp_path: Path, monkeypatch) -> None:
    root, sample = _run(tmp_path)
    monkeypatch.setattr(posthoc_record, "MAX_HASH_BYTES", 4)
    args = _args(sample)
    with posthoc_record.attempt("rerun-programs", args, argv=[]):
        _fake_verb(args)
    (entry,) = _records(root)[0]["outputs"]["files"]
    assert entry["sha256"] is None and entry["bytes"] > 4


def test_index_is_rebuilt_from_records_and_skips_unreadable_ones(tmp_path: Path) -> None:
    root, sample = _run(tmp_path)
    for _ in range(2):
        args = _args(sample)
        with posthoc_record.attempt("rerun-programs", args, argv=[]):
            _fake_verb(args)
    directory = root / "posthoc" / "records"
    (directory / "index.json").unlink()
    (directory / "zz-broken.json").write_text("{not json")
    index = posthoc_record.rebuild_index(directory)
    assert len(index["attempts"]) == 2
    assert {a["state"] for a in index["attempts"]} == {"succeeded"}


def test_concurrent_attempts_get_distinct_records(tmp_path: Path) -> None:
    root, sample = _run(tmp_path)
    script = textwrap.dedent(
        f"""
        import argparse
        from viralscan import posthoc_record
        args = argparse.Namespace(_subcommand="host-loss", run_dir={str(sample)!r}, output=None)
        with posthoc_record.attempt("host-loss", args, argv=[]):
            pass
        """
    )
    env = {**os.environ, "PYTHONPATH": str(Path(posthoc_record.__file__).parents[1])}
    procs = [subprocess.Popen([sys.executable, "-c", script], env=env) for _ in range(4)]
    assert [p.wait() for p in procs] == [0, 0, 0, 0]
    records = _records(root)
    assert len({r["attempt_id"] for r in records}) == 4
    assert len(posthoc_record.rebuild_index(root / "posthoc" / "records")["attempts"]) == 4


def test_dispatch_records_every_posthoc_verb_and_no_other(tmp_path: Path, monkeypatch) -> None:
    root, sample = _run(tmp_path)
    seen = []
    monkeypatch.setattr(
        menu, "_dispatch", lambda args: (seen.append(args._subcommand), _fake_verb(args))
    )
    monkeypatch.setattr(menu, "create_help", lambda: _args(sample))
    menu._main()
    assert seen == ["rerun-programs"] and len(_records(root)) == 1

    monkeypatch.setattr(
        menu, "create_help", lambda: argparse.Namespace(_subcommand="doctor", run_dir=str(sample))
    )
    monkeypatch.setattr(menu, "_dispatch", lambda args: seen.append(args._subcommand))
    menu._main()
    assert seen[-1] == "doctor" and len(_records(root)) == 1
    assert set(posthoc_record.VERBS) >= {
        "evidence",
        "discover-collapse",
        "evidence-plot",
        "translated",
    }


def test_status_lists_attempts(tmp_path: Path) -> None:
    root, sample = _run(tmp_path)
    assert "No post-hoc records" in posthoc_record.status_table(root)
    args = _args(sample)
    with posthoc_record.attempt("rerun-programs", args, argv=[]):
        _fake_verb(args)
    table = posthoc_record.status_table(sample)  # a sample dir finds its run root
    assert "rerun-programs" in table and "succeeded" in table and "SAMPLE" in table


def test_list_valued_source_uses_its_first_directory(tmp_path: Path) -> None:
    root, sample = _run(tmp_path)
    evidence = sample / "posthoc" / "evidence" / "HSV-1"
    evidence.mkdir(parents=True)
    args = argparse.Namespace(
        _subcommand="evidence-plot", evidence_dir=[str(evidence)], output=None
    )
    with posthoc_record.attempt("evidence-plot", args, argv=[]):
        pass
    (record,) = _records(root)
    assert record["source_path"] == str(evidence.resolve()) and record["sample"] == "SAMPLE"


def test_no_record_when_the_source_is_not_a_run(tmp_path: Path) -> None:
    out = tmp_path / "out"
    args = argparse.Namespace(
        _subcommand="evidence", run_dir=str(tmp_path / "nope"), output=str(out)
    )
    with posthoc_record.attempt("evidence", args, argv=[]):
        pass
    assert not list(tmp_path.iterdir())
