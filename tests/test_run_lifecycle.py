"""REVIEW-FIX-D: lock, resume marker, signals, crash litter, summary runtime, resume provenance."""

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from viralscan import menu
from viralscan.run_safety import (
    LOCK_NAME,
    RunSafetyError,
    acquire_run_lock,
    prepare_output_directory,
    record_resume_attempt,
)

MANIFEST = {"run_fingerprint": "f", "options": {}}


def _prepare(path: Path, **kw):
    kw = {"resume": False, "overwrite": False, "yes": True, **kw}
    return prepare_output_directory(path, MANIFEST, **kw)


def test_lock_refuses_live_holder_and_replaces_dead_one(tmp_path: Path, caplog) -> None:
    host = os.uname().nodename
    (tmp_path / LOCK_NAME).write_text(f"{os.getppid()}\n{host}\n")  # live, not us
    with pytest.raises(RunSafetyError, match=LOCK_NAME):
        acquire_run_lock(tmp_path)
    dead = subprocess.Popen([sys.executable, "-c", "pass"])
    dead.wait()
    (tmp_path / LOCK_NAME).write_text(f"{dead.pid}\n{host}\n")
    acquire_run_lock(tmp_path)
    assert "stale lock" in caplog.text
    assert (tmp_path / LOCK_NAME).read_text().split()[0] == str(os.getpid())


def test_lock_blocks_second_run_and_is_released_at_exit(tmp_path: Path) -> None:
    code = (
        "import sys;from pathlib import Path;from viralscan.run_safety import *;"
        "d=Path(sys.argv[1]);acquire_run_lock(d);"
        "prepare_output_directory(d,{'run_fingerprint':'f','options':{}},"
        "resume=False,overwrite=False,yes=True);print('held',flush=True);sys.stdin.read()"
    )
    holder = subprocess.Popen(
        [sys.executable, "-c", code, str(tmp_path)], stdin=subprocess.PIPE, stdout=subprocess.PIPE
    )
    try:
        assert holder.stdout.readline().strip() == b"held"
        # menu takes the lock before prepare, so a refused lock means nothing is wiped.
        with pytest.raises(RunSafetyError, match="Another viralscan run"):
            acquire_run_lock(tmp_path)
        assert (tmp_path / "run_manifest.json").is_file()
    finally:
        holder.stdin.close()
        holder.wait()
    assert not (tmp_path / LOCK_NAME).exists()  # atexit released it


def test_lock_file_is_not_content_and_survives_overwrite(tmp_path: Path) -> None:
    acquire_run_lock(tmp_path)
    assert _prepare(tmp_path) == "new"
    (tmp_path / "junk").write_text("x")
    _prepare(tmp_path, overwrite=True)
    assert (tmp_path / LOCK_NAME).is_file() and not (tmp_path / "junk").exists()


@pytest.mark.parametrize("resume", [False, True])
def test_tmp_litter_is_not_content(tmp_path: Path, resume: bool) -> None:
    (tmp_path / ".run_manifest.json.tmp").write_text("{")
    assert _prepare(tmp_path, resume=resume) == "new"


def test_bad_output_paths_are_run_safety_errors(tmp_path: Path) -> None:
    f = tmp_path / "afile"
    f.write_text("x")
    with pytest.raises(RunSafetyError, match="is a file"):
        _prepare(f)
    with pytest.raises(RunSafetyError, match="Cannot"):
        acquire_run_lock(f / "sub")
    dangling = tmp_path / "dangling"
    dangling.symlink_to(tmp_path / "nowhere" / "deeper")
    try:
        acquire_run_lock(dangling)
        _prepare(dangling)
    except RunSafetyError:
        pass  # a clear error is fine; a FileExistsError traceback is not
    if os.geteuid() != 0:
        ro = tmp_path / "ro"
        ro.mkdir(mode=0o500)
        with pytest.raises(RunSafetyError, match="Cannot"):
            acquire_run_lock(ro / "out")


def test_inventory_ignores_tmp_staging() -> None:
    from viralscan.artifact_inventory import is_ignored

    assert is_ignored("results/artifact_inventory.json.tmp")


def test_sigterm_is_forwarded_to_the_snakemake_group(tmp_path: Path) -> None:
    pidfile = tmp_path / "child.pid"
    child = f"import os,time;open({str(pidfile)!r},'w').write(str(os.getpid()));time.sleep(60)"
    code = (
        "import sys;from viralscan import menu;"
        f"menu._run_snakemake([sys.executable,'-c',{child!r}],what='x',hint='h')"
    )
    parent = subprocess.Popen([sys.executable, "-c", code])
    for _ in range(100):
        if pidfile.exists() and pidfile.read_text():
            break
        time.sleep(0.1)
    pid = int(pidfile.read_text())
    parent.send_signal(signal.SIGTERM)
    assert parent.wait(timeout=20) == 143
    time.sleep(0.2)
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)


def test_summary_runtime_is_cumulative_so_a_noop_resume_keeps_it(tmp_path: Path) -> None:
    (tmp_path / "summary.txt").write_text(
        "Headline\n\nRuntime: 100.0000 seconds.\n\nAmount of genes in data: 2\n"
    )
    menu._write_sample_summary(str(tmp_path), 0.2, 5, 2)
    assert "Runtime: 100.2000 seconds." in (tmp_path / "summary.txt").read_text()


def test_resume_attempt_keeps_first_blocks(tmp_path: Path) -> None:
    (tmp_path / "run_manifest.json").write_text(
        json.dumps({"tool_binaries": {"run": {"kb": "first"}}})
    )
    record_resume_attempt(
        tmp_path,
        {"tool_binaries": {"run": {"kb": "second"}}, "chemistry_detection": {"S": {"c": 1}}},
    )
    m = json.loads((tmp_path / "run_manifest.json").read_text())
    assert m["tool_binaries"]["run"] == {"kb": "first"}
    assert m["chemistry_detection"]["S"] == {"c": 1}  # missing before: filled
    assert m["resume_attempts"][0]["tool_binaries"]["run"] == {"kb": "second"}


# -- end to end through _main with snakemake patched out ---------------------------------


@pytest.fixture
def run_main(tmp_path: Path, monkeypatch):
    from unittest.mock import patch

    for name in ("a_R1.fastq", "a_R2.fastq", "index.idx", "t2g.txt"):
        (tmp_path / name).write_text("@r\nACGT\n+\nIIII\n" if "fastq" in name else "t\tg\n")
    out = tmp_path / "out"
    calls: list[list[str]] = []

    def go(*extra: str) -> list[list[str]]:
        calls.clear()
        argv = [
            "viralscan", "-i", str(tmp_path / "index.idx"), "-t", str(tmp_path / "t2g.txt"),
            "--r1", str(tmp_path / "a_R1.fastq"), "--r2", str(tmp_path / "a_R2.fastq"),
            "-o", str(out), "-x", "10xv3", "--yes", *extra,
        ]  # fmt: skip
        monkeypatch.setattr(sys, "argv", argv)
        with (
            patch.object(menu, "_check_required_tools"),
            patch.object(menu, "_check_cell_caller_tools"),
            patch.object(menu, "_resolve_chemistry", return_value={}),
            patch.object(menu, "_run_snakemake", side_effect=lambda cmd, **k: calls.append(cmd)),
        ):
            menu._main()
        return list(calls)

    go.out = out
    return go


def test_resume_restamps_readmes_before_any_work_that_can_fail(run_main, monkeypatch) -> None:
    run_main()
    out = run_main.out
    assert (out / "run_complete.json").is_file()
    assert "finished" in (out / "README.txt").read_text()

    def boom(*_args):
        raise RuntimeError("boom")

    monkeypatch.setattr(menu, "_count_lines", boom)  # runs after the marker is cleared
    with pytest.raises(RuntimeError):
        run_main("--resume")
    assert not (out / "run_complete.json").exists()
    assert "NOT FINISHED" in (out / "README.txt").read_text()
    assert "NOT FINISHED" in (out / "a" / "README.txt").read_text()


def test_resume_unlocks_snakemake_first_and_logs_the_attempt(run_main) -> None:
    first = run_main()
    assert "--unlock" not in first[0]  # a fresh run only unlocks after success
    second = run_main("--resume")
    assert "--unlock" in second[0] and "all" in second[1]
    manifest = json.loads((run_main.out / "run_manifest.json").read_text())
    assert len(manifest["resume_attempts"]) == 1


def test_run_in_group_captures_and_raises_like_subprocess_run() -> None:
    from viralscan.run_safety import run_in_group

    code = "import sys;print('out');sys.exit(3)"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_in_group([sys.executable, "-c", code], capture_output=True)
    assert err.value.returncode == 3
    assert err.value.stdout == b"out\n"
    assert run_in_group([sys.executable, "-c", code], check=False).returncode == 3
