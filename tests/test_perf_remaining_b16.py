"""Sequential native spool: identical argv/errors, bounded memory, no leaked file."""

import errno
import importlib.util
import subprocess
from pathlib import Path

import pytest

from viralscan import evidence as candidate

spec = importlib.util.spec_from_file_location(
    "original_align", Path(__file__).parent / "_frozen_oracles/frozen_align_d056986.py"
)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)


@pytest.mark.parametrize("failing", [None, "minimap2", "sort", "index"])
def test_native_order_stderr_and_partial_files_match(tmp_path, monkeypatch, failing):
    def run(fn):
        calls = []
        bam = tmp_path / "out.bam"
        bam.unlink(missing_ok=True)
        Path(str(bam) + ".bai").unlink(missing_ok=True)

        def native(cmd, **kwargs):
            calls.append(cmd)
            tool = "minimap2" if cmd[0] == "minimap2" else cmd[1]
            if tool == "sort":
                data = kwargs.get("input")
                if data is None:
                    data = kwargs["stdin"].read()
                assert data == b"complete SAM output\n"
                bam.write_bytes(b"partial sort" if failing == tool else b"complete BAM")
            if failing == tool:
                return subprocess.CompletedProcess(cmd, 17, b"", b"complete stderr\xff\n" * 2000)
            if tool == "minimap2" and hasattr(kwargs.get("stdout"), "write"):
                kwargs["stdout"].write(b"complete SAM output\n")
            if tool == "index":
                Path(str(bam) + ".bai").write_bytes(b"BAI")
            return subprocess.CompletedProcess(cmd, 0, b"complete SAM output\n", b"")

        monkeypatch.setattr(subprocess, "run", native)
        try:
            result, error = fn("reads.fa", "refs.fa", str(bam)), None
        except Exception as exc:
            result, error = None, (type(exc), str(exc))
        return result, error, calls, {p.name: p.read_bytes() for p in tmp_path.iterdir()}

    assert run(original.align_reads_to_viral) == run(candidate.align_reads_to_viral)


def test_spool_creation_failure_keeps_outputs_and_starts_no_child(tmp_path, monkeypatch):
    bam = tmp_path / "preexisting.bam"
    bam.write_bytes(b"old")

    def fail():
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(candidate.tempfile, "TemporaryFile", fail)
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: pytest.fail("started child"))
    with pytest.raises(OSError) as exc:
        candidate.align_reads_to_viral("reads", "refs", str(bam))
    assert exc.value.errno == errno.ENOSPC
    assert bam.read_bytes() == b"old"


def test_spool_write_failure_stops_before_sort_and_closes_handle(tmp_path, monkeypatch):
    real = candidate.tempfile.TemporaryFile
    handles = []

    def spool():
        handle = real()
        handles.append(handle)
        return handle

    monkeypatch.setattr(candidate.tempfile, "TemporaryFile", spool)
    calls = []

    def native(cmd, **kwargs):
        calls.append(cmd)
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(subprocess, "run", native)
    with pytest.raises(OSError):
        candidate.align_reads_to_viral("reads", "refs", str(tmp_path / "out.bam"))
    assert len(calls) == 1 and handles[0].closed


@pytest.mark.parametrize("failing", ["minimap2", "sort", "index"])
def test_real_child_failure_preserves_full_stderr_and_order(tmp_path, monkeypatch, failing):
    import sys

    native_run = subprocess.run
    script = """
import pathlib,sys
stage,bam,fail=sys.argv[1:]
with open(str(pathlib.Path(bam).parent/'calls.txt'),'a') as f:f.write(stage+'\\n')
if stage=='minimap2':sys.stdout.buffer.write(b'complete SAM output\\n')
if stage=='sort':
 assert sys.stdin.buffer.read()==b'complete SAM output\\n'
 pathlib.Path(bam).write_bytes(b'partial BAM' if stage==fail else b'BAM')
if stage==fail:
 sys.stderr.write('complete stderr:'+('E'*262144)+'\\n')
 sys.exit(17)
if stage=='index':pathlib.Path(bam+'.bai').write_bytes(b'BAI')
"""

    def run(fn):
        for p in tmp_path.iterdir():
            p.unlink()
        bam = str(tmp_path / "out.bam")

        def child(cmd, **kwargs):
            stage = "minimap2" if cmd[0] == "minimap2" else cmd[1]
            return native_run([sys.executable, "-c", script, stage, bam, failing], **kwargs)

        monkeypatch.setattr(subprocess, "run", child)
        with pytest.raises(RuntimeError) as caught:
            fn("reads", "refs", bam)
        assert "E" * 262144 in str(caught.value)
        return str(caught.value), {p.name: p.read_bytes() for p in tmp_path.iterdir()}

    assert run(original.align_reads_to_viral) == run(candidate.align_reads_to_viral)


def test_real_spool_disk_failure_stops_before_sort(tmp_path, monkeypatch):
    import sys

    if not Path("/dev/full").exists():
        pytest.skip("requires the kernel ENOSPC test device")
    native_run = subprocess.run
    handles, calls = [], []
    bam = tmp_path / "old.bam"
    bam.write_bytes(b"preexisting")

    def full_spool():
        handle = open("/dev/full", "w+b")  # noqa: SIM115 - caller owns the context
        handles.append(handle)
        return handle

    def child(cmd, **kwargs):
        calls.append(cmd)
        return native_run([sys.executable, "-c", "import os; os.write(1,b'SAM')"], **kwargs)

    monkeypatch.setattr(candidate.tempfile, "TemporaryFile", full_spool)
    monkeypatch.setattr(subprocess, "run", child)
    with pytest.raises(RuntimeError, match="No space left on device"):
        candidate.align_reads_to_viral("reads", "refs", str(bam))
    assert len(calls) == 1 and handles[0].closed
    assert bam.read_bytes() == b"preexisting"
