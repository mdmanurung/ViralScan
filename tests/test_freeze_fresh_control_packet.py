from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path

import pytest

from scripts.freeze_fresh_control_packet import (
    FreshPacketFreezeError,
    freeze_packet,
)


def test_freeze_packet_copies_runner_and_hashes_every_executable_input(
    tmp_path: Path,
) -> None:
    packet = tmp_path / "packet"
    packet.mkdir()
    required = (
        "run_task.sh",
        "tasks.tsv",
        "tasks.small.tsv",
        "tasks.large.tsv",
        "control_inputs.raw.tsv",
    )
    for name in required:
        (packet / name).write_text(f"{name}\n", encoding="utf-8")
    runner = tmp_path / "run_fresh_control.py"
    runner.write_text("print('frozen')\n", encoding="utf-8")

    manifest = freeze_packet(packet_root=packet, runner_source=runner)

    frozen_runner = packet / "source/run_fresh_control.py"
    assert frozen_runner.read_bytes() == runner.read_bytes()
    checksum_rows = {
        relative: digest
        for digest, relative in (
            line.split("  ", maxsplit=1)
            for line in manifest.read_text(encoding="utf-8").splitlines()
        )
    }
    assert set(checksum_rows) == {
        *required,
        "source/run_fresh_control.py",
    }
    assert (
        checksum_rows["source/run_fresh_control.py"]
        == hashlib.sha256(runner.read_bytes()).hexdigest()
    )


def test_freeze_packet_refuses_to_replace_an_existing_snapshot(
    tmp_path: Path,
) -> None:
    packet = tmp_path / "packet"
    packet.mkdir()
    for name in (
        "run_task.sh",
        "tasks.tsv",
        "tasks.small.tsv",
        "tasks.large.tsv",
        "control_inputs.raw.tsv",
    ):
        (packet / name).write_text(f"{name}\n", encoding="utf-8")
    runner = tmp_path / "run_fresh_control.py"
    runner.write_text("print('frozen')\n", encoding="utf-8")
    (packet / "packet.sha256").write_text("prior snapshot\n", encoding="utf-8")

    with pytest.raises(FreshPacketFreezeError, match="already exists"):
        freeze_packet(packet_root=packet, runner_source=runner)

    assert (packet / "packet.sha256").read_text(encoding="utf-8") == "prior snapshot\n"


def test_slurm_runner_uses_explicit_packet_root_after_script_copy(
    tmp_path: Path,
) -> None:
    repository_root = Path(__file__).resolve().parents[1]
    runner = repository_root / "scripts/run_fresh_control_task.sh"
    packet = tmp_path / "packet"
    packet.mkdir()
    for name in (
        "tasks.tsv",
        "tasks.small.tsv",
        "tasks.large.tsv",
        "control_inputs.raw.tsv",
    ):
        (packet / name).write_text("header\n", encoding="utf-8")
    (packet / "run_task.sh").write_bytes(runner.read_bytes())
    (packet / "source").mkdir()
    (packet / "source/run_fresh_control.py").write_text(
        "raise SystemExit('must not execute')\n",
        encoding="utf-8",
    )
    checksum_paths = (
        "control_inputs.raw.tsv",
        "run_task.sh",
        "source/run_fresh_control.py",
        "tasks.large.tsv",
        "tasks.small.tsv",
        "tasks.tsv",
    )
    (packet / "packet.sha256").write_text(
        "".join(
            f"{hashlib.sha256((packet / relative).read_bytes()).hexdigest()}  {relative}\n"
            for relative in checksum_paths
        ),
        encoding="utf-8",
    )
    copied_slurm_script = tmp_path / "slurm-spool/slurm_script"
    copied_slurm_script.parent.mkdir()
    copied_slurm_script.write_bytes(runner.read_bytes())
    copied_slurm_script.chmod(0o755)
    env = os.environ.copy()
    env.update(
        {
            "FRESH_PACKET_ROOT": str(packet),
            "FRESH_TASK_MANIFEST": str(packet / "tasks.small.tsv"),
            "SLURM_ARRAY_TASK_ID": "999",
        }
    )

    completed = subprocess.run(
        [str(copied_slurm_script)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert completed.returncode == 1
    assert "packet.sha256: No such file" not in completed.stderr
    assert "run_task.sh: OK" in completed.stdout
