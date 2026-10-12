"""Live kallisto check of the strand/technology assumptions the evidence replay relies on (M-01).

``replay_exact_target_bus`` re-runs ``kallisto bus`` and must match ``kb count``. Two facts decide
that: no strand flag means *unstranded* for ``-x 10xv3`` (kallisto 0.52.0), so an unset
``config.yaml`` strand is replayed with no flag; and kallisto rejects ``-x 10xv4``, which kb count
rewrites to ``10XV3``.
"""

from __future__ import annotations

import json
import random
import subprocess
from pathlib import Path

import pytest

from viralscan.evidence import have_tools, kallisto_technology
from viralscan.validation import tool_path

_N = 200


def _rc(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _library(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Two transcripts; half the 10x-style pairs carry the sense cDNA, half its reverse complement."""
    rng = random.Random(1)

    def dna(n: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(n))

    transcripts = [dna(1500), dna(1500)]
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">tx{i}\n{s}\n" for i, s in enumerate(transcripts)))
    r1, r2 = tmp_path / "R1.fq", tmp_path / "R2.fq"
    with r1.open("w") as f1, r2.open("w") as f2:
        for i in range(_N):
            start = rng.randint(0, 1300)
            cdna = transcripts[i % 2][start : start + 90]
            f1.write(f"@r{i}\n{dna(28)}\n+\n{'I' * 28}\n")
            f2.write(f"@r{i}\n{cdna if i < _N // 2 else _rc(cdna)}\n+\n{'I' * 90}\n")
    return fasta, r1, r2


def _bus(
    kallisto: str, index: Path, out: Path, x: str, flag: str | None, r1: Path, r2: Path
) -> subprocess.CompletedProcess[bytes]:
    cmd = [kallisto, "bus", "-i", str(index), "-o", str(out), "-x", x]
    return subprocess.run([*cmd, *([flag] if flag else []), str(r1), str(r2)], capture_output=True)


def _aligned(out: Path) -> int:
    return int(json.loads((out / "run_info.json").read_text())["n_pseudoaligned"])


@pytest.mark.integration
def test_no_flag_is_unstranded_and_10xv4_needs_the_kb_rewrite(tmp_path: Path) -> None:
    if have_tools(["kallisto"]):
        pytest.skip("kallisto not on PATH")
    kallisto = tool_path("kallisto")
    assert kallisto is not None
    fasta, r1, r2 = _library(tmp_path)
    index = tmp_path / "idx"
    subprocess.run(
        [kallisto, "index", "-i", str(index), str(fasta)], check=True, capture_output=True
    )

    counts = {}
    for label, flag in (
        ("none", None),
        ("unstranded", "--unstranded"),
        ("forward", "--fr-stranded"),
    ):
        done = _bus(kallisto, index, tmp_path / label, "10xv3", flag, r1, r2)
        assert done.returncode == 0, done.stderr
        counts[label] = _aligned(tmp_path / label)
    assert counts["unstranded"] == _N
    assert counts["none"] == counts["unstranded"]  # the replay's "no flag" == kb count's "no flag"
    assert counts["forward"] == _N // 2  # a stranded run drops the antisense half

    assert _bus(kallisto, index, tmp_path / "raw", "10xv4", None, r1, r2).returncode != 0
    rewritten = kallisto_technology("10xv4")
    assert _bus(kallisto, index, tmp_path / "v4", rewritten, None, r1, r2).returncode == 0
