"""P3-02: the competitive alignment must emit secondaries so exact host/virus ties are ties."""

import random
from collections import Counter

import pytest

from viralscan import evidence as ev
from viralscan.molecule_verdict import molecule_rows

pytestmark = pytest.mark.integration


def test_200_perfect_host_virus_ties_are_all_tie_verdicts(tmp_path):
    missing = ev.have_tools(["minimap2", "samtools"])
    if missing:
        pytest.skip(f"missing tools: {missing}")
    rng = random.Random(7)
    rs = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    reads = [rs(100) for _ in range(200)]
    ref = tmp_path / "ref.fa"
    ref.write_text(
        f">VIRUS|v\n{''.join(rs(50) + r + rs(50) for r in reads)}\n"
        f">HOST|chr\n{''.join(rs(70) + r + rs(70) for r in reads)}\n"
    )
    fasta = tmp_path / "r.fa"
    fasta.write_text("".join(f">CB1_U{i}_{i + 1}\n{r}\n" for i, r in enumerate(reads)))
    bam = ev.align_reads_to_viral(str(fasta), str(ref), str(tmp_path / "raw.bam"), threads=1)
    sam = ev._run(["samtools", "view", bam], capture=True).decode().splitlines()
    verdicts = Counter(r["verdict"] for r in molecule_rows(sam, counted_reads=set(range(1, 201))))
    assert verdicts == {"tie": 200}
