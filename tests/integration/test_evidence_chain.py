"""Integration test: live minimap2 → BAM → samtools coverage → BLAST chain.

Validates the ``viralscan.evidence`` subprocess wrappers against the real
binaries (minimap2, samtools, blastn, makeblastdb) using a fully synthetic,
deterministic viral genome and reads — **no public data, no network access**.

Skipped automatically when any required binary is absent from PATH, so the
default unit suite stays binary-free.  Run explicitly with::

    EVBIN=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/evtools/bin
    PATH="$EVBIN:$PATH" PYTHONPATH=src \\
        python -m pytest -m integration tests/integration/test_evidence_chain.py -v

Why synthetic fixtures
----------------------
The evidence module is designed for extracted viral reads (FASTA) aligned to
a viral reference FASTA.  These inputs are trivial to construct in-process:
a pseudo-random (fixed-seed) genome plus exact substring reads guarantee 100%
alignment identity, exercising every wrapper in one self-contained test.

PLAN S5 status: this test closes the "live wrappers not exercised" gap.
End-to-end ``viralscan evidence`` on a real run-dir is an operational step
(needs kb-python output), not a missing code path.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from viralscan.evidence import (
    align_reads_to_viral,
    blast_identity,
    coverage_table,
    have_tools,
)

# ── Fixture helpers ───────────────────────────────────────────────────────────

_TOOLS = ["minimap2", "samtools", "blastn", "makeblastdb"]

_GENOME_LEN = 2_200  # bp — enough for 40× read coverage at 150 bp reads
_READ_LEN = 150
_N_READS = 40
_SEED = 0  # fixed seed → reproducible, never changes


def _make_genome(length: int, seed: int = _SEED) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _make_reads(genome: str, n: int = _N_READS, read_len: int = _READ_LEN) -> list[str]:
    """Return *n* exact substrings of *genome*, evenly spaced."""
    max_start = len(genome) - read_len
    step = max(1, max_start // n)
    return [genome[i * step : i * step + read_len] for i in range(n)]


def _write_fasta(path: Path, name: str, seq: str) -> None:
    path.write_text(f">{name}\n{seq}\n")


def _write_reads_fasta(path: Path, seqs: list[str]) -> None:
    with open(path, "w") as fh:
        for i, seq in enumerate(seqs):
            fh.write(f">read_{i:04d}\n{seq}\n")


# ── The test ─────────────────────────────────────────────────────────────────


@pytest.mark.integration
def test_live_evidence_chain(tmp_path: Path) -> None:
    """Align synthetic reads → samtools coverage → BLAST; assert live wrappers work."""

    # Skip when any binary is absent (e.g. CI without evtools on PATH).
    missing = have_tools(_TOOLS)
    if missing:
        pytest.skip(f"Required binaries not on PATH: {', '.join(missing)}")

    # ── Build synthetic inputs ────────────────────────────────────────────────
    genome = _make_genome(_GENOME_LEN)
    reads = _make_reads(genome)

    viral_fasta = tmp_path / "virus_A.fasta"
    reads_fasta = tmp_path / "reads.fasta"
    _write_fasta(viral_fasta, "virus_A", genome)
    _write_reads_fasta(reads_fasta, reads)

    # ── 1. Align (minimap2 -ax sr → samtools sort → index) ───────────────────
    out_bam = str(tmp_path / "aligned.bam")
    returned = align_reads_to_viral(str(reads_fasta), str(viral_fasta), out_bam, threads=1)
    assert returned == out_bam, "align_reads_to_viral must return the out_bam path"

    bam_path = Path(out_bam)
    assert bam_path.exists(), "BAM file must exist after alignment"
    assert bam_path.stat().st_size > 0, "BAM file must not be empty"
    assert Path(out_bam + ".bai").exists(), "BAM index (.bai) must be created by samtools index"

    # ── 2. Coverage (samtools coverage → parsed rows) ─────────────────────────
    rows = coverage_table(out_bam)
    assert len(rows) > 0, "coverage_table must return at least one row for covered references"
    rnames = {r["rname"] for r in rows}
    assert "virus_A" in rnames, f"Expected 'virus_A' in coverage rows; got: {rnames}"

    virus_row = next(r for r in rows if r["rname"] == "virus_A")
    assert int(virus_row["numreads"]) > 0, "virus_A must have > 0 mapped reads"
    assert float(virus_row["coverage"]) > 0.0, "virus_A must have non-zero breadth coverage"

    # ── 3. BLAST identity (makeblastdb + blastn) ──────────────────────────────
    blast_work = tmp_path / "blast"
    hits = blast_identity(str(reads_fasta), str(viral_fasta), str(blast_work), threads=1)
    assert len(hits) > 0, "blast_identity must return at least one hit for synthetic reads"

    # Exact-substring reads must align to virus_A at high identity.
    pidents = [float(h["pident"]) for h in hits]
    subjects = [h["subject"] for h in hits]
    assert all(s == "virus_A" for s in subjects), (
        f"All hits should be against virus_A; got: {set(subjects)}"
    )
    assert min(pidents) >= 95.0, (
        f"Exact-substring reads must BLAST at ≥95% identity; min was {min(pidents):.1f}%"
    )
