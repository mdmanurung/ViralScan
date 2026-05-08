"""Direct unit tests for menu.errorhandler().

These tests call ``errorhandler(args)`` with an ``argparse.Namespace``
constructed directly (no sys.argv patching) so each logical branch is
exercised in isolation.  ``os.path.exists`` is mocked globally for each test
so the test suite does not require real files on disk.
"""

from __future__ import annotations

import argparse
from unittest.mock import patch

import pytest

from viralscan.menu import errorhandler


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _args(**kwargs) -> argparse.Namespace:
    """
    Build a Namespace with all required fields set to safe defaults.
    All fake paths end in a valid FASTQ suffix so FASTQ checks pass by default.
    """
    defaults = {
        # Mode selectors
        "ncbi_accession": None,
        "reference": False,
        "gtf": None,
        "fasta": None,
        # Index-mode paths
        "index": "/fake/index.idx",
        "transcripts": "/fake/t2g.txt",
        "f1": None,
        # Samples
        "sample1": "/fake/R1.fastq.gz",
        "sample2": "/fake/R2.fastq.gz",
    }
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def _always_exists(path: str) -> bool:  # noqa: ARG001
    return True


# ---------------------------------------------------------------------------
# NCBI accession mutual exclusion
# ---------------------------------------------------------------------------


class TestNcbiAccessionMutualExclusion:
    def test_ncbi_plus_reference_flag_dies(self) -> None:
        args = _args(ncbi_accession="NC_002021.3", reference=True)
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_ncbi_plus_gtf_dies(self) -> None:
        args = _args(ncbi_accession="NC_002021.3", gtf="/some/file.gtf")
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_ncbi_plus_fasta_dies(self) -> None:
        args = _args(ncbi_accession="NC_002021.3", fasta="/some/file.fasta")
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_ncbi_alone_passes(self) -> None:
        """When only ncbi_accession is set (no reference/gtf/fasta), skip sample check too."""
        # Note: errorhandler still checks samples even in ncbi_accession mode.
        args = _args(ncbi_accession="NC_002021.3")
        with patch("os.path.exists", side_effect=_always_exists):
            errorhandler(args)  # should not raise


# ---------------------------------------------------------------------------
# Reference mode (--reference flag)
# ---------------------------------------------------------------------------


class TestReferenceModeValidation:
    def test_reference_without_gtf_dies(self) -> None:
        args = _args(reference=True, gtf=None, fasta="/fake/virus.fasta")
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_reference_without_fasta_dies(self) -> None:
        args = _args(reference=True, gtf="/fake/virus.gtf", fasta=None)
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_reference_with_nonexistent_fasta_dies(self) -> None:
        args = _args(reference=True, gtf="/fake/virus.gtf", fasta="/no/such.fasta")

        def exists(p: str) -> bool:
            return p != "/no/such.fasta"

        with patch("os.path.exists", side_effect=exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_reference_with_nonexistent_gtf_dies(self) -> None:
        args = _args(reference=True, gtf="/no/such.gtf", fasta="/fake/virus.fasta")

        def exists(p: str) -> bool:
            return p != "/no/such.gtf"

        with patch("os.path.exists", side_effect=exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_reference_valid_paths_passes(self) -> None:
        # index and transcripts must be None to avoid triggering the
        # "you provided an index but also --reference" interactive prompt.
        args = _args(
            reference=True,
            gtf="/fake/virus.gtf",
            fasta="/fake/virus.fasta",
            index=None,
            transcripts=None,
            f1=None,
        )
        with patch("os.path.exists", side_effect=_always_exists):
            errorhandler(args)  # should not raise


# ---------------------------------------------------------------------------
# Index mode (default — no --reference, no --ncbi-accession)
# ---------------------------------------------------------------------------


class TestIndexModeValidation:
    def test_missing_index_dies(self) -> None:
        args = _args(index=None)
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_nonexistent_index_dies(self) -> None:
        args = _args(index="/no/such.idx")

        def exists(p: str) -> bool:
            return p != "/no/such.idx"

        with patch("os.path.exists", side_effect=exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_missing_transcripts_dies(self) -> None:
        args = _args(transcripts=None)
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_nonexistent_transcripts_dies(self) -> None:
        args = _args(transcripts="/no/such_t2g.txt")

        def exists(p: str) -> bool:
            return p != "/no/such_t2g.txt"

        with patch("os.path.exists", side_effect=exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_valid_index_mode_passes(self) -> None:
        args = _args()  # all defaults point to existing fake paths
        with patch("os.path.exists", side_effect=_always_exists):
            errorhandler(args)  # should not raise


# ---------------------------------------------------------------------------
# Sample validation (applies to all modes)
# ---------------------------------------------------------------------------


class TestSampleValidation:
    def test_mismatched_sample_counts_dies(self) -> None:
        args = _args(
            sample1="/fake/A_R1.fastq.gz,/fake/B_R1.fastq.gz",
            sample2="/fake/A_R2.fastq.gz",
        )
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_nonexistent_sample_dies(self) -> None:
        args = _args(sample1="/no/such_R1.fastq.gz")

        def exists(p: str) -> bool:
            return p != "/no/such_R1.fastq.gz"

        with patch("os.path.exists", side_effect=exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    @pytest.mark.parametrize("bad_suffix", [".bam", ".txt", ".csv", ".h5ad"])
    def test_invalid_forward_suffix_dies(self, bad_suffix: str) -> None:
        args = _args(sample1=f"/fake/sample_R1{bad_suffix}")
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    @pytest.mark.parametrize("bad_suffix", [".bam", ".txt", ".csv", ".h5ad"])
    def test_invalid_reverse_suffix_dies(self, bad_suffix: str) -> None:
        args = _args(sample2=f"/fake/sample_R2{bad_suffix}")
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    @pytest.mark.parametrize(
        "s1,s2",
        [
            ("/fake/A_R1.fastq", "/fake/A_R2.fastq"),
            ("/fake/A_R1.fq", "/fake/A_R2.fq"),
            ("/fake/A_R1.fastq.gz", "/fake/A_R2.fastq.gz"),
            ("/fake/A_R1.fq.gz", "/fake/A_R2.fq.gz"),
        ],
    )
    def test_valid_suffixes_pass(self, s1: str, s2: str) -> None:
        args = _args(sample1=s1, sample2=s2)
        with patch("os.path.exists", side_effect=_always_exists):
            errorhandler(args)  # should not raise

    def test_multiple_samples_valid(self) -> None:
        args = _args(
            sample1="/fake/A_R1.fastq.gz,/fake/B_R1.fastq.gz",
            sample2="/fake/A_R2.fastq.gz,/fake/B_R2.fastq.gz",
        )
        with patch("os.path.exists", side_effect=_always_exists):
            errorhandler(args)  # should not raise
