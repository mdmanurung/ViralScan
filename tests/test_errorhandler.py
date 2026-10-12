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
        "host_filter": None,
        "host_index": None,
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

        with patch("os.path.exists", side_effect=exists), pytest.raises(SystemExit) as exc:
            errorhandler(args)
        assert exc.value.code != 0

    def test_reference_with_nonexistent_gtf_dies(self) -> None:
        args = _args(reference=True, gtf="/no/such.gtf", fasta="/fake/virus.fasta")

        def exists(p: str) -> bool:
            return p != "/no/such.gtf"

        with patch("os.path.exists", side_effect=exists), pytest.raises(SystemExit) as exc:
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

        with patch("os.path.exists", side_effect=exists), pytest.raises(SystemExit) as exc:
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

        with patch("os.path.exists", side_effect=exists), pytest.raises(SystemExit) as exc:
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

        with patch("os.path.exists", side_effect=exists), pytest.raises(SystemExit) as exc:
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


class TestHostFilterValidation:
    def test_host_filter_requires_host_index(self) -> None:
        args = _args(host_filter="starsolo", host_index=None)
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_host_index_requires_host_filter(self) -> None:
        args = _args(host_filter=None, host_index="/fake/star-index")
        with patch("os.path.exists", side_effect=_always_exists):
            with pytest.raises(SystemExit) as exc:
                errorhandler(args)
        assert exc.value.code != 0

    def test_missing_host_index_dies(self) -> None:
        args = _args(host_filter="starsolo", host_index="/no/such/star-index")

        def exists(p: str) -> bool:
            return p != "/no/such/star-index"

        with patch("os.path.exists", side_effect=exists), pytest.raises(SystemExit) as exc:
            errorhandler(args)
        assert exc.value.code != 0

    def test_valid_host_filter_pair_passes(self) -> None:
        args = _args(host_filter="starsolo", host_index="/fake/star-index")
        with patch("os.path.exists", side_effect=_always_exists):
            errorhandler(args)


# ---------------------------------------------------------------------------
# AUD C-07 / C-09: interactive prompt without a TTY, pre-run path checks
# ---------------------------------------------------------------------------


def _missing(*bad: str):
    return lambda p: p not in bad


class TestReferenceWithIndexPrompt:
    """C-07: `--reference -i ...` prompted via input() and crashed on EOF."""

    def _args(self, **kw) -> argparse.Namespace:
        return _args(reference=True, gtf="/fake/v.gtf", fasta="/fake/v.fa", **kw)

    def test_eof_dies_with_a_message_not_a_traceback(self, caplog) -> None:
        with (
            patch("os.path.exists", side_effect=_always_exists),
            patch("builtins.input", side_effect=EOFError),
            pytest.raises(SystemExit) as exc,
        ):
            errorhandler(self._args())
        assert exc.value.code != 0
        assert "--yes" in caplog.text

    def test_yes_skips_the_prompt(self) -> None:
        with (
            patch("os.path.exists", side_effect=_always_exists),
            patch("builtins.input", side_effect=AssertionError("must not prompt")),
        ):
            errorhandler(self._args(yes=True))


class TestPreRunPathChecks:
    """C-09: inputs that only failed after kb count are now checked up front."""

    @pytest.mark.parametrize(
        ("attr", "path"),
        [
            ("whitelist", "/no/wl.txt"),
            ("host_h5ad", "/no/host.h5ad"),
            ("cell_types", "/no/types.csv"),
            ("called_cells_file", "/no/called.tsv"),
            ("positive_control_receipt", "/no/receipt.json"),
        ],
    )
    def test_missing_file_dies(self, attr: str, path: str, caplog) -> None:
        args = _args(**{attr: path})
        with patch("os.path.exists", side_effect=_missing(path)), pytest.raises(SystemExit) as exc:
            errorhandler(args)
        assert exc.value.code != 0
        assert path in caplog.text

    def test_existing_files_pass(self) -> None:
        args = _args(whitelist="/f/wl.txt", host_h5ad="/f/h.h5ad", cell_types="/f/c.csv")
        with patch("os.path.exists", side_effect=_always_exists):
            errorhandler(args)

    def test_external_cell_calling_requires_the_file(self, caplog) -> None:
        with patch("os.path.exists", side_effect=_always_exists), pytest.raises(SystemExit):
            errorhandler(_args(cell_calling="external"))
        assert "--called-cells-file" in caplog.text

    def test_called_cells_file_ignored_with_none_warns(self, caplog) -> None:
        with patch("os.path.exists", side_effect=_always_exists):
            errorhandler(_args(cell_calling="none", called_cells_file="/f/c.tsv"))
        assert "ignored" in caplog.text


class TestNumericOptionRanges:
    """C-10: a bad number is a clean error before any output directory is created."""

    @pytest.mark.parametrize(
        ("attr", "bad"),
        [("cores", 0), ("se_threshold", -3), ("umap_n_neighbors", 0), ("min_counts", -1)],
    )
    def test_bad_number_dies_naming_the_option(self, attr: str, bad: int, caplog) -> None:
        with patch("os.path.exists", side_effect=_always_exists), pytest.raises(SystemExit):
            errorhandler(_args(**{attr: bad}))
        assert attr in caplog.text
