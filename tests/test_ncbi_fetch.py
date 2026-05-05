"""Unit tests for viralscan.scripts.ncbi_fetch.

These tests do not hit the network; live integration tests should be marked
with ``@pytest.mark.network`` (see pyproject.toml).
"""

from __future__ import annotations

import textwrap

import pytest

from viralscan.scripts.ncbi_fetch import (
    NCBIFetchError,
    _genbank_to_gtf,
    _parse_location,
    _validate_accession,
    fetch_reference,
)


class TestValidateAccession:
    @pytest.mark.parametrize("acc", ["NC_002021.3", "NC_001512.1", "KX020937.1", "U00096"])
    def test_accepts_valid(self, acc: str) -> None:
        assert _validate_accession(acc) == acc

    @pytest.mark.parametrize(
        "acc", ["", "not-an-accession", "../etc/passwd", "NC_002021.3; rm -rf /"]
    )
    def test_rejects_invalid(self, acc: str) -> None:
        with pytest.raises(NCBIFetchError):
            _validate_accession(acc)


class TestParseLocation:
    def test_simple_range(self) -> None:
        assert _parse_location("1..1024") == [(1, 1024, "+")]

    def test_complement(self) -> None:
        assert _parse_location("complement(1..1024)") == [(1, 1024, "-")]

    def test_join(self) -> None:
        assert _parse_location("join(1..100,200..300)") == [(1, 100, "+"), (200, 300, "+")]

    def test_complement_join(self) -> None:
        assert _parse_location("complement(join(1..100,200..300))") == [
            (1, 100, "-"),
            (200, 300, "-"),
        ]

    def test_partial_markers_stripped(self) -> None:
        assert _parse_location("<1..>1024") == [(1, 1024, "+")]


class TestGenbankToGtf:
    def _record(self, features: str) -> str:
        return (
            textwrap.dedent(
                """\
            LOCUS       NC_TEST                 1024 bp    DNA     linear   VRL
            VERSION     NC_TEST.1
            FEATURES             Location/Qualifiers
            """
            )
            + features
            + "ORIGIN\n//\n"
        )

    def test_extracts_simple_cds(self) -> None:
        features = (
            "     CDS             1..900\n"
            '                     /gene="GAG"\n'
            '                     /product="capsid"\n'
            '                     /protein_id="ABC12345.1"\n'
        )
        gtf = _genbank_to_gtf(self._record(features), "NC_TEST.1")
        assert "NC_TEST.1\tNCBI\texon\t1\t900\t.\t+\t0" in gtf
        assert 'gene_id "GAG"' in gtf
        assert 'transcript_id "ABC12345.1"' in gtf

    def test_complement_strand(self) -> None:
        features = '     CDS             complement(1..500)\n                     /gene="POL"\n'
        gtf = _genbank_to_gtf(self._record(features), "NC_TEST.1")
        assert "\t-\t" in gtf

    def test_no_cds_raises(self) -> None:
        with pytest.raises(NCBIFetchError):
            _genbank_to_gtf(self._record(""), "NC_TEST.1")


class TestFetchReferenceArgValidation:
    def test_no_accessions_raises(self, tmp_path) -> None:
        with pytest.raises(NCBIFetchError):
            fetch_reference([], out_dir=tmp_path, email="me@example.org")

    def test_missing_email_raises(self, tmp_path, monkeypatch) -> None:
        monkeypatch.delenv("NCBI_EMAIL", raising=False)
        with pytest.raises(NCBIFetchError):
            fetch_reference(["NC_002021.3"], out_dir=tmp_path, email=None)
