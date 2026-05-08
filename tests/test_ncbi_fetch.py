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


@pytest.mark.network
class TestFetchReferenceNetworkIntegration:
    """Live integration tests — require internet access.

    NC_002021.3 is Influenza A segment 8 (1027 nt): small, stable RefSeq
    entry unlikely to change or be removed.
    """

    def test_fetch_influenza_a_seg8(self, tmp_path) -> None:
        fasta_path, gtf_path = fetch_reference(
            accessions=["NC_002021.3"],
            out_dir=tmp_path / "ncbi",
            email="viralscan-test@example.org",
        )
        assert fasta_path.exists(), "FASTA file was not created"
        assert fasta_path.stat().st_size > 0, "FASTA file is empty"
        assert gtf_path.exists(), "GTF file was not created"
        assert gtf_path.stat().st_size > 0, "GTF file is empty"
        # Sanity-check FASTA format
        assert fasta_path.read_text().startswith(">"), "FASTA does not start with '>'"
        # Sanity-check GTF has at least one exon record
        assert "exon" in gtf_path.read_text(), "GTF contains no exon records"


# ---------------------------------------------------------------------------
# Audit §3.2 — cache content validation
# ---------------------------------------------------------------------------


class TestCacheValidation:
    """Audit §3.2: truncated/corrupt cached files must be detected and re-downloaded.

    The original _fetch_one() only checks path.exists() and st_size == 0.
    A non-empty but truncated file from a prior interrupted download is
    silently reused, feeding corrupt data to kb ref.

    The fix writes a .sha256 sidecar alongside each cached file and re-downloads
    if the sidecar is missing or the checksum does not match.

    Regression for: audits/2026-05-08-full-pipeline.md §3.2
    """

    VALID_FASTA = ">NC_FAKE1\nATGCATGC\n"
    VALID_GTF = 'NC_FAKE1\tNCBI\texon\t1\t8\t.\t+\t0\tgene_id "X"; transcript_id "X";\n'

    def _write_with_sidecar(self, path, content: str) -> None:
        """Write content and store its SHA-256 in a .sha256 sidecar (as the fix does)."""
        import hashlib
        path.write_text(content)
        sha = hashlib.sha256(content.encode()).hexdigest()
        path.with_suffix(path.suffix + ".sha256").write_text(sha)

    def _make_cache_dir(self, tmp_path, acc: str):
        acc_dir = tmp_path / acc
        acc_dir.mkdir(parents=True, exist_ok=True)
        return acc_dir

    def _patch_fetch(self, monkeypatch, fasta_content: str, efetch_calls: list) -> None:
        """Monkeypatch _efetch and _validate_accession for offline testing."""
        from viralscan.scripts import ncbi_fetch

        def mock_efetch(acc, rettype, email, api_key):
            efetch_calls.append(rettype)
            if rettype == "fasta":
                return fasta_content
            if rettype == "gb":
                return "LOCUS NC_FAKE1\nFEATURES\n     CDS             1..8\n                     /gene=\"X\"\n//\n"
            raise AssertionError(f"Unexpected rettype: {rettype}")

        monkeypatch.setattr(ncbi_fetch, "_efetch", mock_efetch)
        # Bypass accession regex validation for synthetic IDs
        monkeypatch.setattr(ncbi_fetch, "_validate_accession", lambda x: x)

    def test_no_sidecar_triggers_redownload(self, tmp_path, monkeypatch) -> None:
        """
        GIVEN: cached FASTA exists (non-empty) but has NO .sha256 sidecar
        WHEN:  _fetch_one is called
        THEN:  _efetch is called again and the file is refreshed

        Regression for: audits/2026-05-08-full-pipeline.md §3.2
        """
        from viralscan.scripts.ncbi_fetch import _fetch_one

        acc = "NC_FAKE1"
        acc_dir = self._make_cache_dir(tmp_path, acc)
        fasta_path = acc_dir / f"{acc}.fasta"
        gtf_path = acc_dir / f"{acc}.gtf"

        # Corrupt FASTA — non-empty but no sidecar (simulates partial download)
        fasta_path.write_text("partial content — no sidecar")
        # Valid GTF pre-seeded so only FASTA triggers re-fetch
        self._write_with_sidecar(gtf_path, self.VALID_GTF)

        efetch_calls: list = []
        self._patch_fetch(monkeypatch, self.VALID_FASTA, efetch_calls)

        fasta_out, _ = _fetch_one(acc, tmp_path, "test@test.com", None)

        assert "fasta" in efetch_calls, (
            "Expected _efetch(rettype='fasta') to be called when sidecar is missing, "
            f"but efetch calls were: {efetch_calls}"
        )
        assert fasta_out.read_text() == self.VALID_FASTA, (
            "Re-downloaded FASTA content does not match expected valid content."
        )

    def test_mismatched_sidecar_triggers_redownload(self, tmp_path, monkeypatch) -> None:
        """
        GIVEN: cached FASTA exists with a .sha256 sidecar that does NOT match
               (simulates file corruption after download)
        WHEN:  _fetch_one is called
        THEN:  _efetch is called again and the file is refreshed
        """
        from viralscan.scripts.ncbi_fetch import _fetch_one

        acc = "NC_FAKE1"
        acc_dir = self._make_cache_dir(tmp_path, acc)
        fasta_path = acc_dir / f"{acc}.fasta"
        gtf_path = acc_dir / f"{acc}.gtf"

        # FASTA with WRONG sidecar (hash of different content)
        fasta_path.write_text("corrupted content")
        fasta_path.with_suffix(".fasta.sha256").write_text("0" * 64)  # wrong hash
        # Valid GTF pre-seeded
        self._write_with_sidecar(gtf_path, self.VALID_GTF)

        efetch_calls: list = []
        self._patch_fetch(monkeypatch, self.VALID_FASTA, efetch_calls)

        fasta_out, _ = _fetch_one(acc, tmp_path, "test@test.com", None)

        assert "fasta" in efetch_calls, (
            "Expected re-download when sidecar checksum is wrong, "
            f"but efetch calls were: {efetch_calls}"
        )
        assert fasta_out.read_text() == self.VALID_FASTA

    def test_valid_sidecar_skips_redownload(self, tmp_path, monkeypatch) -> None:
        """
        GIVEN: cached FASTA exists with a matching .sha256 sidecar
        WHEN:  _fetch_one is called
        THEN:  _efetch is NOT called (cache hit)
        """
        from viralscan.scripts.ncbi_fetch import _fetch_one

        acc = "NC_FAKE1"
        acc_dir = self._make_cache_dir(tmp_path, acc)
        fasta_path = acc_dir / f"{acc}.fasta"
        gtf_path = acc_dir / f"{acc}.gtf"

        # Both files fully valid with correct sidecars
        self._write_with_sidecar(fasta_path, self.VALID_FASTA)
        self._write_with_sidecar(gtf_path, self.VALID_GTF)

        efetch_calls: list = []
        self._patch_fetch(monkeypatch, self.VALID_FASTA, efetch_calls)

        fasta_out, gtf_out = _fetch_one(acc, tmp_path, "test@test.com", None)

        assert efetch_calls == [], (
            f"Expected no _efetch calls on valid cache hit, but got: {efetch_calls}"
        )
        assert fasta_out.read_text() == self.VALID_FASTA
        assert gtf_out.read_text() == self.VALID_GTF
