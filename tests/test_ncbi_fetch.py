"""Unit tests for viralscan.scripts.ncbi_fetch.

These tests do not hit the network; live integration tests should be marked
with ``@pytest.mark.network`` (see pyproject.toml).
"""

from __future__ import annotations

import re
import textwrap

import pytest

from viralscan.scripts.ncbi_fetch import (
    NCBIFetchError,
    _genbank_to_gtf,
    _locus_fields,
    _parse_location,
    _source_qualifiers,
    _validate_accession,
    catalogue_rows,
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
        assert 'gene_id "NC_TEST.1_GAG"' in gtf
        assert 'gene_name "GAG"' in gtf
        assert 'transcript_id "ABC12345.1"' in gtf

    def test_complement_strand(self) -> None:
        features = '     CDS             complement(1..500)\n                     /gene="POL"\n'
        gtf = _genbank_to_gtf(self._record(features), "NC_TEST.1")
        assert "\t-\t" in gtf

    def test_no_cds_raises(self) -> None:
        with pytest.raises(NCBIFetchError):
            _genbank_to_gtf(self._record(""), "NC_TEST.1")

    def test_gene_id_is_genome_scoped_not_bare_ncbi_symbol(self) -> None:
        """Two genomes may both carry ``/gene="ORF1"``; their gene IDs must differ.

        Across the packaged Anelloviridae panel the bare symbol ``ORF1`` is the
        /product of 150 different genomes. Emitting it as gene_id collapses them
        onto one counting-matrix column.
        """
        features_a = '     CDS             1..900\n                     /gene="ORF1"\n'
        features_b = '     CDS             10..800\n                     /product="ORF1"\n'
        gtf_a = _genbank_to_gtf(self._record(features_a), "ACC_A.1")
        gtf_b = _genbank_to_gtf(self._record(features_b), "ACC_B.1")
        assert 'gene_id "ACC_A.1_ORF1"' in gtf_a
        assert 'gene_id "ACC_B.1_cds1"' in gtf_b
        assert 'gene_name "ORF1"' in gtf_b
        assert re.findall(r'gene_id "([^"]+)"', gtf_a) != re.findall(r'gene_id "([^"]+)"', gtf_b)

    def test_locus_tag_preferred_over_gene_and_protein_id(self) -> None:
        features = (
            "     CDS             1..900\n"
            '                     /locus_tag="TTVgp1"\n'
            '                     /gene="orf2/5"\n'
            '                     /protein_id="NP_817120.1"\n'
        )
        gtf = _genbank_to_gtf(self._record(features), "NC_TEST.1")
        assert 'gene_id "NC_TEST.1_TTVgp1"' in gtf
        assert 'gene_name "orf2/5"' in gtf
        assert 'locus_tag "TTVgp1"' in gtf
        assert 'gene "orf2/5"' in gtf

    def test_slash_in_gene_symbol_is_sanitised(self) -> None:
        features = '     CDS             1..900\n                     /gene="orf2/5"\n'
        gtf = _genbank_to_gtf(self._record(features), "NC_TEST.1")
        assert 'gene_id "NC_TEST.1_orf2_5"' in gtf

    def test_duplicate_symbols_within_one_genome_are_disambiguated(self) -> None:
        features = (
            "     CDS             1..900\n"
            '                     /gene="ORF1"\n'
            "     CDS             100..800\n"
            '                     /gene="ORF1"\n'
        )
        gtf = _genbank_to_gtf(self._record(features), "NC_TEST.1")
        ids = re.findall(r'gene_id "([^"]+)"', gtf)
        assert ids == ["NC_TEST.1_ORF1", "NC_TEST.1_ORF1_dup2"]

    def test_protein_id_fallback_when_no_symbol(self) -> None:
        features = '     CDS             1..900\n                     /protein_id="BAA86944.1"\n'
        gtf = _genbank_to_gtf(self._record(features), "NC_TEST.1")
        assert 'gene_id "NC_TEST.1_BAA86944.1"' in gtf

    def test_ordinal_fallback_when_no_identifier_at_all(self) -> None:
        features = "     CDS             1..900\n"
        gtf = _genbank_to_gtf(self._record(features), "NC_TEST.1")
        assert 'gene_id "NC_TEST.1_cds1"' in gtf

    def test_spliced_gene_emits_one_exon_per_interval(self) -> None:
        features = (
            '     CDS             join(353..711,2564..3077)\n                     /product="VP3"\n'
        )
        gtf = _genbank_to_gtf(self._record(features), "NC_002076.2")
        assert gtf.count("\texon\t") == 2
        assert 'n_exons "2"' in gtf
        assert "\texon\t353\t711\t.\t+\t0" in gtf
        assert "\texon\t2564\t3077\t.\t+\t0" in gtf

    def test_origin_spanning_join_is_flagged_and_kept_in_transcript_order(self) -> None:
        """A circular-genome feature crossing the origin must not be re-sorted.

        Anelloviridae are circular ssDNA, but a record annotates in a linear
        representation: the origin-spanning join is written high-coordinate-first
        and that order *is* the transcript order. Sorting ascending would emit a
        scrambled transcript.
        """
        record = self._record(
            '     CDS             join(4800..5000,1..200)\n                     /product="Rep"\n',
        ).replace("1024 bp", "5000 bp")
        gtf = _genbank_to_gtf(record, "NC_CIRC.1")
        rows = [line.split("\t") for line in gtf.strip().splitlines()]
        assert (rows[0][3], rows[0][4]) == ("4800", "5000")
        assert (rows[1][3], rows[1][4]) == ("1", "200")
        assert 'origin_spanning "true"' in gtf

    def test_internal_join_is_not_flagged_origin_spanning(self) -> None:
        gtf = _genbank_to_gtf(
            self._record("     CDS             join(100..200,300..400)\n"), "NC_TEST.1"
        )
        assert "origin_spanning" not in gtf

    def test_minus_strand_join_is_reversed_into_transcript_order(self) -> None:
        gtf = _genbank_to_gtf(
            self._record("     CDS             complement(join(1..100,300..400))\n"),
            "NC_TEST.1",
        )
        rows = [line.split("\t") for line in gtf.strip().splitlines()]
        assert [(row[3], row[4], row[6]) for row in rows] == [
            ("300", "400", "-"),
            ("1", "100", "-"),
        ]


class TestLocusFields:
    def test_length_topology_and_molecule(self) -> None:
        record = (
            "LOCUS       NC_TEST                 4848 bp    DNA     circular  VRL\n"
            "VERSION     NC_TEST.2\n"
            "FEATURES             Location/Qualifiers\n"
            "     source          1..4848\n"
            '                     /organism="Torque teno virus"\n'
            '                     /isolate="S85"\n'
            "//\n"
        )
        fields = _locus_fields(record)
        assert fields["version"] == "NC_TEST.2"
        assert fields["genome_length"] == 4848
        assert fields["topology"] == "circular"
        assert fields["molecule"] == "DNA"

    def test_version_line_wins_over_locus_token(self) -> None:
        record = "LOCUS       NC_TEST 100 bp DNA linear VRL\nVERSION     NC_TEST.9\n//\n"
        assert _locus_fields(record)["version"] == "NC_TEST.9"


class TestSourceQualifiers:
    def test_isolate_comes_from_the_source_feature(self) -> None:
        record = (
            "LOCUS       NC_TEST 2359 bp DNA linear VRL\n"
            "VERSION     NC_TEST.1\n"
            "FEATURES             Location/Qualifiers\n"
            "     source          1..2359\n"
            '                     /organism="Torque teno virus"\n'
            '                     /isolate="S85"\n'
            "     CDS             1..900\n"
            '                     /product="ORF1"\n'
            "//\n"
        )
        assert _source_qualifiers(record)["isolate"] == "S85"

    def test_absent_source_feature_returns_empty(self) -> None:
        assert _source_qualifiers("LOCUS NC_TEST 10 bp DNA linear VRL\n//\n") == {}


class TestCatalogueRows:
    _RECORD = (
        "LOCUS       NC_002076                4848 bp    DNA     linear   VRL\n"
        "VERSION     NC_002076.2\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..4848\n"
        '                     /organism="Torque teno virus"\n'
        '                     /strain="Tu243"\n'
        "     CDS             589..2901\n"
        '                     /gene="orf1"\n'
        '                     /product="Rep protein"\n'
        '                     /protein_id="NP_817122.1"\n'
        "     CDS             join(353..711,2564..3077)\n"
        '                     /product="VP3"\n'
        '                     /protein_id="NP_817120.1"\n'
        "ORIGIN\n//\n"
    )

    def test_one_row_per_cds_with_exon_structure(self) -> None:
        rows = catalogue_rows("NC_002076.2", self._RECORD)
        assert [row["gene_id"] for row in rows] == [
            "NC_002076.2_orf1",
            "NC_002076.2_NP_817120.1",
        ]
        assert rows[0]["exons"] == "589:2901"
        assert rows[0]["n_exons"] == 1
        assert rows[1]["exons"] == "353:711,2564:3077"
        assert rows[1]["n_exons"] == 2
        assert rows[0]["genome_length"] == 4848

    def test_genotype_is_empty_and_strain_is_carried(self) -> None:
        row = catalogue_rows("NC_002076.2", self._RECORD)[0]
        assert row["source_genotype"] == ""
        assert row["source_strain"] == "Tu243"

    def test_record_without_cds_yields_no_rows(self) -> None:
        record = (
            "LOCUS       KP343852 2359 bp DNA linear VRL\n"
            "VERSION     KP343852.1\n"
            "FEATURES             Location/Qualifiers\n"
            "     source          1..2359\n"
            '                     /isolate="S85"\n'
            "ORIGIN\n//\n"
        )
        assert catalogue_rows("KP343852.1", record) == []


class TestGenbankCache:
    def test_genbank_flatfile_is_cached_with_sidecar(self, tmp_path, monkeypatch) -> None:
        """fetch_genbank writes ``<acc>.gb`` so re-derivation costs no network."""
        from viralscan.scripts import ncbi_fetch

        calls: list[str] = []

        def mock_efetch(acc, rettype, email, api_key):
            calls.append(rettype)
            return (
                "LOCUS       NC_FAKE1 8 bp DNA linear VRL\nVERSION     NC_FAKE1.1\n"
                "FEATURES             Location/Qualifiers\n"
                '     CDS             1..8\n                     /gene="X"\n//\n'
            )

        monkeypatch.setattr(ncbi_fetch, "_efetch", mock_efetch)
        monkeypatch.setattr(ncbi_fetch, "_validate_accession", lambda x: x)

        path, text = ncbi_fetch.fetch_genbank("NC_FAKE1", "me@example.org", None, tmp_path)
        assert path == tmp_path / "NC_FAKE1" / "NC_FAKE1.gb"
        assert "FEATURES" in text
        assert path.with_suffix(".gb.sha256").exists()
        assert calls == ["gb"]

        again = ncbi_fetch.fetch_genbank("NC_FAKE1", "me@example.org", None, tmp_path)
        assert again[1] == text
        assert calls == ["gb"], "second call must be served from cache"

    def test_fetch_one_populates_the_genbank_cache(self, tmp_path, monkeypatch) -> None:
        from viralscan.scripts import ncbi_fetch

        def mock_efetch(acc, rettype, email, api_key):
            if rettype == "fasta":
                return ">NC_FAKE1\nATGCATGC\n"
            return (
                "LOCUS       NC_FAKE1 8 bp DNA linear VRL\nVERSION     NC_FAKE1.1\n"
                "FEATURES             Location/Qualifiers\n"
                '     CDS             1..8\n                     /gene="X"\n//\n'
            )

        monkeypatch.setattr(ncbi_fetch, "_efetch", mock_efetch)
        monkeypatch.setattr(ncbi_fetch, "_validate_accession", lambda x: x)

        _fasta, gtf = ncbi_fetch._fetch_one("NC_FAKE1", tmp_path, "me@example.org", None)
        assert (tmp_path / "NC_FAKE1" / "NC_FAKE1.gb").exists()
        assert 'gene_id "NC_FAKE1_X"' in gtf.read_text()


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
                return 'LOCUS NC_FAKE1\nFEATURES\n     CDS             1..8\n                     /gene="X"\n//\n'
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
