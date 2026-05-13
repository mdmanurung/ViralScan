"""Tests for src/viralscan/scripts/build_reference.py.

All tests here run without network access (no HTTP, no NCBI calls).
Network-dependent integration tests are marked with @pytest.mark.network.
"""

import gzip
import textwrap
from pathlib import Path
from unittest.mock import patch

import pytest

from viralscan.constants import ENSEMBL_SPECIES
from viralscan.scripts.build_reference import (
    _ensembl_species_key,
    _genome_as_transcript_gtf,
)


# ---------------------------------------------------------------------------
# Species lookup
# ---------------------------------------------------------------------------


class TestEnsemblSpeciesKey:
    def test_known_short_name(self):
        assert _ensembl_species_key("human") == "human"

    def test_case_insensitive(self):
        assert _ensembl_species_key("Human") == "human"
        assert _ensembl_species_key("MOUSE") == "mouse"

    def test_spaces_converted_to_underscores(self):
        # Some callers may type "mus musculus" — should resolve
        # to the ensembl name lookup path
        ens_name = ENSEMBL_SPECIES["mouse"][0]  # "mus_musculus"
        assert _ensembl_species_key(ens_name) == "mouse"

    def test_ensembl_name_resolves(self):
        assert _ensembl_species_key("homo_sapiens") == "human"

    def test_unknown_species_raises(self):
        with pytest.raises(ValueError, match="Unknown host species"):
            _ensembl_species_key("unicorn")

    def test_all_registry_entries_round_trip(self):
        for short in ENSEMBL_SPECIES:
            assert _ensembl_species_key(short) == short


# ---------------------------------------------------------------------------
# ENSEMBL_SPECIES registry sanity
# ---------------------------------------------------------------------------


class TestEnsemblSpeciesRegistry:
    def test_human_present(self):
        assert "human" in ENSEMBL_SPECIES

    def test_mouse_present(self):
        assert "mouse" in ENSEMBL_SPECIES

    def test_all_entries_have_two_tuple(self):
        for key, val in ENSEMBL_SPECIES.items():
            assert isinstance(val, tuple) and len(val) == 2, key
            assert all(isinstance(s, str) and s for s in val), key

    def test_ensembl_names_are_lowercase(self):
        for key, (ens, _) in ENSEMBL_SPECIES.items():
            assert ens == ens.lower(), f"{key}: Ensembl name should be lowercase"


# ---------------------------------------------------------------------------
# _genome_as_transcript_gtf
# ---------------------------------------------------------------------------

SIMPLE_FASTA = textwrap.dedent("""\
    >NC_045512.2 Severe acute respiratory syndrome coronavirus 2
    ATTTATTTTCTTATTTAAGAC
    CCAGGTGATGTTTTGGATTTGTCT
    >NC_001477.1 Dengue virus 1
    AGTTGTTAGTCTACGTGGACC
""")


class TestGenomeAsTranscriptGtf:
    def test_returns_string(self):
        result = _genome_as_transcript_gtf(SIMPLE_FASTA, "NC_045512.2")
        assert isinstance(result, str)

    def test_contains_gene_transcript_exon(self):
        result = _genome_as_transcript_gtf(SIMPLE_FASTA, "NC_045512.2")
        features = {line.split("\t")[2] for line in result.splitlines() if line.strip()}
        assert {"gene", "transcript", "exon"}.issubset(features)

    def test_gene_biotype_whole_genome(self):
        result = _genome_as_transcript_gtf(SIMPLE_FASTA, "NC_045512.2")
        assert 'gene_biotype "whole_genome"' in result

    def test_sequence_length_in_coords(self):
        # First seq: 21+24 = 45 bases, second: 21 bases
        result = _genome_as_transcript_gtf(SIMPLE_FASTA, "NC_045512.2")
        lines = result.splitlines()
        # All starts should be 1
        for line in lines:
            parts = line.split("\t")
            assert parts[3] == "1", f"Expected start=1, got {parts[3]}"
        # end for first seq rows should be 45
        first_seq_ends = {int(line.split("\t")[4]) for line in lines if "NC_045512.2_gene1" in line}
        assert first_seq_ends == {45}

    def test_accession_in_ids(self):
        accession = "TEST_ACC"
        result = _genome_as_transcript_gtf(SIMPLE_FASTA, accession)
        assert f'gene_id "{accession}_gene1"' in result

    def test_empty_fasta_returns_empty_string(self):
        result = _genome_as_transcript_gtf("", "ACC")
        assert result == ""

    def test_single_sequence(self):
        fasta = ">SEQ1\nATCGATCGATCG\n"
        result = _genome_as_transcript_gtf(fasta, "ACC")
        lines = [ln for ln in result.splitlines() if ln.strip()]
        assert len(lines) == 3  # gene + transcript + exon

    def test_seqname_matches_fasta_header_first_token(self):
        fasta = ">chr1 some description\nATCG\n"
        result = _genome_as_transcript_gtf(fasta, "ACC")
        for line in result.splitlines():
            assert line.startswith("chr1\t"), line

    def test_multiple_sequences_numbered(self):
        fasta = ">seq1\nAAAA\n>seq2\nCCCC\n>seq3\nGGGG\n"
        result = _genome_as_transcript_gtf(fasta, "VIR")
        assert 'gene_id "VIR_gene1"' in result
        assert 'gene_id "VIR_gene2"' in result
        assert 'gene_id "VIR_gene3"' in result


# ---------------------------------------------------------------------------
# fetch_host_cdna — mock the download layer
# ---------------------------------------------------------------------------


class TestFetchHostCdna:
    def _make_fake_gz(self, content: bytes, path: Path) -> None:
        with gzip.open(path, "wb") as fh:
            fh.write(content)

    def test_raises_on_unknown_species(self, tmp_path):
        from viralscan.scripts.build_reference import fetch_host_cdna

        with pytest.raises(ValueError, match="Unknown host species"):
            fetch_host_cdna("unicorn", tmp_path)

    def test_uses_cache_when_present(self, tmp_path):
        """If cache files already exist, no download should happen."""
        from viralscan.scripts.build_reference import fetch_host_cdna

        cache_dir = tmp_path / "cache"
        out_dir = tmp_path / "out"
        species_cache = cache_dir / "mouse"
        species_cache.mkdir(parents=True)

        # Pre-populate cache with fake files
        fake_cdna = species_cache / "Mus_musculus.GRCm39.cdna.all.fa.gz"
        fake_gtf = species_cache / "Mus_musculus.GRCm39.109.gtf.gz"
        self._make_fake_gz(b">tx1\nATCG\n", fake_cdna)
        self._make_fake_gz(b"# fake gtf\n", fake_gtf)

        with (
            patch(
                "viralscan.scripts.build_reference._list_ensembl_files",
                side_effect=[
                    ["Mus_musculus.GRCm39.cdna.all.fa.gz"],
                    ["Mus_musculus.GRCm39.109.gtf.gz"],
                ],
            ),
            patch("viralscan.scripts.build_reference._download") as mock_dl,
        ):
            result_fasta, result_gtf = fetch_host_cdna("mouse", out_dir, cache_dir)
            # _download should NOT have been called because cache files exist
            mock_dl.assert_not_called()

        assert result_fasta.exists()
        assert result_gtf.exists()


# ---------------------------------------------------------------------------
# build_combined_reference — integration (mocked)
# ---------------------------------------------------------------------------


class TestBuildCombinedReference:
    def test_combines_mocked_host_and_viral_reference_without_kb_ref(self, tmp_path):
        from viralscan.scripts.build_reference import build_combined_reference

        # Prepare fake host cDNA FASTA (.gz) and GTF (.gz)
        host_dir = tmp_path / "host"
        host_dir.mkdir()
        fake_cdna_gz = host_dir / "fake.cdna.all.fa.gz"
        fake_gtf_gz = host_dir / "fake.109.gtf.gz"
        with gzip.open(fake_cdna_gz, "wt") as fh:
            fh.write(">ENST000001\nATCGATCG\n")
        with gzip.open(fake_gtf_gz, "wt") as fh:
            fh.write('chr1\tEnsembl\texon\t1\t8\t.\t+\t.\tgene_id "HOST1";\n')

        # Fake viral FASTA and GTF
        viral_dir = tmp_path / "viral"
        viral_dir.mkdir()
        fake_viral_fasta = viral_dir / "viral.fasta"
        fake_viral_fasta.write_text(">NC_045512.2\nATTTTGGG\n")
        fake_viral_gtf = viral_dir / "viral.gtf"
        fake_viral_gtf.write_text('NC_045512.2\tNCBI\texon\t1\t8\t.\t+\t0\tgene_id "V";\n')

        with (
            patch(
                "viralscan.scripts.build_reference.fetch_host_cdna",
                return_value=(fake_cdna_gz, fake_gtf_gz),
            ),
            patch(
                "viralscan.scripts.ncbi_fetch.fetch_reference",
                return_value=(fake_viral_fasta, fake_viral_gtf),
            ),
        ):
            result = build_combined_reference(
                host_species="human",
                virus_accessions=["NC_045512.2"],
                out_dir=tmp_path / "ref",
                run_kb_ref=False,
            )

        assert result["fasta"] == tmp_path / "ref" / "combined.fa"
        assert result["gtf"] == tmp_path / "ref" / "combined.gtf"
        assert result["index"] is None
        assert result["t2g"] is None
        assert result["fasta"].exists()
        assert result["gtf"].exists()
        assert ">ENST000001" in result["fasta"].read_text()
        assert ">NC_045512.2" in result["fasta"].read_text()
        combined_gtf = result["gtf"].read_text()
        assert 'gene_id "HOST1"' in combined_gtf
        assert 'gene_id "NC_045512.2_gene1"' in combined_gtf

    @pytest.mark.network
    def test_network_build_sars_cov2(self, tmp_path):
        """Full integration test: download SARS-CoV-2 from NCBI + human cDNA subset."""
        # This test is slow and requires network. Only run with -m network.
        from viralscan.scripts.build_reference import build_combined_reference

        result = build_combined_reference(
            host_species="human",
            virus_accessions=["NC_045512.2"],
            out_dir=tmp_path / "ref",
            run_kb_ref=False,
        )
        assert result["fasta"].exists()
        assert result["gtf"].exists()
        assert result["index"] is None
        assert result["t2g"] is None
