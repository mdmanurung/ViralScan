"""Tests for src/viralscan/scripts/build_reference.py.

All tests here run without network access (no HTTP, no NCBI calls).
Network-dependent integration tests are marked with @pytest.mark.network.
"""

import gzip
import json
import textwrap
from pathlib import Path
from unittest.mock import patch

import pytest

from viralscan.anellovirus import load_gene_table
from viralscan.constants import ENSEMBL_SPECIES
from viralscan.scripts.build_reference import (
    _ensembl_species_key,
    _genome_as_transcript_gtf,
    _parse_host_homology_paf,
    build_anellovirus_reference,
    host_cdna_as_gtf,
    index_gtf_by_seqname,
    validate_reference_records,
    viral_gtf_block,
    write_reference_manifest,
)

# ---------------------------------------------------------------------------
# CAT-01: the fetched NCBI GTF must survive into the combined reference
# ---------------------------------------------------------------------------

_HPV_GTF = (
    'NC_001526.4\tNCBI\texon\t1\t1950\t.\t+\t0\tgene_id "NC_001526.4_HpV16gp3"; '
    'transcript_id "NP_041327.2"; gene_name "E1";\n'
    'NC_001526.4\tNCBI\texon\t7604\t7900\t.\t+\t0\tgene_id "NC_001526.4_HpV16gp7"; '
    'transcript_id "NP_041326.1"; gene_name "E7";\n'
    'NC_045512.2\tNCBI\texon\t266\t21555\t.\t+\t0\tgene_id "NC_045512.2_orf1ab"; '
    'transcript_id "YP_009724389.1"; gene_name "ORF1ab";\n'
)
_HPV_FASTA = ">NC_001526.4 Human papillomavirus type 16\nACGTACGTAC\n"


class TestIndexGtfBySeqname:
    """The merged NCBI GTF splits per accession so each block can be reused."""

    def test_groups_lines_under_their_seqname(self) -> None:
        blocks = index_gtf_by_seqname(_HPV_GTF)
        assert set(blocks) == {"NC_001526.4", "NC_045512.2"}
        assert len(blocks["NC_001526.4"]) == 2
        assert len(blocks["NC_045512.2"]) == 1

    def test_skips_comments_and_blank_lines(self) -> None:
        blocks = index_gtf_by_seqname("# header\n\n" + _HPV_GTF)
        assert set(blocks) == {"NC_001526.4", "NC_045512.2"}

    def test_empty_text_yields_no_blocks(self) -> None:
        assert index_gtf_by_seqname("") == {}


class TestViralGtfBlock:
    """CAT-01: real CDS structure beats the whole-genome placeholder."""

    def test_uses_real_ncbi_annotation_when_available(self) -> None:
        block, source = viral_gtf_block(
            _HPV_FASTA,
            "NC_001526.4",
            real_gtf_blocks=index_gtf_by_seqname(_HPV_GTF),
        )
        assert source == "ncbi"
        # The real gene IDs survive; the placeholder's would be NC_001526.4_gene1.
        assert "NC_001526.4_HpV16gp3" in block
        assert "NC_001526.4_HpV16gp7" in block
        assert "_gene1" not in block
        assert "whole_genome" not in block

    def test_falls_back_to_placeholder_only_without_annotation(self) -> None:
        block, source = viral_gtf_block(_HPV_FASTA, "NC_001526.4", real_gtf_blocks={})
        assert source == "placeholder"
        assert 'gene_id "NC_001526.4_gene1"' in block
        assert "whole_genome" in block

    def test_unversioned_accession_still_matches(self) -> None:
        """A FASTA header carrying the bare accession must still find its block."""
        block, source = viral_gtf_block(
            _HPV_FASTA,
            "NC_001526",
            real_gtf_blocks=index_gtf_by_seqname(_HPV_GTF),
        )
        assert source == "ncbi"
        assert "NC_001526.4_HpV16gp3" in block

    def test_anellovirus_catalogue_wins_over_ncbi(self) -> None:
        """The packaged catalogue is preferred; it is the curated CDS structure."""
        accession = "NC_002076.2"
        fasta = f">{accession} Torque teno virus\n{'ACGT' * 40}\n"
        gtf = f'{accession}\tNCBI\texon\t1\t10\t.\t+\t0\tgene_id "{accession}_wrong";\n'
        block, source = viral_gtf_block(
            fasta,
            accession,
            anello_accessions={accession},
            real_gtf_blocks=index_gtf_by_seqname(gtf),
        )
        assert source == "catalogue"
        assert "_wrong" not in block
        assert "TTVgp" in block

    def test_non_anellovirus_accession_ignores_the_catalogue_set(self) -> None:
        block, source = viral_gtf_block(
            _HPV_FASTA,
            "NC_001526.4",
            anello_accessions={"NC_002076.2"},
            real_gtf_blocks=index_gtf_by_seqname(_HPV_GTF),
        )
        assert source == "ncbi"
        assert "NC_001526.4_HpV16gp3" in block


class TestReferenceManifest:
    def test_duplicate_identifier_fails_closed(self, tmp_path):
        fasta = tmp_path / "duplicate-id.fa"
        fasta.write_text(">A\nAAAA\n>A\nCCCC\n")
        with pytest.raises(ValueError, match="Duplicate FASTA identifier"):
            validate_reference_records(fasta)

    def test_duplicate_sequence_fails_closed(self, tmp_path):
        fasta = tmp_path / "duplicate-sequence.fa"
        fasta.write_text(">A\nAAAA\n>B\nAAAA\n")
        with pytest.raises(ValueError, match="Exact duplicate sequences"):
            validate_reference_records(fasta)

    def test_manifest_records_sequence_provenance(self, tmp_path):
        fasta = tmp_path / "reference.fa"
        fasta.write_text(">ENST1\nAAAA\n>NC_1.1\nCCCC\n")
        output = write_reference_manifest(
            fasta,
            tmp_path / "reference_manifest.json",
            profile="curated",
            host_species="human",
            viral_identifiers={"NC_1.1"},
        )

        manifest = json.loads(output.read_text())
        assert manifest["schema_version"] == "3.0.0"
        assert manifest["profile"] == "curated"
        assert len(manifest["fasta_sha256"]) == 64
        records = {record["accession_version"]: record for record in manifest["sequences"]}
        assert records["ENST1"]["taxonomy"] == "human"
        assert records["NC_1.1"]["taxonomy"] == "virus"
        assert records["NC_1.1"]["length"] == 4
        assert len(records["NC_1.1"]["sha256"]) == 64
        assert records["NC_1.1"]["low_complexity_flag"] is True

    def test_host_homology_paf_keeps_raw_best_alignment(self):
        paf = (
            "V1\t100\t0\t50\t+\tchr1\t1000\t10\t60\t45\t50\t60\n"
            "V1\t100\t0\t80\t+\tchr2\t1000\t10\t90\t60\t80\t40\n"
        )
        annotation = _parse_host_homology_paf(paf, {"V1": 100, "V2": 50})
        assert annotation["V1"]["host_homology_best_target"] == "chr2"
        assert annotation["V1"]["host_homology_max_identity"] == pytest.approx(0.75)
        assert annotation["V1"]["host_homology_max_query_coverage"] == pytest.approx(0.8)
        assert annotation["V2"]["host_homology_max_aligned_bases"] == 0


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
# host_cdna_as_gtf — cDNA-level host GTF (regression for the kb-ref hang)
# ---------------------------------------------------------------------------


# Two Ensembl-style cDNA records: a chromosome one and a scaffold one. The
# scaffold header (seqname-looking "KI270728.1") lives only in the description,
# never as the GTF seqname — that is the whole point of the fix.
_ENSEMBL_CDNA = (
    ">ENST00000632684.1 cdna chromosome:GRCh38:14:22438547:22438554:1 "
    "gene:ENSG00000282431.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene\n"
    "ACGTACGT\n"
    ">ENST00000448914.1 cdna scaffold:GRCh38:KI270728.1:100:112:-1 "
    "gene:ENSG00000228985.1 gene_biotype:TR_D_gene\n"
    "TTTTTGGGGG\n"
)


class TestHostCdnaAsGtf:
    def _write_and_read(self, tmp_path: Path) -> str:
        fasta_gz = tmp_path / "cdna.fa.gz"
        with gzip.open(fasta_gz, "wt") as fh:
            fh.write(_ENSEMBL_CDNA)
        out = tmp_path / "host_cdna.gtf"
        n = host_cdna_as_gtf(fasta_gz, out)
        assert n == 2
        return out.read_text()

    def test_seqname_is_transcript_id_not_chromosomal(self, tmp_path):
        gtf = self._write_and_read(tmp_path)
        seqnames = {line.split("\t", 1)[0] for line in gtf.splitlines() if line}
        # Every seqname must be an ENST transcript ID — the FASTA headers.
        assert seqnames == {"ENST00000632684.1", "ENST00000448914.1"}
        # And crucially: no chromosomal / scaffold seqname leaks in.
        for bad in ("14", "KI270728.1", "1"):
            assert bad not in seqnames

    def test_gene_id_from_gene_field(self, tmp_path):
        gtf = self._write_and_read(tmp_path)
        assert 'gene_id "ENSG00000282431.1"' in gtf
        assert 'gene_id "ENSG00000228985.1"' in gtf

    def test_coords_span_full_transcript_length(self, tmp_path):
        gtf = self._write_and_read(tmp_path)
        rows = [ln.split("\t") for ln in gtf.splitlines() if ln.startswith("ENST00000632684.1")]
        # 3 features (gene/transcript/exon), each spanning 1..8 (len("ACGTACGT")).
        assert len(rows) == 3
        for r in rows:
            assert r[3] == "1" and r[4] == "8", r

    def test_every_seqname_matches_a_fasta_header(self, tmp_path):
        """The core kb-ref contract: no GTF seqname without a FASTA record."""
        gtf = self._write_and_read(tmp_path)
        fasta_headers = {
            ln[1:].split()[0] for ln in _ENSEMBL_CDNA.splitlines() if ln.startswith(">")
        }
        gtf_seqnames = {line.split("\t", 1)[0] for line in gtf.splitlines() if line}
        assert gtf_seqnames <= fasta_headers


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
            fh.write(
                ">ENST000001.1 cdna chromosome:GRCh38:1:1:8:1 gene:ENSG000001.1 "
                "gene_biotype:protein_coding\nATCGATCG\n"
            )
        # The chromosomal GTF is intentionally IGNORED by build_combined_reference
        # (its seqname 'chr1' would not match the cDNA header) — kept only to exercise
        # that the combined GTF is NOT built from it.
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
        assert result["manifest"] is not None and result["manifest"].exists()
        assert result["fasta"].exists()
        assert result["gtf"].exists()
        assert ">ENST000001.1" in result["fasta"].read_text()
        assert ">NC_045512.2" in result["fasta"].read_text()
        combined_gtf = result["gtf"].read_text()
        # Host GTF is now cDNA-level: seqname = ENST, gene_id from the gene: field.
        assert 'gene_id "ENSG000001.1"' in combined_gtf
        assert any(ln.startswith("ENST000001.1\t") for ln in combined_gtf.splitlines()), (
            "combined GTF must carry the cDNA-level host seqname"
        )
        # The chromosomal host GTF must NOT leak into the combined GTF (the bug).
        assert 'gene_id "HOST1"' not in combined_gtf
        assert not any(ln.startswith("chr1\t") for ln in combined_gtf.splitlines())
        # CAT-01: the viral GTF fetched from NCBI is carried through, so the real
        # gene survives and the whole-genome placeholder is NOT used. Before the
        # fix this asserted `NC_045512.2_gene1`, the placeholder that replaced it.
        assert 'gene_id "V"' in combined_gtf
        assert 'gene_id "NC_045512.2_gene1"' not in combined_gtf

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


# ---------------------------------------------------------------------------
# build_anellovirus_reference — B.6 tests (network-free)
# ---------------------------------------------------------------------------

_ANELLO_FASTA = textwrap.dedent("""\
    >AB026929.1 Torque teno virus clone
    ATCGATCGATCGATCGATCG
    >NC_002076.2 Torque teno virus 1
    TTTTAAAACCCCGGGG
""")


class TestBuildAnellovirusReference:
    """B.6 — build_anellovirus_reference: FASTA+GTF output and graceful tool-absence handling."""

    def _setup_fake_ncbi(self, tmp_path: Path) -> tuple[Path, Path]:
        fasta = tmp_path / "ncbi" / "merged.fasta"
        fasta.parent.mkdir(parents=True, exist_ok=True)
        fasta.write_text(_ANELLO_FASTA)
        gtf = tmp_path / "ncbi" / "merged.gtf"
        gtf.write_text("")
        return fasta, gtf

    def test_produces_fasta_and_gtf_with_expected_gene_ids(self, tmp_path):
        fasta, gtf = self._setup_fake_ncbi(tmp_path)

        with patch("viralscan.scripts.ncbi_fetch.fetch_reference", return_value=(fasta, gtf)):
            result = build_anellovirus_reference(
                out_dir=tmp_path / "out",
                accessions=["AB026929.1", "NC_002076.2"],
                mask=False,
                cluster=False,
                run_kb_ref=False,
            )

        assert result["fasta"] is not None and result["fasta"].exists()
        assert result["gtf"] is not None and result["gtf"].exists()
        assert result["index"] is None
        assert result["t2g"] is None

        # Real NCBI CDS structure, not a whole-genome placeholder: the panel's
        # anellovirus_accessions.tsv covers both accessions, and a whole-genome
        # single-exon gene is a conservation bucket rather than a measurement.
        gtf_text = result["gtf"].read_text()
        assert "_gene1" not in gtf_text
        assert 'gene_biotype "whole_genome"' not in gtf_text
        assert 'gene_id "AB026929.1_BAA86944.1"' in gtf_text
        assert {
            row["gene_id"] for row in load_gene_table() if row["accession"] == "NC_002076.2"
        } == {
            "NC_002076.2_TTVgp1",
            "NC_002076.2_TTVgp2",
            "NC_002076.2_TTVgp3",
        }
        assert 'gene_id "NC_002076.2_TTVgp3"' in gtf_text

    def test_uncovered_accession_still_gets_a_whole_genome_placeholder(self, tmp_path):
        """kb ref silently drops a sequence with no GTF row, so coverage is total.

        An accession absent from the packaged gene catalogue must still be
        annotated, or the genome becomes neither quantified nor detectable.
        """
        uncovered = "ZZ999999.1"
        fasta = tmp_path / "ncbi" / "merged.fasta"
        fasta.parent.mkdir(parents=True, exist_ok=True)
        fasta.write_text(f">{uncovered} synthetic record\n" + "ACGTACGTAC" * 6 + "\n")
        gtf = tmp_path / "ncbi" / "merged.gtf"
        gtf.write_text("")

        with patch("viralscan.scripts.ncbi_fetch.fetch_reference", return_value=(fasta, gtf)):
            result = build_anellovirus_reference(
                out_dir=tmp_path / "out",
                accessions=[uncovered],
                mask=False,
                cluster=False,
                run_kb_ref=False,
            )

        gtf_text = result["gtf"].read_text()
        assert f'gene_id "{uncovered}_gene1"' in gtf_text
        assert 'gene_biotype "whole_genome"' in gtf_text

    def test_requested_mask_fails_when_dustmasker_absent(self, tmp_path):
        fasta, gtf = self._setup_fake_ncbi(tmp_path)

        with (
            patch("viralscan.scripts.ncbi_fetch.fetch_reference", return_value=(fasta, gtf)),
            patch(
                "viralscan.scripts.build_reference._run_dustmasker", return_value=False
            ) as mock_mask,
        ):
            with pytest.raises(RuntimeError, match="masking was requested"):
                build_anellovirus_reference(
                    out_dir=tmp_path / "out",
                    accessions=["AB026929.1"],
                    mask=True,
                    cluster=False,
                    run_kb_ref=False,
                )

        mock_mask.assert_called_once()

    def test_requested_cluster_fails_when_cdhit_absent(self, tmp_path):
        fasta, gtf = self._setup_fake_ncbi(tmp_path)

        with (
            patch("viralscan.scripts.ncbi_fetch.fetch_reference", return_value=(fasta, gtf)),
            patch(
                "viralscan.scripts.build_reference._run_cdhit_est", return_value=False
            ) as mock_clust,
        ):
            with pytest.raises(RuntimeError, match="clustering was requested"):
                build_anellovirus_reference(
                    out_dir=tmp_path / "out",
                    accessions=["AB026929.1"],
                    mask=False,
                    cluster=True,
                    run_kb_ref=False,
                )

        mock_clust.assert_called_once()

    def test_default_accessions_from_packaged_table(self, tmp_path):
        """When accessions=None, the packaged TSV is loaded and NCBI fetch is called."""
        fasta, gtf = self._setup_fake_ncbi(tmp_path)

        with patch(
            "viralscan.scripts.ncbi_fetch.fetch_reference", return_value=(fasta, gtf)
        ) as mock_fetch:
            result = build_anellovirus_reference(
                out_dir=tmp_path / "out",
                accessions=None,
                mask=False,
                cluster=False,
                run_kb_ref=False,
            )

        # The packaged TSV has >1000 accessions; fetch must have been called with a non-empty list.
        called_accessions = mock_fetch.call_args[0][0]
        assert isinstance(called_accessions, list) and len(called_accessions) > 100
        assert result["fasta"] is not None and result["fasta"].exists()
        assert result["gtf"] is not None and result["gtf"].exists()

    def test_empty_fasta_fails_closed(self, tmp_path):
        empty_fasta = tmp_path / "ncbi" / "merged.fasta"
        empty_fasta.parent.mkdir(parents=True, exist_ok=True)
        empty_fasta.write_text("")
        empty_gtf = tmp_path / "ncbi" / "merged.gtf"
        empty_gtf.write_text("")

        with patch(
            "viralscan.scripts.ncbi_fetch.fetch_reference", return_value=(empty_fasta, empty_gtf)
        ):
            with pytest.raises(ValueError, match="no sequences"):
                build_anellovirus_reference(
                    out_dir=tmp_path / "out",
                    accessions=["AB026929.1"],
                    mask=False,
                    cluster=False,
                    run_kb_ref=False,
                )
