"""Tests for src/viralscan/scripts/build_reference.py.

All tests here run without network access (no HTTP, no NCBI calls).
Network-dependent integration tests are marked with @pytest.mark.network.
"""

import gzip
import json
import random
import shutil
import subprocess
import textwrap
from pathlib import Path
from unittest.mock import patch

import pytest

from viralscan.anellovirus import load_gene_table
from viralscan.constants import ENSEMBL_SPECIES
from viralscan.scripts.build_reference import (
    _ensembl_species_key,
    _fasta_records,
    _genome_as_transcript_gtf,
    _max_tandem_period,
    _parse_host_homology_paf,
    _run_kb_ref,
    build_anellovirus_reference,
    host_cdna_as_gtf,
    index_gtf_by_seqname,
    low_complexity_kmer_counts,
    low_complexity_kmer_fraction,
    low_complexity_report,
    measure_host_homology,
    validate_reference_records,
    viral_gtf_block,
    write_reference_manifest,
)

# ---------------------------------------------------------------------------
# CAT-17: low-complexity k-mers, not N-masking, are what match poly-A reads
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fails", [False, True])
def test_kb_ref_scratch_is_unique_and_cleaned_even_on_failure(tmp_path, monkeypatch, fails):
    monkeypatch.chdir(tmp_path)
    existing = tmp_path / "tmp"
    existing.mkdir()
    (existing / "user-data").write_text("preserve me")
    out = tmp_path / "out"
    out.mkdir()
    scratch = []
    command = ["kb", "ref", "-i", str(out / "index.idx"), "reference.fa", "reference.gtf"]

    def fake_run(cmd, check):
        path = Path(cmd[cmd.index("--tmp") + 1])
        assert check is True
        assert path.parent.parent == out and not path.exists()
        assert path not in scratch
        scratch.append(path)
        path.mkdir()
        (path / "partial").write_text("scratch")
        assert cmd[-2:] == command[-2:]
        if fails:
            raise subprocess.CalledProcessError(2, cmd)

    monkeypatch.setattr("viralscan.scripts.build_reference.subprocess.run", fake_run)
    for _ in range(2):
        if fails:
            with pytest.raises(subprocess.CalledProcessError):
                _run_kb_ref(command, out)
        else:
            _run_kb_ref(command, out)
        assert not scratch[-1].parent.exists()
    assert (existing / "user-data").read_text() == "preserve me"
    assert "--tmp" not in command


def _complex_sequence(length: int) -> str:
    """Deterministic sequence with no short-period structure.

    Do not use ``"ACGT" * n`` as a stand-in for clean sequence: that is a perfect
    4-base tandem repeat, which the gate is designed to flag.
    """
    return "".join("ACGT"[(i * 7 + (i // 3) + (i % 5)) % 4] for i in range(length))


def _homopolymer(sequence: str, run: int) -> str:
    """Sequence carrying a homopolymer of `run` bases at each end."""
    return sequence + "A" * run


class TestLowComplexityKmers:
    """A panel can be almost fully unmasked and still match poly-A tails.

    Measured on the EBV LCL `SRR12682296`: 1.44 % of R2 reads matched the
    anellovirus panel, and 6,437 of 6,437 captured hit reads had *zero* genuine
    anellovirus k-mers once low-complexity k-mers were excluded. The matched
    k-mers were literally `A`*31 and near-neighbours of it.
    """

    def test_pure_homopolymer_is_low_complexity(self):
        low, total = low_complexity_kmer_fraction("A" * 31)
        assert total == 1
        assert low == 1

    def test_long_tail_pushes_fraction_up(self):
        backbone = _complex_sequence(40)  # no runs > 1, no short-period structure
        base_low, base_total = low_complexity_kmer_fraction(backbone)
        assert base_low == 0
        assert base_total == 10
        tailed_low, tailed_total = low_complexity_kmer_fraction(_homopolymer(backbone, 31))
        # 31 A's on the end add 31 windows, all of which span or sit in the run
        assert tailed_total == 41
        assert tailed_low > base_low
        assert tailed_low / tailed_total > base_low / base_total

    def test_dinucleotide_repeat_is_low_complexity(self):
        low, _ = low_complexity_kmer_fraction("AT" * 20)
        assert low > 0

    def test_sub_tiling_trinucleotide_repeat_is_caught(self):
        """A repeat whose unit does not divide k=31 must still be caught.

        `AB303556.1` carries a 28 bp CAG trinucleotide repeat at 2318-2345. It
        produced 1,485 of the 1,852 raw anellovirus reads in the SFL tonsil
        screen (F-010). An exact-tiling test cannot see it: 31 is prime, so a
        3-base unit never tiles a 31-mer. This is the regression that an
        exact-tiling implementation silently misses.
        """
        repeat = "CAG" * 12
        window = repeat[:31]
        assert 31 % 3 != 0  # the reason an exact-tiling test fails
        assert _max_tandem_period(window, 5) == 3
        assert low_complexity_kmer_counts(window)["tandem"] == 1

    def test_periodic_detector_ignores_non_periodic_sequence(self):
        pseudo = "".join("ACGT"[(i * 7 + i // 3 + (i % 5)) % 4] for i in range(31))
        assert _max_tandem_period(pseudo, 5) == 0

    def test_periodic_detector_reports_the_period(self):
        # a homopolymer is periodic under every period, so it reports the max
        assert _max_tandem_period("A" * 31, 5) == 5
        assert _max_tandem_period("CAG" * 11, 5) == 3
        # a pure 3-unit repeat is periodic under 3, and under nothing shorter
        assert _max_tandem_period("CAG" * 11, 2) == 0
        assert _max_tandem_period("GATA" * 8, 4) == 4
        # a strictly alternating string agrees under EVEN shifts only — p=1,3,5 are
        # anti-phase — so the largest agreeing period below 5 is 4, not 5
        assert _max_tandem_period("AT" * 16, 5) == 4

    def test_complex_sequence_is_clean(self):
        # A deterministic non-repetitive sequence: no runs, no short tandem repeat
        seq = "".join("ACGT"[(i * 7 + i // 3) % 4] for i in range(200))
        low, total = low_complexity_kmer_fraction(seq)
        assert total > 0
        assert low == 0

    def test_non_acgt_windows_are_excluded(self):
        # kallisto replaces non-ACGT, so those k-mers cannot reach the index as written.
        # Of the 33 windows here, only the two pure runs survive as candidates.
        low, total = low_complexity_kmer_fraction("A" * 31 + "N" + "C" * 31)
        assert total == 2
        assert low == 2  # both surviving candidates are themselves low-complexity

    def test_poly_a_read_would_have_matched_the_bad_panel(self):
        """The regression itself: a poly-A read's k-mer is in an unmasked panel."""
        panel = "ACGTACGTACGTTTTTACGTACGTACGTACGTACGTACGTACGT"
        read_kmer = "A" * 31
        panel_kmers = {panel[i : i + 31] for i in range(len(panel) - 30)}
        assert read_kmer not in panel_kmers  # sanity: not in a short clean panel
        bad_panel = panel + "A" * 40
        bad_kmers = {bad_panel[i : i + 31] for i in range(len(bad_panel) - 30)}
        assert read_kmer in bad_kmers
        low, _ = low_complexity_kmer_fraction(bad_panel)
        assert low > 0

    def test_report_is_per_record(self, tmp_path):
        fasta = tmp_path / "panel.fa"
        fasta.write_text(
            ">good\n" + _complex_sequence(40) + "\n>bad\n" + _complex_sequence(20) + "A" * 40 + "\n"
        )
        report = low_complexity_report(fasta)
        assert set(report) == {"good", "bad"}
        assert report["good"][0] == 0
        assert report["bad"][0] > 0

    def test_gate_rejects_an_unmasked_panel(self, tmp_path):
        fasta = tmp_path / "panel.fa"
        fasta.write_text(">anello\n" + _complex_sequence(55) + "A" * 45 + "\n")
        with pytest.raises(ValueError, match="k-mers are low-complexity"):
            validate_reference_records(fasta, max_low_complexity_fraction=0.0)

    def test_pure_homopolymer_gate_is_absolute_not_fractional(self, tmp_path):
        """A handful of pure-homopolymer k-mers is enough, so gate on a count.

        Measured: 170 pure 31-mers across 9 records of the deployed panel
        manufactured >1 % of R2 reads, because poly-A reads are common even
        though the k-mers are not.
        """
        # 3000 bp of clean sequence plus 32 A's: a negligible fraction of k-mers,
        # but 2 of them are pure homopolymers
        clean = "ACGTTGCAAGTCAG" * 220 + "A" * 32
        low, total = low_complexity_kmer_fraction(clean)
        counts = low_complexity_kmer_counts(clean)
        assert counts["pure_homopolymer"] == 2
        assert low / total < 0.02  # a strict fraction gate would let this through

        fasta = tmp_path / "panel.fa"
        fasta.write_text(f">anello\n{clean}\n")
        with pytest.raises(ValueError, match="pure-homopolymer"):
            validate_reference_records(fasta, max_pure_homopolymer_kmers=0)

    def test_pure_homopolymer_breakdown(self):
        assert low_complexity_kmer_counts("A" * 31)["pure_homopolymer"] == 1
        assert low_complexity_kmer_counts("A" * 62)["pure_homopolymer"] == 32
        # 30 A's after a clean backbone: long run, but no full pure 31-mer
        counts = low_complexity_kmer_counts("ACGTTGCAAGTCAG" * 10 + "A" * 30)
        assert counts["pure_homopolymer"] == 0
        assert counts["long_run"] > 0

    def test_masked_reference_passes_both_gates(self, tmp_path):
        fasta = tmp_path / "panel.fa"
        fasta.write_text(">anello\n" + _complex_sequence(60) + "\n")
        validate_reference_records(
            fasta, max_low_complexity_fraction=0.0, max_pure_homopolymer_kmers=0
        )

    def test_gate_default_is_backward_compatible(self, tmp_path):
        fasta = tmp_path / "panel.fa"
        fasta.write_text(">anello\n" + _complex_sequence(25) + "A" * 45 + "\n")
        validate_reference_records(fasta)  # no limits -> no raise

    def test_gate_error_names_dlist_cannot_help(self, tmp_path):
        fasta = tmp_path / "panel.fa"
        fasta.write_text(">anello\n" + "A" * 45 + "\n")
        with pytest.raises(ValueError, match="D-list cannot fix this"):
            validate_reference_records(fasta, max_pure_homopolymer_kmers=0)

    def test_manifest_records_kmer_counts(self, tmp_path):
        fasta = tmp_path / "panel.fa"
        fasta.write_text(">anello\n" + _complex_sequence(40) + "A" * 40 + "\n")
        write_reference_manifest(
            fasta,
            tmp_path / "manifest.json",
            profile="curated",
            host_species="homo_sapiens",
            viral_identifiers={"anello"},
        )
        row = json.loads((tmp_path / "manifest.json").read_text())["sequences"][0]
        assert row["low_complexity_kmers"] > 0
        assert 0.0 < row["low_complexity_kmer_fraction"] < 1.0

    def test_masking_with_n_removes_the_kmer(self, tmp_path):
        """The upstream fix: N-masking does reduce the k-mer count."""
        clean = "ACGTACGTACG" * 5 + "A" * 45
        masked = "ACGTACGTACG" * 5 + "N" * 45
        assert low_complexity_kmer_fraction(clean)[0] > 0
        assert low_complexity_kmer_fraction(masked)[0] == 0


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
        # REF-02: a local FASTA proves neither NCBI retrieval nor viral taxonomy.
        assert records["ENST1"]["taxonomy"] == {"organism": "human", "taxid": None}
        assert records["NC_1.1"]["taxonomy"] == {"organism": None, "taxid": None}
        assert records["NC_1.1"]["retrieved_at"] is None
        assert records["NC_1.1"]["length"] == 4
        assert len(records["NC_1.1"]["sha256"]) == 64
        assert records["NC_1.1"]["low_complexity_flag"] is True

    def test_host_homology_paf_keeps_raw_best_alignment(self):
        # Viral panel is the minimap2 reference: query = host contig, target = viral record.
        paf = (
            "chr1\t5000\t10\t60\t+\tV1\t100\t0\t50\t45\t50\t60\n"
            "chr2\t9000\t10\t90\t+\tV1\t100\t0\t80\t60\t80\t40\n"
            "chr3\t9000\t10\t90\t+\tUNKNOWN\t100\t0\t80\t60\t80\t40\n"
        )
        annotation = _parse_host_homology_paf(paf, {"V1": 100, "V2": 50})
        assert annotation["V1"]["host_homology_best_target"] == "chr2"
        assert annotation["V1"]["host_homology_max_identity"] == pytest.approx(0.75)
        # coverage is of the VIRAL genome (80 of 100 bases), not of the host contig
        assert annotation["V1"]["host_homology_max_query_coverage"] == pytest.approx(0.8)
        assert annotation["V1"]["host_homology_max_aligned_bases"] == 80
        assert annotation["V2"]["host_homology_max_aligned_bases"] == 0

    @pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not installed")
    def test_measure_host_homology_finds_a_viral_fragment_in_a_host_genome(self, tmp_path):
        rng = random.Random(1)

        def seq(n):
            return "".join(rng.choice("ACGT") for _ in range(n))

        fragment = seq(600)
        mutated = list(fragment)
        for i in range(0, 600, 25):  # 4 % substitutions: diverged, not identical
            mutated[i] = "ACGT"[("ACGT".index(mutated[i]) + 1) % 4]
        viral_homologous = seq(700) + fragment + seq(700)
        viral_unrelated = seq(2000)
        host = seq(30000) + "".join(mutated) + seq(30000)
        viral_fasta = tmp_path / "viral.fa"
        viral_fasta.write_text(f">V_HOM\n{viral_homologous}\n>V_NONE\n{viral_unrelated}\n")
        host_fasta = tmp_path / "host.fa"
        host_fasta.write_text(f">chrT\n{host}\n")

        result = measure_host_homology(viral_fasta, host_fasta, tmp_path / "out.tsv")

        hom = result["V_HOM"]
        assert hom["host_homology_best_target"] == "chrT"
        assert 500 <= int(hom["host_homology_max_aligned_bases"]) <= 620
        assert float(hom["host_homology_max_identity"]) > 0.9
        # per VIRAL genome: ~600 of 2,000 bases, not 600 of the 60 kb host contig
        assert 0.25 <= float(hom["host_homology_max_query_coverage"]) <= 0.31
        assert result["V_NONE"]["host_homology_max_aligned_bases"] == 0
        assert (tmp_path / "out.tsv").read_text().startswith("accession_version\t")


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
            # Assembly unit test: do not invoke BLAST+ on the tiny mocked input.
            patch(
                "viralscan.scripts.build_reference.mask_low_complexity",
                side_effect=lambda source, target: bool(shutil.copyfile(source, target)),
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


@pytest.mark.skipif(shutil.which("dustmasker") is None, reason="needs BLAST+ dustmasker")
def test_masking_a_poly_a_stretch_passes_the_gate_and_writes_n_not_lowercase(tmp_path):
    """PLAN MASK-01: dustmasker soft-masks; the gate and kallisto upper-case, so only N counts."""
    import random

    from viralscan.scripts.build_reference import (
        build_anellovirus_reference,
        low_complexity_kmer_counts,
    )

    rng = random.Random(3)

    def rnd(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    fasta = tmp_path / "in.fa"
    # One record carries a poly-A stretch; the other arrives already soft-masked (lowercase).
    fasta.write_text(
        f">AB000001.1\n{rnd(400)}{'A' * 60}{rnd(400)}\n>AB000002.1\n{rnd(300)}{'t' * 5}{rnd(300)}\n"
    )

    result = build_anellovirus_reference(
        out_dir=tmp_path / "out", mask=True, run_kb_ref=False, fasta_path=fasta
    )

    text = "".join(
        line.strip()
        for line in result["fasta"].read_text().splitlines()
        if not line.startswith(">")
    )
    assert "N" in text
    assert text == text.upper()
    assert low_complexity_kmer_counts(text)["pure_homopolymer"] == 0


def _combined_build_with_viral_sequence(tmp_path, viral_sequence, **kwargs):
    """Run build_combined_reference on a mocked host and one mocked viral record."""
    from viralscan.scripts.build_reference import build_combined_reference

    host_dir = tmp_path / "host"
    host_dir.mkdir()
    cdna_gz = host_dir / "fake.cdna.all.fa.gz"
    gtf_gz = host_dir / "fake.gtf.gz"
    with gzip.open(cdna_gz, "wt") as fh:
        # A poly-A host record: the gate must look at the viral FASTA only.
        fh.write(
            ">ENST000001.1 cdna chromosome:GRCh38:1:1:8:1 gene:ENSG000001.1\n" + "A" * 80 + "\n"
        )
    with gzip.open(gtf_gz, "wt") as fh:
        fh.write('chr1\tEnsembl\texon\t1\t8\t.\t+\t.\tgene_id "HOST1";\n')
    viral_dir = tmp_path / "viral"
    viral_dir.mkdir()
    fasta = viral_dir / "viral.fasta"
    fasta.write_text(f">NC_045512.2\n{viral_sequence}\n")
    gtf = viral_dir / "viral.gtf"
    gtf.write_text('NC_045512.2\tNCBI\texon\t1\t8\t.\t+\t0\tgene_id "V";\n')
    with (
        patch("viralscan.scripts.build_reference.fetch_host_cdna", return_value=(cdna_gz, gtf_gz)),
        patch("viralscan.scripts.ncbi_fetch.fetch_reference", return_value=(fasta, gtf)),
    ):
        return build_combined_reference(
            host_species="human",
            virus_accessions=["NC_045512.2"],
            out_dir=tmp_path / "ref",
            run_kb_ref=False,
            **kwargs,
        )


def _random_dna(n, seed=5):
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.mark.skipif(shutil.which("dustmasker") is None, reason="needs BLAST+ dustmasker")
def test_combined_build_masks_the_viral_panel_like_the_dedicated_builder(tmp_path):
    """PLAN REF-03: same N-masking and k-mer gate; host cDNA is not touched."""
    seq = _random_dna(400) + "A" * 60 + _random_dna(400, seed=6)

    result = _combined_build_with_viral_sequence(tmp_path, seq)

    records = dict(_fasta_records(result["fasta"]))
    assert "N" in records["NC_045512.2"]
    assert len(records["NC_045512.2"]) == len(seq)  # lengths unchanged: GTF coordinates stay valid
    assert records["ENST000001.1"] == "A" * 80  # host cDNA untouched
    assert "a" not in "".join(
        ln for ln in result["fasta"].read_text().splitlines() if not ln.startswith(">")
    )


def test_combined_build_fails_when_requested_mask_cannot_run(tmp_path):
    with patch("viralscan.scripts.build_reference._run_dustmasker", return_value=False):
        with pytest.raises(RuntimeError, match="Masking was requested"):
            _combined_build_with_viral_sequence(tmp_path, _random_dna(500))


def test_combined_build_with_no_mask_still_gates_low_complexity(tmp_path):
    seq = _random_dna(300) + "A" * 60 + _random_dna(300, seed=6)

    with pytest.raises(RuntimeError, match="low-complexity k-mer gate"):
        _combined_build_with_viral_sequence(tmp_path, seq, mask=False)


def test_build_ref_main_preflights_dustmasker_before_any_download(tmp_path):
    import argparse

    from viralscan.scripts.build_reference import build_ref_main

    args = argparse.Namespace(
        list_species=False, no_kb_ref=True, no_mask=False, genome_dlist=None, host="human"
    )
    with (
        patch("viralscan.scripts.build_reference.shutil.which", return_value=None),
        patch("viralscan.scripts.build_reference.fetch_host_cdna") as fetch,
    ):
        with pytest.raises(SystemExit) as exc:
            build_ref_main(args)

    assert exc.value.code == 2
    fetch.assert_not_called()
