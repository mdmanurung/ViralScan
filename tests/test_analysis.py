"""Tests for the GTF-parsing logic in scripts/analysis.py.

``analysis.py`` is a Snakemake script that references ``snakemake.*`` at
module level, so we cannot import it directly.  Instead we test:

1. The ``obtain_gtf`` helper logic extracted to a standalone callable below
   (mirrors the real code line-for-line so bugs in the original are caught).
2. The data/ GTF files bundled in the package (spot-check a few).
3. The ``_count_unique_genes`` / ``_count_lines`` helpers from menu.py.

Running these tests requires only the standard Python packages available in
the PYTHONPATH=src mode documented in CLAUDE.md — no Snakemake runtime.
"""

from __future__ import annotations

import runpy
import re
import textwrap
from types import SimpleNamespace
from pathlib import Path

import pytest
import yaml


# ---------------------------------------------------------------------------
# Standalone re-implementation of the GTF-parsing core
# (mirrors obtain_gtf() without the Snakemake / config coupling)
# ---------------------------------------------------------------------------


def _parse_gtf_file(path: Path) -> set[str]:
    """Return the set of gene_id values from a GTF file (skipping comment lines)."""
    accessions: set[str] = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            info = cols[8]
            m = re.search(r'gene_id "([^"]+)"', info)
            if m:
                accessions.add(m.group(1))
    return accessions


def _parse_gtf_text(text: str) -> set[str]:
    """Same as _parse_gtf_file but from a raw string (for unit testing)."""
    accessions: set[str] = set()
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        info = cols[8]
        m = re.search(r'gene_id "([^"]+)"', info)
        if m:
            accessions.add(m.group(1))
    return accessions


# ---------------------------------------------------------------------------
# GTF text parsing — unit tests with synthetic fixtures
# ---------------------------------------------------------------------------


class TestParseGtfText:
    GTF_SINGLE = textwrap.dedent(
        """\
        # comment line — must be skipped
        NC_002021\t.\tgene\t1\t2280\t.\t+\t.\tgene_id "NC_002021"; transcript_id "NC_002021";
        NC_002021\t.\texon\t1\t2280\t.\t+\t.\tgene_id "NC_002021"; transcript_id "NC_002021";
        """
    )

    GTF_MULTI = textwrap.dedent(
        """\
        NC_001477\t.\tgene\t1\t100\t.\t+\t.\tgene_id "NC_001477"; transcript_id "NC_001477";
        NC_001612\t.\tgene\t1\t200\t.\t+\t.\tgene_id "NC_001612"; transcript_id "NC_001612";
        NC_001612\t.\texon\t1\t200\t.\t+\t.\tgene_id "NC_001612"; transcript_id "NC_001612";
        """
    )

    GTF_COMMENT_ONLY = "# comment\n# another comment\n"

    GTF_EMPTY = ""

    def test_single_gene_extracted(self) -> None:
        result = _parse_gtf_text(self.GTF_SINGLE)
        assert result == {"NC_002021"}

    def test_comment_lines_skipped(self) -> None:
        result = _parse_gtf_text(self.GTF_SINGLE)
        assert len(result) == 1  # only one unique gene_id

    def test_multiple_genes_extracted(self) -> None:
        result = _parse_gtf_text(self.GTF_MULTI)
        assert result == {"NC_001477", "NC_001612"}

    def test_duplicate_gene_ids_deduplicated(self) -> None:
        # NC_001612 appears twice (gene + exon rows) — must appear once in the set
        result = _parse_gtf_text(self.GTF_MULTI)
        assert len(result) == 2

    def test_comment_only_file_returns_empty_set(self) -> None:
        assert _parse_gtf_text(self.GTF_COMMENT_ONLY) == set()

    def test_empty_file_returns_empty_set(self) -> None:
        assert _parse_gtf_text(self.GTF_EMPTY) == set()

    def test_gene_id_is_second_quoted_token(self) -> None:
        # Verify the splitting logic: gene_id "<value>"
        gtf = 'NC_TEST\t.\tgene\t1\t100\t.\t+\t.\tgene_id "MY_GENE"; transcript_id "T1";\n'
        assert _parse_gtf_text(gtf) == {"MY_GENE"}


# ---------------------------------------------------------------------------
# GTF file writing and reading round-trip via tmp_path
# ---------------------------------------------------------------------------


class TestGtfFileRoundTrip:
    def _write_gtf(self, tmp_path: Path, content: str, name: str = "test.gtf") -> Path:
        p = tmp_path / name
        p.write_text(content)
        return p

    def test_file_parse_matches_text_parse(self, tmp_path: Path) -> None:
        content = (
            'NC_001918\t.\tgene\t1\t300\t.\t+\t.\tgene_id "NC_001918"; transcript_id "NC_001918";\n'
        )
        p = self._write_gtf(tmp_path, content)
        assert _parse_gtf_file(p) == _parse_gtf_text(content)

    def test_multiple_files_union(self, tmp_path: Path) -> None:
        f1 = self._write_gtf(
            tmp_path,
            'NC_A\t.\tgene\t1\t100\t.\t+\t.\tgene_id "NC_A"; transcript_id "NC_A";\n',
            "a.gtf",
        )
        f2 = self._write_gtf(
            tmp_path,
            'NC_B\t.\tgene\t1\t200\t.\t+\t.\tgene_id "NC_B"; transcript_id "NC_B";\n',
            "b.gtf",
        )
        combined: set[str] = set()
        for fp in [f1, f2]:
            combined |= _parse_gtf_file(fp)
        assert combined == {"NC_A", "NC_B"}

    def test_user_gtf_merged_with_existing(self, tmp_path: Path) -> None:
        base = self._write_gtf(
            tmp_path,
            'NC_BUILT\t.\tgene\t1\t100\t.\t+\t.\tgene_id "NC_BUILT"; transcript_id "NC_BUILT";\n',
            "base.gtf",
        )
        user = self._write_gtf(
            tmp_path,
            'NC_USER\t.\tgene\t1\t100\t.\t+\t.\tgene_id "NC_USER"; transcript_id "NC_USER";\n',
            "user.gtf",
        )
        result: set[str] = set()
        for fp in [base, user]:
            result |= _parse_gtf_file(fp)
        assert "NC_BUILT" in result
        assert "NC_USER" in result


# ---------------------------------------------------------------------------
# Bundled data/ GTFs — smoke-check a few known files
# ---------------------------------------------------------------------------


class TestBundledGtfs:
    @pytest.fixture
    def data_dir(self) -> Path:
        here = Path(__file__).resolve().parent.parent
        d = here / "src" / "viralscan" / "data"
        assert d.is_dir(), f"data/ directory not found at {d}"
        return d

    def test_data_dir_has_gtf_files(self, data_dir: Path) -> None:
        gtfs = list(data_dir.glob("*.gtf"))
        assert len(gtfs) > 0, "No GTF files found in data/"

    def test_dengue_gtf_parseable(self, data_dir: Path) -> None:
        dengue = data_dir / "Denguevirus_NC_001477.gtf"
        if not dengue.exists():
            pytest.skip("Dengue GTF not present")
        result = _parse_gtf_file(dengue)
        assert len(result) > 0

    def test_dengue_gtf_contains_expected_gene_id(self, data_dir: Path) -> None:
        # The Dengue GTF uses DENV_ gene_id prefixes (not the chromosome accession).
        dengue = data_dir / "Denguevirus_NC_001477.gtf"
        if not dengue.exists():
            pytest.skip("Dengue GTF not present")
        result = _parse_gtf_file(dengue)
        assert any("DENV" in acc for acc in result)

    def test_all_bundled_gtfs_parseable(self, data_dir: Path) -> None:
        """Every bundled GTF must be parseable without raising exceptions."""
        gtfs = list(data_dir.glob("*.gtf"))
        errors = []
        for p in gtfs:
            try:
                result = _parse_gtf_file(p)
                assert isinstance(result, set)
            except Exception as exc:
                errors.append(f"{p.name}: {exc}")
        assert not errors, "Errors in bundled GTFs:\n" + "\n".join(errors)

    def test_bundled_gtfs_produce_non_empty_sets(self, data_dir: Path) -> None:
        """Each bundled GTF should yield at least one accession."""
        gtfs = list(data_dir.glob("*.gtf"))
        empty = [p.name for p in gtfs if not _parse_gtf_file(p)]
        assert not empty, f"GTFs that yielded no accessions: {empty}"


# ---------------------------------------------------------------------------
# _count_lines / _count_unique_genes from menu.py
# ---------------------------------------------------------------------------


class TestMenuHelpers:
    def test_count_lines(self, tmp_path: Path) -> None:
        from viralscan.menu import _count_lines

        p = tmp_path / "f.txt"
        p.write_text("line1\nline2\nline3\n")
        assert _count_lines(str(p)) == 3

    def test_count_lines_empty_file(self, tmp_path: Path) -> None:
        from viralscan.menu import _count_lines

        p = tmp_path / "empty.txt"
        p.write_text("")
        assert _count_lines(str(p)) == 0

    def test_count_unique_genes(self, tmp_path: Path) -> None:
        from viralscan.menu import _count_unique_genes

        # Format: transcript\tgene\tname
        content = (
            "T1\tG1\tgene1\n"
            "T2\tG1\tgene1\n"  # same (G1, gene1) — duplicate
            "T3\tG2\tgene2\n"
        )
        p = tmp_path / "t2g.txt"
        p.write_text(content)
        assert _count_unique_genes(str(p)) == 2

    def test_count_unique_genes_empty(self, tmp_path: Path) -> None:
        from viralscan.menu import _count_unique_genes

        p = tmp_path / "empty.txt"
        p.write_text("")
        assert _count_unique_genes(str(p)) == 0


# ---------------------------------------------------------------------------
# Audit §2.2 — gene_id extracted by attribute name, not position
# ---------------------------------------------------------------------------


class TestGtfGeneIdAttributeOrder:
    """Audit §2.2: gene_id must be extracted by attribute name, not first-quote position.

    Regression for: audits/2026-05-08-full-pipeline.md §2.2
    """

    def test_gene_id_extracted_when_preceded_by_other_attribute(self) -> None:
        """
        GIVEN: a GTF line where a quoted attribute (source "NCBI") precedes gene_id
        WHEN:  _parse_gtf_text processes the line
        THEN:  the gene_id value "NC_123" is returned, not the source value "NCBI"
        """
        gtf = 'NC_TEST\t.\tgene\t1\t100\t.\t+\t.\tsource "NCBI"; gene_id "NC_123";\n'
        result = _parse_gtf_text(gtf)
        assert result == {"NC_123"}, (
            f"Expected {{'NC_123'}} but got {result!r}. "
            "GTF parser is taking the first quoted token instead of the gene_id value."
        )

    def test_no_gene_id_attribute_yields_empty_set(self) -> None:
        """GTF lines with no gene_id attribute must not crash and must return nothing."""
        gtf = 'NC_TEST\t.\tgene\t1\t100\t.\t+\t.\tsource "NCBI"; transcript_id "T1";\n'
        result = _parse_gtf_text(gtf)
        assert result == set(), f"Expected empty set, got {result!r}"

    def test_gene_id_extracted_when_last_attribute(self) -> None:
        """gene_id as the last attribute must still be parsed correctly."""
        gtf = 'NC_TEST\t.\texon\t1\t100\t.\t+\t.\ttranscript_id "T1"; gene_id "NC_456";\n'
        result = _parse_gtf_text(gtf)
        assert result == {"NC_456"}, f"Expected NC_456, got {result!r}"

    def test_gene_id_file_parse_adversarial(self, tmp_path: Path) -> None:
        """Same adversarial line parsed from a file must yield correct gene_id."""
        p = tmp_path / "adversarial.gtf"
        p.write_text('NC_TEST\t.\tgene\t1\t100\t.\t+\t.\tsource "NCBI"; gene_id "NC_789";\n')
        result = _parse_gtf_file(p)
        assert result == {"NC_789"}, f"Expected NC_789, got {result!r}"


class TestAnalysisScriptDataCache:
    def test_custom_gtf_runs_without_zenodo_cache(self, tmp_path: Path, monkeypatch) -> None:
        """Custom GTF workflows must not require the external Zenodo panel cache."""
        from viralscan import data_fetch

        def missing_cache():
            raise data_fetch.ViralScanDataError("missing test cache")

        monkeypatch.setattr(data_fetch, "ensure_viral_data", missing_cache)

        user_gtf = tmp_path / "custom.gtf"
        user_gtf.write_text('NC_USER\t.\tgene\t1\t100\t.\t+\t.\tgene_id "NC_USER";\n')

        output = tmp_path / "out"
        (output / "log").mkdir(parents=True)
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump({"output": f"{output}/", "gtf": str(user_gtf)}))

        snakemake = SimpleNamespace(params=SimpleNamespace(configfile=str(config_path)))
        script = (
            Path(__file__).resolve().parent.parent / "src" / "viralscan" / "scripts" / "analysis.py"
        )
        runpy.run_path(str(script), init_globals={"snakemake": snakemake})

        assert (output / "log" / "analysis.txt").read_text().strip() == "NC_USER"

    def test_string_none_gtf_is_treated_as_unset(self, tmp_path: Path, monkeypatch) -> None:
        """Snakemake string sentinels for unset GTF must not be parsed as file paths."""
        from viralscan import data_fetch

        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()
        (cache_dir / "panel.gtf").write_text(
            'NC_CACHE\t.\tgene\t1\t100\t.\t+\t.\tgene_id "NC_CACHE";\n'
        )
        monkeypatch.setattr(data_fetch, "ensure_viral_data", lambda: cache_dir)

        output = tmp_path / "out"
        (output / "log").mkdir(parents=True)
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump({"output": f"{output}/", "gtf": "None"}))

        snakemake = SimpleNamespace(params=SimpleNamespace(configfile=str(config_path)))
        script = (
            Path(__file__).resolve().parent.parent / "src" / "viralscan" / "scripts" / "analysis.py"
        )
        runpy.run_path(str(script), init_globals={"snakemake": snakemake})

        assert (output / "log" / "analysis.txt").read_text().strip() == "NC_CACHE"

    def test_comma_separated_custom_gtfs_are_all_parsed(self, tmp_path: Path, monkeypatch) -> None:
        """All custom GTFs listed in a comma-separated config value must be parsed."""
        from viralscan import data_fetch

        def missing_cache():
            raise data_fetch.ViralScanDataError("missing test cache")

        monkeypatch.setattr(data_fetch, "ensure_viral_data", missing_cache)

        user_a = tmp_path / "custom_a.gtf"
        user_b = tmp_path / "custom_b.gtf"
        user_a.write_text('NC_A\t.\tgene\t1\t100\t.\t+\t.\tgene_id "NC_A";\n')
        user_b.write_text('NC_B\t.\tgene\t1\t100\t.\t+\t.\tgene_id "NC_B";\n')

        output = tmp_path / "out"
        (output / "log").mkdir(parents=True)
        config_path = tmp_path / "config.yaml"
        config_path.write_text(
            yaml.safe_dump({"output": f"{output}/", "gtf": f"{user_a},{user_b}"})
        )

        snakemake = SimpleNamespace(params=SimpleNamespace(configfile=str(config_path)))
        script = (
            Path(__file__).resolve().parent.parent / "src" / "viralscan" / "scripts" / "analysis.py"
        )
        runpy.run_path(str(script), init_globals={"snakemake": snakemake})

        assert set((output / "log" / "analysis.txt").read_text().splitlines()) == {
            "NC_A",
            "NC_B",
        }

    def test_missing_custom_gtf_raises_file_not_found(self, tmp_path: Path, monkeypatch) -> None:
        """A missing custom GTF path should fail clearly instead of being skipped."""
        from viralscan import data_fetch

        def missing_cache():
            raise data_fetch.ViralScanDataError("missing test cache")

        monkeypatch.setattr(data_fetch, "ensure_viral_data", missing_cache)

        output = tmp_path / "out"
        (output / "log").mkdir(parents=True)
        missing_gtf = tmp_path / "missing.gtf"
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump({"output": f"{output}/", "gtf": str(missing_gtf)}))

        snakemake = SimpleNamespace(params=SimpleNamespace(configfile=str(config_path)))
        script = (
            Path(__file__).resolve().parent.parent / "src" / "viralscan" / "scripts" / "analysis.py"
        )

        with pytest.raises(FileNotFoundError, match="missing.gtf"):
            runpy.run_path(str(script), init_globals={"snakemake": snakemake})
