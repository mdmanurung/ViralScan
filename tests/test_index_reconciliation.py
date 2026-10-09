"""Tests for the catalogue↔index reconciliation guard (PLAN `CAT-31`, F-015).

A catalogued virus that is not indexed is *undetectable*: no read threshold
recovers a read with no k-mer in the index. Before this guard, nothing compared
``virus_catalog.tsv`` (2,249 accessions) against what the build actually emits,
so 13 of 16 catalogued HPV and all 5 retroviruses were invisible for as long as
the panel existed. The index content is decided by the 195 bundled GTFs plus a
separate anellovirus fetch path, and neither ever looked at the catalogue.

The guard is a build-time pre-flight in ``scripts/build_bundled_panel_ref.py``;
the logic under test here lives in ``viralscan.scripts.build_reference`` so it
can be exercised without NCBI, Ensembl or ``kb``.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

from viralscan.scripts.build_reference import (
    INDEX_EXCLUSIONS_NAME,
    RECONCILIATION_REPORT_COLUMNS,
    RECONCILIATION_REPORT_NAME,
    catalogue_detection_targets,
    format_reconciliation_summary,
    load_index_exclusions,
    normalise_accession,
    reconcile_catalogue_against_panel,
    reconcile_reference_panel,
    reconciliation_failure,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "build_bundled_panel_ref.py"
SPEC = importlib.util.spec_from_file_location("build_bundled_panel_ref", SCRIPT_PATH)
assert SPEC is not None
assert SPEC.loader is not None
build_bundled_panel_ref = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(build_bundled_panel_ref)


def _write_catalogue(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    """Write a minimal catalogue: accession, species, family."""
    text = "accession\taccession_version\tspecies\tfamily\n"
    text += "".join(f"{a}\t{a}.1\t{s}\t{f}\n" for a, s, f in rows)
    path.write_text(text)
    return path


def _write_exclusions(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    """Write an index-exclusion allowlist, with a comment above the header."""
    text = "# deliberate omissions\naccession\treason\tdecided_by\n"
    text += "".join(f"{a}\t{r}\t{d}\n" for a, r, d in rows)
    path.write_text(text)
    return path


def _write_panel(path: Path, headers: list[str]) -> Path:
    """Write a FASTA whose record identifiers are the accessions the panel holds."""
    body = "".join(f">{h}\n{'ACGTGGTACCTGATCGTAG' * 3}\n" for h in headers)
    path.write_text(body)
    return path


def _rows(report: Path) -> list[dict[str, str]]:
    with report.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ---------------------------------------------------------------------------
# Accession normalisation — the step that produced a wrong count in F-015
# ---------------------------------------------------------------------------


class TestNormaliseAccession:
    def test_space_and_underscore_and_version_all_reconcile(self):
        # VIRTUS2's list file spells this "NC 000883.2" while its FASTA uses
        # "NC_000883.2", and the catalogue holds the bare base accession.
        assert normalise_accession("NC 000883.2") == "NC_000883"
        assert normalise_accession("NC_000883.2") == "NC_000883"
        assert normalise_accession("NC_000883") == "NC_000883"

    def test_version_bump_is_not_a_miss(self):
        # A catalogue pinned to .1 and a panel built from .2 are one genome.
        assert normalise_accession("NC_009334.1") == normalise_accession("NC_009334.2")

    def test_bare_unversioned_accessions_survive(self):
        # GenBank-only records have no version and no underscore to fold.
        assert normalise_accession("M74117") == "M74117"
        assert normalise_accession("U31794") == "U31794"

    def test_case_and_surrounding_whitespace_are_ignored(self):
        assert normalise_accession("  nc_001526.4 \n") == "NC_001526"

    def test_fasta_header_space_spelling_is_handled(self):
        # "NC 009334.1" as a whole header token still folds to the base accession.
        assert normalise_accession("NC 009334.1") == normalise_accession("NC_009334.1")


# ---------------------------------------------------------------------------
# Catalogue and allowlist loading
# ---------------------------------------------------------------------------


class TestCatalogueDetectionTargets:
    def test_every_row_is_a_detection_target(self, tmp_path):
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv",
            [("NC_001526", "Human papillomavirus 16", "Papillomaviridae")],
        )
        targets = catalogue_detection_targets(catalogue)
        assert set(targets) == {"NC_001526"}
        assert targets["NC_001526"]["family"] == "Papillomaviridae"
        assert targets["NC_001526"]["species"] == "Human papillomavirus 16"

    def test_absent_catalogue_raises_a_clear_error(self, tmp_path):
        with pytest.raises(ValueError, match="does not exist"):
            catalogue_detection_targets(tmp_path / "absent.tsv")

    def test_malformed_catalogue_without_accession_column_raises(self, tmp_path):
        bad = tmp_path / "catalog.tsv"
        bad.write_text("name\tspecies\nHPV16\tHuman papillomavirus 16\n")
        with pytest.raises(ValueError, match="no 'accession' column"):
            catalogue_detection_targets(bad)

    def test_blank_accession_is_malformed_not_an_exclusion(self, tmp_path):
        bad = tmp_path / "catalog.tsv"
        bad.write_text(
            "accession\tspecies\tfamily\nNC_001526\tHPV16\tPapillomaviridae\n\tHPV18\tPapillomaviridae\n"
        )
        with pytest.raises(ValueError, match="malformed catalogue"):
            catalogue_detection_targets(bad)

    def test_header_only_catalogue_raises(self, tmp_path):
        empty = tmp_path / "catalog.tsv"
        empty.write_text("accession\tspecies\tfamily\n")
        with pytest.raises(ValueError, match="no catalogue rows"):
            catalogue_detection_targets(empty)

    def test_conflicting_duplicate_rows_raise(self, tmp_path):
        bad = tmp_path / "catalog.tsv"
        bad.write_text(
            "accession\tspecies\tfamily\n"
            "NC_001526\tHuman papillomavirus 16\tPapillomaviridae\n"
            "NC_001526.1\tHuman papillomavirus 16\tPapilloviridae\n"
        )
        with pytest.raises(ValueError, match="conflicting family"):
            catalogue_detection_targets(bad)

    def test_identical_duplicate_rows_collapse(self, tmp_path):
        dup = tmp_path / "catalog.tsv"
        row = "NC_001526\tHuman papillomavirus 16\tPapillomaviridae\n"
        dup.write_text("accession\tspecies\tfamily\n" + row + row)
        assert set(catalogue_detection_targets(dup)) == {"NC_001526"}

    def test_panel_scope_limits_targets_to_the_shipped_panel(self, tmp_path):
        """MECH-A: max-panel rows are catalogued for naming, not claimed as shipped."""
        scoped = tmp_path / "catalog.tsv"
        scoped.write_text(
            "accession\tspecies\tfamily\tpanel\n"
            "NC_001526\tHuman papillomavirus 16\tPapillomaviridae\tshipped\n"
            "NC_009334\tHuman herpesvirus 4 type 2\tOrthoherpesviridae\tmax\n"
            "NC_001357\thuman papillomavirus 18\tPapillomaviridae\t\n"
        )
        assert set(catalogue_detection_targets(scoped)) == {"NC_001526", "NC_001357"}
        assert set(catalogue_detection_targets(scoped, panel="max")) == {
            "NC_009334",
            "NC_001357",
        }
        assert len(catalogue_detection_targets(scoped, panel=None)) == 3


class TestLoadIndexExclusions:
    def test_shipped_allowlist_is_present_and_parseable(self):
        # The guard is useless if the file it reads does not ship.
        assert (REPO_ROOT / "src" / "viralscan" / "data" / INDEX_EXCLUSIONS_NAME).is_file()
        assert isinstance(load_index_exclusions(), dict)

    def test_comments_above_the_header_are_ignored(self, tmp_path):
        path = _write_exclusions(
            tmp_path / INDEX_EXCLUSIONS_NAME,
            [("NC_001802", "segment-only record", "M. Manurung 2026-09-27")],
        )
        exclusions = load_index_exclusions(path)
        assert set(exclusions) == {"NC_001802"}
        assert exclusions["NC_001802"]["decided_by"] == "M. Manurung 2026-09-27"

    def test_versioned_exclusion_key_normalises(self, tmp_path):
        path = _write_exclusions(
            tmp_path / INDEX_EXCLUSIONS_NAME,
            [("NC_001802.1", "no CDS", "M. Manurung 2026-09-27")],
        )
        assert set(load_index_exclusions(path)) == {"NC_001802"}

    def test_absent_allowlist_raises_a_clear_error(self, tmp_path):
        with pytest.raises(ValueError, match="Index-exclusion allowlist not found"):
            load_index_exclusions(tmp_path / INDEX_EXCLUSIONS_NAME)

    def test_missing_columns_raise(self, tmp_path):
        bad = tmp_path / INDEX_EXCLUSIONS_NAME
        bad.write_text("accession\treason\nNC_001802\tno CDS\n")
        with pytest.raises(ValueError, match=r"decided_by"):
            load_index_exclusions(bad)

    def test_exclusion_without_a_decider_raises(self, tmp_path):
        # An exclusion with nobody's name against it is the bug wearing a label.
        bad = tmp_path / INDEX_EXCLUSIONS_NAME
        bad.write_text("accession\treason\tdecided_by\nNC_001802\tno CDS\t\n")
        with pytest.raises(ValueError, match="without a decided_by"):
            load_index_exclusions(bad)

    def test_exclusion_without_a_reason_raises(self, tmp_path):
        bad = tmp_path / INDEX_EXCLUSIONS_NAME
        bad.write_text("accession\treason\tdecided_by\nNC_001802\t\tM. Manurung 2026-09-27\n")
        with pytest.raises(ValueError, match="without a reason"):
            load_index_exclusions(bad)


# ---------------------------------------------------------------------------
# The reconciliation itself
# ---------------------------------------------------------------------------


class TestReconciliation:
    def test_clean_panel_produces_an_empty_report_and_no_failure(self, tmp_path):
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv",
            [
                ("NC_001526", "Human papillomavirus 16", "Papillomaviridae"),
                ("NC_007605", "Epstein-Barr virus", "Orthoherpesviridae"),
            ],
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4", "NC_007605.1"])
        report = tmp_path / RECONCILIATION_REPORT_NAME

        result = reconcile_reference_panel(
            panel,
            catalogue=catalogue,
            exclusions=_write_exclusions(tmp_path / "e.tsv", []),
            report=report,
        )

        assert result.rows == []
        assert result.unexplained == []
        assert result.intentional == []
        assert reconciliation_failure(result, strict=True, report=report) is None
        # An empty report is a header alone, not a missing file.
        assert report.read_text() == "\t".join(RECONCILIATION_REPORT_COLUMNS) + "\n"

    def test_unexplained_miss_is_reported_with_a_reason(self, tmp_path):
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv",
            [("NC_001802", "Human immunodeficiency virus 1", "Retroviridae")],
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])
        report = tmp_path / RECONCILIATION_REPORT_NAME

        result = reconcile_reference_panel(
            panel,
            catalogue=catalogue,
            exclusions=_write_exclusions(tmp_path / "e.tsv", []),
            report=report,
        )

        assert result.unexplained == ["NC_001802"]
        assert result.intentional == []
        (row,) = _rows(report)
        assert row["accession"] == "NC_001802"
        assert row["family"] == "Retroviridae"
        assert row["species"] == "Human immunodeficiency virus 1"
        assert row["status"] == "unexplained"
        assert INDEX_EXCLUSIONS_NAME in row["reason"]

    def test_allowlisted_miss_is_reported_as_intentional(self, tmp_path):
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv",
            [("NC_001802", "Human immunodeficiency virus 1", "Retroviridae")],
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])
        report = tmp_path / RECONCILIATION_REPORT_NAME

        result = reconcile_reference_panel(
            panel,
            catalogue=catalogue,
            exclusions=_write_exclusions(
                tmp_path / "e.tsv",
                [
                    (
                        "NC_001802",
                        "deliberately excluded: retroviral reads smear",
                        "M. M. 2026-09-27",
                    )
                ],
            ),
            report=report,
        )

        assert result.intentional == ["NC_001802"]
        assert result.unexplained == []
        # A live exclusion is not also a stale one.
        assert result.stale_exclusions == []
        (row,) = _rows(report)
        assert row["status"] == "intentional"
        assert row["reason"] == "deliberately excluded: retroviral reads smear"
        assert reconciliation_failure(result, strict=True, report=report) is None

    def test_strict_flag_is_the_only_thing_that_turns_a_miss_into_a_failure(self, tmp_path):
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv", [("NC_001802", "HIV-1", "Retroviridae")]
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])
        report = tmp_path / RECONCILIATION_REPORT_NAME
        result = reconcile_reference_panel(
            panel,
            catalogue=catalogue,
            exclusions=_write_exclusions(tmp_path / "e.tsv", []),
            report=report,
        )

        assert reconciliation_failure(result, strict=False, report=report) is None
        failure = reconciliation_failure(result, strict=True, report=report)
        assert failure is not None
        assert "NC_001802" in failure
        assert str(report) in failure

    def test_report_sorts_by_family_then_accession(self, tmp_path):
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv",
            [
                ("NC_001802", "HIV-1", "Retroviridae"),
                ("U31794", "HPV66", "Papillomaviridae"),
                ("M74117", "HPV35", "Papillomaviridae"),
            ],
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])
        report = tmp_path / RECONCILIATION_REPORT_NAME
        reconcile_reference_panel(
            panel,
            catalogue=catalogue,
            exclusions=_write_exclusions(tmp_path / "e.tsv", []),
            report=report,
        )

        assert [row["accession"] for row in _rows(report)] == ["M74117", "U31794", "NC_001802"]

    def test_version_and_spelling_drift_does_not_read_as_a_miss(self, tmp_path):
        # Catalogue says "NC 009334.1", panel says "NC_009334.2", allowlist would
        # say the versioned form. All three are EBV type 2; none is a miss.
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv",
            [("NC 009334", "Epstein-Barr virus type 2", "Orthoherpesviridae")],
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_009334.2"])
        report = tmp_path / RECONCILIATION_REPORT_NAME
        result = reconcile_reference_panel(
            panel,
            catalogue=catalogue,
            exclusions=_write_exclusions(tmp_path / "e.tsv", []),
            report=report,
        )
        assert result.rows == []

    def test_stale_exclusion_is_reported(self, tmp_path):
        # The panel indexes it anyway, so the allowlist is asserting something untrue.
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv", [("NC_001526", "HPV16", "Papillomaviridae")]
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])
        result = reconcile_reference_panel(
            panel,
            catalogue=catalogue,
            exclusions=_write_exclusions(
                tmp_path / "e.tsv", [("NC_001526", "once excluded", "M. M. 2026-09-27")]
            ),
        )
        assert result.stale_exclusions == ["NC_001526"]

    def test_exclusion_of_an_uncatalogued_accession_is_stale(self, tmp_path):
        # The allowlist names something the catalogue does not have at all.
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv", [("NC_001526", "HPV16", "Papillomaviridae")]
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])
        result = reconcile_reference_panel(
            panel,
            catalogue=catalogue,
            exclusions=_write_exclusions(
                tmp_path / "e.tsv",
                [("NC_009334", "EBV type 2, since dropped from the catalogue", "M. M. 2026-09-27")],
            ),
        )
        assert result.stale_exclusions == ["NC_009334"]
        # A stale entry must not also be counted as a live intentional miss.
        assert result.intentional == []

    def test_uncatalogued_emitted_records_are_surfaced(self, tmp_path):
        catalogue = _write_catalogue(
            tmp_path / "catalog.tsv", [("NC_001526", "HPV16", "Papillomaviridae")]
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4", "AF157706"])
        result = reconcile_reference_panel(
            panel, catalogue=catalogue, exclusions=_write_exclusions(tmp_path / "e.tsv", [])
        )
        assert result.uncatalogued == ["AF157706"]

    def test_absent_panel_fasta_raises_a_clear_error(self, tmp_path):
        with pytest.raises(ValueError, match="Panel FASTA to reconcile does not exist"):
            reconcile_reference_panel(
                tmp_path / "viral.fa",
                catalogue=_write_catalogue(tmp_path / "c.tsv", [("NC_1", "x", "y")]),
                exclusions=_write_exclusions(tmp_path / "e.tsv", []),
            )

    def test_reconcile_catalogue_against_panel_accepts_raw_iterables(self):
        result = reconcile_catalogue_against_panel(
            iter(["NC_001526.4", "NC_007605.1"]),
            {"NC_001526": {"family": "Papillomaviridae", "species": "HPV16"}},
            {},
        )
        assert result.rows == []
        assert result.uncatalogued == ["NC_007605"]


# ---------------------------------------------------------------------------
# The operator-facing report text
# ---------------------------------------------------------------------------


class TestSummaryFormatting:
    def test_clean_panel_summary_is_quiet(self, tmp_path):
        result = reconcile_catalogue_against_panel(["NC_001526.4"], {"NC_001526": {}}, {})
        text = format_reconciliation_summary(result, tmp_path / RECONCILIATION_REPORT_NAME)
        assert "every catalogued accession is in the panel" in text
        assert "UNEXPLAINED" not in text

    def test_unexplained_summary_is_loud_and_names_the_accessions(self, tmp_path):
        result = reconcile_catalogue_against_panel(
            ["NC_001526.4"],
            {
                "NC_001802": {
                    "family": "Retroviridae",
                    "species": "Human immunodeficiency virus 1",
                }
            },
            {},
        )
        text = format_reconciliation_summary(result, tmp_path / RECONCILIATION_REPORT_NAME)
        assert "UNEXPLAINED" in text
        assert "UNEXPLAINED : 1" in text
        assert "NC_001802" in text
        assert "Retroviridae" in text
        assert "UNDETECTABLE" in text
        assert "--strict-reconciliation" in text
        assert text.count("=" * 78) >= 2

    def test_summary_truncates_a_long_uncatalogued_list(self, tmp_path):
        emitted = [f"NC_{n:06d}" for n in range(40)]
        result = reconcile_catalogue_against_panel(emitted, {"NC_999999": {}}, {})
        text = format_reconciliation_summary(result, tmp_path / RECONCILIATION_REPORT_NAME)
        assert "(+28 more)" in text

    def test_stale_and_uncatalogued_are_both_flagged(self, tmp_path):
        result = reconcile_catalogue_against_panel(
            ["NC_001526.4", "AF157706"],
            {"NC_001526": {}},
            {"NC_001526": {"reason": "old", "decided_by": "M. M."}},
        )
        text = format_reconciliation_summary(result, tmp_path / RECONCILIATION_REPORT_NAME)
        assert "STALE" in text
        assert "AF157706" in text


# ---------------------------------------------------------------------------
# The build script's wiring — flag contract and exit behaviour
# ---------------------------------------------------------------------------


class TestBuildScriptWiring:
    def _args(
        self,
        tmp_path: Path,
        *extra: str,
        catalogue_rows: list[tuple[str, str, str]] | None = None,
    ):
        rows = catalogue_rows or [("NC_001802", "HIV-1", "Retroviridae")]
        argv = [
            "--out",
            str(tmp_path),
            "--ncbi-email",
            "build@example.org",
            "--catalogue",
            str(_write_catalogue(tmp_path / "catalog.tsv", rows)),
            "--index-exclusions",
            str(_write_exclusions(tmp_path / INDEX_EXCLUSIONS_NAME, [])),
            *extra,
        ]
        return build_bundled_panel_ref._build_arg_parser().parse_args(argv)

    def test_strict_reconciliation_is_opt_in(self, tmp_path):
        # Default must stay non-fatal: the catalogue is 11x the bundled GTF set,
        # so a default-on gate would make the build unusable rather than honest.
        assert self._args(tmp_path).strict_reconciliation is False
        assert self._args(tmp_path, "--strict-reconciliation").strict_reconciliation is True

    def test_clean_panel_exits_cleanly(self, tmp_path, capsys):
        args = self._args(
            tmp_path,
            "--strict-reconciliation",
            catalogue_rows=[("NC_001526", "HPV16", "Papillomaviridae")],
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])

        build_bundled_panel_ref._run_reconciliation(panel, args)

        assert (tmp_path / RECONCILIATION_REPORT_NAME).is_file()
        assert "every catalogued accession is in the panel" in capsys.readouterr().out

    def test_unexplained_miss_exits_nonzero_under_strict(self, tmp_path, capsys):
        args = self._args(tmp_path, "--strict-reconciliation")
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])

        with pytest.raises(SystemExit) as excinfo:
            build_bundled_panel_ref._run_reconciliation(panel, args)

        assert excinfo.value.code != 0
        assert "NC_001802" in str(excinfo.value.code)
        # The report and the loud banner are still produced, so the operator can
        # act on the failure rather than just read it.
        assert (tmp_path / RECONCILIATION_REPORT_NAME).is_file()
        assert "UNEXPLAINED" in capsys.readouterr().out

    def test_unexplained_miss_warns_but_succeeds_without_strict(self, tmp_path, capsys):
        args = self._args(tmp_path)
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])

        build_bundled_panel_ref._run_reconciliation(panel, args)

        assert "UNEXPLAINED : 1" in capsys.readouterr().out
        assert (tmp_path / RECONCILIATION_REPORT_NAME).is_file()

    def test_allowlisted_miss_survives_strict(self, tmp_path):
        args = self._args(tmp_path, "--strict-reconciliation")
        _write_exclusions(
            tmp_path / INDEX_EXCLUSIONS_NAME,
            [("NC_001802", "retroviral reads smear across the panel", "M. M. 2026-09-27")],
        )
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])

        build_bundled_panel_ref._run_reconciliation(panel, args)

        (row,) = _rows(tmp_path / RECONCILIATION_REPORT_NAME)
        assert row["status"] == "intentional"

    def test_absent_catalogue_exits_with_a_message_not_a_traceback(self, tmp_path):
        args = self._args(tmp_path)
        args.catalogue = tmp_path / "absent.tsv"
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])

        with pytest.raises(SystemExit) as excinfo:
            build_bundled_panel_ref._run_reconciliation(panel, args)

        assert "reconciliation could not run" in str(excinfo.value.code)
        assert "does not exist" in str(excinfo.value.code)

    def test_malformed_catalogue_exits_with_a_message_not_a_traceback(self, tmp_path):
        args = self._args(tmp_path)
        (tmp_path / "catalog.tsv").write_text("name\tspecies\nHPV16\tHuman papillomavirus 16\n")
        panel = _write_panel(tmp_path / "viral.fa", ["NC_001526.4"])

        with pytest.raises(SystemExit) as excinfo:
            build_bundled_panel_ref._run_reconciliation(panel, args)

        assert "no 'accession' column" in str(excinfo.value.code)


# ---------------------------------------------------------------------------
# The shipped panel really is short of the catalogue (F-015, reproduced)
# ---------------------------------------------------------------------------


class TestShippedCatalogueIsReconciled:
    """The guard is only worth shipping if it is currently red.

    If this test ever goes green on its own, someone added a genome and the
    allowlist, which is the intended path — the test then documents the new
    floor rather than the old defect.
    """

    def test_every_catalogued_anellovirus_is_reachable_but_breadth_is_not(self):
        from viralscan.anellovirus import load_accession_table
        from viralscan.virus_catalog import catalogue_path

        targets = catalogue_detection_targets(catalogue_path())

        # `emitted` is built from tracked inputs only: the anellovirus accession
        # table plus the frozen broad-discovery accession list. The bundled GTFs
        # are git-ignored (6bb5c64), so reading them made this fail on a clean tree.
        emitted = {normalise_accession(row["accession"]) for row in load_accession_table()}
        broad = catalogue_path().with_name("broad_discovery_accessions.tsv")
        with open(broad, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                emitted.add(normalise_accession(row["accession_version"]))

        missing = set(targets) - emitted
        # Rewritten 2026-09-28. This test used to assert that the Retroviridae and
        # Papillomaviridae gaps were *present* -- i.e. it encoded F-015's measured
        # sensitivity loss as the expected state. The CAT-31/32/33 additions closed
        # those gaps, so the assertion inverted: the curated breadth that was
        # missing must now be reachable.
        families = {targets[key]["family"] for key in missing}
        assert "Retroviridae" not in families, (
            "all five Retroviridae are indexed again; update this test deliberately"
        )
        assert "Papillomaviridae" not in families, (
            "catalogued HPV breadth is indexed again; update this test deliberately"
        )

        # The only remaining anellovirus misses are the CAT-05 duplicate pair:
        # NC_038359.1 and AB303562.1 are byte-identical, so at most one may be
        # indexed. Both are catalogued, so this asserts the pair rather than
        # requiring the family to be clean.
        anello_missing = {key for key in missing if targets[key]["family"] == "Anelloviridae"}
        assert anello_missing <= {"NC_038359", "AB303562"}, (
            "the only expected anellovirus reconciliation gap is the CAT-05 duplicate "
            f"pair; got {sorted(anello_missing)}"
        )

    def test_exclusions_file_parses_and_has_the_documented_columns(self):
        import importlib.resources

        path = Path(
            str(importlib.resources.files("viralscan.data").joinpath(INDEX_EXCLUSIONS_NAME))
        )
        text = path.read_text()
        header = next(line for line in text.splitlines() if not line.startswith("#"))
        assert header.split("\t") == ["accession", "reason", "decided_by"]
        # The only decision on file: the CAT-05 byte-identical duplicate of NC_038359.1.
        assert set(load_index_exclusions(path)) == {"AB303562"}


def test_reconciliation_report_filename_matches_plan():
    # PLAN.md `CAT-31` names the file; renaming it silently breaks the contract.
    assert RECONCILIATION_REPORT_NAME == "catalogued_not_indexed.tsv"
    plan = (REPO_ROOT / "PLAN.md").read_text()
    assert RECONCILIATION_REPORT_NAME in plan


def test_report_is_tab_separated_with_the_documented_columns():
    assert RECONCILIATION_REPORT_COLUMNS == ("accession", "family", "species", "status", "reason")
