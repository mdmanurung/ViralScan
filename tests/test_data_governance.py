from pathlib import Path

from scripts.check_data_governance import violations

ROOT = Path(__file__).resolve().parents[1]


def test_rejects_restricted_artifacts_and_ship_scope_paths(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "src").mkdir()
    (tmp_path / "src" / "module.py").write_text('PATH = "/exports/private/data"')
    problems = violations([Path("patient.fastq.gz"), Path("src/module.py")])
    assert any("restricted biological artifact" in item for item in problems)
    assert any("institutional absolute path" in item for item in problems)


def test_allows_synthetic_fixture_names_and_nonshipping_docs(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "fixture.tsv").write_text("synthetic")
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "plan.md").write_text("/exports/example")
    assert violations([Path("tests/fixture.tsv"), Path("docs/plan.md")]) == []


def test_allows_only_checksum_pinned_synthetic_fastq_fixtures(monkeypatch) -> None:
    monkeypatch.chdir(ROOT)
    fixtures = [
        Path("tests/data/evidence_tiny/R1.fastq"),
        Path("tests/data/evidence_tiny/R2.fastq"),
    ]
    assert violations(fixtures) == []


def test_rejects_tampered_synthetic_fastq_fixture(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    fixture = tmp_path / "tests/data/evidence_tiny/R1.fastq"
    fixture.parent.mkdir(parents=True)
    fixture.write_text("not the pinned synthetic fixture\n")
    problems = violations([Path("tests/data/evidence_tiny/R1.fastq")])
    assert problems == [
        "restricted biological artifact is tracked: tests/data/evidence_tiny/R1.fastq"
    ]
