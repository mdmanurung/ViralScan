"""CAT-14 (host cross-talk), CAT-15 (programme markers) and CAT-16 (scale check) gates."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from viralscan import reference_gates as gates
from viralscan.gene_programs import load_catalogue


def _homology_table(path: Path, rows: list[dict]) -> Path:
    cols = [
        "accession_version",
        "host_homology_status",
        "host_homology_max_identity",
        "host_homology_max_query_coverage",
        "host_homology_max_aligned_bases",
        "host_homology_best_target",
        "host_homology_loci",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, cols, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in cols})
    return path


def _loci(*spans: tuple[int, int]) -> str:
    return json.dumps([{"viral_start": a, "viral_end": b} for a, b in spans])


def test_host_fraction_is_the_union_of_loci(tmp_path: Path) -> None:
    table = _homology_table(
        tmp_path / "t.tsv",
        [
            {
                "accession_version": "A.1",
                "host_homology_status": "measured",
                "host_homology_loci": _loci((0, 100), (50, 150), (400, 500)),
            },
            {
                "accession_version": "B.1",
                "host_homology_status": "measured",
                "host_homology_loci": "[]",
            },
        ],
    )
    report, problems = gates.host_crosstalk_report(table, {"A.1": 1000, "B.1": 500})
    assert problems == []
    by = {r["accession_version"]: r for r in report}
    assert by["A.1"]["host_homologous_bases"] == 250  # [0,150) + [400,500)
    assert by["A.1"]["host_homologous_fraction"] == 0.25
    assert by["B.1"]["host_homologous_fraction"] == 0.0


def test_unmeasured_or_missing_accessions_fail_and_are_never_zero(tmp_path: Path) -> None:
    table = _homology_table(
        tmp_path / "t.tsv",
        [{"accession_version": "A.1", "host_homology_status": "not_measured"}],
    )
    report, problems = gates.host_crosstalk_report(table, {"A.1": 100, "MISSING.1": 100})
    assert len(problems) == 2
    assert report[0]["host_homologous_fraction"] == ""  # not 0


def test_cli_writes_report_and_exit_code(tmp_path: Path) -> None:
    table = _homology_table(
        tmp_path / "t.tsv",
        [
            {
                "accession_version": "A.1",
                "host_homology_status": "measured",
                "host_homology_loci": "[]",
            }
        ],
    )
    fasta = tmp_path / "v.fa"
    fasta.write_text(">A.1 x\nACGT\nACGT\n")
    out = tmp_path / "report.tsv"
    assert (
        gates.main(
            [
                "host-crosstalk",
                "--table",
                str(table),
                "--viral-fasta",
                str(fasta),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    row = next(csv.DictReader(out.open(), delimiter="\t"))
    assert row["genome_length"] == "8"
    fasta.write_text(">A.1\nACGT\n>B.1\nACGT\n")
    assert (
        gates.main(
            [
                "host-crosstalk",
                "--table",
                str(table),
                "--viral-fasta",
                str(fasta),
                "--out",
                str(out),
            ]
        )
        == 1
    )


def _summary(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    path.write_text(
        "virus_name\tdetected\tinfected_called\n" + "".join("\t".join(r) + "\n" for r in rows)
    )
    return path


def test_host_only_negative_reports_no_call(tmp_path: Path) -> None:
    clean = _summary(tmp_path / "a.tsv", [("EBV", "False", "0"), ("HSV-1", "False", "0")])
    assert gates.host_only_calls(clean) == []
    assert gates.main(["host-only", str(clean)]) == 0
    dirty = _summary(tmp_path / "b.tsv", [("EBV", "False", "0"), ("HHV-6B", "True", "3")])
    assert gates.host_only_calls(dirty) == ["HHV-6B"]
    assert gates.main(["host-only", str(dirty)]) == 1
    assert gates.host_only_calls(_summary(tmp_path / "c.tsv", [("X", "False", "2")])) == ["X"]


# --- CAT-15 ---


def _all_marker_ids(catalogue: list[dict[str, str]]) -> set[str]:
    return {r["gene_id_bundled"] for r in catalogue}


def test_programmes_pass_when_every_marker_is_in_the_index() -> None:
    cat = load_catalogue()
    report, problems = gates.check_programme_markers(_all_marker_ids(cat) | {"ENSG1"})
    assert problems == []
    assert {r["virus"] for r in report} == {r["virus"] for r in cat}


def test_programmes_fail_when_a_marker_does_not_resolve() -> None:
    cat = load_catalogue()
    ids = _all_marker_ids(cat)
    victim = next(r for r in cat if r["virus"] == "Epstein-Barr virus")
    ids.discard(victim["gene_id_bundled"])
    ids.discard(victim["gene_id_merged"])
    _, problems = gates.check_programme_markers(ids)
    assert len(problems) == 1 and "Epstein-Barr virus" in problems[0]
    # an index without the virus at all fails too, and an unknown virus name is flagged
    _, none = gates.check_programme_markers({"ENSG1"}, ["Human herpesvirus 8"])
    assert none and "Human herpesvirus 8" in none[0]
    assert gates.check_programme_markers(ids, ["No such virus"])[1]


def test_programmes_cli_reads_a_t2g(tmp_path: Path) -> None:
    cat = load_catalogue()
    t2g = tmp_path / "t.t2g"
    t2g.write_text("".join(f"tx{i}\t{g}\n" for i, g in enumerate(sorted(_all_marker_ids(cat)))))
    assert gates.main(["programmes", "--t2g", str(t2g)]) == 0
    t2g.write_text("tx0\tENSG1\n")
    assert gates.main(["programmes", "--t2g", str(t2g)]) == 1


# --- CAT-16 ---


def _gtf(path: Path, genes: dict[str, int]) -> Path:
    lines = []
    for gene, length in genes.items():
        lines.append(
            f'{gene}_contig\tt\texon\t1\t{length}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{gene}_t";\n'
        )
    path.write_text("".join(lines))
    return path


def test_scale_report_growth_memory_and_budgets(tmp_path: Path) -> None:
    old = _gtf(tmp_path / "old.gtf", {"A": 100})
    new = _gtf(tmp_path / "new.gtf", {f"G{i}": 100 for i in range(120)})
    report, problems = gates.scale_report(new, old, n_cells=1000)
    assert problems == []
    assert report["n_genes"] == 120 and report["gene_growth_factor"] == 120.0
    assert report["warnings"]  # beyond the row's ~10-50x
    assert report["dense_gib_per_copy"] == round(
        1000 * report["max_genes_per_virus"] * 8 / 2**30, 3
    )
    _, over = gates.scale_report(new, old, n_cells=10**9, max_genes=100, max_dense_gib=0.001)
    assert len(over) == 2


def test_scale_cli(tmp_path: Path) -> None:
    gtf = _gtf(tmp_path / "g.gtf", {"A": 10})
    out = tmp_path / "scale.json"
    assert gates.main(["scale", "--gtf", str(gtf), "--out", str(out)]) == 0
    assert json.loads(out.read_text())["n_genes"] == 1
    assert gates.main(["scale", "--gtf", str(gtf), "--max-genes", "0"]) == 1
