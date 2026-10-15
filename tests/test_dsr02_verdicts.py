"""scripts/dsr02_verdicts.py: verdicts from synthetic evidence directories."""

import gzip
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "dsr02_verdicts", Path(__file__).parent.parent / "scripts" / "dsr02_verdicts.py"
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)

QC = "count_layer\treference\treference_class\treads\tcomplex_body_fraction\n"


def _ev(tmp_path, name, qc_rows, flags=None, reads=("ACGT",)):
    d = tmp_path / name
    d.mkdir()
    (d / "alignment_qc.tsv").write_text(
        QC + "".join(f"deduplicated\t{ref}\t{cls}\t{n}\t{c}\n" for ref, cls, n, c in qc_rows)
    )
    flags = {"host_homology": "not_flagged", "low_complexity": "not_flagged", **(flags or {})}
    (d / "interpretation_flags.tsv").write_text(
        "flag\tstatus\n" + "".join(f"{k}\t{v}\n" for k, v in flags.items())
    )
    (d / "viral_reads.fasta").write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(reads)))
    return d


def test_verdict_rules(tmp_path):
    assert mod.verdict(_ev(tmp_path, "none", []))["verdict"] == "no_support"
    # every read host, flags unflagged because BLAST masked all reads (GSM5725695)
    assert (
        mod.verdict(_ev(tmp_path, "allhost", [("HOST|chr1", "host", 40, 1.0)]))["verdict"]
        == "host_best"
    )
    assert (
        mod.verdict(
            _ev(
                tmp_path,
                "hflag",
                [("HOST|chr1", "host", 1, 0), ("VIRUS|x", "virus", 5, 0.1)],
                {"host_homology": "flagged"},
            )
        )["verdict"]
        == "host_best"
    )
    assert (
        mod.verdict(_ev(tmp_path, "lc", [("VIRUS|x", "virus", 5, 0.1)]))["verdict"]
        == "low_complexity"
    )
    got = mod.verdict(
        _ev(tmp_path, "ok", [("VIRUS|x", "virus", 5, 0.9)], reads=(mod.TSO + "AC", "ACGT"))
    )
    assert got["verdict"] == "viral_best" and got["tso_reads"] == 1


def test_clean_viral_reads_are_not_low_complexity(tmp_path):
    """VERDICT-01: EBV (1.52 M of 1.55 M reads on virus, templated bodies) was labelled low_complexity."""
    rows = [("VIRUS|EBV", "virus", 1_522_203, 0.97), ("HOST|chr1", "host", 2326, 1.0)]
    assert mod.verdict(_ev(tmp_path, "ebv", rows))["verdict"] == "viral_best"


def test_complex_fraction_is_read_weighted_and_blank_is_not_zero(tmp_path):
    # a 1-read low-complexity reference must not outvote 1000 clean reads
    rows = [("VIRUS|a", "virus", 1000, 0.95), ("VIRUS|b", "virus", 1, 0.0)]
    assert mod.verdict(_ev(tmp_path, "mix", rows))["verdict"] == "viral_best"
    # a missing value is "not measured", not "no templated body"
    assert (
        mod.verdict(_ev(tmp_path, "blank", [("VIRUS|x", "virus", 5, "")]))["verdict"]
        == "viral_best"
    )


def test_raw_duplicate_weight_does_not_decide_complexity(tmp_path):
    d = _ev(tmp_path, "duplicates", [("VIRUS|x", "virus", 2, 0.5)])
    with (d / "alignment_qc.tsv").open("a") as fh:
        fh.write("raw\tVIRUS|x\tvirus\t100\t0.01\n")
    assert mod.verdict(d)["verdict"] == "viral_best"


def _molecules(d, virus, host, tie=0):
    rows = "".join(
        f"counted\t{k}\t{n}\n"
        for k, n in (("virus_best", virus), ("host_best", host), ("tie", tie), ("unaligned", 0))
    )
    (d / "molecule_verdict_summary.tsv").write_text("scope\tverdict\tmolecules\n" + rows)
    return d


def test_molecule_summary_decides_over_the_run_level_host_flag(tmp_path):
    """VERDICT-01: host_homology is one number for the run; the counted molecules decide the call."""
    flagged = {"host_homology": "flagged"}
    rows = [("VIRUS|x", "virus", 5, 0.9)]
    got = mod.verdict(_molecules(_ev(tmp_path, "pos", rows, flagged), virus=30, host=2))
    assert got["verdict"] == "viral_best" and got["virus_best_molecules"] == 30
    assert (
        mod.verdict(_molecules(_ev(tmp_path, "neg", rows, flagged), virus=3, host=40))["verdict"]
        == "host_best"
    )
    # ties are not evidence for the virus
    assert (
        mod.verdict(_molecules(_ev(tmp_path, "tie", rows), virus=1, host=0, tie=5))["verdict"]
        == "host_best"
    )
    assert (
        mod.verdict(_molecules(_ev(tmp_path, "none", rows), virus=0, host=0))["verdict"]
        == "no_support"
    )


def _table(d, rows):
    """molecule_verdicts.tsv.gz as `viralscan evidence` writes it: (cb, verdict, counted)."""
    with gzip.open(d / "molecule_verdicts.tsv.gz", "wt") as fh:
        fh.write("cb\tub\treads\tcounted\tvirus_as\thost_as\tverdict\n")
        for i, (cb, v, counted) in enumerate(rows):
            fh.write(f"{cb}\tU{i}\t1\t{counted}\t\t\t{v}\n")
    return d


def test_reference_cell_restriction_drops_off_list_molecules(tmp_path):
    """VERDICT-01: molecules in barcodes outside the frozen cells must not decide the call."""
    qc = [("VIRUS|x", "virus", 5, 0.9)]
    rows = [("A", "virus_best", True), ("B", "virus_best", True)] + [
        ("EMPTY", "host_best", True) for _ in range(5)
    ]
    d = _table(_molecules(_ev(tmp_path, "m", qc), virus=2, host=5), rows)
    assert mod.verdict(d)["verdict"] == "host_best"  # all cells: 2 of 7 molecules on the virus
    got = mod.verdict(d, {"A", "B"})
    assert got["verdict"] == "viral_best" and got["molecule_scope"] == "restricted"
    assert (got["virus_best_molecules"], got["host_best_molecules"]) == (2, 0)
    # the reverse: virus molecules in empty droplets must not flatter a host-best call in the cells
    flip = [("A", "host_best", True)] * 3 + [("EMPTY", "virus_best", True)] * 9
    d2 = _table(_molecules(_ev(tmp_path, "f", qc), virus=9, host=3), flip)
    assert mod.verdict(d2)["verdict"] == "viral_best"
    assert mod.verdict(d2, {"A"})["verdict"] == "host_best"
    # uncounted molecules never count
    d3 = _table(_molecules(_ev(tmp_path, "u", qc), virus=0, host=0), [("A", "virus_best", False)])
    assert mod.verdict(d3, {"A"})["verdict"] == "no_support"


def test_unrestricted_evidence_is_labelled_not_passed_off_as_restricted(tmp_path):
    qc = [("VIRUS|x", "virus", 5, 0.9)]
    old = _molecules(_ev(tmp_path, "old", qc), virus=30, host=2)  # no molecule_verdicts.tsv.gz
    assert mod.verdict(old, {"A"})["molecule_scope"] == "no_molecule_table"
    assert mod.verdict(old)["molecule_scope"] == "all_cells"


def test_cli_reads_frozen_cells_by_sample_prefix(tmp_path, capsys):
    root = tmp_path / "ev"
    root.mkdir()
    qc = [("VIRUS|x", "virus", 5, 0.9)]
    _table(
        _molecules(_ev(root, "ds1__S1__combined_off__HHV_6B", qc), virus=1, host=5),
        [("A", "virus_best", True)] + [("EMPTY", "host_best", True)] * 5,
    )
    cells = tmp_path / "reference_cells"
    cells.mkdir()
    (cells / "ds1__S1.tsv").write_text("barcode\nA\n")
    mod.main([str(root), "--reference-cells", str(cells)])
    row = capsys.readouterr().out.splitlines()[1].split("\t")
    assert row[-2:] == ["restricted", "viral_best"]
    (cells / "ds1__S1.tsv").unlink()
    mod.main([str(root), "--reference-cells", str(cells)])
    assert capsys.readouterr().out.splitlines()[1].split("\t")[-2] == "no_frozen_cells"
