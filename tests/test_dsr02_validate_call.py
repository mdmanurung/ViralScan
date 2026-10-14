"""DSR-02: scripts/dsr02_validate_call.py (per-read verdicts with a repeat class, parameterised driver)."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import shlex
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location(
    "dsr02_validate_call",
    Path(__file__).resolve().parents[1] / "scripts" / "dsr02_validate_call.py",
)
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)

TELO = "TTAGGG" * 8
COMPLEX = "ACGTTGCAAGCTTCGATCGGATCCATGCAGTCAAGTCGA"


def sam(name, flag, ref, seq, score):
    tags = [] if score is None else [f"AS:i:{score}"]
    return "\t".join([name, str(flag), ref, "1", "60", "*", "*", "0", "0", seq, "*", *tags])


def rows_of(lines, counted=None):
    return {r["read"]: r for r in mod.read_verdict_rows(lines, counted)}


def test_repeat_classes_and_precedence():
    assert mod.repeat_class(TELO) == "telomeric"
    assert mod.repeat_class("CCCTAA" * 8) == "telomeric"  # the other strand
    assert mod.repeat_class("GG" + TELO[1:]) == "telomeric"  # any frame
    assert mod.repeat_class(mod.TSO + COMPLEX) == "tso"
    assert mod.repeat_class(COMPLEX + mod.TSO_RC) == "tso"
    assert mod.repeat_class("TTAGGG" + COMPLEX * 3) == "none"  # one hexamer is not a repeat array
    assert mod.repeat_class("A" * 40) == "low_complexity"
    assert mod.repeat_class(COMPLEX) == "none"
    assert mod.repeat_class(mod.TSO + TELO) == "tso"  # first match wins


def test_per_read_verdict_and_telomere_class():
    """F-025: telomere-repeat reads that tie or lose to host must be labelled, not just counted."""
    lines = [
        "@SQ\tSN:VIRUS|x\tLN:100",
        sam("AAAC_GGGT_1|a", 0, "VIRUS|x", COMPLEX, 60),
        sam("AAAC_GGGT_1|a", 256, "HOST|chr1", "*", 40),
        sam("AAAC_GGGT_2|b", 0, "HOST|chr1", TELO, 90),
        sam("AAAC_GGGT_2|b", 256, "VIRUS|x", "*", 90),
        sam("CCCA_TTTT_3|c", 0, "HOST|chr2", COMPLEX, 70),
        sam("CCCA_TTTT_4|d", 4, "*", COMPLEX, None),
    ]
    got = rows_of(lines, counted={1, 2})
    assert got["AAAC_GGGT_1|a"]["verdict"] == "virus_best"
    assert got["AAAC_GGGT_1|a"]["host_as"] == 40 and got["AAAC_GGGT_1|a"]["virus_as"] == 60
    tie = got["AAAC_GGGT_2|b"]
    assert tie["verdict"] == "tie" and tie["repeat_class"] == "telomeric"
    assert got["CCCA_TTTT_3|c"]["verdict"] == "host_best"
    assert got["CCCA_TTTT_4|d"]["verdict"] == "unaligned"
    assert [got[k]["counted"] for k in sorted(got)] == [1, 1, 0, 0]
    assert got["AAAC_GGGT_1|a"]["read_number"] == 1 and got["AAAC_GGGT_1|a"]["cb"] == "AAAC"


def test_write_read_verdicts_files(tmp_path):
    with gzip.open(tmp_path / "read_lineage.tsv.gz", "wt") as fh:
        fh.write("read_number\tassigned_weight\n1\t1.0\n2\t0\n")
    lines = [
        sam("AAAC_GGGT_1|a", 0, "VIRUS|x", COMPLEX, 60),
        sam("AAAC_GGGT_2|b", 0, "HOST|c", TELO, 80),
    ]
    summary = mod.write_read_verdicts(lines, tmp_path, tmp_path / "read_lineage.tsv.gz")
    assert {(r["verdict"], r["repeat_class"]): r["reads"] for r in summary} == {
        ("virus_best", "none"): 1,
        ("host_best", "telomeric"): 1,
    }
    with gzip.open(tmp_path / "read_verdicts.tsv.gz", "rt") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["counted"] for r in table] == ["1", "0"]
    assert (tmp_path / "read_verdict_summary.tsv").is_file()


def test_molecular_summary_is_not_inflated_by_five_duplicate_reads(tmp_path):
    names = [f"CB_U1_{i}|r" for i in range(1, 6)] + ["CB_U2_6|r"]
    lines = [sam(n, 0, "VIRUS|x", COMPLEX, 80) for n in names]
    lines += [sam(n, 256, "HOST|c", "*", 60) for n in names]
    with gzip.open(tmp_path / "read_lineage.tsv.gz", "wt") as fh:
        fh.write("read_number\tassigned_weight\n2\t1\n6\t1\n")
    raw = mod.write_read_verdicts(lines, tmp_path, tmp_path / "read_lineage.tsv.gz")
    assert sum(r["reads"] for r in raw) == 6
    with (tmp_path / "molecule_read_verdict_summary.tsv").open() as fh:
        summary = list(csv.DictReader(fh, delimiter="\t"))
    assert sum(int(r["molecules"]) for r in summary) == 2
    with gzip.open(tmp_path / "molecule_read_verdicts.tsv.gz", "rt") as fh:
        molecules = list(csv.DictReader(fh, delimiter="\t"))
    assert len(molecules) == 2 and all(
        r["counted"] == "1" and r["host_as"] == "60" for r in molecules
    )


def _calls(tmp_path):
    p = tmp_path / "evidence" / "calls.tsv"
    p.parent.mkdir()
    p.write_text(
        "dataset\tarm\tsample\trun_dir\tvirus\tmolecules\tanellovirus\texpected\n"
        "ebv\tcombined_off\tS1\t/runs/S1/S1\tHuman betaherpesvirus 6B\t8.0\tfalse\tfalse\n"
        "hpv16\tcombined_off\tS2\t/runs/S2/S2\tHPV 118\t5.0\tfalse\tfalse\n"
    )
    return p


def test_driver_takes_dataset_specifics_from_the_row_and_arguments(tmp_path, capsys):
    calls = _calls(tmp_path)
    base = [
        "--calls",
        str(calls),
        "--host-fasta",
        "host.fa",
        "--panel-fasta",
        "panel.fa",
        "--dry-run",
    ]
    assert mod.main([*base, "--row", "0"]) == 0
    first = shlex.split(capsys.readouterr().out)
    assert first[first.index("--run-dir") + 1] == "/runs/S1/S1"
    assert first[first.index("--virus") + 1] == "Human betaherpesvirus 6B"
    assert first[first.index("-o") + 1] == str(
        calls.parent / "ebv__S1__combined_off__Human_betaherpesvirus_6B"
    )
    assert first[first.index("--viral-fasta") + 1] == "panel.fa"
    mod.main([*base, "--row", "1", "--panel-for", "hpv16=hpv_panel.fa"])
    second = shlex.split(capsys.readouterr().out)
    assert second[second.index("--viral-fasta") + 1] == "hpv_panel.fa"
    assert second[second.index("--run-dir") + 1] == "/runs/S2/S2"


def test_evidence_command_is_an_argv_list_with_the_virus_intact():
    cmd = mod.evidence_command("/r", "Human betaherpesvirus 6B", "/o", "p.fa", "h.fa", 4)
    assert isinstance(cmd, list) and "Human betaherpesvirus 6B" in cmd
    assert cmd[cmd.index("--cores") + 1] == "4"


def test_missing_inputs_are_errors(tmp_path):
    calls = _calls(tmp_path)
    with pytest.raises(SystemExit):
        mod.main(["--calls", str(calls), "--row", "0", "--host-fasta", "h.fa", "--dry-run"])
    with pytest.raises(SystemExit):
        mod.main(
            [
                "--calls",
                str(calls),
                "--row",
                "7",
                "--host-fasta",
                "h.fa",
                "--panel-fasta",
                "p",
                "--dry-run",
            ]
        )


def test_done_directory_is_not_rerun_and_postprocess_runs(tmp_path, monkeypatch):
    calls = _calls(tmp_path)
    out = calls.parent / "ebv__S1__combined_off__Human_betaherpesvirus_6B"
    out.mkdir()
    (out / "interpretation_flags.tsv").write_text("flag\tstatus\n")
    ran, post = [], []
    monkeypatch.setattr(mod.subprocess, "run", lambda *a, **k: ran.append(a))
    monkeypatch.setattr(mod, "postprocess", lambda d, **k: post.append((d, k)))
    argv = ["--calls", str(calls), "--row", "0", "--host-fasta", "h", "--panel-fasta", "p"]
    mod.main(argv)
    assert not ran and post == [(out, {"keep_intermediates": False})]
    (out / "interpretation_flags.tsv").unlink()
    mod.main(argv)
    assert len(ran) == 1 and ran[0][0][0:4] == [
        mod.sys.executable,
        "-m",
        "viralscan.menu",
        "evidence",
    ]


def test_postprocess_stamps_code_and_drops_bulky_files(tmp_path, monkeypatch):
    (tmp_path / "competitive_reads.raw.bam").write_bytes(b"")
    (tmp_path / "competitive_host_target.fasta").write_text(">x\nA\n")
    (tmp_path / "blast").mkdir()
    (tmp_path / "blast" / "competitive_db.nsq").write_text("x")
    (tmp_path / "blast" / "keep.tsv").write_text("x")

    class FakePopen:
        returncode = 0
        stdout = [sam("AAAC_GGGT_1|a", 0, "VIRUS|x", COMPLEX, 60)]

        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    real_popen = mod.subprocess.Popen  # software_identity() runs git through the same module

    def popen(cmd, *a, **k):
        return FakePopen() if cmd[0] == "samtools" else real_popen(cmd, *a, **k)

    monkeypatch.setattr(mod.subprocess, "Popen", popen)
    summary = mod.postprocess(tmp_path)
    assert summary == [{"verdict": "virus_best", "repeat_class": "none", "reads": 1}]
    assert (tmp_path / "code_sha.txt").read_text().strip()
    assert not (tmp_path / "competitive_host_target.fasta").exists()
    assert not (tmp_path / "blast" / "competitive_db.nsq").exists()
    assert (tmp_path / "blast" / "keep.tsv").exists()
