"""CLI behaviour: exit codes, flag aliases, fast --version (PLAN CLI-MATURE-01)."""

import argparse
import gzip
import subprocess
import sys

import pytest

from viralscan import menu
from viralscan.cli.parsers import build_parser


def _sub(name):
    return build_parser()._subparsers._group_actions[0].choices[name]


def test_keyboard_interrupt_exits_130(monkeypatch):
    def boom():
        raise KeyboardInterrupt

    monkeypatch.setattr(menu, "_main", boom)
    with pytest.raises(SystemExit) as exc:
        menu.main()
    assert exc.value.code == 130


def test_declining_reference_prompt_exits_1(monkeypatch):
    args = build_parser().parse_args(
        ["run", "--reference", "--fasta", __file__, "--gtf", __file__, "-i", __file__]
        + ["--r1", "a", "--r2", "b", "-o", "o"]
    )
    monkeypatch.setattr("builtins.input", lambda _prompt: "n")
    with pytest.raises(SystemExit) as exc:
        menu.errorhandler(args)
    assert exc.value.code == 1


def test_flag_combination_error_exits_2():
    args = build_parser().parse_args(
        ["run", "--reference", "--fasta", __file__, "--r1", "a", "--r2", "b", "-o", "o"]
    )
    with pytest.raises(SystemExit) as exc:
        menu.errorhandler(args)
    assert exc.value.code == 2


@pytest.mark.parametrize("cmd", ["star-index", "translated-index", "translated", "discover"])
def test_threads_is_alias_of_cores(cmd):
    flags = {f for a in _sub(cmd)._actions for f in a.option_strings}
    assert {"--cores", "-c", "--threads"} <= flags


def test_translated_index_output_and_cores_aliases():
    sub = _sub("translated-index")
    ns = sub.parse_args(
        ["--proteins-dir", "p", "--host-mask", "h", "--output", "d", "--cores", "3"]
    )
    assert (ns.out_dir, ns.threads) == ("d", 3)
    ns = sub.parse_args(
        ["--proteins-dir", "p", "--host-mask", "h", "--out-dir", "d", "--threads", "3"]
    )
    assert (ns.out_dir, ns.threads) == ("d", 3)


def test_verbose_quiet_exclusive():
    with pytest.raises(SystemExit) as exc:
        build_parser().parse_args(["doctor", "--verbose", "--quiet"])
    assert exc.value.code == 2


def test_version_without_numpy_or_yaml():
    code = (
        "import sys; sys.modules['numpy']=None; sys.modules['yaml']=None; "
        "sys.argv=['viralscan','--version']; "
        "from viralscan.__main__ import main; main()"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.startswith("viralscan ")


def test_bare_viralscan_prints_help_to_stderr_exit_2(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["viralscan"])
    with pytest.raises(SystemExit) as exc:
        menu.main()
    assert exc.value.code == 2
    assert "subcommand" in capsys.readouterr().err


# ---- CLI-EXIT-02: usage errors exit 2, no tracebacks ---------------------------------------


@pytest.fixture
def run_argv(tmp_path):
    """A valid minimal `run` command line; tests append the flags under test."""
    paths = {n: tmp_path / n for n in ("idx", "t2g", "a_R1.fastq", "a_R2.fastq")}
    for p in paths.values():
        p.write_text("")
    return [
        "run", "-i", str(paths["idx"]), "-t", str(paths["t2g"]),
        "--r1", str(paths["a_R1.fastq"]), "--r2", str(paths["a_R2.fastq"]),
        "-o", str(tmp_path / "out"),
    ]  # fmt: skip


def _exit_code(fn, *a):
    with pytest.raises(SystemExit) as exc:
        fn(*a)
    return exc.value.code


def test_valid_run_argv_passes_errorhandler(run_argv):
    menu.errorhandler(build_parser().parse_args(run_argv))  # guards the other tests


@pytest.mark.parametrize(
    "flags",
    [
        ["--cores", "0"],
        ["--se-threshold", "-5"],
        ["--cell-calling", "external"],
        ["--star-confirm"],
        ["--auto-evidence"],
    ],
)
def test_run_flag_errors_exit_2(run_argv, flags):
    args = build_parser().parse_args(run_argv + flags)
    assert _exit_code(menu.errorhandler, args) == 2


def test_r1_r2_count_mismatch_exits_2(run_argv, tmp_path):
    argv = run_argv + ["--r1", f"{tmp_path / 'a_R1.fastq'},{tmp_path / 'a_R1.fastq'}"]
    assert _exit_code(menu.errorhandler, build_parser().parse_args(argv)) == 2


def test_ncbi_accession_with_index_or_transcripts_exits_2(run_argv, caplog):
    args = build_parser().parse_args(run_argv + ["--ncbi-accession", "NC_000001"])
    assert _exit_code(menu.errorhandler, args) == 2
    assert "-i/-t" in caplog.text


def test_bad_run_config_combination_fails_before_network(run_argv):
    """RunConfig-only rule (positive-control scope without a gene) exits 2 up front."""
    args = build_parser().parse_args(run_argv + ["--positive-control-scope", "panel_mechanics"])
    menu.errorhandler(args)  # passes: only RunConfig knows this rule
    assert _exit_code(menu._precheck_run_config, args) == 2


def test_precheck_accepts_reference_build(run_argv, tmp_path):
    args = build_parser().parse_args(run_argv + ["--anello-align"])
    menu._precheck_run_config(args)  # needs no built index


def test_ignored_flags_warn(run_argv, caplog):
    argv = run_argv + [
        "--ncbi-email",
        "a@b.c",
        "--emptydrops-seed",
        "7",
        "--cell-calling",
        "none",
        "--no-multimapping",
    ]
    argv += ["--multimap-em-tol", "0.1"]
    with caplog.at_level("WARNING"):
        menu._warn_ignored_flags(build_parser().parse_args(argv))
    text = caplog.text
    assert "--ncbi-email ignored" in text
    assert "--emptydrops-seed ignored" in text
    assert "--multimap-em-tol ignored" in text


def test_default_flags_do_not_warn(run_argv, caplog):
    with caplog.at_level("WARNING"):
        menu._warn_ignored_flags(build_parser().parse_args(run_argv))
    assert caplog.text == ""


def test_build_ref_missing_host_exits_2(tmp_path):
    from viralscan.scripts.build_reference import build_ref_main

    args = build_parser().parse_args(
        ["build-ref", "--output", str(tmp_path), "--virus-accessions", "NC_1"]
    )
    assert _exit_code(build_ref_main, args) == 2


@pytest.mark.parametrize(
    "manifest",
    [
        "[1, 2]",
        '{"options": [1], "input_fingerprints": {"nocolon": "x"}}',
        '{"options": {}, "input_fingerprints": {"a:zz": "x"}}',
        '{"options": "s", "input_fingerprints": ["a"]}',
        '{"options": {}, "input_fingerprints": {"a:0": "x"}, "input_paths": 5}',
    ],
)
def test_validate_run_malformed_manifest_reports_issues(tmp_path, manifest):
    from viralscan.validation import validate_run

    (tmp_path / "run_manifest.json").write_text(manifest)
    assert validate_run(tmp_path)["ok"] is False


def test_json_report_unwritable_path_exits_2(tmp_path, monkeypatch, capsys):
    argv = ["viralscan", "validate-run", "--run-dir", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", argv + ["--json", str(tmp_path / "no" / "x.json")])
    assert _exit_code(menu.main) == 2
    assert "Cannot write --json" in capsys.readouterr().err


def test_evidence_plot_missing_dir_exits_2(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(
        sys, "argv", ["viralscan", "evidence-plot", "--evidence-dir", str(tmp_path / "nope")]
    )
    assert _exit_code(menu.main) == 2
    assert "Traceback" not in capsys.readouterr().err


@pytest.mark.parametrize("case", ["outfile", "dir_as_t", "dir_whitelist", "trunc_gz", "dir_r1"])
def test_bad_inputs_give_one_error_line(tmp_path, case):
    good = tmp_path / "R1.fastq"
    good.write_text("@r\nACGT\n+\nIIII\n")
    trunc = tmp_path / "t.fastq.gz"
    trunc.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n" * 5000)[:-20])
    run = ["run", "--r1", str(good), "--r2", str(good), "-i", str(good)]
    argv = {
        "outfile": [*run, "-o", str(good)],
        "dir_as_t": [*run, "-t", str(tmp_path), "-o", str(tmp_path / "o")],
        "dir_whitelist": ["check-whitelist", "--r1", str(good), "--whitelist", str(tmp_path)],
        "trunc_gz": ["check-chemistry", "--r1", str(trunc)],
        "dir_r1": ["check-chemistry", "--r1", str(tmp_path)],
    }[case]
    r = subprocess.run([sys.executable, "-m", "viralscan", *argv], capture_output=True, text=True)
    assert r.returncode in (1, 2), r.stderr
    assert "Traceback" not in r.stderr, r.stderr
    assert "ERROR" in r.stderr or r.returncode == 2, r.stderr


def test_translated_index_checks_inputs_before_download(monkeypatch, tmp_path):
    from viralscan import translated

    monkeypatch.setattr(translated, "fetch_palmdb", lambda *a: pytest.fail("downloaded"))
    ns = argparse.Namespace(
        out_dir=str(tmp_path / "x"),
        proteins_dir="/nonexist",
        host_mask="/nonexist",
        palmdb=True,
        threads=1,
    )
    with pytest.raises(translated.TranslatedError):
        translated.run_translated_index(ns)
