"""One CLI grammar (PLAN CLI-UNIFY-01): `viralscan run`, removed spellings named, no prefix matching."""

import pytest

from viralscan.cli.parsers import REMOVED_OPTIONS, parse_args

RUN = ["run", "-o", "out", "--r1", "a_R1.fq.gz", "--r2", "a_R2.fq.gz"]


def test_run_keeps_the_config_dests():
    args = parse_args(RUN + ["--no-visual", "--cache-dir", "c", "--cell-type-labels", "l.csv"])
    assert args._subcommand == "run"
    assert (args.sample1, args.sample2, args.output) == ("a_R1.fq.gz", "a_R2.fq.gz", "out")
    assert (args.visual, args.reference, args.multimapping) == (False, False, True)
    assert (args.data_cache_dir, args.cell_types) == ("c", "l.csv")


def test_bare_invocation_warns_and_runs(capsys):
    assert parse_args(RUN[1:])._subcommand == "run"
    assert "deprecated" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("sub", "old"), [(sub, old) for sub, table in REMOVED_OPTIONS.items() for old in table]
)
def test_removed_option_names_its_replacement(sub, old, capsys):
    with pytest.raises(SystemExit) as exc:
        parse_args([sub, old, "x"])
    assert exc.value.code == 2
    assert f"{old} was removed; use {REMOVED_OPTIONS[sub][old]}" in capsys.readouterr().err


def test_boolean_value_form_is_refused(capsys):
    with pytest.raises(SystemExit):
        parse_args(RUN + ["--visual", "False"])
    assert "--no-visual" in capsys.readouterr().err


@pytest.mark.parametrize(
    "argv",
    [
        RUN + ["--multi"],  # was a prefix of --multimapping
        ["validate-run", "--run", "r"],
        ["hostresponse", "--run-dir", "s", "--host-h5ad", "h", "--min-u", "3"],
    ],
)
def test_no_prefix_matching(argv):
    with pytest.raises(SystemExit):
        parse_args(argv)


@pytest.mark.parametrize(
    "argv",
    [
        # The upstream (mahfouzlab/ViralScan) README forms, kept working until CLI-UNIFY-01.
        ["-t", "t.txt", "-i", "i.idx", "-o", "o/", "-s1", "a_1.fq.gz", "-s2", "a_2.fq.gz"],
        ["-reference", "True", "-fasta", "g.fa", "-gtf", "g.gtf", "-o", "o/", "-s1", "a_1.fq.gz"],
    ],
)
def test_upstream_readme_forms_fail_with_the_replacement(argv, capsys):
    with pytest.raises(SystemExit) as exc:
        parse_args(argv)
    assert exc.value.code == 2
    assert "was removed; use" in capsys.readouterr().err


def test_hostresponse_run_dir_input_and_optional_output(capsys):
    base = ["hostresponse", "--run-dir", "run", "--host-h5ad", "h"]
    args = parse_args(base)
    assert (args.run_dir, args.output) == ("run", None)
    assert parse_args(base + ["--output", "hr"]).output == "hr"
    with pytest.raises(SystemExit):
        parse_args(["hostresponse", "--host-h5ad", "h"])  # --run-dir is required
    capsys.readouterr()
    with pytest.raises(SystemExit) as exc:
        parse_args(["hostresponse", "-o", "run", "--host-h5ad", "h"])
    assert exc.value.code == 2
    assert "-o was removed; use --run-dir" in capsys.readouterr().err
