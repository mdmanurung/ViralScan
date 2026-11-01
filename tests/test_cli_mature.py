"""CLI behaviour: exit codes, flag aliases, fast --version (PLAN CLI-MATURE-01)."""

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
