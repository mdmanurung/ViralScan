"""``viralscan demo`` (PLAN DEMO-01): argument handling, bundled data and the digest comparison.

The end-to-end run is ``tests/integration/test_demo.py``.
"""

import hashlib
import json
from pathlib import Path

import pytest

from viralscan import demo
from viralscan.cli.parsers import build_parser

GOLDEN = Path(__file__).parent / "integration" / "golden" / "tiny_run.json"


def _fake_run(root: Path, tables: dict[str, str]) -> Path:
    """A run root with one sample directory, the layout ``output_digests`` reads."""
    sample = root / "S1"
    (sample / "results").mkdir(parents=True)
    (sample / "config.yaml").write_text("x: 1\n")
    for name, text in tables.items():
        (sample / "results" / name).write_text(text)
    (sample / "summary.txt").write_text("Samples: 1\nRuntime: 12 s\n")
    (root / "combined_summary.tsv").write_text("a\tb\n")
    return root


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def test_parser_defaults_and_flags() -> None:
    args = build_parser().parse_args(["demo"])
    assert (args._subcommand, args.output, args.cores, args.keep, args.force) == (
        "demo",
        "viralscan_demo",
        2,
        False,
        False,
    )
    args = build_parser().parse_args(["demo", "--output", "d", "--cores", "4", "--keep", "--force"])
    assert (args.output, args.cores, args.keep, args.force) == ("d", 4, True, True)


def test_bundled_inputs_are_small_and_match_the_golden_digests() -> None:
    files = [demo.demo_data() / name for name in (*demo.DEMO_INPUTS, "expected_digests.json")]
    assert all(path.is_file() for path in files)
    assert sum(path.stat().st_size for path in files) < 100_000
    # One pin: the packaged copy is the golden file the integration test compares with.
    assert demo.expected_digests() == json.loads(GOLDEN.read_text(encoding="utf-8"))


def test_output_digests_of_a_fake_run(tmp_path: Path) -> None:
    run = _fake_run(tmp_path, {"viral_summary.tsv": "v\n", "notes.json": "{}"})
    assert demo.output_digests(run) == {
        "combined_summary.tsv": _sha("a\tb\n"),
        "results/viral_summary.tsv": _sha("v\n"),
        "summary.txt": _sha("Samples: 1\n"),  # the Runtime line is dropped; .json is ignored
    }


def test_compare_digests_names_every_kind_of_difference() -> None:
    got = {"a": "1", "b": "2", "d": "4"}
    want = {"a": "1", "b": "9", "c": "3"}
    assert demo.compare_digests(got, want) == {
        "a": "PASS",
        "b": "FAIL",
        "c": "MISSING",
        "d": "UNEXPECTED",
    }


def test_missing_tool_exits_2_with_the_env_hint(monkeypatch, capsys, tmp_path: Path) -> None:
    monkeypatch.setattr(demo, "missing_tools", lambda tools: ["kb"])
    assert demo.run_demo(tmp_path / "d") == 2
    err = capsys.readouterr().err
    assert "kb" in err and "viralscan check env --profile full" in err
    assert not (tmp_path / "d").exists()


def test_existing_output_is_refused_without_force(monkeypatch, capsys, tmp_path: Path) -> None:
    monkeypatch.setattr(demo, "missing_tools", lambda tools: [])
    monkeypatch.setattr(demo, "tool_path", lambda name: "/bin/true")
    out = tmp_path / "d"
    out.mkdir()
    (out / "mine.txt").write_text("keep me")
    assert demo.run_demo(out) == 2
    assert "already holds files" in capsys.readouterr().err
    assert (out / "mine.txt").read_text() == "keep me"


@pytest.fixture
def fake_pipeline(monkeypatch, tmp_path: Path):
    """Tools present, ``kallisto``/``run`` replaced by a step that writes a fake run."""
    from viralscan import validation

    tables = {"viral_summary.tsv": "v\n"}
    state = {"tables": tables, "check_ok": True}

    def step(label, cmd, cwd, log):
        if label == "viralscan run":
            _fake_run(Path(cmd[cmd.index("-o") + 1]), state["tables"])
        return True

    monkeypatch.setattr(demo, "missing_tools", lambda tools: [])
    monkeypatch.setattr(demo, "tool_path", lambda name: "/bin/true")
    monkeypatch.setattr(demo, "_step", step)
    monkeypatch.setattr(
        validation, "validate_run", lambda run, verify_inputs=True: {"ok": state["check_ok"]}
    )
    monkeypatch.setattr(validation, "tool_provenance", lambda: {})
    monkeypatch.setattr(
        demo,
        "expected_digests",
        lambda: {
            "combined_summary.tsv": _sha("a\tb\n"),
            "results/viral_summary.tsv": _sha("v\n"),
            "summary.txt": _sha("Samples: 1\n"),
        },
    )
    return state, tmp_path / "d"


def test_pass_exits_0_and_leaves_only_the_report(fake_pipeline, capsys) -> None:
    _, out = fake_pipeline
    assert demo.run_demo(out) == 0
    text = capsys.readouterr().out
    assert "PASS      results/viral_summary.tsv" in text and "demo: PASS" in text
    assert sorted(p.name for p in out.iterdir()) == [demo.REPORT]
    assert json.loads((out / demo.REPORT).read_text())["ok"] is True


def test_keep_leaves_the_run(fake_pipeline) -> None:
    _, out = fake_pipeline
    assert demo.run_demo(out, keep=True) == 0
    assert (out / "run" / "S1" / "config.yaml").is_file()


def test_a_changed_output_exits_1_and_keeps_the_run(fake_pipeline, capsys) -> None:
    state, out = fake_pipeline
    state["tables"] = {"viral_summary.tsv": "changed\n"}
    assert demo.run_demo(out) == 1
    assert "FAIL      results/viral_summary.tsv" in capsys.readouterr().out
    assert (out / "run").is_dir()
    assert json.loads((out / demo.REPORT).read_text())["ok"] is False


def test_a_failed_check_run_exits_1(fake_pipeline, capsys) -> None:
    state, out = fake_pipeline
    state["check_ok"] = False
    assert demo.run_demo(out) == 1
    assert "FAIL      check run" in capsys.readouterr().out


def test_force_replaces_a_previous_demo(fake_pipeline) -> None:
    _, out = fake_pipeline
    assert demo.run_demo(out, keep=True) == 0
    assert demo.run_demo(out) == 2  # refused without --force
    assert demo.run_demo(out, force=True) == 0


def test_cores_below_one_exits_2(monkeypatch) -> None:
    from viralscan import menu

    monkeypatch.setattr("sys.argv", ["viralscan", "demo", "--cores", "0"])
    with pytest.raises(SystemExit) as raised:
        menu.main()
    assert raised.value.code == 2
