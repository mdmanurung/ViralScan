"""``run`` needs ``--strand`` unless tchem confirms a 3' chemistry (PLAN STRAND-EXPLICIT-01).

No kb and no tchem run here: the chemistry blocks are built from the same dataclasses the
run records in ``run_manifest.json``.
"""

import json
import logging

import pytest

from viralscan import kbparams, menu
from viralscan.chemistry import Detection
from viralscan.tchem import Opinion, TchemCall


def _blocks(layout="10x 3' v2", status="agree", strand=None, sample="s"):
    det = Detection("10xv2", "bundled on-list", "test", 26, None, 10, 1000).as_block()
    call = TchemCall("match", "high", layout, "ok", None, None, strand, "test")
    det["tchem"] = Opinion(status, f"tchem {status}", call).as_block()
    return {sample: det}


def _args(strand=None, resume=False):
    return type(
        "A",
        (),
        {
            "strand": strand,
            "resume": resume,
            "sample1": "/x/s_R1.fastq.gz",
            "sample2": "/x/s_R2.fastq.gz",
            "technology": "10xv2",
            "tchem": "auto",
        },
    )()


@pytest.fixture(autouse=True)
def _id(monkeypatch):
    monkeypatch.setattr(menu, "_sample_id", lambda path: "s")
    monkeypatch.setattr(menu, "split_comma_paths", lambda text: text.split(","))


def _stops(args, blocks, out, caplog):
    with caplog.at_level(logging.ERROR), pytest.raises(SystemExit) as exc:
        menu._require_explicit_strand(args, out, blocks)
    return exc.value.code, caplog.text


@pytest.mark.parametrize("layout", ["10x 3' v2", "10x 3' v3/v3.1", "Drop-seq"])
def test_confirmed_three_prime_runs_without_the_flag(tmp_path, layout):
    menu._require_explicit_strand(_args(), tmp_path, _blocks(layout))


def test_five_prime_without_the_flag_stops_with_advice(tmp_path, caplog):
    code, text = _stops(_args(), _blocks("10x 5' v1/v2", strand="reverse"), tmp_path, caplog)
    assert code == 2
    assert "--strand is required" in text
    assert "10x 5' v1/v2" in text
    assert "advice: --strand reverse" in text


def test_five_prime_with_an_explicit_strand_runs(tmp_path):
    menu._require_explicit_strand(_args("reverse"), tmp_path, _blocks("10x 5' v1/v2"))


def test_tchem_unavailable_without_the_flag_stops(tmp_path, caplog):
    blocks = _blocks()
    blocks["s"]["tchem"] = Opinion("unavailable", "tchem is not installed").as_block()
    code, text = _stops(_args(), blocks, tmp_path, caplog)
    assert code == 2 and "not confirmed" in text


def test_no_block_without_the_flag_stops(tmp_path, caplog, monkeypatch):
    monkeypatch.setattr(
        kbparams, "recommend_kb_parameters", lambda *a, **k: (_ for _ in ()).throw(OSError("x"))
    )
    code, text = _stops(_args(), {}, tmp_path, caplog)
    assert code == 2 and "no advice available" in text


def test_auto_is_unaffected(tmp_path):
    menu._require_explicit_strand(_args("auto"), tmp_path, _blocks("10x 5' v1/v2"))


def _manifest(out, strand):
    out.mkdir(exist_ok=True)
    (out / "run_manifest.json").write_text(json.dumps({"options": {"strand": strand}}))


def test_resume_with_a_recorded_strand_skips_the_gate(tmp_path):
    _manifest(tmp_path, "reverse")
    menu._require_explicit_strand(_args(resume=True), tmp_path, _blocks("10x 5' v1/v2"))


def test_resume_without_a_strand_stops_when_kb_would_rerun(tmp_path, caplog):
    _manifest(tmp_path, None)
    code, _ = _stops(_args(resume=True), _blocks("10x 5' v1/v2"), tmp_path, caplog)
    assert code == 2


def test_resume_without_a_strand_passes_when_kb_is_done(tmp_path):
    _manifest(tmp_path, None)
    (tmp_path / "s" / "log").mkdir(parents=True)
    (tmp_path / "s" / "log" / "kb.done").write_text("")
    menu._require_explicit_strand(_args(resume=True), tmp_path, _blocks("10x 5' v1/v2"))


def test_only_agreeing_three_prime_counts():
    assert kbparams.is_confirmed_three_prime(_blocks()["s"]["tchem"])
    assert not kbparams.is_confirmed_three_prime(_blocks(status="disagree")["s"]["tchem"])
    assert not kbparams.is_confirmed_three_prime(None)
