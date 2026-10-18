"""TusslChemistry second opinion (PLAN CHEM-02): parsing, comparison and the run policy.

No real ``tchem`` is needed: a small fake executable answers ``--format json`` and
``--emit kb`` like the documented tool does (report schema 1.x; the kb command as a
``bc:umi:seq`` technology string). The exact text of the real ``--emit kb`` output
is not pinned by these tests; ``parse_emit`` only needs ``-x``, ``-w`` and
``--strand`` as tokens of a shell command.
"""

from __future__ import annotations

import argparse
import json
import os
import stat
import sys
from pathlib import Path

import pytest

from viralscan import chemistry, menu, tchem

V2 = "0,0,16:0,16,26:1,0,0"
V3 = "0,0,16:0,16,28:1,0,0"


def _report(verdict="match", layout="10x 3' v2", confidence="high", schema="1.1", **extra):
    rep = {
        "tool": "tussl-chemistry",
        "version": "0.1.0",
        "schema_version": schema,
        "verdict": verdict,
        "confidence": confidence if verdict == "match" else None,
        "layout": {"name": layout, "validation": {"level": "validated"}} if verdict == "match" else None,
        "candidates": ["10x 3' v2", "10x 3' v3"] if verdict == "ambiguous" else [],
        "cdna_strand": "forward" if verdict == "match" else None,
        "hint": "no list matched" if verdict == "no_match" else None,
        "unknown_future_key": {"x": 1},
    }
    rep.update(extra)
    return rep


def _fake(tmp_path: Path, report: dict, emit: str = "", json_rc=0, emit_rc=0, junk=False) -> str:
    script = tmp_path / "tchem"
    script.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "a = sys.argv[1:]\n"
        "if '--emit' in a:\n"
        f"    sys.stdout.write({emit!r}); sys.exit({emit_rc})\n"
        f"sys.stdout.write({'not json' if junk else json.dumps(report)!r}); sys.exit({json_rc})\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return str(script)


def _emit(tech=V2, wl="/lists/v2.txt", strand="forward") -> str:
    return (
        f"kb count \\\n  -i $KB_INDEX -g $KB_T2G \\\n  -x {tech} \\\n  -w {wl} \\\n"
        f"  --strand {strand} \\\n  R1.fastq.gz R2.fastq.gz\n"
    )


# Captured from tchem 0.1.0 (release v0.1.0, sha256 9112e501...488e) on a 10x 3' v2 library (100,000 reads,
# 2026-10-18): the exact text of `--emit kb` and the keys of `--format json`.
REAL_EMIT = """# tchem 0.1.0: 10x 3' v2 (confidence HIGH)
# validation: error_under_3pct: 626 held-out HIGH calls from 335 studies, none wrong (error under 0.9% at 95%); 27 more agree only with labels read off the reads (not counted); 85 disputed (label names another kit)
# set KB_INDEX: a kallisto index (kb ref -i)
# set KB_T2G: the transcript-to-gene map (kb ref -g)
kb count \\
  -i "${KB_INDEX:?set KB_INDEX to a kallisto index (kb ref -i)}" \\
  -g "${KB_T2G:?set KB_T2G to the transcript-to-gene map (kb ref -g)}" \\
  -x 0,0,16:0,16,26:1,0,0 \\
  -w /lists/10x_v2_permit.txt \\
  --strand forward \\
  -o kb_count /data/S_R1.fastq.gz /data/S_R2.fastq.gz
"""
REAL_REPORT = {
    "tool": "tussl-chemistry", "version": "0.1.0", "schema_version": "1.1", "verdict": "match",
    "confidence": "high",
    "layout": {"name": "10x 3' v2", "geometry": "1{b[16]u[10]x:}2{r:}", "barcode_lists": ["/lists/10x_v2_permit.txt"],
               "note": "", "long_read_kit": "3prime:v2",
               "validation": {"applies_to": "high", "level": "error_under_3pct", "held_out_high_calls": 626,
                              "studies": 335, "family_only_calls": 0, "disputed": 85, "wrong": 0,
                              "error_bound_95": 0.008955223880597015}, "family_validation": None},
    "candidates": [], "files": {"barcode_read": "R1", "cdna_read": "R2"}, "cdna_strand": "forward",
    "library": {"kind": "gene expression"}, "swapped": False, "strand_evidence": {}, "contradictions": [],
    "notes": [], "hint": None, "not_checked": ["10x Flex v1"], "reads_sampled": 100000, "read_kind": "short",
    "instrument": {}, "long_read_evidence": None, "fingerprint": None,
}


class TestParse:
    def test_real_tchem_0_1_0_emit_and_report(self):
        assert tchem.parse_emit(REAL_EMIT) == (V2, "/lists/10x_v2_permit.txt", "forward")
        c = tchem.parse_report(REAL_REPORT, REAL_EMIT)
        assert (c.verdict, c.layout, c.technology, c.strand, c.validation, c.version) == (
            "match", "10x 3' v2", V2, "forward", "error_under_3pct", "0.1.0")
        assert tchem.compare("10xv2", c).status == "agree"

    def test_multiome_gex_geometry_agrees_with_a_10xv3_run(self):
        emit = _emit(V3, "/lists/10x_arc_v1_gex_permit.txt")
        c = tchem.parse_report(_report(layout="10x Multiome GEX"), emit)
        assert tchem.compare("10xv3", c).status == "agree"

    def test_emit_with_continuations_and_placeholders(self):
        assert tchem.parse_emit(_emit()) == (V2, "/lists/v2.txt", "forward")

    def test_emit_without_a_list_or_strand(self):
        assert tchem.parse_emit(f"kb count -x {V2} R1 R2") == (V2, None, None)

    def test_emit_equals_form(self):
        assert tchem.parse_emit(f"kb count -x={V2} --strand=reverse a b") == (V2, None, "reverse")

    def test_match_report(self):
        c = tchem.parse_report(_report(), _emit())
        assert (c.verdict, c.layout, c.technology, c.strand, c.validation) == (
            "match", "10x 3' v2", V2, "forward", "validated")
        assert "high" in c.reason

    def test_ambiguous_and_no_match(self):
        a = tchem.parse_report(_report("ambiguous"))
        assert a.verdict == "ambiguous" and "10x 3' v3" in a.reason
        n = tchem.parse_report(_report("no_match"))
        assert n.verdict == "no_match" and n.reason == "no list matched"

    def test_other_tool_and_major_schema_refused(self):
        with pytest.raises(tchem.TchemError, match="not from tussl-chemistry"):
            tchem.parse_report({**_report(), "tool": "other"})
        with pytest.raises(tchem.TchemError, match="not 1.x"):
            tchem.parse_report(_report(schema="2.0"))
        tchem.parse_report(_report(schema="1.9"))  # a minor step is fine

    def test_unknown_verdict_refused(self):
        with pytest.raises(tchem.TchemError, match="unknown tchem verdict"):
            tchem.parse_report(_report(verdict="maybe"))


class TestCompare:
    def _call(self, tech=V2, verdict="match"):
        return tchem.parse_report(_report(verdict), _emit(tech) if verdict == "match" else None)

    def test_same_geometry_agrees_for_a_named_or_string_technology(self):
        assert tchem.compare("10xv2", self._call()).status == "agree"
        assert tchem.compare(V2, self._call()).status == "agree"

    def test_different_umi_length_disagrees(self):
        op = tchem.compare("10xv2", self._call(V3))
        assert op.status == "disagree" and "10xv2" in op.reason and "16 + 12" in op.reason

    def test_no_match_is_no_call(self):
        assert tchem.compare("10xv2", self._call(verdict="no_match")).status == "no_call"

    def test_match_without_kb_command_is_no_call(self):
        c = tchem.parse_report(_report(), None)
        assert tchem.compare("10xv2", c).status == "no_call"

    def test_no_technology_reports_info(self):
        op = tchem.compare(None, self._call())
        assert op.status == "info" and f"-x {V2}" in op.reason and "--strand forward" in op.reason

    def test_unparseable_technology_is_no_call(self):
        assert tchem.compare("not-a-chem", self._call()).status == "no_call"


class TestBinary:
    def test_explicit_beats_env_beats_path(self, tmp_path, monkeypatch):
        a = _fake(tmp_path, _report())
        monkeypatch.setenv(tchem.ENV_BIN, "/nonexistent/tchem")
        assert tchem.find_binary(a) == a
        assert tchem.find_binary() is None  # the env path is set but not executable: no PATH fallback
        monkeypatch.delenv(tchem.ENV_BIN)
        monkeypatch.setenv("PATH", str(tmp_path))
        assert tchem.find_binary() == str(tmp_path / "tchem")

    def test_whitelists_from_env(self, monkeypatch):
        monkeypatch.delenv(tchem.ENV_TCHEM_LISTS, raising=False)
        monkeypatch.setenv(tchem.ENV_LISTS, "/lists")
        assert tchem.find_whitelists() == "/lists"
        assert tchem.find_whitelists("/x") == "/x"
        monkeypatch.delenv(tchem.ENV_LISTS)
        monkeypatch.setenv(tchem.ENV_TCHEM_LISTS, "/tchem-lists")  # tchem's own variable is the last fallback
        assert tchem.find_whitelists() == "/tchem-lists"


class TestSecondOpinion:
    def test_agree(self, tmp_path):
        b = _fake(tmp_path, _report(), _emit())
        op = tchem.second_opinion("10xv2", "R1.fq.gz", "R2.fq.gz", binary=b)
        assert op.status == "agree" and op.call is not None and op.call.technology == V2

    def test_disagree_with_x(self, tmp_path):
        b = _fake(tmp_path, _report(layout="10x 3' v3"), _emit(V3))
        assert tchem.second_opinion("10xv2", "R1", None, binary=b).status == "disagree"

    def test_no_match_and_emit_exit_3(self, tmp_path):
        b = _fake(tmp_path, _report("no_match"))
        assert tchem.second_opinion("10xv2", "R1", None, binary=b).status == "no_call"
        b = _fake(tmp_path, _report(), "", emit_rc=3)  # a match tchem has no kb command for
        assert tchem.second_opinion("10xv2", "R1", None, binary=b).status == "no_call"

    def test_missing_binary_by_mode(self, monkeypatch):
        monkeypatch.delenv(tchem.ENV_BIN, raising=False)
        monkeypatch.setenv("PATH", "")
        assert tchem.second_opinion("10xv2", "R1", None, mode="auto").status == "unavailable"
        assert tchem.second_opinion("10xv2", "R1", None, mode="off").status == "unavailable"
        with pytest.raises(tchem.TchemError, match="no tchem binary"):
            tchem.second_opinion("10xv2", "R1", None, mode="require")

    @pytest.mark.parametrize("kw", [{"json_rc": 1}, {"junk": True}])
    def test_failure_is_unavailable_in_auto_and_an_error_in_require(self, tmp_path, kw):
        b = _fake(tmp_path, _report(), _emit(), **kw)
        assert tchem.second_opinion("10xv2", "R1", None, binary=b).status == "unavailable"
        with pytest.raises(tchem.TchemError):
            tchem.second_opinion("10xv2", "R1", None, mode="require", binary=b)

    def test_unknown_mode(self):
        with pytest.raises(ValueError):
            tchem.second_opinion("10xv2", "R1", None, mode="always")


class TestEnforce:
    def test_disagree_and_no_call_stop(self):
        for status in ("disagree", "no_call"):
            with pytest.raises(chemistry.ChemistryError, match="tchem second opinion failed"):
                tchem.enforce([tchem.Opinion("agree", "ok"), tchem.Opinion(status, "why")])

    def test_force_and_non_blocking_statuses_pass(self):
        tchem.enforce([tchem.Opinion("disagree", "why")], force=True)
        tchem.enforce([tchem.Opinion(s, "x") for s in ("agree", "unavailable", "info")])


class TestMenuWiring:
    """``_resolve_chemistry`` runs the second opinion after the built-in resolve."""

    def _args(self, tmp_path, tech="10xv2", **kw):
        r1 = tmp_path / "s_R1.fastq.gz"
        r1.write_bytes(b"")
        base = dict(sample1=str(r1), sample2=None, whitelist=None, technology=tech,
                    force_technology=False, tchem="auto", tchem_bin=None, tchem_whitelists=None)
        base.update(kw)
        return argparse.Namespace(**base)

    @pytest.fixture(autouse=True)
    def _builtin_says_v2(self, monkeypatch):
        det = chemistry.Detection("10xv2", "bundled on-list", "97.3% on the 10xv2 list", 26, None, 10, 100)
        monkeypatch.setattr(chemistry, "detect", lambda samples, wl: [det] * len(samples))

    def test_agreement_is_recorded_in_the_manifest_block(self, tmp_path):
        args = self._args(tmp_path, tchem_bin=_fake(tmp_path, _report(), _emit()))
        blocks = menu._resolve_chemistry(args)
        (block,) = blocks.values()
        assert block["tchem"]["status"] == "agree" and args.technology == "10xv2"

    def test_disagreement_stops_the_run_unless_forced(self, tmp_path):
        b = _fake(tmp_path, _report(layout="10x 3' v3"), _emit(V3))
        with pytest.raises(SystemExit):
            menu._resolve_chemistry(self._args(tmp_path, tchem_bin=b))
        blocks = menu._resolve_chemistry(self._args(tmp_path, tchem_bin=b, force_technology=True))
        assert next(iter(blocks.values()))["tchem"]["status"] == "disagree"

    def test_off_and_missing_never_stop_or_add_a_block(self, tmp_path, monkeypatch):
        monkeypatch.delenv(tchem.ENV_BIN, raising=False)
        monkeypatch.setenv("PATH", "")
        blocks = menu._resolve_chemistry(self._args(tmp_path, tchem="off"))
        assert "tchem" not in next(iter(blocks.values()))
        blocks = menu._resolve_chemistry(self._args(tmp_path))
        assert next(iter(blocks.values()))["tchem"]["status"] == "unavailable"

    def test_require_without_a_binary_stops(self, tmp_path, monkeypatch):
        monkeypatch.delenv(tchem.ENV_BIN, raising=False)
        monkeypatch.setenv("PATH", "")
        with pytest.raises(SystemExit):
            menu._resolve_chemistry(self._args(tmp_path, tchem="require"))


@pytest.mark.skipif(
    not (os.environ.get("VIRALSCAN_TCHEM") and os.environ.get("VIRALSCAN_TCHEM_TEST_R1")),
    reason="needs the real tchem binary: set VIRALSCAN_TCHEM, VIRALSCAN_TCHEM_WHITELISTS and "
    "VIRALSCAN_TCHEM_TEST_R1/_R2/_X (a real library and the -x it should agree with)",
)
def test_real_binary_agrees_on_a_real_library():
    op = tchem.second_opinion(
        os.environ.get("VIRALSCAN_TCHEM_TEST_X"),
        os.environ["VIRALSCAN_TCHEM_TEST_R1"],
        os.environ.get("VIRALSCAN_TCHEM_TEST_R2"),
        mode="require",
    )
    assert op.status == "agree", op.reason
