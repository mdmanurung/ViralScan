# ruff: noqa: B023  (closures below are called within the loop iteration that binds them)
"""B24 strict oracle: candidate gate == frozen baseline gate, classifier runs once per record.

The frozen pre-B24 module is loaded by path (read-only) and compared field by field:
exact log records, exact error class/message/cause, exact report dict (order + values).
The oracle is a verbatim copy of ``scripts/build_reference.py`` at commit 5e551c9 (tests/_frozen_oracles).
"""

import contextlib
import importlib.util
import logging
import math
import random
import sys
from pathlib import Path

import pytest

from viralscan.scripts import build_reference as cand

BASELINE = Path(__file__).resolve().parent / "_frozen_oracles" / "frozen_build_reference_5e551c9.py"
assert BASELINE.is_file(), BASELINE


def _load_baseline():
    spec = importlib.util.spec_from_file_location(
        "viralscan_b24_baseline_build_reference", BASELINE
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


base = _load_baseline()
assert Path(cand.__file__).resolve() != BASELINE


def rnd(seed, n):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def run_len(seed, total, run, b="A"):
    n = total + 30
    left = (n - run) // 2
    return rnd(seed, left) + b * run + rnd(seed + 1, n - run - left)


def exact_low(total, low=20, run=12):
    """Seed-search a record with exactly *low* low-complexity k-mers of *total* (5% boundary)."""
    for seed in range(1000, 5000, 7):
        seq = run_len(seed, total, run)
        c = base.low_complexity_kmer_counts(seq)
        if (
            c["total"] == total
            and sum(c[k] for k in ("pure_homopolymer", "long_run", "few_bases", "tandem")) == low
        ):
            return seq
    raise AssertionError


CLEAN = rnd(1, 300)
CASES = {
    "empty": "",
    "header_only": ">x\n",
    "clean": f">a\n{CLEAN}\n",
    "dup_id_distinct": f">d\n{CLEAN}\n>d\n{rnd(2, 300)}\n",
    "dup_id_last_value": f">d\n{'A' * 80}\n>x\n{rnd(3, 300)}\n>d\n{CLEAN}\n",
    "dup_seq": f">a\n{CLEAN}\n>b\n{CLEAN.lower()}\n",
    "empty_record": f">a\n{CLEAN}\n>b\n",
    "n_and_iupac": f">a\n{rnd(4, 40)}{'N' * 50}{rnd(5, 40)}\n>b\n{rnd(6, 60)}RYKM{rnd(7, 60)}\n",
    "short": f">s\nACGT\n>t\n{rnd(8, 30)}\n>u\n{rnd(9, 31)}\n",
    "pure_limit_2": f">h\n{run_len(11, 2000, 32)}\n",
    "pure_limit_3": f">h\n{run_len(11, 2000, 33)}\n",
    "frac_exact": f">b\n{exact_low(400)}\n",
    "frac_above": f">b\n{exact_low(399)}\n",
    "frac_below": f">b\n{exact_low(401)}\n",
    "ties": f">t1\n{run_len(61, 400, 12)}\n>t2\n{run_len(62, 400, 12)}\n>c\n{CLEAN}\n",
    "eleven_tied": "".join(f">o{i:02d}\n{run_len(100 + 2 * i, 300, 40)}\n" for i in range(11)),
    "twelve_tied": "".join(f">o{i:02d}\n{run_len(100 + 2 * i, 300, 40)}\n" for i in range(12)),
    "twelve_distinct": "".join(
        f">u{i:02d}\n{run_len(200 + 2 * i, 250 + 17 * i, 12 + i)}\n" for i in range(12)
    ),
    "tandem": f">t\n{('CAG' * 30)}\n>c\n{CLEAN}\n",
    "few_bases": f">f\n{('AACC' * 20)}\n>c\n{CLEAN}\n",
}


def capture(module, fasta, mask, caplog):
    caplog.clear()
    with caplog.at_level(logging.DEBUG, logger="viralscan"):
        try:
            module._enforce_low_complexity_gate(fasta, mask=mask, label="T")
            err = None
        except Exception as exc:  # noqa: BLE001
            cause = exc.__cause__
            err = (
                type(exc),
                str(exc),
                None if cause is None else (type(cause), str(cause)),
            )
    return [(r.levelname, r.getMessage()) for r in caplog.records], err


@pytest.mark.parametrize("mask", [True, False])
@pytest.mark.parametrize("name", sorted(CASES))
def test_gate_matches_baseline_exactly(tmp_path, caplog, name, mask):
    fasta = tmp_path / "in.fa"
    fasta.write_text(CASES[name])
    b_logs, b_err = capture(base, fasta, mask, caplog)
    c_logs, c_err = capture(cand, fasta, mask, caplog)
    # error classes differ by module object only for RuntimeError/ValueError builtins: equal
    assert c_logs == b_logs
    assert c_err == b_err


@pytest.mark.parametrize("name", sorted(CASES))
def test_report_and_validate_unchanged(tmp_path, name):
    fasta = tmp_path / "in.fa"
    fasta.write_text(CASES[name])
    assert list(cand.low_complexity_report(fasta).items()) == list(
        base.low_complexity_report(fasta).items()
    )
    recs = base._fasta_records(fasta)
    grid = [None, 0, 0.0, 0.05, 1.0]
    if recs and recs[0][1]:
        c = base.low_complexity_kmer_counts(recs[0][1])
        low = sum(c[k] for k in ("pure_homopolymer", "long_run", "few_bases", "tandem"))
        if c["total"]:
            f = low / c["total"]
            grid += [f, math.nextafter(f, 0.0), math.nextafter(f, 2.0)]
    for frac in grid:
        for pure in (None, 0, 1, 2, 3):

            def go(m):
                try:
                    return (
                        "ok",
                        m.validate_reference_records(
                            fasta, max_low_complexity_fraction=frac, max_pure_homopolymer_kmers=pure
                        ),
                    )
                except Exception as exc:  # noqa: BLE001
                    return (type(exc), str(exc))

            assert go(cand) == go(base)


@pytest.mark.parametrize("name", ["clean", "ties", "twelve_tied", "dup_id_distinct", "n_and_iupac"])
def test_classifier_runs_once_per_record_in_candidate(tmp_path, monkeypatch, name):
    fasta = tmp_path / "in.fa"
    fasta.write_text(CASES[name])
    n_records = len(cand._fasta_records(fasta))
    calls = {"cand": 0, "base": 0}
    for key, mod in (("cand", cand), ("base", base)):
        orig = mod.low_complexity_kmer_counts

        def counted(*a, _k=key, _o=orig, **kw):
            calls[_k] += 1
            return _o(*a, **kw)

        monkeypatch.setattr(mod, "low_complexity_kmer_counts", counted)
        with contextlib.suppress(ValueError, RuntimeError):
            mod._enforce_low_complexity_gate(fasta, mask=False, label="T")
    assert calls["cand"] == n_records
    # baseline: report + pure scan always; verdict scan only when structure passes
    assert calls["base"] in (2 * n_records, 3 * n_records)
