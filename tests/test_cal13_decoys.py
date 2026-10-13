"""CAL-13: scripts/cal13_decoys.py (k-let shuffle decoys, stratified empirical p, target-decoy q)."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import itertools
import random
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

_spec = importlib.util.spec_from_file_location(
    "cal13_decoys", Path(__file__).resolve().parents[1] / "scripts" / "cal13_decoys.py"
)
cal = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cal)

_R = random.Random(1)
SEQ = "".join(_R.choice("ACGT") for _ in range(600))


def _kmers(seq: str, k: int) -> Counter:
    return Counter(seq[i : i + k] for i in range(len(seq) - k + 1))


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------ shuffle


@pytest.mark.parametrize("k", [1, 2, 3, 6])
def test_klet_composition_and_ends_preserved_exactly(k):
    for seed in range(5):
        out = cal.klet_shuffle(SEQ, k, random.Random(seed))
        assert len(out) == len(SEQ)
        assert _kmers(out, k) == _kmers(SEQ, k)
        assert out[: k - 1] == SEQ[: k - 1] and out[len(out) - k + 1 :] == SEQ[len(SEQ) - k + 1 :]


def test_klet_shuffle_really_shuffles_and_keeps_non_acgt():
    seq = "ACGNNTRACGTTGCAAGGCTNAACCGGTATACGC"
    outs = {cal.klet_shuffle(seq, 3, random.Random(s)) for s in range(10)}
    assert len(outs) > 1 and seq not in outs
    assert all(_kmers(o, 3) == _kmers(seq, 3) for o in outs)


def test_klet_shuffle_is_uniform_over_admissible_sequences():
    seq, k = "AABABBAAB", 2
    admissible = {
        "".join(t)
        for t in itertools.product("AB", repeat=len(seq))
        if t[0] == seq[0] and t[-1] == seq[-1] and _kmers("".join(t), k) == _kmers(seq, k)
    }
    assert len(admissible) > 3
    rng = random.Random(0)
    n = 6000
    seen = Counter(cal.klet_shuffle(seq, k, rng) for _ in range(n))
    assert set(seen) == admissible
    expected = n / len(admissible)
    assert all(abs(c - expected) < 6 * expected**0.5 for c in seen.values())


def test_klet_shuffle_rejects_short_sequences():
    with pytest.raises(ValueError):
        cal.klet_shuffle("ACG", 6, random.Random(0))


# ------------------------------------------------------------------ make


def _fa(path: Path, records: dict[str, str]) -> Path:
    path.write_text("".join(f">{k} description\n{v}\n" for k, v in records.items()))
    return path


def _make(tmp_path, name, panel, extra=(), seed=7, k=6, copies=1):
    out = tmp_path / name
    argv = ["make", "--panel-fasta", str(panel), "--out-dir", str(out), "--seed", str(seed)]
    argv += ["--k", str(k), "--copies", str(copies)]
    if extra:
        argv += ["--unrelated-host-fasta", *map(str, extra)]
    cal.main(argv)
    return out


def test_make_writes_fasta_gtf_manifest(tmp_path):
    panel = _fa(tmp_path / "p.fa", {"NC_000001.1": SEQ, "NC_000002.1": SEQ[::-1]})
    unrel = _fa(tmp_path / "u.fa", {"NC_PLANT.1": "ACGT" * 30})
    out = _make(tmp_path, "d", panel, [unrel], copies=2)
    manifest = _tsv(out / "decoys_manifest.tsv")
    assert [r["decoy_id"] for r in manifest] == [
        "DECOY_NC_000001.1_1",
        "DECOY_NC_000001.1_2",
        "DECOY_NC_000002.1_1",
        "DECOY_NC_000002.1_2",
        "DECOY_NC_PLANT.1",
    ]
    fasta = {line[1:]: next_ for line, next_ in _pairs((out / "decoys.fa").read_text())}
    for row in manifest:
        seq = fasta[row["decoy_id"]]
        assert int(row["length"]) == len(seq)
        assert row["sha256"] == hashlib.sha256(seq.encode()).hexdigest()
    first = manifest[0]
    assert (first["class"], first["source_id"], first["seed"], first["k"]) == (
        "kmer_shuffle",
        "NC_000001.1",
        "7",
        "6",
    )
    assert _kmers(fasta["DECOY_NC_000001.1_1"], 6) == _kmers(SEQ, 6)
    assert fasta["DECOY_NC_000001.1_1"] != fasta["DECOY_NC_000001.1_2"]
    plant = manifest[-1]
    assert plant["class"] == "unrelated_host" and plant["seed"] == plant["k"] == ""
    assert fasta["DECOY_NC_PLANT.1"] == "ACGT" * 30
    gtf = [ln.split("\t") for ln in (out / "decoys.gtf").read_text().splitlines()]
    assert len(gtf) == 3 * len(manifest)
    assert {g[0] for g in gtf} == set(fasta)
    assert all(g[0].startswith("DECOY_") and f'gene_name "{g[0]}"' in g[8] for g in gtf)
    assert {(g[0], g[3], g[4]) for g in gtf} == {(k, "1", str(len(v))) for k, v in fasta.items()}


def _pairs(text: str):
    lines = text.splitlines()
    cur, seq = None, []
    for ln in lines:
        if ln.startswith(">"):
            if cur:
                yield cur, "".join(seq)
            cur, seq = ln, []
        else:
            seq.append(ln)
    yield cur, "".join(seq)


def test_make_is_deterministic_and_independent_of_record_order(tmp_path):
    a = _fa(tmp_path / "a.fa", {"X1": SEQ, "X2": SEQ[::-1]})
    b = _fa(tmp_path / "b.fa", {"X2": SEQ[::-1], "X1": SEQ})
    o1, o2 = _make(tmp_path, "o1", a), _make(tmp_path, "o2", a)
    assert (o1 / "decoys.fa").read_bytes() == (o2 / "decoys.fa").read_bytes()
    assert (o1 / "decoys_manifest.tsv").read_bytes() == (o2 / "decoys_manifest.tsv").read_bytes()
    o3 = _make(tmp_path, "o3", b)
    by_id = lambda d: {r["decoy_id"]: r["sha256"] for r in _tsv(d / "decoys_manifest.tsv")}  # noqa: E731
    assert by_id(o1) == by_id(o3)
    o4 = _make(tmp_path, "o4", a, seed=8)
    assert by_id(o1) != by_id(o4)


def test_make_refuses_prefixed_and_duplicate_ids(tmp_path):
    with pytest.raises(SystemExit):
        _make(tmp_path, "x", _fa(tmp_path / "p.fa", {"DECOY_NC_1": SEQ}))
    panel = _fa(tmp_path / "q.fa", {"NC_1": SEQ})
    unrel = _fa(tmp_path / "u.fa", {"NC_1_1": SEQ})  # DECOY_NC_1_1 = shuffle id of NC_1
    with pytest.raises(SystemExit):
        _make(tmp_path, "y", panel, [unrel])


# ------------------------------------------------------------------ scoring maths


def _rows(spec):
    """(id, is_decoy, stratum, score) tuples -> score_library rows."""
    return [{"target_id": i, "is_decoy": d, "stratum": s, "score": float(x)} for i, d, s, x in spec]


A_DECOYS = [(f"d{i}", True, "A", x) for i, x in enumerate([0, 0, 1, 0])]
A_TARGETS = [
    ("t5", False, "A", 5),
    ("t3", False, "A", 3),
    ("t0a", False, "A", 0),
    ("t0b", False, "A", 0),
]  # noqa: E501


def test_empirical_p_and_td_q_hand_computed():
    out = {r["target_id"]: r for r in cal.score_library(_rows(A_DECOYS + A_TARGETS))}
    # m = 4 decoys; b = decoys >= score
    assert (out["t5"]["b"], out["t5"]["p_empirical"]) == (0, pytest.approx(1 / 5))
    assert (out["t3"]["b"], out["t3"]["p_empirical"]) == (0, pytest.approx(1 / 5))
    assert (out["t0a"]["b"], out["t0a"]["p_empirical"]) == (4, pytest.approx(1.0))
    # ratio r = 4/4 = 1. FDR(0) = (4+1)/(1*4) -> 1; FDR(3) = (0+1)/(1*2) = .5; FDR(5) = 1/1 -> 1,
    # then running minimum over lower thresholds gives q(5) = .5
    assert out["t0a"]["q_td"] == pytest.approx(1.0)
    assert out["t3"]["q_td"] == pytest.approx(0.5)
    assert out["t5"]["q_td"] == pytest.approx(0.5)
    # BH on p = [.2, .2, 1, 1]: .8, .4, 1, 1 -> monotone .4, .4, 1, 1
    assert out["t5"]["q_bh"] == pytest.approx(0.4) and out["t0b"]["q_bh"] == pytest.approx(1.0)


def test_td_q_uses_decoy_to_target_ratio():
    # 2 decoys per target: expected false targets above s is D/2.
    q = cal.td_qvalues([9.0, 0.0], [0.0, 0.0, 0.0, 0.0])
    # s=0: D=4, T=2, r=2 -> (5)/(2*2) -> capped 1; s=9: D=0, T=1 -> 1/(2*1) = .5
    assert q == [pytest.approx(0.5), pytest.approx(1.0)]


def test_stratification_is_honoured():
    spec = A_DECOYS + [("db1", True, "B", 2), ("db2", True, "B", 2)]
    spec += [("ta", False, "A", 2), ("tb", False, "B", 2)]
    out = {r["target_id"]: r for r in cal.score_library(_rows(spec))}
    assert out["ta"]["n_decoys"] == 4 and out["ta"]["b"] == 0
    assert out["ta"]["p_empirical"] == pytest.approx(1 / 5)
    assert out["tb"]["n_decoys"] == 2 and out["tb"]["b"] == 2  # only B's decoys count
    assert out["tb"]["p_empirical"] == pytest.approx(1.0)
    # a stratum with targets but no decoys is never significant
    lone = cal.score_library(_rows([("t", False, "Z", 9)]))[0]
    assert (lone["n_decoys"], lone["p_empirical"], lone["q_td"]) == (0, 1.0, 1.0)


def test_expected_false_calls_hand_computed():
    spec = (
        A_DECOYS + A_TARGETS + [("db1", True, "B", 2), ("db2", True, "B", 2), ("tb", False, "B", 2)]
    )  # noqa: E501
    res = cal.expected_false_calls(
        _rows(spec), 1.0, rng=np.random.default_rng(0), n_draws=4000, level=0.95
    )
    # A: 4 targets * 1/4 decoys over threshold = 1; B: 1 target * 2/2 = 1
    assert res["expected_false_calls"] == pytest.approx(2.0)
    assert (res["n_targets"], res["n_decoys"], res["n_called"]) == (5, 6, 3)
    assert res["ci_lo"] < 2.0 < res["ci_hi"]
    # no decoy over the threshold: point estimate 0 but the Jeffreys interval is not zero-width
    zero = cal.expected_false_calls(
        _rows(spec), 99.0, rng=np.random.default_rng(0), n_draws=4000, level=0.95
    )
    assert zero["expected_false_calls"] == 0 and zero["ci_hi"] > 0 and zero["n_called"] == 0
    # strata with targets but no decoys are excluded and counted
    only = cal.expected_false_calls(
        _rows([("t", False, "Z", 9)]), 1.0, rng=np.random.default_rng(0), n_draws=10, level=0.9
    )
    assert only["n_targets_in_strata_without_decoys"] == 1 and only["expected_false_calls"] == 0


# ------------------------------------------------------------------ fdr CLI


def _calls(path: Path, libs: dict[str, list[tuple]]) -> Path:
    lines = ["library_id\ttarget_id\tis_decoy\tstratum\tscore"]
    for lib, spec in libs.items():
        lines += [f"{lib}\t{i}\t{str(d).lower()}\t{s}\t{x}" for i, d, s, x in spec]
    path.write_text("\n".join(lines) + "\n")
    return path


def _fdr(tmp_path, calls, *extra):
    out = tmp_path / "s"
    cal.main(
        [
            "fdr",
            "--calls",
            str(calls),
            "--out-dir",
            str(out),
            "--threshold",
            "1",
            "--seed",
            "3",
            *extra,
        ]
    )
    return _tsv(out / "target_pvalues.tsv"), _tsv(out / "library_summary.tsv")


def test_fdr_cli_per_library_and_protocol(tmp_path):
    calls = _calls(
        tmp_path / "c.tsv",
        {"L1": A_DECOYS + A_TARGETS, "L2": [("d", True, "A", 0), ("t", False, "A", 4)]},
    )  # noqa: E501
    proto = tmp_path / "protocol.yaml"
    proto.write_text(
        "false_discovery:\n  expected_false_L2_calls_per_library: 1.5\n  interval: [0, 3]\n"
    )
    targets, libs = _fdr(tmp_path, calls, "--protocol", str(proto))
    assert {r["library_id"] for r in targets} == {"L1", "L2"}
    t = {(r["library_id"], r["target_id"]): r for r in targets}
    assert float(t["L1", "t5"]["p_empirical"]) == pytest.approx(0.2)
    assert float(t["L2", "t"]["p_empirical"]) == pytest.approx(0.5)  # m = 1, b = 0
    s = {r["library_id"]: r for r in libs}
    assert float(s["L1"]["expected_false_calls"]) == pytest.approx(1.0)
    assert (
        s["L1"]["prereg_target"] == "1.5" and s["L1"]["expected_within_prereg_interval"] == "True"
    )
    assert s["L1"]["observed_within_prereg_interval"] == "True"  # 2 calls in [0, 3]
    assert float(s["L2"]["expected_false_calls"]) == 0.0


def test_fdr_cli_tbd_protocol_is_unset_and_nothing_is_baked_in(tmp_path):
    calls = _calls(tmp_path / "c.tsv", {"L1": A_DECOYS + A_TARGETS})
    proto = tmp_path / "protocol.yaml"
    proto.write_text(
        "false_discovery:\n  expected_false_L2_calls_per_library: TBD_PREREGISTERED\n"
        "  interval: TBD_PREREGISTERED\n"
    )
    _, libs = _fdr(tmp_path, calls, "--protocol", str(proto))
    assert libs[0]["prereg_target"] == "NA" and libs[0]["expected_within_prereg_interval"] == "NA"
    with pytest.raises(SystemExit):  # the threshold is not defaulted
        cal.main(["fdr", "--calls", str(calls), "--out-dir", str(tmp_path / "z"), "--seed", "1"])


def test_fdr_cli_rejects_bad_input_and_missing_manifest_decoys(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("library_id\ttarget_id\tis_decoy\tstratum\tscore\nL\tt\tmaybe\tA\t1\n")
    with pytest.raises(SystemExit):
        _fdr(tmp_path, bad)
    dup = _calls(tmp_path / "dup.tsv", {"L": [("t", False, "A", 1), ("t", False, "A", 2)]})
    with pytest.raises(SystemExit):
        _fdr(tmp_path, dup)
    calls = _calls(tmp_path / "c.tsv", {"L1": A_DECOYS + A_TARGETS})
    manifest = tmp_path / "m.tsv"
    manifest.write_text("decoy_id\nd0\nd1\nd2\nd3\nd_missing\n")
    with pytest.raises(SystemExit):
        _fdr(tmp_path, calls, "--manifest", str(manifest))
    manifest.write_text("decoy_id\nd0\nd1\nd2\nd3\n")
    _fdr(tmp_path, calls, "--manifest", str(manifest))
