"""Hand-computed fixtures for the M04-SCORE additions to the v3 scorer (PLAN M04-SCORE).

Hopped/ambient molecules in a negative are scored three ways, truth tiers are validated and never
pooled, D1 manifests are exempt from requested == realised, and a pre-M04 manifest scores exactly as
before. No acceptance threshold appears anywhere: the tests assert counts and ratios only.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import v3_score as score  # noqa: E402

from tests._truth_panel_real import write_background, write_real_inputs  # noqa: E402
from tests.test_generate_truth_panel import read_tsv  # noqa: E402
from tests.test_generate_truth_panel_m04 import d1_make, hop_make  # noqa: E402


def row(origin, cb, umi, virus="HSV-1", **extra):
    """A truth-manifest row with only the fields the scorer reads."""
    planted = origin in ("cell", "hopped", "ambient")
    base = dict.fromkeys(score.required_columns("truth_manifest"), "x") | {
        "origin": origin,
        "cell_barcode": cb,
        "umi": umi,
        "planted_virus_id": virus if planted else "",
        "planted_gene_id": "G1" if planted else "",
        "is_planted_molecule": "true" if planted else "false",
        "is_infected_cell": "false",
        "partition": "training",
    }
    return base | extra


HOP = dict(hop_donor_sample_id="sD", hop_cb_preserved="true", hop_landing="cell")


def lib(sid, called, rows, tier="T1", **extra):
    return dict(sample_id=sid, tier=tier, called=called, truth_rows=rows) | extra


# ── three numbers ────────────────────────────────────────────────────────────


def hop_fixture():
    return [
        # L1 called HSV-1 on two molecules, both hopped with provenance: explained, excluded
        lib(
            "L1",
            {"HSV-1": {("cbA", "u1"), ("cbA", "u2")}},
            [row("hopped", "cbA", "u1", **HOP), row("hopped", "cbA", "u2", **HOP)],
        ),
        # L2 called EBV on a molecule nothing in the manifest explains: a false call that stays
        lib("L2", {"EBV": {("cbB", "u3")}}, [row("unknown", "cbB", "u3")]),
        # L3 clean and uncalled
        lib("L3", {}, [row("unknown", "cbC", "u4")]),
        # L4 holds one hopped molecule (landed in an empty droplet) and was not called: excluded
        lib("L4", {}, [row("hopped", "cbD", "u5", **(HOP | {"hop_landing": "empty"}))]),
    ]


def test_negative_with_hops_reports_three_numbers_with_the_unconditional_as_headline():
    r = score.score_negative_libraries(hop_fixture(), viruses=["EBV", "HSV-1", "KSHV"])
    assert r["headline"] == "unconditional" and r["tier"] == "T1"
    u, e, c = r["unconditional"], r["excluded_with_evidence"], r["contamination_attributable"]
    # (i) L1 and L2 called, 4 complete libraries
    assert (u["n"], u["false_calls"], u["rate"]) == (4, 2, 0.5)
    assert u["per_virus"] == {
        "EBV": dict(false_calls=1, rate=0.25),
        "HSV-1": dict(false_calls=1, rate=0.25),
        "KSHV": dict(false_calls=0, rate=0.0),
    }
    # (ii) L1 and L4 leave; L2 and L3 remain, one false call
    assert (e["n_excluded"], e["n"], e["false_calls"], e["rate"]) == (2, 2, 1, 0.5)
    assert [x["sample_id"] for x in e["excluded"]] == ["L1", "L4"]
    assert e["excluded"][0]["evidence"] == dict(
        n_molecules=2,
        n_hopped=2,
        n_ambient=0,
        n_with_hop_provenance=2,
        n_cb_preserved=2,
        donors=["sD"],
        landing={"cell": 2},
    )
    assert e["excluded"][1]["evidence"]["landing"] == {"empty": 1}
    assert e["per_virus"]["EBV"] == dict(false_calls=1, rate=0.5)
    assert e["per_virus"]["HSV-1"] == dict(false_calls=0, rate=0.0)
    # (iii) one of the two false calls is contamination, over the same 4 libraries as (i)
    assert (c["n"], c["false_calls"], c["rate"]) == (4, 1, 0.25)
    assert (c["residual_false_calls"], c["residual_rate"]) == (1, 0.25)
    assert u["false_calls"] == c["false_calls"] + c["residual_false_calls"]


def test_negative_with_ambient_only():
    libraries = [
        lib(
            "A1",
            {"HSV-1": {("cbA", "u1")}},
            [row("ambient", "cbA", "u1"), row("ambient", "cbA", "u2")],
            tier="T4",
        ),
        lib("A2", {}, [row("unknown", "cbB", "u3")], tier="T4"),
    ]
    r = score.score_negative_libraries(libraries)
    assert r["unconditional"]["rate"] == 0.5
    e = r["excluded_with_evidence"]
    assert (e["n_excluded"], e["n"], e["false_calls"], e["rate"]) == (1, 1, 0, 0.0)
    ev = e["excluded"][0]["evidence"]
    assert (ev["n_ambient"], ev["n_hopped"], ev["n_with_hop_provenance"], ev["landing"]) == (
        2,
        0,
        0,
        {},
    )
    c = r["contamination_attributable"]
    assert (c["false_calls"], c["rate"], c["residual_false_calls"]) == (1, 0.5, 0)


def test_a_call_with_one_unexplained_molecule_stays_a_false_call():
    libraries = [
        lib(
            "L",
            {"HSV-1": {("cbA", "u1"), ("cbA", "u9")}},  # u9 is background, not contamination
            [row("hopped", "cbA", "u1", **HOP), row("unknown", "cbA", "u9")],
        )
    ]
    r = score.score_negative_libraries(libraries)
    assert r["excluded_with_evidence"]["n_excluded"] == 0
    assert r["excluded_with_evidence"]["false_calls"] == 1
    assert r["contamination_attributable"]["false_calls"] == 0
    assert r["contamination_attributable"]["residual_false_calls"] == 1


def test_a_call_with_no_supporting_molecules_is_never_attributed():
    r = score.score_negative_libraries(
        [lib("L", {"HSV-1": set()}, [row("hopped", "c", "u", **HOP)])]
    )
    assert r["unconditional"]["false_calls"] == 1
    assert r["contamination_attributable"]["false_calls"] == 0


def test_incomplete_libraries_are_left_out_and_counted():
    libraries = hop_fixture() + [lib("L5", {"EBV": {("c", "u")}}, [], complete=False)]
    r = score.score_negative_libraries(libraries)
    assert r["n_incomplete"] == 1 and r["unconditional"]["n"] == 4


def test_a_negative_library_with_own_cell_planted_molecules_is_refused():
    with pytest.raises(ValueError, match="cell-origin planted"):
        score.score_negative_libraries([lib("L", {}, [row("cell", "c", "u")])])


def test_the_text_report_prints_all_three_with_the_headline_labelled():
    text = score.format_negative_report(score.score_negative_libraries(hop_fixture()))
    assert text.splitlines()[1] == "(i) HEADLINE unconditional any-virus false call: 2/4 = 0.5000"
    assert "(ii) excluded with evidence (2 libraries excluded): 1/2 = 0.5000" in text
    assert "(iii) contamination-attributable: 1/4 = 0.2500" in text
    assert "residual (not attributable): 1/4" in text


def test_report_has_no_threshold_or_verdict_field():
    r = score.score_negative_libraries(hop_fixture())
    blob = json.dumps(r).lower()
    assert not any(word in blob for word in ("pass", "fail", "threshold", "accept"))


# ── truth tiers ──────────────────────────────────────────────────────────────


def test_tier_pooling_is_refused_and_by_tier_reports_each_separately():
    mixed = [lib("L1", {}, [], tier="T1"), lib("L2", {}, [], tier="T4")]
    with pytest.raises(ValueError, match="refusing to pool"):
        score.score_negative_libraries(mixed)
    by_tier = score.score_negative_libraries_by_tier(mixed)
    assert list(by_tier) == ["T1", "T4"]
    assert all(rep["unconditional"]["n"] == 1 for rep in by_tier.values())


@pytest.mark.parametrize("bad", ["T0", "T6", "t1", "", "1", "tier1", None])
def test_bad_tier_string_is_refused(bad):
    with pytest.raises((ValueError, TypeError)):
        score.check_tier(bad)
    with pytest.raises((ValueError, TypeError)):
        score.score_negative_libraries([lib("L", {}, [], tier=bad)])


def test_vocabulary_is_t1_to_t5_and_t5_is_not_scored():
    assert sorted(score.TIERS) == ["T1", "T2", "T3", "T4", "T5"]
    assert [score.check_tier(t) for t in ("T1", "T2", "T3", "T4")] == ["T1", "T2", "T3", "T4"]
    with pytest.raises(ValueError, match="not scored"):
        score.check_tier("T5")


# ── D1 exemption, unknown origin ─────────────────────────────────────────────

REALISED = dict(
    requested_collision=2,
    collision_same=1,
    collision_disjoint=1,
    requested_ambient=3,
    ambient=3,
    requested_hopped=3,
    hopped=0,
)


def test_d1_is_exempt_from_requested_equals_realised_and_others_are_not():
    assert score.realised_shortfalls({"realised": REALISED}) == [("hopped", 3, 0)]
    d1 = {"realised": REALISED, "background_library_id": "bgabc", "twin_sample_id": "s1"}
    assert score.realised_shortfalls(d1) == []
    assert score.realised_shortfalls({"realised": REALISED | {"hopped": 3}}) == []


def test_unknown_origin_is_accepted_and_never_a_truth_positive(tmp_path):
    columns = score.required_columns("truth_manifest")
    rows = [
        row("unknown", "cbA", "u1"),
        row("cell", "cbB", "u2"),
        row("unknown", "cbC", "u3", is_infected_cell="true"),
    ]
    path = tmp_path / "truth_manifest.tsv"
    header = "\t".join(columns)
    path.write_text(
        header + "\n" + "\n".join("\t".join(r[c] for c in columns) for r in rows) + "\n"
    )
    read = score.read_truth_manifest(path)
    assert score.truth_positive_molecules(read) == {("cbB", "u2", "HSV-1")}
    bad = [
        row("unknown", "cbA", "u1")
        | {"is_planted_molecule": "true", "planted_virus_id": "V", "planted_gene_id": "G"}
    ]
    path.write_text(header + "\n" + "\t".join(bad[0][c] for c in columns) + "\n")
    with pytest.raises(ValueError, match="unknown-origin background"):
        score.read_truth_manifest(path)


# ── real generator output: hop provenance and D1 ─────────────────────────────


@pytest.fixture(scope="module")
def hop_panel(tmp_path_factory):
    root = tmp_path_factory.mktemp("hop")
    _, _, onlist = write_real_inputs(root / "in")
    hop_make(root / "panel", "--hop-preserve-cb", "--onlist", str(onlist))
    return root / "panel"


def test_generator_hop_manifest_is_read_and_its_hopped_rows_are_contamination(hop_panel):
    path = hop_panel / "truth" / "truth_manifest.tsv"
    rows = score.read_truth_manifest(path)
    assert set(score.HOP_COLUMNS) <= set(rows[0])
    found = score._contamination(rows)
    hopped = [r for r in rows if r["origin"] == "hopped"]
    assert hopped
    assert {k for k, r in found.items() if r["origin"] == "hopped"} == {
        (r["cell_barcode"], r["umi"], r["planted_virus_id"]) for r in hopped
    }
    assert all(r["hop_cb_preserved"] == "true" and r["hop_landing"] for r in hopped)
    # a library holding such molecules is evidence-excluded only when its call is all contamination
    sid = hopped[0]["biological_sample_id"]
    mine = [r for r in rows if r["biological_sample_id"] == sid and r["origin"] == "hopped"]
    called = {r["planted_virus_id"]: set() for r in mine}
    for r in mine:
        called[r["planted_virus_id"]].add((r["cell_barcode"], r["umi"]))
    own = [r for r in rows if r["biological_sample_id"] == sid and r["origin"] != "cell"]
    report = score.score_negative_libraries([lib(sid, called, own, tier="T1")])
    assert report["excluded_with_evidence"]["n_excluded"] == 1
    assert report["contamination_attributable"]["false_calls"] == 1


def test_generator_d1_manifest_scores_with_unknown_background_and_is_exempt(tmp_path):
    background = write_background(tmp_path / "in")
    manifest = d1_make(tmp_path / "panel", background)
    rows = score.read_truth_manifest(tmp_path / "panel" / "truth" / "truth_manifest.tsv")
    assert any(r["origin"] == "unknown" for r in rows)
    truth = score.truth_positive_molecules(rows)
    assert truth and not {
        (r["cell_barcode"], r["umi"]) for r in rows if r["origin"] == "unknown"
    } & {(cb, umi) for cb, umi, _ in truth}
    for detail in manifest["samples"].values():
        assert score.realised_shortfalls(detail) == []
    assert any(
        d["realised"]["requested_hopped"] > 0
        for d in manifest["samples"].values()
        if d.get("realised")
    )
    # the twin (all background, nothing planted) is a clean negative library
    twin = next(d for d in manifest["samples"].values() if d.get("twin_of"))
    twin_rows = [
        r
        for r in rows
        if r["biological_sample_id"] == next(s for s, d in manifest["samples"].items() if d is twin)
    ]
    assert twin_rows and not score.truth_positive_molecules(twin_rows)
    assert read_tsv(tmp_path / "panel" / "truth" / "truth_manifest.tsv")


# ── a pre-M04 manifest scores exactly as before ──────────────────────────────

PRE_M04_ROWS = [
    '{"origin": "cell", "cell_barcode": "cbA", "umi": "u1", "planted_virus_id": "HSV-1"}',
    '{"origin": "hopped", "cell_barcode": "cbB", "umi": "u2", "planted_virus_id": "HSV-1"}',
    '{"origin": "ambient", "cell_barcode": "cbC", "umi": "u3", "planted_virus_id": "EBV"}',
]


def test_pre_m04_manifest_reads_and_scores_byte_identically(tmp_path):
    """No hop columns, no unknown origin: the reader returns the rows untouched and a negative
    scored from them gives one fixed JSON string (the expected value was written by hand)."""
    columns = score.required_columns("truth_manifest")
    assert not set(columns) & set(score.HOP_COLUMNS)
    lines = ["\t".join(columns)]
    for text in PRE_M04_ROWS:
        extra = json.loads(text)
        planted = {
            "is_planted_molecule": "true",
            "planted_gene_id": "G1",
            "is_infected_cell": "false",
            "partition": "training",
        }
        full = dict.fromkeys(columns, "x") | planted | extra
        lines.append("\t".join(full[c] for c in columns))
    path = tmp_path / "truth_manifest.tsv"
    path.write_text("\n".join(lines) + "\n")
    rows = score.read_truth_manifest(path)
    assert all(list(r) == columns for r in rows)  # exactly the contract columns, nothing added
    assert [r["origin"] for r in rows] == ["cell", "hopped", "ambient"]
    assert score.truth_positive_molecules(rows) == {
        ("cbA", "u1", "HSV-1"),
        ("cbB", "u2", "HSV-1"),
        ("cbC", "u3", "EBV"),
    }
    negative = [
        lib("N1", {"HSV-1": {("cbB", "u2")}}, rows[1:2], tier="T1"),
        lib("N2", {}, [], tier="T1"),
    ]
    out = json.dumps(score.score_negative_libraries(negative), sort_keys=True)
    expected = (
        '{"contamination_attributable": {"false_calls": 1, "n": 2, "rate": 0.5, "residual_false_calls": 0, "residual_rate": 0.0},'
        ' "excluded_with_evidence": {"excluded": [{"evidence": {"donors": [], "landing": {}, "n_ambient": 0, "n_cb_preserved": 0, "n_hopped": 1, "n_molecules": 1, "n_with_hop_provenance": 0}, "sample_id": "N1"}],'
        ' "false_calls": 0, "n": 1, "n_excluded": 1, "per_virus": {"HSV-1": {"false_calls": 0, "rate": 0.0}}, "rate": 0.0},'
        ' "headline": "unconditional", "n_incomplete": 0, "tier": "T1",'
        ' "unconditional": {"false_calls": 1, "n": 2, "per_virus": {"HSV-1": {"false_calls": 1, "rate": 0.5}}, "rate": 0.5}}'
    )
    assert out == expected


def test_default_generator_manifest_reads_as_it_did_before_m04(tmp_path):
    """Digest of the reader's rows on a default (no hop, no D1) panel, taken with the pre-M04 reader."""
    import hashlib

    hop_make(tmp_path / "panel")
    rows = score.read_truth_manifest(tmp_path / "panel" / "truth" / "truth_manifest.tsv")
    assert len(rows) == 15200
    digest = hashlib.sha256(json.dumps(rows).encode()).hexdigest()
    assert digest == "d9dee8e83a7f4e0bca171377e909dc2b327ba9edf4c2fe1f642a75acef26dd8d"
