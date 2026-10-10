"""Hand-computed fixtures for the v3 scorer and partition (PLAN VAL-06/VAL-07). Code only, no outcomes."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import v3_partition as part  # noqa: E402
import v3_score as score  # noqa: E402

# ── precision / recall / F1 ──────────────────────────────────────────────────


def test_prf_hand_counted():
    # truth {1..5}, called {4,5,6,7}: tp=2 fp=2 fn=3 -> P=2/4, R=2/5, F1=2*2/(4+5)
    r = score.precision_recall_f1({1, 2, 3, 4, 5}, [4, 5, 6, 7])
    assert (r["tp"], r["fp"], r["fn"]) == (2, 2, 3)
    assert r["precision"] == pytest.approx(0.5) and r["recall"] == pytest.approx(0.4)
    assert r["f1"] == pytest.approx(4 / 9)
    assert r["f1"] == pytest.approx(2 * 0.5 * 0.4 / (0.5 + 0.4))


def test_prf_undefined_is_none_not_zero():
    nothing_called = score.precision_recall_f1({1, 2}, [])
    assert nothing_called["precision"] is None and nothing_called["recall"] == 0
    assert nothing_called["f1"] == 0
    nothing_planted = score.precision_recall_f1([], {9})
    assert nothing_planted["recall"] is None and nothing_planted["precision"] == 0
    both_empty = score.precision_recall_f1([], [])
    assert both_empty["precision"] is both_empty["recall"] is both_empty["f1"] is None


def test_prf_perfect_and_duplicates():
    assert score.precision_recall_f1({1, 2}, [1, 1, 2, 2])["f1"] == 1.0  # a molecule counts once


# ── AUPRC ────────────────────────────────────────────────────────────────────


def test_auprc_hand_computed():
    # thresholds .9: P=1,R=.5 ; .8: P=.5,R=.5 ; .7: P=2/3,R=1 -> 0.5*1 + 0*.5 + 0.5*2/3
    assert score.auprc([1, 0, 1, 0], [0.9, 0.8, 0.7, 0.1]) == pytest.approx(5 / 6)


def test_auprc_perfect_worst_and_ties():
    assert score.auprc([1, 1, 0, 0], [0.9, 0.8, 0.2, 0.1]) == pytest.approx(1.0)
    # both positives ranked last: R jumps at P = 1/3 then 2/4
    assert score.auprc([0, 0, 1, 1], [0.9, 0.8, 0.2, 0.1]) == pytest.approx(
        0.5 * (1 / 3) + 0.5 * 0.5
    )
    # a tied pair enters as one threshold: TP=1, seen=2 -> P=.5, R=1 -> AP=.5, in any input order
    assert score.auprc([1, 0], [0.5, 0.5]) == pytest.approx(0.5)
    assert score.auprc([0, 1], [0.5, 0.5]) == pytest.approx(0.5)


def test_auprc_no_positive_is_undefined_and_length_checked():
    assert score.auprc([0, 0], [0.4, 0.2]) is None
    with pytest.raises(ValueError):
        score.auprc([1], [0.1, 0.2])


# ── negative false calls, sibling confusion, host-virus allocation ───────────


def test_negative_false_call_rate_excludes_incomplete_rows():
    called = {"a": False, "b": True, "c": False, "failed": True}
    r = score.negative_false_call_rate(called, ["a", "b", "c"])  # 'failed' did not complete
    assert (r["n"], r["false_calls"]) == (3, 1) and r["rate"] == pytest.approx(1 / 3)
    assert score.negative_false_call_rate(called, [])["rate"] is None
    with pytest.raises(ValueError, match="without a call status"):
        score.negative_false_call_rate({"a": False}, ["a", "zz"])


def test_sibling_confusion_denominator_is_planted_target_molecules():
    planted = {1: "EBV1", 2: "EBV1", 3: "EBV1", 4: "EBV1", 5: "EBV2"}
    assigned = {1: "EBV1", 2: "EBV2", 3: "EBV2", 4: None, 5: "EBV1"}
    r = score.sibling_confusion(planted, assigned, target="EBV1", sibling="EBV2")
    assert (r["n_target"], r["n_to_sibling"]) == (4, 2) and r["rate"] == pytest.approx(0.5)
    assert score.sibling_confusion(planted, assigned, "HSV1", "HSV2")["rate"] is None


def test_host_virus_allocation_error_reports_both_directions():
    r = score.host_virus_allocation_error(
        truth_host={1, 2, 3, 4},
        truth_viral={10, 11, 12},
        assigned_viral={1, 10, 11},
        assigned_host={2, 3, 12},
    )
    assert r["host_to_viral"] == pytest.approx(1 / 4)
    assert r["viral_to_host"] == pytest.approx(1 / 3)
    assert (r["n_host"], r["n_viral"]) == (4, 3)


# ── partition ────────────────────────────────────────────────────────────────


def test_apportion_largest_remainder_hand_computed():
    # N=12, seats=round(3.6)=4, quotas 4/3 each: floors 1,1,1; the spare seat goes to 'a' (tie, ASCII)
    assert part.apportion({"a": 4, "b": 4, "c": 4}, 0.3) == {"a": 2, "b": 1, "c": 1}
    # N=13, seats=round(3.9)=4, quotas 3.077 / 0.923: floors 3,0; spare seat to the larger remainder 'b'
    assert part.apportion({"a": 10, "b": 3}, 0.3) == {"a": 3, "b": 1}
    # round half up, not banker's: 0.5 * 5 = 2.5 -> 3 seats, not 2
    assert sum(part.apportion({"a": 5}, 0.5).values()) == 3


@pytest.mark.parametrize("fraction", [0.1, 0.3, 0.5, 0.77])
def test_apportion_seats_sum_to_the_global_total_and_fit_strata(fraction):
    sizes = {"a": 4, "b": 7, "c": 4, "d": 9, "e": 1}
    seats = part.apportion(sizes, fraction)
    assert sum(seats.values()) == int(fraction * sum(sizes.values()) + 0.5)
    assert all(0 <= seats[k] <= sizes[k] for k in sizes)


def test_split_takes_holdout_in_ascii_order_and_is_order_independent():
    members = {"x": ["s3", "s1", "s2", "s4"], "y": ["s9", "s8"]}
    labels = part.split(members, 0.3)  # N=6, seats=round(1.8)=2: quotas x=1.33 y=0.67 -> x 1, y 1
    assert labels["s1"] == part.HOLDOUT and labels["s2"] == part.TRAINING
    assert labels["s8"] == part.HOLDOUT and labels["s9"] == part.TRAINING
    assert part.split({"y": ["s8", "s9"], "x": ["s4", "s2", "s1", "s3"]}, 0.3) == labels


def test_split_rejects_bad_inputs():
    with pytest.raises(ValueError):
        part.split({"x": ["s1"], "y": ["s1"]}, 0.3)  # one sample in two strata
    for bad in (0, 1, -0.2, 1.5):
        with pytest.raises(ValueError):
            part.apportion({"a": 3}, bad)
    with pytest.raises(ValueError):
        part.apportion({}, 0.3)


def test_sample_id_is_opaque_and_stable():
    a = part.sample_id(7, "ds", "abund03|homol0|10xv3", 1)
    assert a == part.sample_id(7, "ds", "abund03|homol0|10xv3", 1)
    assert a != part.sample_id(7, "ds", "abund03|homol0|10xv3", 2) and a != part.sample_id(
        8, "ds", "abund03|homol0|10xv3", 1
    )
    assert len(a) == 13 and a[0] == "s" and "abund" not in a


# ── truth manifest validation ────────────────────────────────────────────────


def manifest(tmp_path, mutate=lambda row: row, drop=None):
    columns = score.required_columns("truth_manifest")
    row = dict.fromkeys(columns, "x") | {
        "origin": "cell",
        "planted_virus_id": "EBV1",
        "planted_gene_id": "BZLF1",
        "is_infected_cell": "true",
        "is_planted_molecule": "true",
        "partition": "training",
    }
    row = mutate(row)
    header = [c for c in columns if c != drop]
    path = tmp_path / "truth_manifest.tsv"
    path.write_text("\t".join(header) + "\n" + "\t".join(row[c] for c in header) + "\n")
    return path


def test_valid_manifest_is_read(tmp_path):
    rows = score.read_truth_manifest(manifest(tmp_path))
    assert len(rows) == 1 and rows[0]["planted_virus_id"] == "EBV1"


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda r: r | {"planted_virus_id": ""}, "without planted_virus_id"),
        (lambda r: r | {"origin": "artefact:polyG"}, "artefact row labelled as planted"),
        (lambda r: r | {"partition": "validation"}, "not training/holdout"),
        (lambda r: r | {"is_planted_molecule": "maybe"}, "not a boolean"),
    ],
)
def test_malformed_manifest_rows_are_rejected(tmp_path, mutation, message):
    with pytest.raises(ValueError, match=message):
        score.read_truth_manifest(manifest(tmp_path, mutation))


def test_missing_contract_column_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="lacks columns .*background_library_id"):
        score.read_truth_manifest(manifest(tmp_path, drop="background_library_id"))


def test_artefact_row_that_is_not_planted_is_accepted(tmp_path):
    path = manifest(
        tmp_path,
        lambda r: (
            r
            | {
                "origin": "artefact:polyG",
                "planted_virus_id": "",
                "planted_gene_id": "",
                "is_planted_molecule": "false",
            }
        ),
    )
    assert score.read_truth_manifest(path)[0]["origin"] == "artefact:polyG"
