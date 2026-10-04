"""The planted 3' reads must carry an untemplated poly-A tail.

Without it every planted read is a clean genome slice, so a poly-A or
homopolymer gate scored on this set kills nothing and looks free — while on real
data it removes exactly the genuine 3'-end reads (Biomni review, 2026-10-04).
These tests pin the tail so that failure mode cannot come back silently.
"""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "plant_anello_10x.py"


@pytest.fixture(scope="module")
def plant():
    spec = importlib.util.spec_from_file_location("plant_anello_10x", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def seq(plant):
    # Deterministic, complex, and long enough that no slice wraps by accident.
    rng = random.Random(0)
    return "".join(rng.choice("ACGT") for _ in range(3000))


def test_a_read_at_the_cleavage_site_is_all_tail(plant, seq):
    read, body_len = plant.polya_read(seq, 1000, 1000, random.Random(1))
    assert body_len == 0
    assert read == "A" * plant.READ_LEN


def test_a_read_crossing_the_site_is_part_body_part_tail(plant, seq):
    cleavage, body = 1000, 30
    read, body_len = plant.polya_read(seq, cleavage - body, cleavage, random.Random(1))
    assert body_len == body
    assert read[:body] == seq[cleavage - body : cleavage]
    assert read[body:] == "A" * (plant.READ_LEN - body)
    assert len(read) == plant.READ_LEN


def test_a_read_far_upstream_has_no_tail(plant, seq):
    read, body_len = plant.polya_read(seq, 500, 1000, random.Random(1))
    assert body_len == plant.READ_LEN
    assert read == seq[500 : 500 + plant.READ_LEN]
    assert not read.endswith("AAAAAAAAAA")


def test_the_tail_is_untemplated_not_a_genome_slice(plant):
    """The point of the fixture: the A-run must not be copied from the genome."""
    # A genome whose post-cleavage sequence is deliberately not poly-A.
    seq = "C" * 1000 + "G" * 2000
    read, body_len = plant.polya_read(seq, 990, 1000, random.Random(1))
    assert body_len == 10
    assert read[:10] == "C" * 10
    assert set(read[10:]) == {"A"}, "tail must be untemplated A, not the genome's G run"


def test_every_3p_planted_read_is_recorded_with_its_body_length(plant, seq, tmp_path):
    """``truth.tsv`` must carry body_len/tail_len so recovery can be scored by it."""
    import csv
    import io

    genomes = {"G1": {"seq": seq, "orf1_start": 600, "last_cds_end": 2600}}
    selection = [{
        "genome": "G1", "genus": "Betatorquevirus", "identity_band": "<85",
        "length": len(seq), "orf1_start": 600, "last_cds_end": 2600,
    }]
    r1, r2, tf = io.StringIO(), io.StringIO(), io.StringIO()
    truth = csv.writer(tf, delimiter="\t")
    plant.LEVELS = (3,)  # keep it small; the shape is what matters
    n = plant.plant_reads(selection, genomes, ["A" * 16], random.Random(7), r1, r2, truth)

    rows = [r for r in csv.reader(io.StringIO(tf.getvalue()), delimiter="\t") if r]
    assert len(rows) == n
    by_window = {}
    for row in rows:
        by_window.setdefault(row[4], []).append((int(row[9]), int(row[10])))

    # 5p / uniform are internal fragments: fully templated, no tail.
    for window in ("5p", "uniform"):
        assert all(b == plant.READ_LEN and t == 0 for b, t in by_window[window])
    # 3p spans the cleavage site, so at least one read must carry a real tail.
    assert any(t > 0 for _, t in by_window["3p"]), "no 3p read carries a poly-A tail"
    assert all(b + t == plant.READ_LEN for b, t in by_window["3p"])
