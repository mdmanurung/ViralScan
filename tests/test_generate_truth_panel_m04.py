"""Truth-panel generator, M04-GEN: CB-preserving hops (Arm X) and the D1 real-background mode.

Everything is written under ``tmp_path``. The default-output guard is the AUD-18 hash test in
``tests/test_generate_truth_panel_real.py``; the tests here pin that the new options change nothing
unless asked for, and that what they add is what the scorer needs.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import generate_truth_panel as gt  # noqa: E402

from tests._truth_panel_real import write_real_inputs  # noqa: E402
from tests.test_generate_truth_panel import read_tsv, tree_hashes  # noqa: E402

HOP_ARGS = [
    "--generation-seed",
    "1",
    "--error-rate",
    "0",
    "--chemistries",
    "10xv3",
    "--datasets",
    "synthetic_factorial",
    "synthetic_mixed_host_virus",
    "--abundance",
    "100",
]


def hop_make(out: Path, *extra: str) -> dict:
    return gt.generate(gt.parse_args([*HOP_ARGS, "--out", str(out), *extra]))


@pytest.fixture(scope="module")
def hop_panel(tmp_path_factory):
    root = tmp_path_factory.mktemp("hop")
    _, _, onlist = write_real_inputs(root / "in")
    flags = ["--hop-preserve-cb", "--onlist", str(onlist)]
    return root, hop_make(root / "panel", *flags)


def test_default_and_explicit_no_flag_are_byte_identical_and_add_no_columns(tmp_path):
    hop_make(tmp_path / "a")
    hop_make(tmp_path / "b", "--no-hop-preserve-cb")
    assert tree_hashes(tmp_path / "a") == tree_hashes(tmp_path / "b")
    header = (tmp_path / "a" / "truth" / "truth_manifest.tsv").read_text().split("\n", 1)[0]
    assert not set(gt.HOP_EXTRA) & set(header.split("\t"))
    assert "hop_mode" not in (tmp_path / "a" / "run_manifest.json").read_text()


def test_cb_preserving_hop_keeps_the_donor_barcode_and_umi(hop_panel):
    root, manifest = hop_panel
    assert manifest["hop_mode"] == "cb-preserving"
    rows = read_tsv(root / "panel" / "truth" / "truth_manifest.tsv")
    by_molecule = {
        r["planted_molecule_id"]: r
        for r in rows
        if r["origin"] == "cell" and r["is_planted_molecule"] == "true"
    }
    hopped = [r for r in rows if r["origin"] == "hopped"]
    assert hopped
    for r in hopped:
        donor = by_molecule[r["hop_donor_molecule_id"]]
        assert (r["cell_barcode"], r["umi"]) == (donor["cell_barcode"], donor["umi"])
        assert (
            r["hop_donor_sample_id"] == donor["biological_sample_id"] != r["biological_sample_id"]
        )
        assert r["hop_cb_preserved"] == "true"
        assert (r["planted_virus_id"], r["planted_gene_id"]) == (
            donor["planted_virus_id"],
            donor["planted_gene_id"],
        )
    assert all(not r[c] for r in rows if r["origin"] != "hopped" for c in gt.HOP_EXTRA), (
        "only hopped rows carry hop provenance"
    )


def test_hop_landing_agrees_with_the_recipients_barcode_pool(hop_panel):
    root, _ = hop_panel
    rows = [
        r
        for r in read_tsv(root / "panel" / "truth" / "truth_manifest.tsv")
        if r["origin"] == "hopped"
    ]
    role = {}
    for sid in {r["biological_sample_id"] for r in rows}:
        pool = read_tsv(root / "panel" / sid / "barcode_pool.tsv")
        role[sid] = {p["cell_barcode"]: p["role"] for p in pool}
    seen = set()
    for r in rows:
        expected = role[r["biological_sample_id"]].get(r["cell_barcode"], "other")
        assert r["hop_landing"] == expected
        seen.add(expected)
    assert "other" in seen and seen & {"cell", "empty"}, seen  # the on-list makes coincidence real


def test_a_cell_landing_does_not_make_the_recipient_cell_infected(hop_panel):
    root, _ = hop_panel
    rows = read_tsv(root / "panel" / "truth" / "truth_manifest.tsv")
    infected = {
        (r["biological_sample_id"], r["cell_barcode"])
        for r in rows
        if r["origin"] == "cell" and r["is_planted_molecule"] == "true"
    }
    for r in rows:
        if r["origin"] == "hopped":
            key = (r["biological_sample_id"], r["cell_barcode"])
            assert (r["is_infected_cell"] == "true") == (key in infected)
