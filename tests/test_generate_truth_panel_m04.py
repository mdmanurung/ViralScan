"""Truth-panel generator, M04-GEN: CB-preserving hops (Arm X) and the D1 real-background mode.

Everything is written under ``tmp_path``. The default-output guard is the AUD-18 hash test in
``tests/test_generate_truth_panel_real.py``; the tests here pin that the new options change nothing
unless asked for, and that what they add is what the scorer needs.
"""

from __future__ import annotations

import gzip
import io
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import generate_truth_panel as gt  # noqa: E402
import plant_anello_10x as planter  # noqa: E402

from tests._truth_panel_real import write_background, write_real_inputs  # noqa: E402
from tests.test_generate_truth_panel import fastq_records, read_tsv, tree_hashes  # noqa: E402

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


# ── D1: a real background library, planted arm and unplanted twin ─────────────

D1_ARGS = ["--generation-seed", "1", "--error-rate", "0", "--n-cells", "20", "--n-empty", "8"]
D1_ARGS += ["--abundance", "10", "100"]


def d1_make(out: Path, background, *extra: str) -> dict:
    args = [*D1_ARGS, "--background-r1", str(background[0]), "--background-r2", str(background[1])]
    return gt.generate(gt.parse_args([*args, "--out", str(out), *extra]))


def raw_fastq(path: Path) -> list[tuple[str, str, str]]:
    with gzip.open(path, "rt") as handle:
        lines = handle.read().split("\n")
    return [(lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines) - 1, 4)]


@pytest.fixture(scope="module")
def d1(tmp_path_factory):
    root = tmp_path_factory.mktemp("d1")
    background = write_background(root / "in")
    manifest = d1_make(root / "panel", background)
    twins = {sid: d["twin_sample_id"] for sid, d in manifest["samples"].items()}
    return root, background, manifest, {k: v for k, v in twins.items() if v}


def test_d1_writes_each_planted_sample_with_an_unplanted_twin(d1):
    root, _, manifest, twins = d1
    # one planted sample per abundance x homology cell (Q121), one background library
    assert len(twins) == 4
    assert len(manifest["samples"]) == 8
    for planted, twin in twins.items():
        detail, twin_detail = manifest["samples"][planted], manifest["samples"][twin]
        assert twin_detail["twin_of"] == planted and not twin_detail["twin_sample_id"]
        assert detail["background_library_id"] == twin_detail["background_library_id"]
        assert detail["partition"] == twin_detail["partition"] == "training"
    assert manifest["holdout_samples"] == [] and manifest["unsplittable_strata"] == []


def test_d1_background_is_verbatim_and_identical_in_both_arms(d1):
    root, background, _, twins = d1
    original = [
        (seq1, qual1, seq2, qual2)
        for (_, seq1, qual1), (_, seq2, qual2) in zip(
            raw_fastq(background[0]), raw_fastq(background[1])
        )
    ]
    n = len(original)
    for planted, twin in twins.items():
        arms = {}
        for sid in (planted, twin):
            r1 = raw_fastq(root / "panel" / sid / f"{sid}_R1.fastq.gz")
            r2 = raw_fastq(root / "panel" / sid / f"{sid}_R2.fastq.gz")
            arms[sid] = (r1, r2)
        t1, t2 = arms[twin]
        p1, p2 = arms[planted]
        assert len(t1) == len(t2) == n < len(p1) == len(p2)
        assert [(a[1], a[2], b[1], b[2]) for a, b in zip(t1, t2)] == original
        assert p1[:n] == t1 and p2[:n] == t2  # same reads, same names, same qualities


def test_d1_truth_marks_background_unknown_and_plants_only_into_real_cells(d1):
    root, background, manifest, twins = d1
    rows = read_tsv(root / "panel" / "truth" / "truth_manifest.tsv")
    cells = {
        r["cell_barcode"]
        for r in read_tsv(root / "panel" / next(iter(twins)) / "barcode_pool.tsv")
        if r["role"] == "cell"
    }
    real = {(cb, umi) for cb, umi, *_ in gt.load_background(*background, "10xv3", 0, 20, 8).reads}
    for r in rows:
        if r["biological_sample_id"] in twins.values():
            assert r["origin"] == "unknown" and r["is_planted_molecule"] == "false"
            assert r["viral_abundance_level"] == "none" and not r["twin_sample_id"]
            assert r["is_infected_cell"] == "false"
        elif r["biological_sample_id"] in twins:
            assert r["twin_sample_id"] == twins[r["biological_sample_id"]]
            if r["is_planted_molecule"] == "true":
                assert r["origin"] in ("cell", "ambient")  # no hop: one library is one hop unit
                assert (r["cell_barcode"], r["umi"]) not in real
                assert r["origin"] != "cell" or r["cell_barcode"] in cells
            elif r["host_virus_homology_level"] == "repeat_homology":
                assert r["origin"] in ("unknown", "challenge")
            else:
                assert r["origin"] == "unknown"
    planted = [r for r in rows if r["is_planted_molecule"] == "true"]
    assert planted and {r["background_library_id"] for r in rows} == {
        manifest["inputs"]["background"]["library_id"]
    }
    assert not any(r["origin"].startswith(("artefact", "hopped")) for r in rows)
    assert {r["low_complexity_level"] for r in rows} == {"0"}


def test_d1_crosses_abundance_by_homology_per_library_with_challenge_reads_in_planted_only(d1):
    """Q121 (SCI-05 round 13 M8): every library carries each abundance x homology cell and its twin."""
    root, _, manifest, twins = d1
    samples = manifest["samples"]
    level = lambda sid: tuple(samples[sid]["stratum"].split("|")[:2])  # noqa: E731  (abundance, homology)
    cells = {level(sid) for sid in twins}
    assert cells == {(a, h) for a in ("10", "100") for h in ("none", "repeat_homology")}
    rows = read_tsv(root / "panel" / "truth" / "truth_manifest.tsv")
    challenged = {r["biological_sample_id"] for r in rows if r["origin"] == "challenge"}
    homology_planted = {sid for sid in twins if level(sid)[1] == "repeat_homology"}
    assert (
        challenged == homology_planted
    )  # none-level samples and every twin carry no challenge read
    homology_rows = {
        r["biological_sample_id"]
        for r in read_tsv(root / "panel" / "truth" / "host_homology_manifest.tsv")
    }
    assert homology_rows == homology_planted


def test_d1_cells_are_the_most_frequent_barcodes_and_the_inputs_are_hashed_without_paths(d1):
    _, background, manifest, _ = d1
    entry = manifest["inputs"]["background"]
    assert entry["n_pairs"] == 20 * 40 + 8 * 3 + 30 and entry["min_cell_reads"] == 40
    assert entry["r1"]["sha256"] == gt.sha256_file(background[0])
    assert str(background[0].parent) not in json.dumps(manifest)


def test_d1_is_reproducible_and_leaves_no_trace_without_the_option(d1, tmp_path):
    root, background, _, _ = d1
    d1_make(tmp_path / "again", background)
    assert tree_hashes(tmp_path / "again") == tree_hashes(root / "panel")
    hop_make(tmp_path / "plain")
    assert "inputs" not in json.loads((tmp_path / "plain" / "run_manifest.json").read_text())


def test_d1_with_hop_preserve_cb_has_no_hops(tmp_path):
    background = write_background(tmp_path / "in")
    manifest = d1_make(tmp_path / "panel", background, "--hop-preserve-cb")
    assert all(d.get("realised", {}).get("hopped", 0) == 0 for d in manifest["samples"].values())


def test_d1_defaults_follow_the_background_chemistry(tmp_path):
    ns = gt.parse_args(["--generation-seed", "1", "--out", "x", "--background-r1", "a"])
    assert (ns["datasets"], ns["chemistries"]) == (["synthetic_factorial"], ["10xv3"])
    ns = gt.parse_args(["--generation-seed", "1", "--out", "x"])
    assert len(ns["datasets"]) == 4 and ns["chemistries"] == ["10xv2", "10xv3", "drop-seq"]


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        (["--datasets", "synthetic_host_only"], "synthetic_factorial only"),
        (["--chemistries", "10xv2"], "one chemistry"),
        (["--n-cells", "40", "--n-empty", "20"], "distinct ACGT barcodes"),
        (["--onlist", "x"], "--onlist is for synthetic barcodes"),
    ],
)
def test_d1_refuses_a_mismatched_panel(tmp_path, extra, message):
    background = write_background(tmp_path / "in")
    with pytest.raises(ValueError, match=message):
        d1_make(tmp_path / "panel", background, *extra)


def test_d1_refuses_one_mate_and_a_wrong_r1_length(tmp_path):
    background = write_background(tmp_path / "in")
    with pytest.raises(ValueError, match="given together"):
        gt.generate(
            gt.parse_args(
                [*D1_ARGS, "--background-r1", str(background[0]), "--out", str(tmp_path / "a")]
            )
        )
    short = write_background(tmp_path / "short", r1_len=26)
    with pytest.raises(ValueError, match="needs 28"):
        d1_make(tmp_path / "b", short)


def test_copy_background_still_streams_the_same_block_into_every_handle(tmp_path):
    r1, r2 = write_background(tmp_path / "in")
    first, second = io.StringIO(), io.StringIO()
    n = planter.copy_background(r1, r2, 5, [(first, second), (io.StringIO(), io.StringIO())])
    assert n == 5
    assert first.getvalue() == "".join(
        "".join(rec1) for rec1, _ in planter.iter_background(r1, r2, 5)
    )
    assert len(fastq_records(r1)) > 5
