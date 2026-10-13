"""SENS-07: per-dataset planted-read sensitivity (planter geometry, plant, score, control receipt)."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import random
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from tests._exact_control import ACCESSION, write_exact_control
from viralscan.sensitivity import minimum_molecules_for_detection, verify_exact_control
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import plant_anello_10x as planter  # noqa: E402


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sens = _load("sens07_planted_sensitivity")


def _seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


# ── the planter takes a dataset's own geometry ───────────────────────────────────────────────
def test_umi_and_polya_read_follow_the_requested_geometry():
    assert len(planter.umi(random.Random(0), 8)) == 8
    read, body = planter.polya_read(_seq(300, 1), 100, 120, random.Random(0), read_len=50)
    assert (len(read), body) == (50, 20) and read[20:] == "A" * 30


def _plant(tmp_path, *, circular=True, taken=None, error_rate=0.0, levels=(100,), **geometry):
    genome = _seq(200, 2)
    selection = [
        {"genome": "G", "genus": "fam", "identity_band": "exact", "length": 200,
         "orf1_start": 0, "last_cds_end": 200}
    ]  # fmt: skip
    r1, r2 = io.StringIO(), io.StringIO()
    out = io.StringIO()
    import csv as _csv

    n = planter.plant_reads(
        selection, {"G": {"seq": genome}}, ["AAAACCCCGGGGTTTT"], random.Random(3), r1, r2,
        _csv.writer(out, delimiter="\t"), levels=levels, windows=("uniform",),
        error_rate=error_rate, taken=taken, circular=circular, **geometry,
    )  # fmt: skip
    fq = lambda h: h.getvalue().splitlines()[1::4]  # noqa: E731
    return genome, n, fq(r1), fq(r2)


def test_plant_reads_uses_the_given_read_and_umi_length(tmp_path):
    _, n, r1, r2 = _plant(tmp_path, read_len=60, umi_len=8)
    assert n == 100
    assert {len(x) for x in r1} == {16 + 8}
    assert {len(x) for x in r2} == {60}


def test_planted_pairs_never_reuse_a_taken_barcode_umi(tmp_path):
    # 4-nt UMIs (about 250 usable) with 100 molecules in one cell collide without a guard.
    pre = {("AAAACCCCGGGGTTTT", planter.umi(random.Random(i), 4)) for i in range(60)}
    seed_pairs = set(pre)
    _, _, r1, _ = _plant(tmp_path, umi_len=4, read_len=60, taken=pre)
    pairs = [(x[:16], x[16:]) for x in r1]
    assert len(set(pairs)) == len(pairs), "planted molecules collided with each other"
    assert not set(pairs) & seed_pairs, "a planted molecule reused a background (CB, UMI)"


def test_linear_genome_reads_do_not_wrap(tmp_path):
    genome, _, _, r2 = _plant(tmp_path, circular=False, read_len=60)
    assert all(read in genome for read in r2)
    genome, _, _, r2 = _plant(tmp_path, circular=True, read_len=60, levels=(400,))
    assert any(read not in genome for read in r2), "the circular default wraps the genome end"


# ── plant ────────────────────────────────────────────────────────────────────────────────────
READ_LEN = 60


def _background(path: Path, n_cells=40, per_cell=15, r1_len=28) -> set[tuple[str, str]]:
    rng = random.Random(7)
    cells = ["".join(rng.choice("ACGT") for _ in range(16)) for _ in range(n_cells)]
    seen = set()
    with (
        gzip.open(path / "bgR1.fastq.gz", "wt") as f1,
        gzip.open(path / "bgR2.fastq.gz", "wt") as f2,
    ):
        for i in range(n_cells * per_cell):
            cb = cells[i % n_cells]
            umi = "".join(rng.choice("ACGT") for _ in range(r1_len - 16))
            seen.add((cb, umi))
            f1.write(f"@bg{i} 1\n{cb + umi}\n+\n{'I' * r1_len}\n")
            f2.write(f"@bg{i} 2\n{_seq(READ_LEN, i)}\n+\n{'I' * READ_LEN}\n")
    return seen


def _inputs(tmp_path: Path, **bg):
    bg_pairs = _background(tmp_path, **bg)
    (tmp_path / "exact.fa").write_text(f">{ACCESSION}\n{_seq(900, 11)}\n>E2\n{_seq(700, 12)}\n")
    (tmp_path / "held.fa").write_text(f">H1\n{_seq(800, 13)}\n")
    (tmp_path / "families.tsv").write_text(
        f"accession\tfamily\n{ACCESSION}\tFamA\nE2\tFamB\nH1\tFamA\n"
    )
    return bg_pairs


def _args(tmp_path: Path, **over) -> argparse.Namespace:
    base = dict(
        out_dir=tmp_path / "out", technology="10xv3",
        background_r1=tmp_path / "bgR1.fastq.gz", background_r2=tmp_path / "bgR2.fastq.gz",
        exact_fasta=tmp_path / "exact.fa", heldout_fasta=tmp_path / "held.fa",
        families=tmp_path / "families.tsv", paf=None, per_family=1, levels=[5, 10],
        cells_per_condition=3, n_background=0, error_rate=0.0, seed=1, gene_spans=None,
        control_genes=[], index=None, reference_manifest=None, identity=None,
        count_layer="X", count_strategy="kb_count",
    )  # fmt: skip
    return argparse.Namespace(**{**base, **over})


def test_plant_writes_the_datasets_own_geometry_and_disjoint_conditions(tmp_path):
    bg_pairs = _inputs(tmp_path)
    manifest = sens.plant(_args(tmp_path))
    out = tmp_path / "out"
    assert (manifest["read_length"], manifest["cb_len"], manifest["umi_len"]) == (READ_LEN, 16, 12)
    assert manifest["status"] == "exploratory" and manifest["thresholds"] is None
    with gzip.open(out / "planted" / "R2.fastq.gz", "rt") as fh:
        reads = fh.read().splitlines()[1::4]
    with gzip.open(out / "unplanted" / "R2.fastq.gz", "rt") as fh:
        n_twin = len(fh.read().splitlines()[1::4])
    assert n_twin == 600 and len(reads) == 600 + manifest["n_planted_reads"]
    assert {len(r) for r in reads} == {READ_LEN}  # not the 90 nt ANDET-09e default
    truth = list(csv.DictReader((out / "truth.tsv").read_text().splitlines(), delimiter="\t"))
    assert len(truth) == 3 * (5 + 10)  # two exact genomes + one held-out, two levels each
    assert not {(r["cb"], r["umi"]) for r in truth} & bg_pairs
    blocks: dict[tuple, set] = {}
    for r in truth:
        blocks.setdefault((r["genome"], r["level"]), set()).add(r["cb"])
    cells = [cb for block in blocks.values() for cb in block]
    assert len(cells) == len(set(cells)), "two conditions share a cell"
    assert {r["identity_band"] for r in truth if r["genome"] == "H1"} == {"heldout:unknown"}
    assert {r["identity_band"] for r in truth if r["genome"] == "E2"} == {"exact"}


def test_plant_refuses_what_it_cannot_do_honestly(tmp_path):
    _inputs(tmp_path)
    (tmp_path / "held.fa").write_text(f">{ACCESSION}\n{_seq(800, 13)}\n")  # also in the panel
    with pytest.raises(ValueError, match="absent from the panel"):
        sens.plant(_args(tmp_path))
    _inputs(tmp_path)
    with pytest.raises(ValueError, match="barcodes"):
        sens.plant(_args(tmp_path, cells_per_condition=30))
    _inputs(tmp_path, r1_len=26)
    with pytest.raises(ValueError, match="barcode \\+ UMI"):
        sens.plant(_args(tmp_path))


def _fastq(path: Path, mate: int) -> list[list[str]]:
    with gzip.open(path / f"R{mate}.fastq.gz", "rt") as fh:
        lines = fh.read().splitlines()
    return [lines[i : i + 4] for i in range(0, len(lines), 4)]


def test_trim_r1_to_geometry_trims_both_twins_identically_and_declares_it(tmp_path):
    _inputs(tmp_path, r1_len=40)  # 28 nt of barcode + UMI plus 12 nt of extra R1
    with pytest.raises(ValueError, match="--trim-r1-to-geometry"):
        sens.plant(_args(tmp_path))
    manifest = sens.plant(_args(tmp_path, trim_r1_to_geometry=True))
    out = tmp_path / "out"
    assert manifest["r1_trim"] == {
        "trim_r1_to_geometry": True,
        "original_length": 40,
        "trimmed_length": 28,
        "applied_to": ["planted", "unplanted"],
    }
    assert manifest["read_changes"] == [
        "R1 trimmed from 40 to 28 nt (barcode + UMI), identically in planted and unplanted"
    ]
    on_disk = json.loads((out / "plant_manifest.json").read_text())
    assert on_disk["r1_trim"] == manifest["r1_trim"]
    planted, twin = _fastq(out / "planted", 1), _fastq(out / "unplanted", 1)
    assert planted[: len(twin)] == twin  # the twin is exactly the planted library's background
    assert {len(r[1]) for r in planted} == {len(r[3]) for r in planted} == {28}
    with gzip.open(tmp_path / "bgR1.fastq.gz", "rt") as fh:
        original = [ln for i, ln in enumerate(fh.read().splitlines()) if i % 4 == 1]
    assert [r[1] for r in twin] == [seq[:28] for seq in original]
    # R2 is not touched, and a planted pair never reuses a background (CB, UMI) after trimming
    assert _fastq(out / "planted", 2)[: len(twin)] == _fastq(out / "unplanted", 2)
    bg = {(r[1][:16], r[1][16:]) for r in twin}
    truth = csv.DictReader((out / "truth.tsv").read_text().splitlines(), delimiter="\t")
    assert not {(r["cb"], r["umi"]) for r in truth} & bg


def test_untrimmed_run_records_no_read_change(tmp_path):
    _inputs(tmp_path)
    manifest = sens.plant(_args(tmp_path))
    assert manifest["read_changes"] == [] and manifest["r1_trim"]["trim_r1_to_geometry"] is False
    assert manifest["r1_trim"]["original_length"] == manifest["r1_trim"]["trimmed_length"] == 28


def test_trim_never_rescues_an_r1_shorter_than_barcode_plus_umi(tmp_path):
    _inputs(tmp_path, r1_len=26)  # 10xv3 needs 28
    with pytest.raises(ValueError, match="too short to trim"):
        sens.plant(_args(tmp_path, trim_r1_to_geometry=True))


def test_per_family_selection_does_not_depend_on_the_other_families():
    accs = [f"{f}{i}" for f in "ABC" for i in range(8)]
    fams = {acc: acc[0] for acc in accs}

    def picks(keep, seed=7, family_seed=7):
        sub = {a: f for a, f in fams.items() if f in keep}
        rows = sens.select_genomes(
            {"exact": {a: "ACGT" for a in sub}, "heldout": {}},
            sub, 2, {}, random.Random(seed), family_seed,
        )  # fmt: skip
        return {f: sorted(r["genome"] for r in rows if r["genus"] == f) for f in keep}

    full = picks({"A", "B", "C"})
    assert picks({"A"})["A"] == full["A"] and picks({"B", "C"})["C"] == full["C"]
    assert picks({"A", "B", "C"}, seed=999) == full  # the shared RNG is not consulted
    assert full != picks({"A", "B", "C"}, family_seed=8)  # but the seed matters
    assert all(len(v) == 2 for v in full.values())


def test_default_selection_still_uses_the_shared_rng():
    # round 2b reproducibility: with no family_seed the draw is the original one-RNG sorted draw
    exact = {f"{f}{i}": "ACGT" for f in "AB" for i in range(6)}
    fams = {acc: acc[0] for acc in exact}
    rows = sens.select_genomes({"exact": exact}, fams, 2, {}, random.Random(3))
    rng = random.Random(3)
    want = sorted(rng.sample([f"A{i}" for i in range(6)], 2)) + sorted(
        rng.sample([f"B{i}" for i in range(6)], 2)
    )
    assert [r["genome"] for r in rows] == want


def test_plant_per_family_mode_is_recorded_and_stable_across_family_sets(tmp_path):
    _inputs(tmp_path)
    one = sens.plant(_args(tmp_path, selection_seed_mode="per-family"))
    assert one["selection_seed_mode"] == "per-family"
    sel1 = (tmp_path / "out" / "selection.tsv").read_text().splitlines()
    fam = (tmp_path / "families.tsv").read_text().replace("E2\tFamB\n", "")
    (tmp_path / "families.tsv").write_text(fam)  # FamB no longer planted
    sens.plant(_args(tmp_path, selection_seed_mode="per-family", out_dir=tmp_path / "out2"))
    rows2 = (tmp_path / "out2" / "selection.tsv").read_text().splitlines()
    assert len(rows2) < len(sel1) and all(line in sel1 for line in rows2)


# ── score ────────────────────────────────────────────────────────────────────────────────────
def _truth(genome, family, band, level, cells):
    return [
        {"genome": genome, "genus": family, "identity_band": band, "level": str(level), "cb": cb}
        for cb in cells
        for _ in range(level // len(cells))
    ]


def test_score_recovers_per_condition_capture_and_lod():
    truth = _truth("E", "FamA", "exact", 10, ["c1", "c2"]) + _truth(
        "E", "FamA", "exact", 100, ["c3", "c4"]
    )
    counts = {("c1", "gA"): 4, ("c2", "gA"): 4, ("c3", "gA"): 90, ("c4", "gB"): 7, ("c9", "gA"): 50}
    twin = {("c1", "gA"): 1}
    rows, summary = sens.score(truth, counts, {"gA": "FamA", "gB": "FamB"}, twin, threshold=1)
    low, high = rows
    assert (low["planted"], low["recovered_raw"], low["background_in_twin"], low["recovered"]) == (
        10,
        8,
        1,
        7,
    )
    assert (high["recovered"], high["other_viral"]) == (
        90,
        7,
    )  # the cell outside the block is ignored
    (fam,) = summary
    assert fam["capture"] == pytest.approx(97 / 110) and fam["capture_status"] == "measured"
    assert fam["lod95_true_molecules"] == pytest.approx(
        minimum_molecules_for_detection(0.95, 1) / (97 / 110)
    )
    assert fam["capture_ci_low"] < fam["capture"] < fam["capture_ci_high"]
    assert fam["lowest_level_all_detected"] == 10


def test_score_labels_failed_and_over_recovered_controls():
    truth = _truth("E", "FamA", "exact", 10, ["c1"]) + _truth(
        "H", "FamB", "heldout:90-95", 10, ["c2"]
    )
    rows, summary = sens.score(truth, {("c2", "gB"): 12}, {"gA": "FamA", "gB": "FamB"})
    by = {r["family"]: r for r in summary}
    assert by["FamA"]["capture_status"] == "failed" and by["FamA"]["lod95_true_molecules"] == ""
    assert (
        by["FamB"]["capture_status"] == "over_recovered"
        and by["FamB"]["lod95_true_molecules"] == ""
    )
    assert by["FamB"]["mode"] == "heldout" and by["FamA"]["lowest_level_all_detected"] == ""


def test_load_counts_reads_a_kb_dir_and_an_h5ad_layer(tmp_path):
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    kb = tmp_path / "kb"
    kb.mkdir()
    (kb / "cells_x_genes.barcodes.txt").write_text("c1\nc2\n")
    (kb / "cells_x_genes.genes.txt").write_text("gA\ngB\n")
    (kb / "cells_x_genes.mtx").write_text(
        "%%MatrixMarket matrix coordinate real general\n%\n2 2 2\n1 2 3\n2 1 5\n"
    )
    assert sens.load_counts(kb) == {("c1", "gB"): 3.0, ("c2", "gA"): 5.0}
    adata = ad.AnnData(sp.csr_matrix(np.array([[1.0, 0], [0, 2.0]])))
    adata.obs_names, adata.var_names = ["c1", "c2"], ["gA", "gB"]
    adata.layers["resolved"] = sp.csr_matrix(np.array([[0, 4.0], [0, 0]]))
    adata.write_h5ad(tmp_path / "a.h5ad")
    assert sens.load_counts(tmp_path / "a.h5ad") == {("c1", "gA"): 1.0, ("c2", "gB"): 2.0}
    assert sens.load_counts(tmp_path / "a.h5ad", "resolved") == {("c1", "gB"): 4.0}


# ── controls: the emitted receipt verifies ───────────────────────────────────────────────────
@pytest.mark.parametrize("match_sequence", [True, False])
def test_control_receipt_binds_the_planted_sequence(tmp_path, match_sequence):
    _inputs(tmp_path)
    gene = f"{ACCESSION}_ctl"
    index, _, _ = write_exact_control(tmp_path / "ref", gene)
    manifest = tmp_path / "ref" / "reference_manifest.json"
    if match_sequence:
        sequence = planter.read_fasta(tmp_path / "exact.fa")[ACCESSION]
        reference = json.loads(manifest.read_text())
        reference["sequences"][0].update(
            sha256=hashlib.sha256(sequence.encode()).hexdigest(), length=len(sequence)
        )
        reference["fasta_sha256"] = sens.sha256_file(tmp_path / "exact.fa")
        manifest.write_text(json.dumps(reference))
    table = VirusIdentityTable(
        (
            GeneIdentity(
                gene, ACCESSION, "catalogued", True, virus_key="taxid:123", virus_name="Target"
            ),
        )
    )
    identity = table.write_tsv(tmp_path / "virus_identity.tsv")
    (tmp_path / "spans.tsv").write_text(
        f"virus_id\tgene_id\tstart\tend\tstrand\n{ACCESSION}\t{gene}\t101\t400\t+\n"
    )
    args = _args(
        tmp_path, gene_spans=tmp_path / "spans.tsv", control_genes=[gene], index=index,
        reference_manifest=manifest, identity=identity,
    )  # fmt: skip
    if not match_sequence:
        with pytest.raises(ValueError, match="planted exact sequence .* does not match"):
            sens.plant(args)
        assert not list((tmp_path / "out").rglob("*.receipt.json"))
        return
    manifest_out = sens.plant(args)
    (row,) = manifest_out["controls"]
    assert (
        row["expected_molecules"] == 15 and "--positive-control-scope exact_sequence" in row["args"]
    )
    selection = list(
        csv.DictReader(
            (tmp_path / "out" / "selection.tsv").read_text().splitlines(), delimiter="\t"
        )
    )
    assert [r["genome"] for r in selection if r["genome"].startswith(ACCESSION)] == [
        f"{ACCESSION}|{gene}"
    ]
    config = SimpleNamespace(
        positive_control_receipt=row["receipt"], positive_control_gene=gene,
        positive_control_virus_key="taxid:123", positive_control_expected_molecules=15.0,
        multimapping=False, multimap_method="x", index=str(index),
    )  # fmt: skip
    assert verify_exact_control(config, table, "X")["verification_status"] == "verified"


@pytest.mark.parametrize("group", ["genus", "virus_name", "virus_key"])
def test_score_refuses_groups_not_recorded_in_planting_truth(tmp_path, group):
    argv = [
        "score",
        "--out-dir",
        str(tmp_path),
        "--truth",
        "truth.tsv",
        "--counts-h5ad",
        "counts.h5ad",
        "--identity",
        "identity.tsv",
    ]
    assert sens.parser().parse_args(argv).group_column == "family"
    with pytest.raises(SystemExit):
        sens.parser().parse_args([*argv, "--group-column", group])
    args = sens.parser().parse_args(argv)
    args.group_column = group
    with pytest.raises(ValueError, match="truth records family labels"):
        sens.run_score(args)
