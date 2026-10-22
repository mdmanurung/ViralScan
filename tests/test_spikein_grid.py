import gzip
import importlib
import subprocess
import sys
from pathlib import Path

import pytest

A = Path(__file__).parent.parent / "analysis" / "spikein_detection_limit"
sys.path.insert(0, str(A))
grid = importlib.import_module("grid")
run_condition = importlib.import_module("run_condition")
score_condition = importlib.import_module("score_condition")


def test_grid_row_count_ids_and_replicates():
    rows = grid.build_conditions()
    spiked = [r for r in rows if r["level"] > 0]
    assert len(spiked) == 5 * 3 * 2 * 2 == 60
    assert len(rows) == 60 + 2 * grid.N_NULL_REPS >= 66
    ids = [r["cond_id"] for r in rows]
    assert len(set(ids)) == len(ids)
    assert not any(
        "_" in i for i in ids
    )  # viralscan derives the sample id before the first underscore
    for bg in grid.BACKGROUNDS:
        nulls = [r for r in rows if r["background"] == bg and r["level"] == 0]
        assert len({r["bg_rep"] for r in nulls}) >= 3
    assert {r["cells"] for r in spiked} == {5, 100}


def _tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


def _kb(arm_dir, mols, evidence=()):
    _tsv(
        arm_dir / "virus_identity.tsv",
        ["gene_id", "genome_accession"],
        [("T1_g1", "T1.1"), ("T1_g2", "T1.1"), ("O1_g", "O1.1"), ("O2_g", "O2.1")],
    )
    with gzip.open(arm_dir / "kb_molecules.tsv.gz", "wt") as fh:
        fh.write("barcode\tgene_id\tunique\tambiguous\n")
        fh.writelines("\t".join(map(str, m)) + "\n" for m in mols)
    _tsv(
        arm_dir / "multimap_evidence.tsv",
        ["gene_id", "viral_molecules_total_est", "evidence_tier"],
        evidence,
    )


COND = {"cond_id": "X-L003-s1-c005", "background": "X", "source": "S", "level": "3", "seed": "1",
        "layout": "c005", "cells": "5", "bg_rep": "1"}  # fmt: skip


def _bundle(tmp_path, spiked=True):
    b = tmp_path / "bundle"
    if spiked:
        _tsv(
            b / "truth.tsv",
            ["name", "cb", "class", "contigs"],
            [
                (f"SPIKE:a:{i}", "AAA" if i < 2 else "BBB", "sense_in_gene", "VIRUS|T1.1")
                for i in range(4)
            ],
        )
    _tsv(
        b / "star/star_spike_fates.tsv",
        ["name", "on_truth", "fate"],
        [(f"SPIKE:a:{i}", 1, "source") for i in range(3)],  # spike 3 is absent
    )
    _tsv(b / "star/star_contig_reads.tsv", ["contig", "n_any", "n_primary"],
         [("VIRUS|T1.1", 7, 7), ("VIRUS|O2.1", 2, 2)])  # fmt: skip
    mols = [("AAA", "T1_g1", 2, 0.5), ("BBB", "T1_g2", 0, 1), ("CCC", "T1_g1", 5, 0),
            ("CCC", "O1_g", 4, 0), ("CCC", "O2_g", 3, 0)]  # fmt: skip
    ev = [("T1_g1", 7.5, "candidate_unique"), ("T1_g2", 1, "candidate_virus_ambiguous"),
          ("O1_g", 4, "candidate_unique")]  # fmt: skip
    _kb(b / "old", mols, ev)
    _kb(b / "flank", [("AAA", "T1_g1", 1, 0)], ev[:1])
    return b


def test_scorer_hand_computed(tmp_path):
    b = _bundle(tmp_path)
    base = tmp_path / "base"
    _kb(base / "old", [("AAA", "T1_g1", 1, 0)])
    _kb(base / "flank", [])
    rows = {r["arm"]: r for r in score_condition.score_bundle(b, COND, base)}
    old = rows["old"]
    assert old["n_spiked"] == 4
    assert (old["kb_unique"], old["kb_ambiguous_allocated"], old["kb_total"]) == (2.0, 1.5, 3.5)
    assert old["evidence_tier"] == "candidate_unique"
    assert old["evidence_truth_total_all_barcodes"] == 8.5
    assert old["kb_total_in_baseline"] == 1.0
    assert old["star_n_source_contig_any"] == 3
    assert (old["call_kb"], old["call_star"], old["call_combined"]) == (True, True, True)
    # floor: O1 has 4 molecules and no STAR read (unsupported); O2 has 3 and 2 STAR reads (a false call)
    assert old["kb_unsupported_molecules"] == 4.0 and old["n_false_calls"] == 1
    assert old["kb_total_viral_molecules"] == 15.5
    fl = rows["flank"]
    assert fl["kb_total"] == 1.0 and fl["call_kb"] is False and fl["call_combined"] is False
    st = rows["star"]
    assert (st["star_n_virus"], st["star_n_absent"], st["call_star"]) == (3, 1, True)
    assert (st["star_nonspike_reads_total"], st["star_contigs_with_reads"]) == (9, 2)
    out = tmp_path / "scores.tsv"
    score_condition.write_scores(list(rows.values()), out)
    got = score_condition.read_tsv(out)
    assert [r["arm"] for r in got] == ["old", "flank", "star"] and got[0]["call_kb"] == "True"


def test_scorer_unspiked_has_floor_and_no_calls(tmp_path):
    b = _bundle(tmp_path, spiked=False)
    (b / "star/star_spike_fates.tsv").write_text("name\ton_truth\tfate\n")
    rows = {r["arm"]: r for r in score_condition.score_bundle(b, COND)}
    old = rows["old"]
    assert old["n_spiked"] == 0 and old["call_combined"] is None
    # no truth set: every accession counts; unsupported = O1 (4); T1 has STAR reads, O2 too -> false calls: T1 (7), O2 (3)
    assert old["kb_unsupported_molecules"] == 4.0 and old["n_false_calls"] == 2


def test_extract_star_lines(tmp_path):
    lines = [
        "SPIKE:a:0\t0\tVIRUS|T1.1\t5\t255\n",
        "SPIKE:a:1\t256\tVIRUS|T2.1\t5\t0\n",  # secondary on the other truth contig
        "SPIKE:a:1\t0\tVIRUS|T1.1\t9\t0\n",
        "SPIKE:a:2\t0\tVIRUS|O1.1\t9\t0\n",
        "bg1\t0\tVIRUS|O1.1\t9\t0\n",
        "bg2\t256\tVIRUS|O1.1\t9\t0\n",
    ]
    run_condition.extract_star_lines(lines, "VIRUS|T1.1,VIRUS|T2.1", tmp_path)
    fates = {r["name"]: r for r in score_condition.read_tsv(tmp_path / "star_spike_fates.tsv")}
    assert fates["SPIKE:a:1"]["on_truth"] == "1" and fates["SPIKE:a:2"]["on_truth"] == "0"
    assert fates["SPIKE:a:2"]["fate"] == "other_virus"
    reads = {r["contig"]: r for r in score_condition.read_tsv(tmp_path / "star_contig_reads.tsv")}
    assert reads == {"VIRUS|O1.1": {"contig": "VIRUS|O1.1", "n_any": "2", "n_primary": "1"}}


def test_dry_run_prints_the_real_commands_and_writes_nothing(tmp_path, capsys):
    cond = {k: str(v) for k, v in next(
        r for r in grid.build_conditions() if r["cond_id"] == "D03-L030-s1-c005").items()}  # fmt: skip
    run_condition.run(cond, dry=True, root=tmp_path)
    out = capsys.readouterr().out
    cmds = [ln for ln in out.splitlines() if ln.startswith("+ ")]
    runs = [c for c in cmds if "viralscan.menu run" in c]
    assert len(runs) == 2 and sum("--star-confirm" in c for c in runs) == 1
    assert (
        "star_confirm" in runs[0] or "--star-confirm" in runs[0]
    )  # STAR arm rides on the first (old) run
    for c in runs:
        for want in ("--host-filter starsolo", "-x 10xv2", "--strand reverse", "--cell-calling external",
                     "--gtf", "D03-L030-s1-c005_R1.fastq.gz", "--no-visual", "--yes"):  # fmt: skip
            assert want in c
        assert (
            "-s1 " not in c and " -gtf" not in c
        )  # old single-dash spellings are rejected by the product
    assert "panel.idx" in runs[0] and "panel_flank.idx" in runs[1]
    assert any("spikein.py spike" in c and "--levels 30" in c and "--cells 5" in c for c in cmds)
    assert any("spikein.py bg" in c and "--seed 1" in c for c in cmds)
    assert any("score_condition.py" in c for c in cmds)
    assert list(tmp_path.iterdir()) == []


def test_dry_run_ebv_uses_v3_geometry_and_truth_set(tmp_path, capsys):
    cond = {k: str(v) for k, v in next(
        r for r in grid.build_conditions() if r["cond_id"] == "R01-L300-s2-c100").items()}  # fmt: skip
    run_condition.run(cond, dry=True, root=tmp_path)
    out = capsys.readouterr().out
    assert "-x 10xv3" in out and "--strand forward" in out and "--max-nh 2" in out
    assert (
        "VIRUS|NC_007605.1,VIRUS|NC_009334.1" in out
        and "--cells 100" in out
        and "--umi-len 12" in out
    )


def test_null_condition_has_no_spike_step(tmp_path, capsys):
    cond = {k: str(v) for k, v in next(
        r for r in grid.build_conditions() if r["cond_id"] == "R01-NULL-r3").items()}  # fmt: skip
    run_condition.run(cond, dry=True, root=tmp_path)
    assert "spikein.py spike" not in capsys.readouterr().out


def test_rmtree_inside_refuses_outside(tmp_path):
    root = tmp_path / "cond"
    (root / "a").mkdir(parents=True)
    run_condition.rmtree_inside(root / "a", root)
    assert not (root / "a").exists()
    for bad in (root, tmp_path, tmp_path / "other"):
        with pytest.raises(SystemExit):
            run_condition.rmtree_inside(bad, root)


def test_sbatch_syntax_and_array_row():
    r = subprocess.run(["bash", "-n", str(A / "grid.sbatch")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "SLURM_ARRAY_TASK_ID" in (A / "grid.sbatch").read_text()
