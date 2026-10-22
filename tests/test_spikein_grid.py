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
        ["gene_id", "genome_accession", "viral"],
        [
            ("T1_g1", "T1.1", "true"),
            ("T1_g2", "T1.1", "true"),
            ("O1_g", "O1.1", "true"),
            ("O2_g", "O2.1", "true"),
            ("H_g", "ENST1.1", "false"),
        ],
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


# ------------------------------------------------------------------ resume / idempotence / gates (review fixes)
def _cond(cid):
    return {
        k: str(v)
        for k, v in next(r for r in grid.build_conditions() if r["cond_id"] == cid).items()
    }


class _Calls(list):
    def __init__(self):
        super().__init__()
        self.state = {"star_reads": "5"}


@pytest.fixture
def calls(monkeypatch):
    """Record every ctx.sh command; a score command writes a scores.tsv the way score_condition would."""
    log = _Calls()
    state = log.state

    def fake_sh(self, cmd, log_=None, cwd=None, **kw):
        log.append([str(x) for x in cmd])
        if any("score_condition.py" in str(x) for x in cmd):
            b = self.d / "bundle"
            b.mkdir(parents=True, exist_ok=True)
            cols = "arm\tstar_n_source_contig_any"
            (b / "scores.tsv").write_text(f"{cols}\nold\t\nflank\t\nstar\t{state['star_reads']}\n")

    monkeypatch.setattr(run_condition.Ctx, "sh", fake_sh)
    return log


def test_done_marker_skips_everything(tmp_path, calls):
    d = tmp_path / "cond" / "R01-NULL-r1"
    d.mkdir(parents=True)
    (d / "DONE").write_text("{}")
    run_condition.run(_cond("R01-NULL-r1"), root=tmp_path)
    assert calls == []


def test_stage_markers_skip_arms_and_run_only_scoring(tmp_path, calls):
    d = tmp_path / "cond" / "R01-NULL-r1"
    d.mkdir(parents=True)
    for arm in grid.ARMS:
        (d / f"stage_{arm}.done").write_text("x")
    run_condition.run(_cond("R01-NULL-r1"), root=tmp_path)
    assert len(calls) == 1 and "score_condition.py" in calls[0][2] and (d / "DONE").exists()
    # --baseline-bundle is deliberately absent: the baseline is subtracted in aggregate
    assert "--baseline-bundle" not in calls[0]


def _spiked_ready(tmp_path, level=100):
    cid = f"D03-L{level:03d}-s1-c005"
    d = tmp_path / "cond" / cid
    for arm in grid.ARMS:
        d.mkdir(parents=True, exist_ok=True)
        (d / f"stage_{arm}.done").write_text("x")
    (d / "spikes").mkdir()
    (d / "bundle").mkdir()
    stem = d / "spikes" / f"spike_L{level}_s1"
    stem.with_suffix(".truth.tsv").write_text("name\tcb\nSPIKE:a:0\tAAA\n")
    stem.with_suffix(".json").write_text("{}")
    return cid, d


def test_zero_star_reads_gate_fails_loudly_at_level_100(tmp_path, calls):
    cid, d = _spiked_ready(tmp_path)
    calls.state["star_reads"] = "0"
    with pytest.raises(SystemExit):
        run_condition.run(_cond(cid), root=tmp_path)
    assert not (d / "DONE").exists() and (d / "spikes").exists()


def test_nonzero_star_reads_finishes_and_cleans_spikes(tmp_path, calls):
    cid, d = _spiked_ready(tmp_path)
    run_condition.run(_cond(cid), root=tmp_path)
    assert (d / "DONE").exists() and not (d / "spikes").exists()
    assert (d / "bundle" / "truth.tsv").exists()


def test_zero_star_reads_allowed_below_gate_level(tmp_path, calls):
    cid, d = _spiked_ready(tmp_path, level=30)
    calls.state["star_reads"] = "0"
    run_condition.run(_cond(cid), root=tmp_path)
    assert (d / "DONE").exists()


def _stage_spike_ctx(tmp_path, monkeypatch, called, truth_cb="AAA"):
    ctx = run_condition.Ctx(_cond("D03-L003-s1-c005"), False, tmp_path)
    ctx.bgd.mkdir(parents=True)
    (ctx.bgd / "called_cells.txt").write_text(called)
    monkeypatch.setattr(run_condition.Ctx, "onlist", lambda self: "onlist.txt")

    def fake_sh(self, cmd, log=None, cwd=None):
        t = self.d / "spikes" / "spike_L3_s1.truth.tsv"
        t.parent.mkdir(parents=True, exist_ok=True)
        t.write_text(f"name\tcb\nSPIKE:a:0\t{truth_cb}\n")

    monkeypatch.setattr(run_condition.Ctx, "sh", fake_sh)
    return ctx


def test_spike_barcodes_must_be_in_called_list(tmp_path, monkeypatch):
    ctx = _stage_spike_ctx(tmp_path, monkeypatch, "CCC\nDDD\n")
    with pytest.raises(SystemExit):
        run_condition.stage_spike(ctx)
    assert not (ctx.d / "spikes" / "SPIKES_OK").exists()


def test_stage_spike_ok_then_skipped(tmp_path, monkeypatch):
    ctx = _stage_spike_ctx(tmp_path, monkeypatch, "AAA\nDDD\n")
    run_condition.stage_spike(ctx)
    assert (ctx.d / "spikes" / "SPIKES_OK").exists()
    monkeypatch.setattr(run_condition.Ctx, "sh", lambda *a, **k: pytest.fail("re-ran spike"))
    run_condition.stage_spike(ctx)


def _gz(path, n):
    with gzip.open(path, "wt") as fh:
        fh.writelines(f"@r{i}\nACGT\n+\nIIII\n" for i in range(n))


def test_count_pairs_cached_atomic_and_stage_bg_subsamples(tmp_path, monkeypatch):
    r1 = tmp_path / "r1.fastq.gz"
    _gz(r1, 10)
    ctx = run_condition.Ctx(_cond("D03-NULL-r1"), False, tmp_path)
    assert run_condition.count_pairs(ctx, r1) == 10
    cache = tmp_path / "bg" / "D03.total_pairs"
    assert cache.read_text() == "10" and not list((tmp_path / "bg").glob("*.part"))
    r1.unlink()  # a cache hit must not touch the FASTQ again
    assert run_condition.count_pairs(ctx, r1) == 10

    ctx.bg = {**ctx.bg, "r1": str(r1)}
    monkeypatch.setattr(grid, "BG_TARGET_PAIRS", 5)
    seen = []

    def fake_sh(self, cmd, log=None, cwd=None):
        seen.append([str(x) for x in cmd])
        for m in ("R1", "R2"):
            (self.bgd / f"tmp_{m}.fastq.gz").write_text("x")

    monkeypatch.setattr(run_condition.Ctx, "sh", fake_sh)
    monkeypatch.setattr(run_condition.Ctx, "onlist", lambda self: "onlist.txt")
    spikein = importlib.import_module("spikein")
    monkeypatch.setattr(spikein, "load_whitelist", lambda p: set())
    monkeypatch.setattr(spikein, "top_barcodes", lambda *a, **k: ["AAA", "CCC"])
    run_condition.stage_bg(ctx)
    assert seen[0][seen[0].index("--frac") + 1] == "0.50000000"  # subsampled, not 1.0
    assert (ctx.bgd / "READY").exists() and (
        ctx.bgd / "called_cells.txt"
    ).read_text() == "AAA\nCCC\n"
    run_condition.stage_bg(ctx)
    assert len(seen) == 1  # READY: second call builds nothing


def _h5ad(path, var_names, gene_ids, obs):
    import anndata
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    m = sp.csr_matrix(np.array([[2, 0], [0, 1]], dtype=float)[: len(obs)])
    ad = anndata.AnnData(
        X=m, obs=pd.DataFrame(index=obs), var=pd.DataFrame({"gene_id": gene_ids}, index=var_names)
    )
    ad.layers["counts_unique_viral"] = m
    ad.layers["counts_ambiguous_allocated"] = sp.csr_matrix(m.shape)
    path.parent.mkdir(parents=True, exist_ok=True)
    ad.write_h5ad(path)


def _sample(tmp_path, var_names, obs):
    s = tmp_path / "sample"
    _h5ad(s / "kb-python/counts_unfiltered/adata_multimap.h5ad", var_names, ["V1", "V2"], obs)
    _tsv(
        s / "results/virus_identity.tsv",
        ["gene_id", "genome_accession", "viral"],
        [("V1", "T.1", "true")],
    )
    (tmp_path / "dest").mkdir(exist_ok=True)
    return s, tmp_path / "dest"


def test_extract_kb_matches_gene_id_column_and_strips_barcode_suffix(tmp_path):
    pytest.importorskip("anndata")
    s, dest = _sample(tmp_path, ["geneA", "geneB"], ["AAA-1", "CCC-1"])  # names are NOT gene ids
    run_condition.extract_kb(s, dest, {"AAA"}, 100)
    with gzip.open(dest / "kb_molecules.tsv.gz", "rt") as fh:
        rows = fh.read().splitlines()
    assert rows == ["barcode\tgene_id\tunique\tambiguous", "AAA\tV1\t2\t0"]


def test_extract_kb_fails_on_no_viral_column_or_no_spike_barcode(tmp_path):
    pytest.importorskip("anndata")
    s, dest = _sample(tmp_path, ["geneA", "geneB"], ["AAA", "CCC"])
    _tsv(
        s / "results/virus_identity.tsv",
        ["gene_id", "genome_accession", "viral"],
        [("ZZ", "T.1", "true")],
    )
    with pytest.raises(SystemExit, match="no viral gene"):
        run_condition.extract_kb(s, dest)
    _tsv(
        s / "results/virus_identity.tsv",
        ["gene_id", "genome_accession", "viral"],
        [("V1", "T.1", "true")],
    )
    with pytest.raises(SystemExit, match="spike barcodes"):
        run_condition.extract_kb(s, dest, {"GGG"}, 100)
    run_condition.extract_kb(s, dest, {"GGG"}, 30)  # below the gate level: allowed


def test_run_arm_skips_a_finished_arm(tmp_path, calls):
    ctx = run_condition.Ctx(_cond("R01-NULL-r1"), False, tmp_path)
    ctx.d.mkdir(parents=True)
    (ctx.d / "stage_old.done").write_text("x")
    run_condition.run_arm(ctx, "old")
    assert calls == []


# ------------------------------------------------------------------------------ PLAN SPIKE-03 scorer fix
def _hosty(arm_dir):
    """virus_identity as the product writes it: viral true/false + role, host ENSG rows, a contaminant."""
    _tsv(
        arm_dir / "virus_identity.tsv",
        ["gene_id", "genome_accession", "viral", "role"],
        [("T1_g1", "T1.1", "true", ""), ("O1_g", "O1.1", "True", ""), ("ENSG1", "ENST1.1", "false", ""),
         ("C_g", "C1.1", "true", "contaminant")],
    )  # fmt: skip


def test_viral_identity_filters_host_and_contaminant_case_insensitive(tmp_path):
    _hosty(tmp_path)
    assert score_condition.viral_identity(tmp_path / "virus_identity.tsv") == {
        "T1_g1": "T1.1", "O1_g": "O1.1"
    }  # fmt: skip
    _tsv(tmp_path / "nocol.tsv", ["gene_id", "genome_accession"], [("a", "b")])
    with pytest.raises(ValueError, match="viral"):
        score_condition.viral_identity(tmp_path / "nocol.tsv")


def test_host_genes_do_not_count_as_viral_molecules(tmp_path):
    b = _bundle(tmp_path)
    mols = [
        ("AAA", "T1_g1", 2, 0.5),
        ("CCC", "O1_g", 4, 0),
        ("CCC", "ENSG1", 30000, 0),
        ("CCC", "C_g", 9, 0),
    ]
    _kb(b / "old", mols)
    _hosty(b / "old")
    _kb(b / "flank", [("AAA", "T1_g1", 1, 0)])
    _hosty(b / "flank")
    old = {r["arm"]: r for r in score_condition.score_bundle(b, COND)}["old"]
    assert old["kb_total"] == 2.5  # truth side unchanged
    assert old["kb_unsupported_molecules"] == 4.0  # O1 only: no host ENSG1, no contaminant C_g
    assert old["kb_total_viral_molecules"] == 6.5


def test_extract_kb_writes_viral_rows_only(tmp_path):
    pytest.importorskip("anndata")
    s, dest = _sample(tmp_path, ["geneA", "geneB"], ["AAA", "CCC"])
    _tsv(
        s / "results/virus_identity.tsv",
        ["gene_id", "genome_accession", "viral", "role"],
        [("V1", "T.1", "true", ""), ("V2", "ENST.1", "false", "")],
    )
    run_condition.extract_kb(s, dest, {"AAA"}, 100)
    with gzip.open(dest / "kb_molecules.tsv.gz", "rt") as fh:
        genes = {ln.split("\t")[1] for ln in fh.read().splitlines()[1:]}
    assert genes == {"V1"}
