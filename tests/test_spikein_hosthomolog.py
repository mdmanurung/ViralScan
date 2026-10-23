import collections
import gzip
import importlib
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

A = Path(__file__).parent.parent / "analysis" / "spikein_detection_limit"
sys.path.insert(0, str(A))
spikein = importlib.import_module("spikein")
hosthomolog = importlib.import_module("hosthomolog")

T1, T2, HOST, OTHER = "VIRUS|T1", "VIRUS|T2", "HOST|chr1", "VIRUS|OTHER"
DNA = ("ACGTTGCAAC" * 10)[:98]
needs_samtools = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools missing")


def _bam(tmp_path, alns):
    sq = "".join(f"@SQ\tSN:{c}\tLN:5000\n" for c in (T1, T2, HOST, OTHER))
    body = "".join(
        f"{n}\t{fl}\t{c}\t{pos}\t255\t98M\t*\t0\t0\t{DNA}\t{'F' * 98}\tNH:i:{nh}\n"
        for n, fl, c, pos, nh in alns
    )
    sam = tmp_path / "x.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n" + sq + body)
    bam = str(tmp_path / "x.bam")
    subprocess.run(["samtools", "sort", "-o", bam, str(sam)], check=True)
    subprocess.run(["samtools", "index", bam], check=True)
    return bam


ALNS = (
    [(f"u{i}", 0, T1, 200 + i, 1) for i in range(5)]  # NH=1 inside
    + [("m0", 0, T1, 300, 2), ("m0", 256, T2, 300, 2)]  # NH=2 both inside
    + [("h0", 0, T1, 300, 2), ("h0", 256, HOST, 300, 2)]  # host partner
    + [("o0", 0, T1, 310, 2), ("o0", 256, OTHER, 300, 2)]  # other-virus partner
    + [("n3", 0, T1, 300, 30), ("n3", 256, HOST, 300, 30)]  # NH above max
)


@needs_samtools
def test_pool_mode_keeps_exactly_the_dropped_reads(tmp_path, capsys):
    bam = _bam(tmp_path, ALNS)
    g = {T1: [(100, 4000, "+")], T2: [(100, 4000, "+")]}
    kw = dict(max_nh=20)
    default = spikein.source_pool(bam, f"{T1},{T2}", g, "forward", 90, 1, **kw)
    capsys.readouterr()
    hh = spikein.source_pool(
        bam, f"{T1},{T2}", g, "forward", 90, 1, pool_mode="outside_truth", **kw
    )
    err = capsys.readouterr().err
    assert len(default) == 6 and len(hh) == 2  # u0-4 + m0 vs h0, o0
    assert collections.Counter(r["pos"] for r in hh) == {300: 1, 310: 1}
    assert {r["nh"] for r in hh} == {2} and '"inside_truth_set": 6' in err
    assert '"nh_above_max": 1' in err
    # default mode unchanged (explicit pool_mode == implicit)
    assert default == spikein.source_pool(
        bam, f"{T1},{T2}", g, "forward", 90, 1, pool_mode="default", **kw
    )
    assert (
        len(spikein.source_pool(bam, T1, g, "forward", 90, 1, max_nh=1, pool_mode="outside_truth"))
        == 0
    )


def _spike_args(**kw):
    a = dict(
        bam="x", gtf="g", bg_r1="b", out_dir=None, contig=T1, max_nh=20, pool="default", levels=[3],
        seeds=[1], cells=2, whitelist=None, cb_len=16, umi_len=10, r1_len=26, r2_len=90,
        source_orient="forward", bg_orient="forward",
    )  # fmt: skip
    return type("A", (), {**a, **kw})()


def _run_spike(tmp_path, monkeypatch, pool, size):
    monkeypatch.setattr(
        spikein, "source_pool",
        lambda *a, **k: [{"seq": DNA, "qual": "F" * 98, "cls": "sense_in_gene", "pos": i, "flag": 0,
                          "contig": T1, "nh": 2, "src_contig": T1} for i in range(size)],
    )  # fmt: skip
    monkeypatch.setattr(spikein, "top_barcodes", lambda *a, **k: ["A" * 16, "C" * 16])
    monkeypatch.setattr(spikein, "load_genes", lambda *a, **k: [])
    monkeypatch.setattr(spikein, "background_umis", lambda *a, **k: set())
    out = tmp_path / pool
    spikein.cmd_spike(_spike_args(out_dir=str(out), pool=pool))
    return out


def test_truth_columns_default_identical_and_pool_label_otherwise(tmp_path, monkeypatch):
    d = _run_spike(tmp_path, monkeypatch, "default", 10)
    h = _run_spike(tmp_path, monkeypatch, "outside_truth", 10)
    td = (d / "spike_L3_s1.truth.tsv").read_text().splitlines()
    th = (h / "spike_L3_s1.truth.tsv").read_text().splitlines()
    assert "pool" not in td[0] and th[0] == td[0] + "\tpool"
    assert [x + "\toutside_truth" for x in td[1:]] == th[1:]
    assert '"pool"' not in (d / "spike_L3_s1.json").read_text()


def test_empty_or_small_pool_fails_loudly(tmp_path, monkeypatch):
    with pytest.raises(SystemExit, match="has 0 reads"):
        _run_spike(tmp_path, monkeypatch, "outside_truth", 0)
    with pytest.raises(SystemExit, match="need >= 3"):
        _run_spike(tmp_path / "s", monkeypatch, "outside_truth", 2)


def test_conditions_30_rows_per_pool_unique_ids():
    for pool in ("outside_truth", "default"):
        rows = hosthomolog.build_conditions(pool)
        assert len(rows) == 30 and {r["cells"] for r in rows} == {100}
        ids = [r["cond_id"] for r in rows]
        assert len(set(ids)) == 30 and not any("_" in i for i in ids)
    assert not {r["cond_id"] for r in hosthomolog.build_conditions("default")} & {
        r["cond_id"] for r in hosthomolog.build_conditions("outside_truth")
    }


def test_host_loss_fraction_from_fake_lineage(tmp_path):
    truth = tmp_path / "t.tsv"
    names = [f"SPIKE:c:{i}" for i in range(4)]
    truth.write_text(
        "name\tcb\tumi\n" + "".join(f"{n}\tAAAA\tU{i // 2}\n" for i, n in enumerate(names))
    )  # molecules U0 (reads 0,1), U1 (reads 2,3)
    lin = tmp_path / "lineage.tsv.gz"
    with gzip.open(lin, "wt") as fh:
        fh.write("read_id\tfilter_decision\treason\n")
        for i, n in enumerate(names + ["bg1"]):
            fh.write(f"{n}\t{'removed' if i in (0, 2, 3) else 'retained'}\tx\n")
    r = hosthomolog.host_loss(truth, lin, tmp_path / "w")
    assert (r["d15_truth_fragments"], r["d15_removed_fragments"]) == (4, 3)
    assert r["d15_loss_fraction"] == 0.75
    assert (r["d16_truth_molecules"], r["d16_lost_molecules"], r["d16_loss_fraction"]) == (
        2,
        1,
        0.5,
    )
    assert r["not_in_lineage"] == 0
    lin2 = tmp_path / "l2.tsv.gz"
    with gzip.open(lin2, "wt") as fh:
        fh.write("read_id\tfilter_decision\treason\nSPIKE:c:0\tremoved\tx\n")
    assert hosthomolog.host_loss(truth, lin2, tmp_path / "w")["not_in_lineage"] == 3


def test_run_writes_loss_table_with_stubbed_star(tmp_path, monkeypatch):
    root, bgd = tmp_path / "root", tmp_path / "gridroot" / "bg" / "D03-r1"
    bgd.mkdir(parents=True)
    (bgd / "READY").write_text("{}")
    for m in (1, 2):
        (bgd / f"bg_R{m}.fastq.gz").write_bytes(gzip.compress(b""))
    monkeypatch.setattr(hosthomolog.grid, "ROOT", tmp_path / "gridroot")
    cond = {k: str(v) for k, v in hosthomolog.build_conditions("outside_truth")[0].items()}
    cond["level"] = "4"

    def fake_spike(ctx, pool):
        sp = ctx.d / "spikes"
        sp.mkdir(parents=True)
        (sp / "spike_L4_s1.truth.tsv").write_text(
            "name\tcb\tumi\n" + "".join(f"SPIKE:c:{i}\tAAAA\tU{i}\n" for i in range(4))
        )
        (sp / "spike_L4_s1.json").write_text('{"pool_class_mix": {"sense_in_gene": 40}}')
        for m in (1, 2):
            (sp / f"spike_L4_s1_R{m}.fastq.gz").write_bytes(gzip.compress(b""))

    def fake_star(ctx, param, out_dir):
        out_dir.mkdir(parents=True)
        gone = 2 if param == "pinned" else 1  # removed spikes
        with gzip.open(out_dir / "fragment_lineage.tsv.gz", "wt") as fh:
            fh.write("read_id\tfilter_decision\treason\n")
            for i in range(4):
                fh.write(f"SPIKE:c:{i}\t{'removed' if i < gone else 'retained'}\tx\n")
        (out_dir / "host_filter_audit.tsv").write_text(
            "category\tfragments\ninput\t5\nretained_host_unmapped\t3\n"
        )

    monkeypatch.setattr(hosthomolog, "stage_spike", fake_spike)
    monkeypatch.setattr(hosthomolog, "run_star_filter", fake_star)
    hosthomolog.run(cond, root=root)
    rows = [
        r.split("\t")
        for r in (root / "cond" / cond["cond_id"] / "host_loss.tsv").read_text().splitlines()
    ]
    h = rows[0]
    by = {r[h.index("star_params")]: r for r in rows[1:]}
    assert by["pinned"][h.index("d15_loss_fraction")] == "0.5"
    assert by["star-default"][h.index("d15_loss_fraction")] == "0.25"
    assert by["pinned"][h.index("pool_size")] == "40"
    assert not (root / "cond" / cond["cond_id"] / "work_pinned").exists()
