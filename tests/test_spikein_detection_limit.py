import gzip
import importlib.util
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location(
    "spikein", Path(__file__).parent.parent / "analysis" / "spikein_detection_limit" / "spikein.py"
)
spikein = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(spikein)
_sspec = importlib.util.spec_from_file_location(
    "spike_score", Path(_spec.origin).with_name("score.py")
)
score = importlib.util.module_from_spec(_sspec)
_sspec.loader.exec_module(score)


def test_allocate_sums_to_level_and_keeps_proportions():
    mix = {"sense_in_gene": 68, "flank3_unannotated": 9, "antisense_in_gene": 19, "intergenic": 4}
    for n in (3, 10, 30, 100, 300):
        assert sum(spikein.allocate(mix, n).values()) == n
    assert spikein.allocate(mix, 100) == {
        "sense_in_gene": 68,
        "flank3_unannotated": 9,
        "antisense_in_gene": 19,
        "intergenic": 4,
    }


def test_make_spike_geometry_trim_and_rc():
    read = {"seq": "AAAACCCCGG", "qual": "ABCDEFGHIJ"}
    r1, r2 = spikein.make_spike(read, "SPIKE:x:0", "C" * 16, "G" * 10, 28, 6, rc=False)
    assert r1.split("\n")[1] == "C" * 16 + "G" * 10 + "AA" and len(r1.split("\n")[3]) == 28
    assert r2.split("\n")[1:4:2] == ["AAAACC", "ABCDEF"]  # 5' end kept
    _, r2rc = spikein.make_spike(read, "SPIKE:x:0", "C" * 16, "G" * 10, 28, 6, rc=True)
    assert r2rc.split("\n")[1:4:2] == ["GGTTTT", "FEDCBA"]  # trim first, then reverse-complement


def test_classify_matches_flank_audit_rule():
    genes = [(1000, 2000, "+")]
    # reverse-oriented library: RNA strand is opposite the read; flag 16 -> RNA "+" -> sense
    assert spikein.classify(16, 1500, "50M", genes, "reverse") == "sense_in_gene"
    assert spikein.classify(0, 1500, "50M", genes, "reverse") == "antisense_in_gene"
    assert spikein.classify(16, 2500, "50M", genes, "reverse") == "flank3_unannotated"
    assert spikein.classify(16, 9000, "50M", genes, "reverse") == "intergenic"
    assert spikein.classify(0, 1500, "50M", genes, "forward") == "sense_in_gene"


def test_spike_cli_unique_umis_and_truth(tmp_path, monkeypatch):
    pool = [
        {"seq": "ACGT" * 30, "qual": "F" * 120, "cls": c, "pos": 100 + i, "flag": 16}
        for i, c in enumerate(["sense_in_gene"] * 30 + ["antisense_in_gene"] * 10)
    ]
    monkeypatch.setattr(spikein, "source_pool", lambda *a, **k: pool)
    monkeypatch.setattr(spikein, "top_barcodes", lambda *a, **k: ["A" * 16, "C" * 16])
    monkeypatch.setattr(spikein, "load_genes", lambda *a, **k: [])
    monkeypatch.setattr(spikein, "background_umis", lambda *a, **k: set())
    spikein.main(
        [
            "spike",
            "--bam",
            "x",
            "--gtf",
            "x",
            "--bg-r1",
            "x",
            "--out-dir",
            str(tmp_path),
            "--levels",
            "8",
            "--cb-len",
            "16",
            "--umi-len",
            "10",
            "--r1-len",
            "26",
            "--r2-len",
            "90",
            "--bg-orient",
            "forward",
        ]
    )
    rows = [
        ln.split("\t") for ln in (tmp_path / "spike_L8_s1.truth.tsv").read_text().splitlines()[1:]
    ]
    assert len(rows) == 8 and len({(r[1], r[2]) for r in rows}) == 8
    assert {r[1] for r in rows} == {"A" * 16, "C" * 16}
    with gzip.open(tmp_path / "spike_L8_s1_R2.fastq.gz", "rt") as fh:
        r2 = fh.read().split("\n")
    assert r2[1] == ("ACGT" * 30)[:90][::-1].translate(
        str.maketrans("ACGT", "TGCA")
    )  # source reverse vs bg forward -> RC


def _write_fq(path, n, tag):
    with gzip.open(path, "wt") as fh:
        for i in range(n):
            fh.write(f"@r{i}\n{tag * 8}\n+\n{'F' * 8}\n")


def _bg(tmp_path, name, *extra):
    pre = str(tmp_path / name)
    spikein.main(
        [
            "bg",
            "--r1",
            str(tmp_path / "a_R1.fq.gz"),
            "--r2",
            str(tmp_path / "a_R2.fq.gz"),
            "--out-prefix",
            pre,
            *extra,
        ]
    )

    def ids(r):
        with gzip.open(f"{pre}_{r}.fastq.gz", "rt") as fh:
            return [ln.strip() for ln in fh][0::4]

    return ids("R1"), ids("R2")


def test_bg_frac_lockstep_and_seed(tmp_path):
    _write_fq(tmp_path / "a_R1.fq.gz", 400, "A")
    _write_fq(tmp_path / "a_R2.fq.gz", 400, "C")
    r1a, r2a = _bg(tmp_path, "s1", "--frac", "0.25", "--seed", "1")
    r1a2, _ = _bg(tmp_path, "s1b", "--frac", "0.25", "--seed", "1")
    r1b, r2b = _bg(tmp_path, "s2", "--frac", "0.25", "--seed", "2")
    assert r1a == r2a and r1b == r2b  # mates stay in lockstep
    assert r1a == r1a2  # same seed reproduces
    assert r1a != r1b and 40 < len(r1a) < 160  # different seed differs; ~25% kept
    r1n, r2n = _bg(tmp_path, "n", "--n-pairs", "7")
    assert r1n == r2n == [f"@r{i}" for i in range(7)]


def test_spike_barcodes_onlist_and_umis_avoid_background(tmp_path):
    with gzip.open(tmp_path / "bg_R1.fq.gz", "wt") as fh:
        for cb, n in (
            ("T" * 16, 9),
            ("A" * 16, 5),
            ("C" * 16, 3),
        ):  # T is the top barcode but off-list
            for _ in range(n):
                fh.write(f"@x\n{cb}{'G' * 10}\n+\n{'F' * 26}\n")
    (tmp_path / "wl.txt").write_text(f"{'A' * 16}\n{'C' * 16}\n")
    assert spikein.top_barcodes(
        str(tmp_path / "bg_R1.fq.gz"),
        16,
        2,
        onlist=spikein.load_whitelist(str(tmp_path / "wl.txt")),
    ) == ["A" * 16, "C" * 16]
    with pytest.raises(SystemExit):
        spikein.top_barcodes(str(tmp_path / "bg_R1.fq.gz"), 16, 3, onlist={"A" * 16})
    assert spikein.background_umis(str(tmp_path / "bg_R1.fq.gz"), 16, 10, ["A" * 16]) == {
        ("A" * 16, "G" * 10)
    }


def test_star_fates_headline_counts_secondary_on_source_contig():
    lines = [
        "SPIKE:a:0\t256\tVIRUS|AF157706.1\t1\n",  # secondary on source only
        "SPIKE:a:0\t0\tHOST|chr1\t1\n",
        "SPIKE:a:1\t0\tVIRUS|OTHER.1\t1\n",
        "SPIKE:a:2\t16\tVIRUS|AF157706.1\t1\n",
        "bg\t0\tVIRUS|AF157706.1\t1\n",
    ]
    fate, any_src = score.star_fates(lines, "VIRUS|AF157706.1")
    assert any_src == {"SPIKE:a:0", "SPIKE:a:2"}
    assert fate == {"SPIKE:a:0": "host", "SPIKE:a:1": "other_virus", "SPIKE:a:2": "source"}


def test_synth_pool_origin_span_and_internal_reads_are_substrings_of_the_circle():
    genome = "".join("ACGT"[(i * 7 + i // 3) % 4] for i in range(400))
    circle = genome + genome  # every read of a circular genome is a substring of the doubled record
    pool = spikein.synth_pool(genome, "origin", 90, [], "forward", 30, seed=1)
    assert set(pool) == {"origin_span", "internal"}
    for r in pool["origin_span"]:
        assert (
            len(r["seq"]) == 90 and r["seq"] in circle and r["seq"] not in genome
        )  # crosses the junction
    for r in pool["internal"]:
        assert r["seq"] in genome


def test_synth_pool_short_lengths_and_reverse_orientation():
    genome = "ACGTTGCAAC" * 40
    pool = spikein.synth_pool(genome, "short", 90, [30, 50], "reverse", 10, seed=2)
    assert set(pool) == {"internal", "short_30", "short_50"}
    assert {len(r["seq"]) for r in pool["short_30"]} == {30}
    assert all(spikein.revcomp(r["seq"]) in genome for r in pool["short_50"])


def test_synth_cli_equal_class_split_and_truth(tmp_path, monkeypatch):
    fa = tmp_path / "g.fa"
    fa.write_text(">VIRUS|MZ1.1 anello\n" + "ACGTTGCAAC" * 40 + "\n")
    monkeypatch.setattr(spikein, "top_barcodes", lambda *a, **k: ["A" * 16, "C" * 16])
    monkeypatch.setattr(spikein, "background_umis", lambda *a, **k: set())
    spikein.main(
        ["synth", "--fasta", str(fa), "--contig", "MZ1.1", "--kind", "origin", "--bg-r1", "x",
         "--out-dir", str(tmp_path), "--levels", "10", "--cb-len", "16", "--umi-len", "10",
         "--r1-len", "26", "--r2-len", "90", "--bg-orient", "reverse"]
    )  # fmt: skip
    rows = [
        ln.split("\t") for ln in (tmp_path / "spike_L10_s1.truth.tsv").read_text().splitlines()[1:]
    ]
    assert len(rows) == 10 and {r[3] for r in rows} == {"origin_span", "internal"}
    assert sum(r[3] == "origin_span" for r in rows) == 5
