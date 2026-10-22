import collections
import csv
import gzip
import importlib.util
import json
import shutil
import subprocess
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


# ---- EBV-style multi-contig truth set -------------------------------------------------------------------

T1, T2, HOST, OTHER = "VIRUS|T1", "VIRUS|T2", "HOST|chr1", "VIRUS|OTHER"


def _bam(tmp_path, alns):
    """Sorted+indexed BAM from (name, flag, contig, pos, nh, seq) tuples; all 98M."""
    sq = "".join(f"@SQ\tSN:{c}\tLN:5000\n" for c in (T1, T2, HOST, OTHER))
    body = "".join(
        f"{n}\t{fl}\t{c}\t{pos}\t255\t{len(sq_)}M\t*\t0\t0\t{sq_}\t{'F' * len(sq_)}\tNH:i:{nh}\n"
        for n, fl, c, pos, nh, sq_ in alns
    )
    sam = tmp_path / "x.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n" + sq + body)
    bam = str(tmp_path / "x.bam")
    subprocess.run(["samtools", "sort", "-o", bam, str(sam)], check=True)
    subprocess.run(["samtools", "index", bam], check=True)
    return bam


def _alns():
    dna = ("ACGTTGCAAC" * 10)[:98]
    out = []
    for i in range(6):  # NH=1 on T1, inside gene (T1 gene 100-4000 '+', forward library)
        out.append((f"u{i}", 0, T1, 200 + i, 1, dna))
    for i in range(4):  # NH=2 T1+T2, primary alternates; flank-class read on T1 only for i==0
        pos1 = 4500 if i == 0 else 300 + i  # 4500 is 500 nt past the gene end -> flank3
        out += [(f"m{i}", 0 if i % 2 else 256, T1, pos1, 2, dna)]
        out += [(f"m{i}", 256 if i % 2 else 0, T2, 300 + i, 2, dna)]
    out += [("h0", 0, T1, 300, 2, dna), ("h0", 256, HOST, 300, 2, dna)]  # host partner
    out += [("o0", 0, T1, 300, 2, dna), ("o0", 256, OTHER, 300, 2, dna)]  # other-virus partner
    out += [("n3", 0, T1, 300, 3, dna), ("n3", 256, T2, 300, 3, dna), ("n3", 256, T2, 400, 3, dna)]
    return out


needs_samtools = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools missing")


@needs_samtools
def test_source_pool_nh_policy_and_truth_set_filtering(tmp_path, capsys):
    bam = _bam(tmp_path, _alns())
    g = {T1: [(100, 4000, "+")], T2: [(100, 4000, "+")]}
    one = spikein.source_pool(bam, T1, g, "forward", 90, 1, max_nh=1)
    assert {r["nh"] for r in one} == {1} and len(one) == 6  # NH>1 reads dropped, as before
    pool = spikein.source_pool(bam, f"{T1},{T2}", g, "forward", 90, 1, max_nh=2)
    err = capsys.readouterr().err
    names = collections.Counter(r["nh"] for r in pool)
    assert names == {1: 6, 2: 4}  # h0 (host) and o0 (other virus) excluded; n3 above max_nh
    assert '"alignment_outside_truth_set": 2' in err and '"nh_above_max": 1' in err
    # classified on the first truth contig (T1) even when T2 is primary; truth-set contig is the PRIMARY
    m = sorted((r for r in pool if r["nh"] == 2), key=lambda r: r["pos"])
    assert [r["src_contig"] for r in m] == [T1] * 4
    assert collections.Counter(r["contig"] for r in m) == {T1: 2, T2: 2}
    assert collections.Counter(r["cls"] for r in m) == {"sense_in_gene": 3, "flank3_unannotated": 1}
    assert all(r["seq"] == ("ACGTTGCAAC" * 10)[:98] for r in pool)  # flag-0 reads stay as stored


@needs_samtools
def test_ebv_style_cli_forward_forward_no_rc_geometry_truth_and_class_mix(tmp_path, monkeypatch):
    bam = _bam(tmp_path, _alns())
    (tmp_path / "g.gtf").write_text(
        "".join(f'{c}\tx\tgene\t100\t4000\t.\t+\t.\tgene_id "g";\n' for c in (T1, T2))
    )
    monkeypatch.setattr(spikein, "top_barcodes", lambda *a, **k: ["A" * 16, "C" * 16])
    monkeypatch.setattr(spikein, "background_umis", lambda *a, **k: set())
    out = tmp_path / "o"
    spikein.main(
        ["spike", "--bam", bam, "--gtf", str(tmp_path / "g.gtf"), "--bg-r1", "x", "--out-dir",
         str(out), "--contig", f"{T1},{T2}", "--max-nh", "2", "--levels", "10", "--cb-len", "16",
         "--umi-len", "12", "--r1-len", "28", "--r2-len", "90", "--source-orient", "forward",
         "--bg-orient", "forward"]
    )  # fmt: skip
    meta = json.loads((out / "spike_L10_s1.json").read_text())
    assert meta["rc"] is False and meta["pool_class_mix"] == {
        "sense_in_gene": 9, "flank3_unannotated": 1, "antisense_in_gene": 0, "intergenic": 0
    }  # fmt: skip
    rows = list(csv.DictReader((out / "spike_L10_s1.truth.tsv").open(), delimiter="\t"))
    # largest-remainder split of 10 over 9:1 is exactly 9 sense + 1 flank3: proportions preserved
    assert collections.Counter(r["class"] for r in rows) == {
        "sense_in_gene": 9,
        "flank3_unannotated": 1,
    }
    assert {r["contigs"] for r in rows} == {f"{T1},{T2}"} and {r["nh_stratum"] for r in rows} <= {
        "1",
        "2",
    }
    assert {r["contig"] for r in rows} <= {T1, T2} and {r["src_contig"] for r in rows} == {T1}
    with gzip.open(out / "spike_L10_s1_R1.fastq.gz", "rt") as fh:
        assert {len(x.strip()) for x in list(fh)[1::4]} == {28}
    with gzip.open(out / "spike_L10_s1_R2.fastq.gz", "rt") as fh:
        r2 = list(fh)[1::4]
    assert {len(x.strip()) for x in r2} == {90}
    assert {x.strip() for x in r2} == {
        ("ACGTTGCAAC" * 10)[:90]
    }  # 98 -> 5' 90 nt, not reverse-complemented
    # score.py picks the truth set up from the truth table and counts any alignment on any truth contig
    lines = [
        "SPIKE:a:0\t256\tVIRUS|T2\t1\n", "SPIKE:a:0\t0\tHOST|chr1\t1\n",  # secondary on T2 only
        "SPIKE:a:1\t0\tVIRUS|OTHER\t1\n", "SPIKE:a:2\t0\tVIRUS|T1\t1\n",
    ]  # fmt: skip
    ts = score.truth_contigs(None, rows, tmp_path)
    fate, anyc = score.star_fates(lines, ts)
    assert ts == {T1, T2} and anyc == {"SPIKE:a:0", "SPIKE:a:2"}
    assert fate == {"SPIKE:a:0": "host", "SPIKE:a:1": "other_virus", "SPIKE:a:2": "source"}
    assert score.truth_contigs("A,B", [], "x") == {"A", "B"}


def test_truth_contigs_falls_back_to_meta_json(tmp_path):
    (tmp_path / "t.json").write_text(json.dumps({"contig": "VIRUS|X,VIRUS|Y"}))
    assert score.truth_contigs(None, [{"name": "n"}], tmp_path / "t.truth.tsv") == {
        "VIRUS|X",
        "VIRUS|Y",
    }
