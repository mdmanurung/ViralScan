import gzip
import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "spikein", Path(__file__).parent.parent / "analysis" / "spikein_detection_limit" / "spikein.py"
)
spikein = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(spikein)


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
