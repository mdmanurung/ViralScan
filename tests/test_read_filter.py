"""DEF-01 read-artefact filter (``scripts/read_filter.py``)."""

from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path

import pytest

from tests._fastq import write_fastq
from viralscan import anello_align as aa
from viralscan.runconfig import RunConfig
from viralscan.scripts import read_filter as rf
from viralscan.scripts.host_filter import lost_truth_counts

_RNG = random.Random(7)
BODY = "".join(_RNG.choice("ACGT") for _ in range(60))
CB_UMI = "ACGTACGTACGTACGTAAACCCGGGTTT"  # 28 nt, 10x v3


@pytest.mark.parametrize(
    "r2",
    [
        "G" * 90,  # two-colour no-signal
        "A" * 90,
        "CAG" * 30,
        aa.TSO + "A" * 65,  # TSO then poly-A, the VAL-01 planter's TSO+polyA class
    ],
    ids=["polyG", "polyA", "CAG", "TSO+polyA"],
)
def test_val01_artefact_classes_are_removed(r2: str) -> None:
    """Every artefact class the VAL-01 planter injects (tract >= 40 nt) is dropped."""
    assert rf.classify(CB_UMI, r2, 28)


@pytest.mark.parametrize(
    "r2",
    [
        BODY[:20] + "A" * 70,
        BODY + "A" * 30,
        BODY + "G" * 30,
        BODY + BODY[:30],
        "T" * 30 + BODY,  # 5' R2 is antisense: the mirror of [body][poly-A]
        aa._revcomp(BODY[:20] + "A" * 70),
    ],
    ids=[
        "20nt-body+polyA",
        "body+polyA",
        "body+trailing-polyG",
        "complex",
        "polyT+body",
        "antisense-20nt-body",
    ],
)
def test_genuine_reads_are_retained(r2: str) -> None:
    assert rf.classify(CB_UMI, r2, 28) == ""


def test_r1_tso_counts_only_in_the_barcode_umi_span() -> None:
    assert rf.classify(aa.TSO[:28], BODY, 28) == "r1_tso"
    # Past CB+UMI (a long 5' R1 reading into cDNA) it is not a barcode chimera.
    assert rf.classify(CB_UMI + aa.TSO + BODY, BODY, 28) == ""


def test_reason_order_is_fixed() -> None:
    """A pair failing every check gets the first reason, so audits are deterministic."""
    assert rf.classify(aa.TSO[:28], aa.TSO + "T" * 65, 28) == "r1_tso"
    assert rf.classify(CB_UMI, aa.TSO + "T" * 65, 28) == "r2_reagent"


def test_filter_pairs_writes_retained_pairs_lineage_and_audit(tmp_path: Path) -> None:
    pairs = [
        ("keep-a", CB_UMI, BODY + "A" * 30),
        ("gone-g", CB_UMI, "G" * 90),
        ("keep-b", CB_UMI, BODY),
        ("gone-tso", aa.TSO[:28], BODY),
    ]
    write_fastq(tmp_path / "in_R1.fastq.gz", [(f"{n}/1", r1) for n, r1, _ in pairs])
    write_fastq(tmp_path / "in_R2.fastq.gz", [(f"{n}/2", r2) for n, _, r2 in pairs])
    out = tmp_path / "out"
    out.mkdir()

    counts = rf.filter_pairs(
        str(tmp_path / "in_R1.fastq.gz"), str(tmp_path / "in_R2.fastq.gz"), out, "10xv3"
    )

    assert counts == {
        "input": 4,
        "retained": 2,
        "r1_tso": 1,
        "r2_reagent": 0,
        "r2_no_complex_body": 1,
    }
    with gzip.open(out / "R2.fastq.gz", "rt") as handle:
        assert [line.strip() for line in handle][1::4] == [BODY + "A" * 30, BODY]
    with gzip.open(out / "fragment_lineage.tsv.gz", "rt") as handle:
        lineage = [
            (r["read_id"], r["filter_decision"], r["reason"])
            for r in csv.DictReader(handle, delimiter="\t")
        ]
    assert lineage == [
        ("keep-a", "retained", "complex_body"),
        ("gone-g", "removed", "r2_no_complex_body"),
        ("keep-b", "retained", "complex_body"),
        ("gone-tso", "removed", "r1_tso"),
    ]
    with (out / "read_filter_audit.tsv").open() as handle:
        audit = {r["category"]: r["fragments"] for r in csv.DictReader(handle, delimiter="\t")}
    assert audit["removed_r2_no_complex_body"] == "1"
    assert audit["pct_retained"] == "50.00"
    assert audit["param:homopolymer_run"] == "15"
    assert audit["param:min_body_entropy"] == "2.0"

    # The lineage has the host filter's schema, so VAL-01 loss accounting reuses it.
    truth = tmp_path / "truth.tsv"
    truth.write_text("read_id\tlabel\nkeep-a\tviral\ngone-g\tviral\n")
    assert (
        lost_truth_counts(str(truth), str(out / "fragment_lineage.tsv.gz"))["d15_removed_fragments"]
        == 1
    )


def test_mate_mismatch_fails(tmp_path: Path) -> None:
    write_fastq(tmp_path / "r1.fastq.gz", [("a/1", CB_UMI)])
    write_fastq(tmp_path / "r2.fastq.gz", [("b/2", BODY)])
    with pytest.raises(ValueError, match="mate mismatch"):
        rf.filter_pairs(
            str(tmp_path / "r1.fastq.gz"), str(tmp_path / "r2.fastq.gz"), tmp_path, "10xv3"
        )


@pytest.mark.parametrize(
    ("read_filter", "host_index", "expected_dir"),
    [
        ("artefact", None, "read_filtered"),
        ("artefact", "/host", "read_filtered"),
        ("off", "/host", "host_filtered"),
        ("", None, None),
    ],
)
def test_kb_reads_the_last_filter_output(read_filter, host_index, expected_dir) -> None:
    cfg = {
        "output": "/out/",
        "index": "i",
        "transcripts": "t",
        "sample1": "s1.fastq.gz",
        "sample2": "s2.fastq.gz",
        "gtf": "",
        "fasta": "",
        "visual": "false",
        "f1": "",
        "reference": "false",
        "umap": "false",
        "technology": "10xv3",
        "whitelist": "",
        "multimapping": "true",
        "host_index": host_index or "",
        "read_filter": read_filter,
    }
    config = RunConfig.from_snakemake_config(cfg)
    assert config.read_filter == (read_filter or "off")
    if expected_dir is None:
        assert config.kb_r1 == "s1.fastq.gz"
    else:
        assert config.kb_r1 == f"/out/{expected_dir}/R1.fastq.gz"


TSO30 = aa.TSO + "ATGGG"  # F-028: the 30 nt the HPV77 false positives carry


def _trim_pairs(tmp_path: Path, r2s: list[str]) -> tuple[dict[str, int], list[str], list[str]]:
    write_fastq(tmp_path / "r1.fastq.gz", [(f"p{i}/1", CB_UMI) for i in range(len(r2s))])
    write_fastq(tmp_path / "r2.fastq.gz", [(f"p{i}/2", r2) for i, r2 in enumerate(r2s)])
    out = tmp_path / "out"
    out.mkdir()
    counts = rf.filter_pairs(
        str(tmp_path / "r1.fastq.gz"), str(tmp_path / "r2.fastq.gz"), out, "10xv3", "tso-trim"
    )
    with gzip.open(out / "R2.fastq.gz", "rt") as handle:
        lines = [line.strip() for line in handle]
    with gzip.open(out / "fragment_lineage.tsv.gz", "rt") as handle:
        reasons = [r["reason"] for r in csv.DictReader(handle, delimiter="\t")]
    return counts, lines, reasons


def test_tso_trim_cuts_only_a_leading_tso(tmp_path: Path) -> None:
    inside = BODY[:20] + TSO30 + BODY[20:]  # a TSO inside the read is not touched
    counts, lines, reasons = _trim_pairs(tmp_path, [TSO30 + BODY, BODY, inside, TSO30 + BODY[:10]])
    assert counts == {
        "input": 4,
        "retained": 3,
        "r2_short_after_trim": 1,
        "tso_trimmed": 1,
    }
    assert lines[1::4] == [BODY, BODY, inside]
    assert lines[3::4] == ["I" * len(BODY), "I" * len(BODY), "I" * len(inside)]  # quality trimmed too
    assert reasons == ["tso_trimmed", "untouched", "untouched", "r2_short_after_trim"]


def test_tso_trim_removes_the_f028_junction() -> None:
    """The HPV77 read: TSO + GGG + CAG tract. Trimmed, no ``ACATGGGGCAG`` junction k-mer is left."""
    read = TSO30 + "GCAGCAGCAGCAGCAGCAGCAGCAGAGACCTCTCCACTTTCCCTTAGCCCCTCTGCTG"
    _, rec = rf.trim_tso(["@a\n", read + "\n", "+\n", "I" * len(read) + "\n"])
    assert "ACATGGGG" not in rec[1] and rec[1].startswith("GCAGCAG")


def test_unknown_mode_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unknown read-filter mode"):
        rf.filter_pairs("a", "b", tmp_path, "10xv3", "bogus")
