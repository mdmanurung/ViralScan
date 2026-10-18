"""B18 strict oracle: fused molecule reducer == frozen two-phase baseline.

The frozen pre-B18 module is loaded by path (read-only, a verbatim copy of
``src/viralscan/molecule_verdict.py`` in the wave-2 frozen baseline, commit
``faf5667``) and compared against the shipped module over the cases the card names:
secondary hits, supplementary/unmapped-only reads, malformed names and ``AS`` tags,
score ties, interleaved molecules, and counted-read mappings — including the card's
explicit rule that molecule order follows the first **usable scored alignment**, not
the first raw QNAME occurrence.

The writer is compared on decompressed bytes, summary, and the gzip stream with its
wall-clock ``mtime`` field masked (the header carries the build time, so raw bytes are
only comparable within the same second). Parse-failure precedence and the
parse-before-write boundary are pinned too, and the last test fails on the pre-change
two-phase writer rather than merely describing it.
"""

from __future__ import annotations

import csv
import gc
import gzip
import importlib.util
import sys
import tracemalloc
from pathlib import Path

import pytest

from viralscan import molecule_verdict as cand

BASELINE = (
    Path(__file__).resolve().parent / "_frozen_oracles" / "frozen_molecule_verdict_faf5667.py"
)
assert BASELINE.is_file(), BASELINE


def _load_baseline():
    spec = importlib.util.spec_from_file_location("viralscan_b18_frozen_verdict", BASELINE)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


base = _load_baseline()
assert Path(cand.__file__).resolve() != BASELINE


def _sam(qname: str, flag: int, ref: str, score: int | None) -> str:
    tags = [] if score is None else [f"AS:i:{score}"]
    return "\t".join(
        [qname, str(flag), ref, "1", "60", "50M", "*", "0", "0", "A" * 50, "I" * 50, *tags]
    )


CASES: dict[str, list[str]] = {
    "secondary-hits": [
        _sam("AAA_U1_1|r", 0, "VIRUS|x", 90),
        _sam("AAA_U1_1|r", 256, "HOST|chr1", 60),
        _sam("AAA_U1_1|r", 256, "HOST|chr2", 61),
        _sam("AAA_U2_2|r", 0, "HOST|chr1", 80),
        _sam("AAA_U2_2|r", 256, "VIRUS|x", 70),
    ],
    "tie": [_sam("AAA_U1_1|r", 0, "HOST|chr1", 75), _sam("AAA_U1_1|r", 256, "VIRUS|x", 75)],
    "supplementary-and-unmapped-only": [
        _sam("AAA_U1_1|r", 4, "*", None),
        _sam("AAA_U2_2|r", 0, "VIRUS|x", 50),
        _sam("AAA_U2_2|r", 2048, "HOST|chr1", 99),
        _sam("AAA_U3_3|r", 2048, "VIRUS|x", 99),
        _sam("AAA_U4_4|r", 4, "VIRUS|x", None),
    ],
    "malformed-names": [
        _sam("NOUNDERSCORE", 0, "VIRUS|x", 50),
        _sam("_UB_1|r", 0, "VIRUS|x", 50),
        _sam("CB__1|r", 0, "VIRUS|x", 50),
        _sam("CB_UB|r", 0, "VIRUS|x", 50),
        _sam("CB_UB_notanumber|r", 0, "HOST|chr1", 50),
        _sam("CB_UB_7|x", 0, "HOST|chr1", 60),
    ],
    "no-as-tag": [_sam("AAA_U1_1|r", 0, "VIRUS|x", None), _sam("AAA_U2_2|r", 256, "HOST|chr1", 70)],
    "unmapped-reference-star": [_sam("AAA_U1_1|r", 0, "*", 50)],
    "header-lines-and-blanks": [
        "@HD\tVN:1.6",
        "@SQ\tSN:HOST|chr1",
        "",
        "AAA_U1_1|r",
    ],
    "interleaved-molecules": [
        _sam("AAA_U1_1|r", 0, "VIRUS|x", 90),
        _sam("BBB_U1_2|r", 0, "HOST|chr1", 91),
        _sam("AAA_U1_3|r", 0, "HOST|chr1", 95),
        _sam("CCC_U1_4|r", 0, "VIRUS|x", 40),
        _sam("BBB_U1_5|r", 0, "VIRUS|x", 99),
    ],
    "first-raw-occurrence-unusable": [
        # the first raw QNAME occurrence is supplementary and the second is unmapped,
        # so the molecule order must come from the later usable alignment
        _sam("ZZZ_U9_1|r", 2048, "HOST|chr1", 99),
        _sam("YYY_U9_2|r", 4, "*", None),
        _sam("AAA_U1_3|r", 0, "VIRUS|x", 50),
        _sam("ZZZ_U9_4|r", 0, "HOST|chr1", 60),
    ],
    "order-not-alphabetical": [
        # usable-first order (ZZZ, AAA, MMM) is not the sorted order (AAA, MMM, ZZZ)
        _sam("ZZZ_U9_1|r", 0, "VIRUS|x", 50),
        _sam("AAA_U1_2|r", 0, "HOST|chr1", 60),
        _sam("MMM_U5_3|r", 0, "HOST|chr1", 60),
    ],
    "empty-input": [],
    "only-headers": ["@HD\tVN:1.6", "@PG\tID:samtools"],
    "truncated-fields": ["AAA_U1_1|r\t0\tVIRUS|x", _sam("AAA_U2_2|r", 0, "HOST|chr1", 60)],
    "counted-boundaries": [
        _sam("AAA_U1_1|r", 0, "VIRUS|x", 50),
        _sam("AAA_U2_2|r", 0, "HOST|chr1", 60),
        _sam("AAA_U3_3|r", 0, "HOST|chr1", 60),
    ],
}


def _counted(case: str) -> set[int]:
    if case == "counted-boundaries":
        return {1, 3}
    return {1, 2}


def _gzip_sha_without_mtime(path: Path) -> str:
    """gzip bytes with the 4-byte wall-clock mtime field zeroed.

    The header also carries the output file name (``FNAME``), which both arms write as
    ``molecule_verdicts.tsv.gz``; only the ``mtime`` field is time-dependent.
    """
    import hashlib

    data = bytearray(path.read_bytes())
    assert data[:2] == b"\x1f\x8b", "not a gzip stream"
    data[4:8] = b"\x00\x00\x00\x00"
    return hashlib.sha256(bytes(data)).hexdigest()


def _write_tsv(path: Path, rows: list[tuple[int, str]]) -> Path:
    """A minimal lineage TSV with (read_number, assigned_weight)."""
    with gzip.open(path, "wt", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(["read_number", "read_id", "cb", "ub", "assigned_weight"])
        for read_number, weight in rows:
            w.writerow([read_number, "r", "CB", "UB", weight])
    return path


@pytest.fixture
def lineage(tmp_path: Path) -> Path:
    return _write_tsv(
        tmp_path / "read_lineage.tsv.gz", [(1, "1.0"), (2, "0.0"), (3, "0.5"), (4, "")]
    )


@pytest.mark.parametrize("case", list(CASES), ids=list(CASES))
def test_rows_identical_to_frozen_baseline(case):
    lines = CASES[case]
    counted = _counted(case)
    got = list(cand.molecule_rows(lines, counted))
    want = list(base.molecule_rows(lines, counted))
    assert got == want
    assert [tuple(r.items()) for r in got] == [tuple(r.items()) for r in want]


@pytest.mark.parametrize("case", list(CASES), ids=list(CASES))
def test_declared_molecule_order_is_first_usable_alignment(case):
    """The card's rule: order follows the first usable scored alignment."""
    lines = CASES[case]
    rows = list(cand.molecule_rows(lines, set()))
    order = [(r["cb"], r["ub"]) for r in rows]
    expected: list[tuple[str, str]] = []
    for line in lines:
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if cand.side_score(fields) is None:
            continue
        parts = fields[0].split("_")
        if len(parts) < 3 or not parts[0] or not parts[1]:
            continue
        key = (parts[0], parts[1])
        if key not in expected:
            expected.append(key)
    assert order == expected


@pytest.mark.parametrize("case", list(CASES), ids=list(CASES))
def test_writer_is_byte_identical(case, tmp_path, lineage):
    sam = CASES[case]
    for arm in ("cand", "base"):
        (tmp_path / arm).mkdir()
    got = cand.write_molecule_verdicts(sam, lineage, tmp_path / "cand")
    want = base.write_molecule_verdicts(sam, lineage, tmp_path / "base")
    assert got == want

    got_gz = gzip.decompress((tmp_path / "cand" / "molecule_verdicts.tsv.gz").read_bytes())
    want_gz = gzip.decompress((tmp_path / "base" / "molecule_verdicts.tsv.gz").read_bytes())
    assert got_gz == want_gz
    assert _gzip_sha_without_mtime(
        tmp_path / "cand" / "molecule_verdicts.tsv.gz"
    ) == _gzip_sha_without_mtime(tmp_path / "base" / "molecule_verdicts.tsv.gz")
    assert (tmp_path / "cand" / "molecule_verdict_summary.tsv").read_bytes() == (
        tmp_path / "base" / "molecule_verdict_summary.tsv"
    ).read_bytes()


def test_molecule_groups_and_read_counts(lineage, tmp_path):
    """reads counts unique usable read names per molecule; counted follows the lineage."""
    sam = [
        _sam("AAA_U1_1|r", 0, "VIRUS|x", 50),
        _sam("AAA_U1_1|r", 256, "HOST|chr1", 40),  # same read name: counted once
        _sam("AAA_U1_2|r", 0, "HOST|chr1", 60),  # a second read of the same molecule
        _sam("AAA_U2_3|r", 0, "HOST|chr1", 60),
        _sam("BBB_U1_4|r", 4, "*", None),  # no usable alignment: no molecule
    ]
    rows = list(cand.molecule_rows(sam, {1, 3}))
    assert rows == list(base.molecule_rows(sam, {1, 3}))
    assert [(r["cb"], r["ub"], r["reads"], r["counted"]) for r in rows] == [
        ("AAA", "U1", 2, True),
        ("AAA", "U2", 1, True),
    ]
    summary = cand.write_molecule_verdicts(sam, lineage, tmp_path)
    got = {(r["scope"], r["verdict"]): r["molecules"] for r in summary}
    assert got[("all", "host_best")] == 2
    assert got[("counted", "host_best")] == 2
    assert got[("counted", "virus_best")] == 0


def test_parse_failure_precedence_and_no_output(tmp_path, lineage):
    """SAM faults raise in stream order, before any output file exists."""
    bad_flag = [
        "@HD\tVN:1.6",
        _sam("AAA_U1_1|r", 0, "VIRUS|x", 50),
        "AAA_U2_2|r\tnope\tVIRUS|x\t1\t60\t50M\t*\t0\t0\tA\tI\tAS:i:5",
    ]
    bad_as = [
        _sam("AAA_U1_1|r", 0, "VIRUS|x", 50),
        _sam("AAA_U2_2|r", 0, "HOST|chr1", 60).replace("AS:i:60", "AS:i:zz"),
    ]
    for name, sam in (("bad-flag", bad_flag), ("bad-as", bad_as)):
        for arm in ("cand", "base"):
            (tmp_path / f"{name}-{arm}").mkdir()
        with pytest.raises(ValueError) as cand_exc:
            cand.write_molecule_verdicts(sam, lineage, tmp_path / f"{name}-cand")
        with pytest.raises(ValueError) as base_exc:
            base.write_molecule_verdicts(sam, lineage, tmp_path / f"{name}-base")
        assert str(cand_exc.value) == str(base_exc.value), name
        assert type(cand_exc.value) is type(base_exc.value)
        assert not (tmp_path / f"{name}-cand" / "molecule_verdicts.tsv.gz").exists()
        assert not (tmp_path / f"{name}-base" / "molecule_verdicts.tsv.gz").exists()


def test_lineage_error_precedes_sam_error(tmp_path):
    """A missing lineage file still fails before the SAM is parsed at all."""
    sam = ["not-a-sam-line\tat\tall", _sam("AAA_U1_1|r", 0, "VIRUS|x", 50)]
    for arm, mod in (("cand", cand), ("base", base)):
        (tmp_path / arm).mkdir()
        with pytest.raises(FileNotFoundError) as exc:
            mod.write_molecule_verdicts(sam, tmp_path / "missing.tsv.gz", tmp_path / arm)
        assert "missing.tsv.gz" in str(exc.value)
        assert not (tmp_path / arm / "molecule_verdicts.tsv.gz").exists()


def test_public_best_scores_api_is_unchanged():
    sam = CASES["secondary-hits"] + CASES["malformed-names"]
    assert cand.best_scores(sam) == base.best_scores(sam)
    assert list(cand.best_scores(sam)) == list(base.best_scores(sam))


def _fixture(n_reads: int, reads_per_molecule: int) -> list[str]:
    lines: list[str] = []
    for i in range(n_reads):
        mol = i // reads_per_molecule
        cb, ub = f"CB{mol % 509:04d}", f"U{mol:07d}"
        qname = f"{cb}_{ub}_{i}|r{i}"
        lines.append(_sam(qname, 0, "VIRUS|v1" if i % 3 == 0 else "HOST|chr1", 60 + (i % 40)))
        lines.append(_sam(qname, 256, "HOST|chr1" if i % 3 == 0 else "VIRUS|v1", 40 + (i % 30)))
    return lines


def _peak(fn) -> int:
    gc.collect()
    tracemalloc.start()
    try:
        fn()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


@pytest.mark.parametrize("rpm", [8, 1])
def test_writer_memory_is_below_frozen_baseline(rpm, tmp_path, lineage):
    """The shipped writer must not materialize the row list the frozen one holds.

    The thresholds are the measured ratios with margin: 0.53/0.73 shipped against
    0.42/0.91 for a writer that materializes the row list on top of the fused parser
    (``cards/b18/memory_b182.json``), so the 1-read-per-molecule case fails on that
    variant while the shipped code passes.
    """
    sam = _fixture(40_000, rpm)
    (tmp_path / "cand").mkdir()
    (tmp_path / "base").mkdir()
    cand_peak = _peak(lambda: cand.write_molecule_verdicts(sam, lineage, tmp_path / "cand"))
    base_peak = _peak(lambda: base.write_molecule_verdicts(sam, lineage, tmp_path / "base"))
    assert gzip.decompress(
        (tmp_path / "cand" / "molecule_verdicts.tsv.gz").read_bytes()
    ) == gzip.decompress((tmp_path / "base" / "molecule_verdicts.tsv.gz").read_bytes())
    assert base_peak > 8 * (1 << 20), (base_peak, rpm)
    assert cand_peak < base_peak * {1: 0.8, 8: 0.6}[rpm], (cand_peak, base_peak, rpm)
