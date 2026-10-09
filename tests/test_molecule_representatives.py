"""EVID-CORR-01 contract: one deterministic representative per corrected molecule.

A *molecule* is a corrected (CB, UMI) on one reference. Every diagnostic layer
(BAM deduplication, read-start profile, alignment QC, per-cell QC, samtools
coverage/depth) must count that molecule once, however many alignments it has at
different starts, strands or CIGARs, and must pick the same representative
whatever the input order. Fixtures are hand-computable; expected values are
written down, never derived from the code under test.

Representative rank (best first): highest MAPQ (255 = unavailable ranks lowest),
highest AS (missing ranks after present), lowest NM (missing ranks after present),
longest aligned query span (M/I/=/X; clips excluded), then lexical
(qname, pos, flag, cigar, full record).
"""

from __future__ import annotations

import itertools
import shutil
import subprocess
from pathlib import Path

import pytest

from viralscan import evidence
from viralscan.evidence import (
    _alignment_qc_from_text,
    _parse_sam_read_starts,
    _per_cell_qc_from_text,
    deduplicate_umi_sam,
    select_molecule_representatives,
)

HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:VIRUS|v\tLN:1000\n@SQ\tSN:HOST|h\tLN:1000\n"


def rec(qname, flag, rname, pos, mapq=60, cigar="50M", tags=("AS:i:50", "NM:i:0")):
    return "\t".join(
        [qname, str(flag), rname, str(pos), str(mapq), cigar, "*", "0", "0", "*", "*", *tags]
    )


def body(sam: str) -> list[str]:
    return [ln for ln in sam.splitlines() if ln and not ln.startswith("@")]


# ---- one molecule, two starts -------------------------------------------------

TWO_STARTS = [
    rec("CB1_U1_1", 0, "VIRUS|v", 100),
    rec("CB1_U1_2", 0, "VIRUS|v", 700),  # same molecule, a different start
]


def test_same_molecule_at_two_starts_survives_once_in_bam_dedup():
    out = body(deduplicate_umi_sam(HEADER + "\n".join(TWO_STARTS)))
    assert len(out) == 1


def test_same_molecule_at_two_starts_counts_once_in_start_profile():
    rows = _parse_sam_read_starts("\n".join(TWO_STARTS), dedup="umi")
    assert sum(int(r["n_read_starts"]) for r in rows) == 1
    assert rows[0]["n_reads"] == 1


@pytest.mark.parametrize(
    "other",
    [
        rec("CB1_U1_2", 0, "VIRUS|v", 101),  # different start
        rec("CB1_U1_2", 16, "VIRUS|v", 100),  # different strand
        rec("CB1_U1_2", 0, "VIRUS|v", 100, cigar="25M1D25M"),  # different CIGAR
        rec("CB1_U1_2", 16, "VIRUS|v", 900, cigar="10S40M"),  # all three at once
    ],
)
def test_start_strand_or_cigar_never_splits_a_molecule(other):
    sam = HEADER + rec("CB1_U1_1", 0, "VIRUS|v", 100) + "\n" + other
    assert len(body(deduplicate_umi_sam(sam))) == 1
    assert _parse_sam_read_starts(sam, dedup="umi")[0]["n_reads"] == 1


# ---- deterministic tie-break --------------------------------------------------


@pytest.mark.parametrize(
    "better, worse",
    [
        # MAPQ dominates every later criterion (worse has better AS/NM/span).
        (
            rec("z_U_1", 0, "VIRUS|v", 900, mapq=60, tags=("AS:i:1", "NM:i:9")),
            rec("a_U_1", 0, "VIRUS|v", 100, mapq=40, tags=("AS:i:99", "NM:i:0")),
        ),
        # equal MAPQ: higher AS
        (
            rec("z_U_1", 0, "VIRUS|v", 900, tags=("AS:i:60", "NM:i:9")),
            rec("a_U_1", 0, "VIRUS|v", 100, tags=("AS:i:50", "NM:i:0")),
        ),
        # equal MAPQ/AS: lower NM
        (
            rec("z_U_1", 0, "VIRUS|v", 900, tags=("AS:i:50", "NM:i:0")),
            rec("a_U_1", 0, "VIRUS|v", 100, tags=("AS:i:50", "NM:i:2")),
        ),
        # equal MAPQ/AS/NM: longer aligned query span (45M5S=45 beats 40M10S=40)
        (
            rec("z_U_1", 0, "VIRUS|v", 900, cigar="45M5S", tags=("AS:i:50", "NM:i:0")),
            rec("a_U_1", 0, "VIRUS|v", 100, cigar="40M10S", tags=("AS:i:50", "NM:i:0")),
        ),
        # a present AS/NM beats a missing one
        (
            rec("z_U_1", 0, "VIRUS|v", 900, tags=("AS:i:5", "NM:i:5")),
            rec("a_U_1", 0, "VIRUS|v", 100, tags=()),
        ),
        # MAPQ 255 means "unavailable": ranks below a real MAPQ of 0
        (rec("z_U_1", 0, "VIRUS|v", 900, mapq=0), rec("a_U_1", 0, "VIRUS|v", 100, mapq=255)),
    ],
)
def test_representative_rank_order(better, worse):
    # Fixture uses names "<x>_U_1": CB = "z"/"a", so make both the same molecule.
    better = better.replace("z_U_1", "CB_U_zz")
    worse = worse.replace("a_U_1", "CB_U_aa")
    kept, _ = select_molecule_representatives([worse, better])
    assert kept == [better]
    kept, _ = select_molecule_representatives([better, worse])
    assert kept == [better]


def test_full_tie_falls_back_to_lexical_qname_then_start():
    a = rec("CB_U_a", 0, "VIRUS|v", 900)
    b = rec("CB_U_b", 0, "VIRUS|v", 100)  # lower start but later name
    assert select_molecule_representatives([b, a])[0] == [a]
    same_name_late = rec("CB_U_a", 0, "VIRUS|v", 901)
    assert select_molecule_representatives([same_name_late, a])[0] == [a]


def test_input_order_permutation_gives_identical_artifacts():
    records = [
        rec("CB1_U1_1", 0, "VIRUS|v", 100, mapq=30),
        rec("CB1_U1_2", 16, "VIRUS|v", 700, mapq=60),  # winner for CB1/U1/VIRUS
        rec("CB1_U1_3", 0, "HOST|h", 50),  # separate reference -> own survivor
        rec("CB2_U2_1", 0, "VIRUS|v", 300),
        rec("CB2_U2_2", 0, "VIRUS|v", 300),  # complete tie -> lexical name
    ]
    ref_out = deduplicate_umi_sam(HEADER + "\n".join(records))
    ref_starts = _parse_sam_read_starts("\n".join(records), dedup="umi")
    for perm in itertools.permutations(records):
        assert deduplicate_umi_sam(HEADER + "\n".join(perm)) == ref_out
        assert _parse_sam_read_starts("\n".join(perm), dedup="umi") == ref_starts
    kept = body(ref_out)
    assert len(kept) == 3  # CB1/U1 on VIRUS and on HOST, CB2/U2 on VIRUS
    assert rec("CB1_U1_2", 16, "VIRUS|v", 700, mapq=60) in kept
    assert rec("CB2_U2_1", 0, "VIRUS|v", 300) in kept  # name tie-break


def test_same_representative_in_bam_dedup_and_start_profile():
    records = [
        rec("CB1_U1_1", 0, "VIRUS|v", 100, mapq=30),
        rec("CB1_U1_2", 16, "VIRUS|v", 700, cigar="50M", mapq=60),
    ]
    kept = body(deduplicate_umi_sam(HEADER + "\n".join(records)))
    assert kept == [records[1]]
    # strand-aware 5' start of the winner: 0-based 699 + 50 - 1 = 748.
    rows = _parse_sam_read_starts("\n".join(records), dedup="umi")
    assert [(r["position"], r["n_read_starts"]) for r in rows] == [(748, 1)]
    rows = _parse_sam_read_starts("\n".join(kept), dedup="none")
    assert [(r["position"], r["n_read_starts"]) for r in rows] == [(748, 1)]


# ---- references stay separate -------------------------------------------------


def test_different_references_keep_one_survivor_each():
    records = [
        rec("CB1_U1_1", 0, "VIRUS|v", 100),
        rec("CB1_U1_2", 0, "HOST|h", 100),
        rec("CB1_U1_3", 0, "VIRUS|v", 500),
    ]
    out = body(deduplicate_umi_sam(HEADER + "\n".join(records)))
    assert sorted(ln.split("\t")[2] for ln in out) == ["HOST|h", "VIRUS|v"]
    starts = _parse_sam_read_starts("\n".join(records), dedup="umi")
    assert {(r["reference"], r["n_reads"]) for r in starts} == {("VIRUS|v", 1), ("HOST|h", 1)}


# ---- filters: one policy everywhere -------------------------------------------

EXCLUDED = {
    "unmapped": rec("CB9_U9_1", 4, "*", 0),
    "secondary": rec("CB9_U9_1", 0x100, "VIRUS|v", 100),
    "supplementary": rec("CB9_U9_1", 0x800, "VIRUS|v", 100),
    "qcfail": rec("CB9_U9_1", 0x200, "VIRUS|v", 100),
    "duplicate": rec("CB9_U9_1", 0x400, "VIRUS|v", 100),
}


@pytest.mark.parametrize("name", sorted(EXCLUDED))
def test_excluded_flags_are_ignored_by_every_layer(name):
    keep = rec("CB1_U1_1", 0, "VIRUS|v", 100)
    sam = "\n".join([EXCLUDED[name], keep])
    assert select_molecule_representatives(sam.splitlines())[0] == [keep]
    assert body(deduplicate_umi_sam(HEADER + sam)) == [keep]
    for dedup in ("umi", "none"):
        rows = _parse_sam_read_starts(sam, dedup=dedup)
        assert sum(int(r["n_read_starts"]) for r in rows) == 1
    (qc,) = _alignment_qc_from_text(HEADER, sam, "")
    assert qc["reads"] == 1 and qc["molecules"] == 1
    (cell,) = _per_cell_qc_from_text(sam)
    assert cell["reads"] == 1 and cell["molecules"] == 1


def test_malformed_records_are_dropped_not_fatal():
    good = rec("CB1_U1_1", 0, "VIRUS|v", 100)
    junk = ["short\t0\tVIRUS|v", "x\tNaN\tVIRUS|v\t1\t60\t50M\t*\t0\t0\t*\t*", ""]
    assert select_molecule_representatives([*junk, good])[0] == [good]


# ---- unresolved lineage -------------------------------------------------------


def test_missing_lineage_is_read_level_and_reported_unresolved_not_a_molecule():
    named = rec("CB1_U1_1", 0, "VIRUS|v", 100)
    anon1 = rec("readA", 0, "VIRUS|v", 200)
    anon2 = rec("readB", 0, "VIRUS|v", 200)
    kept, unresolved = select_molecule_representatives([named, anon1, anon2])
    assert sorted(kept) == sorted([named, anon1, anon2])  # never merged together
    assert unresolved == 2
    (qc,) = _alignment_qc_from_text(HEADER, "\n".join(kept), "")
    assert qc["reads"] == 3 and qc["molecules"] == 1  # only CB1/U1 is a molecule


# ---- sorted, indexable output -------------------------------------------------


def test_dedup_output_is_coordinate_sorted_in_header_reference_order():
    records = [
        rec("CB1_U1_1", 0, "HOST|h", 10),
        rec("CB2_U2_1", 0, "VIRUS|v", 900),
        rec("CB3_U3_1", 0, "VIRUS|v", 20),
    ]
    out = deduplicate_umi_sam(HEADER + "\n".join(records))
    assert out.startswith(HEADER)
    keys = [(ln.split("\t")[2], int(ln.split("\t")[3])) for ln in body(out)]
    assert keys == [("VIRUS|v", 20), ("VIRUS|v", 900), ("HOST|h", 10)]  # @SQ order


# ---- QC / starts / coverage agree --------------------------------------------

RAW = [
    rec("CB1_U1_1", 0, "VIRUS|v", 100, mapq=30),
    rec("CB1_U1_2", 16, "VIRUS|v", 400, mapq=60),
    rec("CB1_U1_3", 0, "VIRUS|v", 700, mapq=60, cigar="25M1D25M"),
    rec("CB2_U2_1", 0, "VIRUS|v", 100),
    rec("CB3_U3_1", 0, "HOST|h", 100),
    rec("CB3_U3_2", 0, "HOST|h", 300),
    EXCLUDED["supplementary"],
    EXCLUDED["duplicate"],
]
# Hand count: VIRUS|v has molecules CB1/U1 and CB2/U2 -> 2; HOST|h has CB3/U3 -> 1.
EXPECTED = {"VIRUS|v": 2, "HOST|h": 1}


def test_qc_per_cell_and_starts_agree_on_the_deduplicated_layer():
    dedup = deduplicate_umi_sam(HEADER + "\n".join(RAW))
    qc = {r["reference"]: r for r in _alignment_qc_from_text(HEADER, dedup, "")}
    starts: dict[str, int] = {}
    for r in _parse_sam_read_starts(dedup, dedup="none"):
        starts[str(r["reference"])] = int(r["n_reads"])  # n_reads repeats per reference
    direct = {
        str(r["reference"]): int(r["n_reads"])
        for r in _parse_sam_read_starts("\n".join(RAW), dedup="umi")
    }
    assert {k: v["reads"] for k, v in qc.items()} == EXPECTED
    assert {k: v["molecules"] for k, v in qc.items()} == EXPECTED
    assert starts == direct == EXPECTED
    cells = {(r["cell_barcode"], r["reference_class"]): r for r in _per_cell_qc_from_text(dedup)}
    assert {k: (v["reads"], v["molecules"]) for k, v in cells.items()} == {
        ("CB1", "virus"): (1, 1),
        ("CB2", "virus"): (1, 1),
        ("CB3", "host"): (1, 1),
    }


# ---- samtools: the same policy reaches coverage / depth -----------------------


def test_samtools_invocations_pin_flag_mapq_and_baseq_policy(monkeypatch):
    calls: list[list[str]] = []

    def fake_run(cmd, **_kw):
        calls.append(list(cmd))
        return b""

    monkeypatch.setattr(evidence, "_run", fake_run)
    evidence.coverage_table("x.bam")
    evidence.coverage_depth_points("x.bam")
    evidence._run(["samtools", "view", "x.bam"])
    cov = next(c for c in calls if c[:2] == ["samtools", "coverage"])
    depth = next(c for c in calls if c[:2] == ["samtools", "depth"])
    assert cov[cov.index("--ff") : cov.index("--ff") + 2] == ["--ff", "0xF04"]
    assert "0x800" in depth  # supplementary excluded on top of samtools' defaults
    for c in (cov, depth):
        assert c[c.index("-q") + 1] == "0" and c[c.index("-Q") + 1] == "0"


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")
def test_samtools_coverage_numreads_equals_python_counts(tmp_path: Path):
    sam = tmp_path / "raw.sam"
    seq = "A" * 50
    # `rec` writes SEQ/QUAL as "*"; give real ones so samtools depth sees bases.
    lines = [
        ln.replace("\t*\t0\t0\t*\t*", f"\t*\t0\t0\t{seq}\t{'I' * 50}")
        for ln in RAW
        if ln.split("\t")[2] != "*"
    ]
    sam.write_text(HEADER + "\n".join(lines) + "\n")
    raw_bam = tmp_path / "raw.bam"
    subprocess.run(["samtools", "sort", "-o", str(raw_bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(raw_bam)], check=True)
    dedup_bam = evidence.deduplicate_bam(str(raw_bam), str(tmp_path / "d.bam"), "umi")
    cov = {r["rname"]: int(r["numreads"]) for r in evidence.coverage_table(dedup_bam)}
    qc = {str(r["reference"]): int(r["reads"]) for r in evidence.alignment_qc_table(dedup_bam)}
    starts = {
        str(r["reference"]): int(r["n_reads"])
        for r in evidence.read_start_distribution(dedup_bam, dedup="none")
    }
    assert cov == qc == starts == EXPECTED
    # raw layer follows the same flag policy: supplementary and duplicate excluded.
    raw_cov = {r["rname"]: int(r["numreads"]) for r in evidence.coverage_table(str(raw_bam))}
    raw_qc = {
        str(r["reference"]): int(r["reads"]) for r in evidence.alignment_qc_table(str(raw_bam))
    }
    assert raw_cov == raw_qc == {"VIRUS|v": 4, "HOST|h": 2}


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")
def test_deletion_policy_is_the_same_for_coverage_and_depth(tmp_path: Path):
    # 25M1D25M spans 51 reference bases, 1 of them deleted: both samtools
    # views must agree that the deleted base is not covered (50 covered bases).
    seq = "A" * 50
    line = rec("CB1_U1_1", 0, "VIRUS|v", 100, cigar="25M1D25M").replace(
        "\t*\t0\t0\t*\t*", f"\t*\t0\t0\t{seq}\t{'I' * 50}"
    )
    sam = tmp_path / "d.sam"
    sam.write_text(HEADER + line + "\n")
    bam = tmp_path / "d.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    (cov,) = evidence.coverage_table(str(bam))
    covered_by_depth = [p for p in evidence.coverage_depth_points(str(bam)) if p["depth"] > 0]
    assert int(cov["covbases"]) == len(covered_by_depth) == 50
