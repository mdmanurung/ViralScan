"""Q137 molecule-support flags (src/viralscan/read_support.py) and the UMI-deduplication regression.

The fixture is the PCR-duplicate case: 5 reads of one molecule plus 1 read of another. Every count that feeds a
verdict, flag, threshold, summary number or plot must say 2 molecules, never 6 reads.
"""

from __future__ import annotations

import random
import shutil
import subprocess

import pytest

from viralscan import evidence as ev
from viralscan import molecule_verdict as mv
from viralscan import read_support as rs

SEQ = "".join(random.Random(3).choices("ACGT", k=90))


def _sam(name: str, flag: int, ref: str, pos: int, cigar: str, score: int, seq: str = SEQ) -> str:
    return "\t".join(
        [
            name,
            str(flag),
            ref,
            str(pos),
            "60",
            cigar,
            "*",
            "0",
            "0",
            seq,
            "I" * len(seq),
            f"AS:i:{score}",
        ]
    )


def _duplicates() -> list[str]:
    """Molecule U1 sequenced 5 times (PCR duplicates), molecule U2 once; both fully virus-aligned."""
    names = [f"CB_U1_{i}|r" for i in range(1, 6)] + ["CB_U2_6|r"]
    return [_sam(n, 0, "VIRUS|v", 100, "90M", 150) for n in names]


def test_dust_masks_a_homopolymer_and_leaves_random_sequence() -> None:
    masked = rs.dust_mask("ACGTTGCAGC" + "A" * 30 + "TGCATCGGAT")
    assert masked[:10].isupper() and masked[40:].isupper() and masked[10:40].islower()
    assert rs.dust_mask(SEQ).isupper()


@pytest.mark.skipif(not shutil.which("sdust"), reason="sdust not installed")
def test_dust_matches_sdust() -> None:
    rng = random.Random(5)
    seqs = []
    for _ in range(150):
        unit = rng.choice(["A", "CAG", "AT", "TTTGC", "G", "ACGTT"])
        seqs.append(
            "".join(rng.choices("ACGT", k=rng.randint(5, 40)))
            + unit
            * rng.randint(
                3, 60 // len(unit)
            )  # a repeat longer than the window ends differently in sdust
            + "".join(rng.choices("ACGT", k=rng.randint(5, 40)))
        )
    fasta = "".join(f">{i}\n{s}\n" for i, s in enumerate(seqs))
    out = subprocess.run(
        ["sdust", "-w", "64", "-t", "20", "-"],
        input=fasta,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    expected: list[set[int]] = [set() for _ in seqs]
    for line in out.splitlines():
        i, a, b = line.split("\t")[:3]
        expected[int(i)].update(range(int(a), int(b)))
    for s, want in zip(seqs, expected, strict=True):
        assert {i for i, c in enumerate(rs.dust_mask(s)) if c.islower()} == want


def test_pcr_duplicates_are_two_molecules_in_every_count() -> None:
    lines = _duplicates()
    reps, _by_read, n_reads = mv.select_representatives(lines)
    assert len(reps) == 2 and sum(n_reads.values()) == 6
    assert len(mv.deduplicate_competitive_alignments(lines)) == 2
    assert len(list(mv.molecule_rows(lines, set()))) == 2  # verdicts: one row per molecule
    rows = rs.molecule_support_rows(lines)
    summary = rs.summarise_support(rows)
    assert (summary["n_molecules"], summary["n_reads"], summary["n_pass_all_gates"]) == (2, 6, 2)
    # alignment QC: the raw layer is labelled reads, the deduplicated layer is the molecules
    header = "@SQ\tSN:VIRUS|v\tLN:1000\n"
    raw = ev._alignment_qc_from_text(header, "\n".join(lines), "")
    dedup = ev._alignment_qc_from_text(
        header, "\n".join(mv.deduplicate_competitive_alignments(lines)), ""
    )
    assert (raw[0]["reads"], raw[0]["molecules"]) == (6, 2)
    assert (dedup[0]["reads"], dedup[0]["molecules"]) == (2, 2)
    assert ev._per_cell_qc_from_text("\n".join(lines))[0]["molecules"] == 2
    starts = ev._parse_sam_read_starts("\n".join(lines), dedup="umi")
    assert starts[0]["n_reads"] == 2 and starts[0]["n_read_starts"] == 2


def test_interpretation_flags_ignore_the_raw_layer() -> None:
    def row(layer: str, reads: int, window: float) -> dict[str, object]:
        return {
            "count_layer": layer,
            "reference_class": "virus",
            "reads": reads,
            "max_50bp_window_fraction": window,
            "breadth_1x": 0.5,
            "host_competitive_fraction": 0.0,
        }

    # 600 duplicated reads on one start look like a hotspot; the 2 molecules do not
    flags = {
        r["flag"]: r
        for r in ev.interpretation_flags(
            [row("raw", 600, 1.0), row("deduplicated", 2, 0.5)], [], []
        )
    }
    assert (
        flags["coverage_hotspot"]["status"] == "not_assessed"
    )  # <= 2 molecules: uncalibrated, not flagged
    flags = {
        r["flag"]: r
        for r in ev.interpretation_flags(
            [row("raw", 600, 1.0), row("deduplicated", 3, 0.3)], [], []
        )
    }
    assert (
        flags["coverage_hotspot"]["status"] == "not_flagged"
        and flags["coverage_hotspot"]["metric"] == 0.3
    )


def test_gates_pass_a_clean_molecule_and_flag_short_low_complexity_and_human() -> None:
    clean = _sam("CB_U1_1|r", 0, "VIRUS|v", 100, "90M", 150)
    short = _sam("CB_U2_2|r", 0, "VIRUS|v", 100, "40M50S", 80)
    poly = "ACGTT" + "A" * 80 + "TGCAT"
    lowc = _sam("CB_U3_3|r", 0, "VIRUS|v", 100, "90M", 150, poly)
    human = [
        _sam("CB_U4_4|r", 0, "HOST|chr1", 5, "90M", 160),
        _sam("CB_U4_4|r", 256, "VIRUS|v", 7, "90M", 150, "*"),
    ]
    rows = {r["cb"] + r["umi"]: r for r in rs.molecule_support_rows([clean, short, lowc, *human])}
    assert rows["CBU1"]["passes_all_gates"] and rows["CBU1"]["dust_unmasked_bases"] == 90
    assert not rows["CBU2"]["gate_aligned_length"] and rows["CBU2"]["aligned_bases"] == 40
    assert not rows["CBU3"]["gate_dust"] and rows["CBU3"]["aligned_AT_warning"]
    assert not rows["CBU4"]["gate_competitive"] and rows["CBU4"]["verdict"] == "host_best"
    assert len(rows) == 4  # flags remove no molecule


def test_hotspot_and_dnajb6_from_the_data_file_and_call_level_known_artefact() -> None:
    def hot(i: int, host_pos: int) -> list[str]:
        n = f"CB_U{i}_{i}|r"
        return [
            _sam(n, 0, "HOST|chr7", host_pos, "90M", 170),
            _sam(n, 256, "VIRUS|KP343820.1", 2798, "31M59S", 62, "*"),
        ]

    lines = hot(1, 157386214) + hot(2, 157386240)
    rows = rs.molecule_support_rows(lines)
    assert all(r["host_decoy"] and r["hotspot_only"] for r in rows)
    summary = rs.summarise_support(rows)
    assert summary["known_artefact"] == "human_DNAJB6;low_complexity_hotspot"
    assert summary["n_pass_all_gates"] == 0
    # one molecule that is neither: no reason covers every molecule of the call
    clean = [_sam("CB_U9_9|r", 0, "VIRUS|v", 100, "90M", 150)]
    assert rs.summarise_support(rs.molecule_support_rows(lines + clean))["known_artefact"] == ""
    assert {r["kind"] for r in rs.load_loci()} == {"hotspot", "host_decoy"}


def test_reverse_strand_alignment_uses_original_read_coordinates() -> None:
    rc = SEQ.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    row = rs.molecule_support_rows([_sam("CB_U1_1|r", 16, "VIRUS|v", 100, "20S70M", 150, rc)])[0]
    assert row["aligned_bases"] == 70 and row["trimmed_query_bases"] == 90
