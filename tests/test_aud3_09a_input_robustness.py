"""AUD3-09a: four small input-robustness fixes (P3-06, P4-04, P4-06, P4-03)."""

from __future__ import annotations

import logging
import random
from pathlib import Path

import pytest

from viralscan import chemistry
from viralscan import read_support as rs
from viralscan import whitelist_preflight as wp
from viralscan.scripts.host_filter import starsolo_barcode_args

SEQ = "".join(random.Random(3).choices("ACGT", k=90))


def _sam(name: str, seq: str = SEQ) -> str:
    qual = "I" * len(seq) if seq != "*" else "*"
    fields = [name, "0", "VIRUS|v", "100", "60", "90M", "*", "0", "0", seq, qual, "AS:i:150"]
    return "\t".join(fields)


# P3-06 ---------------------------------------------------------------------
def test_seq_star_viral_primary_is_skipped_with_one_warning(
    caplog: pytest.LogCaptureFixture,
) -> None:
    sam = [_sam("CB_U1_1|r", "*"), _sam("CB_U2_2|r"), _sam("CB_U3_3|r", "*")]
    with caplog.at_level(logging.WARNING, logger="viralscan"):
        rows = rs.molecule_support_rows(sam)
    assert [r["umi"] for r in rows] == ["U2"]
    warnings = [m for m in caplog.messages if "SEQ=*" in m]
    assert len(warnings) == 1 and "2 " in warnings[0]


def test_seq_star_with_queries_fasta_is_kept() -> None:
    rows = rs.molecule_support_rows([_sam("CB_U1_1|r", "*")], {"CB_U1_1|r": SEQ})
    assert len(rows) == 1 and rows[0]["aligned_bases"] == 90


# P4-04 ---------------------------------------------------------------------
def test_suffixed_whitelist_matches_like_plain(tmp_path: Path) -> None:
    rng = random.Random(1)
    bcs = ["".join(rng.choices("ACGT", k=16)) for _ in range(50)]
    (tmp_path / "wl.txt").write_text("\n".join(bcs) + "\n")
    (tmp_path / "wl_suf.txt").write_text("\n".join(b + "-1" for b in bcs) + "\n")
    reads = "".join(f"@r{i}\n{bcs[i % 50]}{'A' * 10}\n+\n{'I' * 26}\n" for i in range(200))
    (tmp_path / "r1.fq").write_text(reads)
    plain = wp.check_whitelist(str(tmp_path / "r1.fq"), str(tmp_path / "wl.txt"), "10xv2")
    suffixed = wp.check_whitelist(str(tmp_path / "r1.fq"), str(tmp_path / "wl_suf.txt"), "10xv2")
    assert plain.ok and suffixed.ok
    assert suffixed.match_rate == plain.match_rate == 1.0


def test_load_whitelist_strips_only_terminal_gem_group(tmp_path: Path) -> None:
    p = tmp_path / "wl.txt"
    p.write_text("AAAA-1\nCCCC-12\nGGGG\nTTTT-0\n")
    assert wp.load_whitelist(str(p)) == {"AAAA", "CCCC", "GGGG", "TTTT-0"}


# P4-06 ---------------------------------------------------------------------
@pytest.mark.parametrize(
    "geometry",
    [
        "1,0,16:1,16,26:0,0,0",  # barcode + UMI on file 1
        "0,8,24:0,24,34:1,0,0",  # barcode at R1 offset 8
        "0,0,16:1,0,10:0,0,0",  # UMI on another file
        "0,0,16:0,20,30:1,0,0",  # gap between barcode and UMI
    ],
)
def test_explicit_geometry_that_downstream_cannot_honour_is_refused(geometry: str) -> None:
    with pytest.raises(ValueError, match="file index|offset|start"):
        chemistry.cb_umi_geometry(geometry)
    with pytest.raises(ValueError):
        starsolo_barcode_args(geometry, None)


def test_supported_geometries_are_unchanged() -> None:
    assert chemistry.cb_umi_geometry("0,0,16:0,16,28:1,0,0") == (16, 12)
    assert chemistry.cb_umi_geometry("0,0,12:0,12,20:1,0,0") == (12, 8)
    for name, c in chemistry.CHEMISTRIES.items():
        assert chemistry.cb_umi_geometry(name) == (c.cb_len, c.umi_len)
    args = starsolo_barcode_args("10xv3", None)
    assert args[args.index("--soloCBstart") + 1] == "1"
    assert args[args.index("--soloUMIstart") + 1] == "17"
    assert args[args.index("--soloUMIlen") + 1] == "12"


def test_resolve_still_compares_lengths_of_explicit_geometry() -> None:
    """Plain kb count takes -x verbatim, so a file-1 geometry must not fail the chemistry check."""
    assert chemistry._same("1,0,16:1,16,26:0,0,0", "10xv2", "user on-list")
