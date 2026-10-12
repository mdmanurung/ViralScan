"""CAT-19 multimapping contract and TONSIL-02 native 5' recipe, pinned by measured numbers."""

from __future__ import annotations

import random
import re
from pathlib import Path

from viralscan import chemistry as ch
from viralscan.strand import infer_strand

SNAKEFILE = Path(__file__).resolve().parents[1] / "src" / "viralscan" / "Snakefile"


def test_kb_count_keeps_the_unique_layer() -> None:
    """CAT-19: multi-gene UMIs are recovered by ViralScan's own multimapping step (default on, test_cli,
    from the BUS file), not by `kb count --mm`/`--tcc`, which would turn kb's matrix into a
    mixed layer, and not by an EM flag on kb."""
    text = SNAKEFILE.read_text()
    kb_lines = [ln for ln in text.splitlines() if re.search(r"\bkb count\b", ln)]
    assert kb_lines
    for line in kb_lines:
        assert not re.search(r"--(mm|tcc|em)\b", line), line


def test_tonsil_x223_geometry_is_10xv3() -> None:
    """The row's `-x 0,0,16:0,16,28:1,0,0` is the 10xv3 geometry: 16 bp barcode + 12 bp UMI."""
    assert ch.cb_umi_geometry("0,0,16:0,16,28:1,0,0") == ch.cb_umi_geometry("10xv3") == (16, 12)


def test_tonsil_pilot_calls_10xv3_with_the_cellranger_list_and_reverse_strand() -> None:
    """Pilot job 25702321/25702343 (2026-10-06): R1 28 bp, cellranger raw list 73.1 % on 100k
    reads, p_pseudoaligned forward 8.2 / reverse 51.4 / unstranded 59.2 %."""
    rng = random.Random(0)
    reads = ["".join(rng.choice("ACGT") for _ in range(28)) for _ in range(200)]
    rates = {"10xv2": 0.0002, "10xv3": 0.0052, "10xv4": 0.0506, "user": 0.7307}
    got = ch.classify(reads, rates, user_list="cellranger_raw_whitelist.txt")
    assert (got.chemistry, got.basis, got.umi_len) == ("10xv3", "user on-list", 12)
    assert infer_strand({"forward": 8.2, "reverse": 51.4, "unstranded": 59.2}) == "reverse"
    assert ch.kb_whitelist_arg("10xv3", "wl.txt") == "wl.txt"
