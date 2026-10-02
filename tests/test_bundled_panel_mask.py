"""CAT-42: flag-gated Anelloviridae low-complexity 31-mer window mask."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_bundled_panel_ref.py"
_spec = importlib.util.spec_from_file_location("build_bundled_panel_ref", SCRIPT)
bp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bp)

_rng = random.Random(0)
RAND = "".join(_rng.choice("ACGT") for _ in range(300))
POLYA_BROKEN = "AAAAAAAAAATAAAAAAAAATAAAAAAAAAATAAAAAAAAATAAAAAAAAAATAAAAAAAA"  # no run >= 20
GC_TAIL = "GGGGGGGGGGGCCGGGGGGGGGGGGCCCCCCCCCCCC"  # G/C only, no run >= 20
CT_REPEAT = "CT" * 40  # max base 16/31


def _fa(tmp_path, records):
    p = tmp_path / "viral.fa"
    p.write_text("".join(f">{a} test\n{s}\n" for a, s in records))
    return p


def _seqs(p):
    return {n.split()[0]: s for n, s in bp._read_fasta(p)}


def test_interrupted_polya_and_gc_tail_masked_in_anello(tmp_path):
    p = _fa(
        tmp_path,
        [("KP343822.1", RAND + POLYA_BROKEN + RAND), ("KP343842.1", RAND + GC_TAIL + RAND)],
    )
    res = bp._mask_lowcomplexity_kmers(p, {"KP343822.1", "KP343842.1"})
    seqs = _seqs(p)
    a = seqs["KP343822.1"]
    assert set(a[300 : 300 + len(POLYA_BROKEN)]) == {"N"}
    assert res["KP343822.1"] >= len(POLYA_BROKEN)
    assert set(seqs["KP343842.1"][300 : 300 + 31]) == {"N"}  # first window: G23 of 31
    assert len(a) == 600 + len(POLYA_BROKEN)


def test_normal_and_ct_repeat_untouched(tmp_path):
    p = _fa(tmp_path, [("KP1.1", RAND), ("KP2.1", RAND + CT_REPEAT + RAND)])
    assert bp._mask_lowcomplexity_kmers(p, {"KP1.1", "KP2.1"}) == {}
    assert _seqs(p)["KP2.1"] == RAND + CT_REPEAT + RAND


def test_non_anellovirus_record_untouched(tmp_path):
    p = _fa(tmp_path, [("NC_001806.2", RAND + POLYA_BROKEN + RAND)])
    assert bp._mask_lowcomplexity_kmers(p, {"KP343822.1"}) == {}
    assert "N" not in _seqs(p)["NC_001806.2"]


def test_anellovirus_accessions_from_catalogue():
    accs = bp._anellovirus_accessions()
    assert {"KP343822.1", "KP343842.1"} <= accs
    assert not accs & {"NC_001526.4", "NC_007605.1", "NC_001806.2"}


def test_flag_is_off_by_default():
    parser = bp._build_arg_parser()
    assert parser.parse_args(["--out", "x"]).lowcomplexity_kmer_mask is False
    assert parser.parse_args(["--out", "x", "--lowcomplexity-kmer-mask"]).lowcomplexity_kmer_mask
