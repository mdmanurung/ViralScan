"""PANEL-01 WP1b: scripts/panel_kmer_sharing.py counts shared k-mers on tiny synthetic genomes."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import pytest

pytest.importorskip("numpy")

SPEC = importlib.util.spec_from_file_location(
    "panel_kmer_sharing", Path(__file__).resolve().parents[1] / "scripts" / "panel_kmer_sharing.py"
)
pks = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pks)

K = 11


def _seq(seed: int, n: int = 200) -> bytes:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n)).encode()


def _revcomp(seq: bytes) -> bytes:
    return seq.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]


def test_canonical_kmers_strand_and_ambiguity() -> None:
    s = _seq(1)
    assert (pks.canonical_kmers(s, K) == pks.canonical_kmers(_revcomp(s), K)).all()
    # the N kills every window that spans it
    assert (
        pks.canonical_kmers(b"A" * 12 + b"N" + b"A" * 12, K).tolist()
        == pks.canonical_kmers(b"A" * 12, K).tolist()
    )


def test_sharing_fractions() -> None:
    a = _seq(1)
    panel = [a]  # panel genome, group 0
    twin_strand = _revcomp(a)  # candidate sharing everything with the panel, same group
    unrelated = _seq(2)  # candidate sharing nothing, own group
    cross = a[:100] + _seq(3, 100)  # half shared with a different-group panel genome
    seqs = panel + [twin_strand, unrelated, cross]
    kmers = [pks.canonical_kmers(s, K) for s in seqs]
    stats = pks.sharing(kmers, [0, 0, 1, 2], [True, False, False, False], [1, 2, 3])
    assert stats[1]["frac_sibling"] == 1.0 and stats[1]["frac_not_in_panel"] == 0.0
    assert stats[2]["frac_other_group"] < 0.1 and stats[2]["frac_unique"] > 0.9
    assert 0.4 < stats[3]["frac_other_group"] < 0.7 and stats[3]["frac_sibling"] == 0.0


def test_partners_rank_the_genome_that_shares_most() -> None:
    a, b = _seq(1), _seq(2)
    near_a = a[:150] + _seq(5, 50)  # shares 150 of 200 bases with a, nothing with b
    kmers = [pks.canonical_kmers(s, K) for s in (a, b, near_a)]
    stats = pks.sharing(kmers, [0, 1, 2], [True, True, False], [2], partner_below=0.9)
    (partner, shared), *rest = stats[2]["partners"]
    assert partner == 0 and shared > 100 and not rest  # b shares nothing, so it is not listed
    assert "partners" not in pks.sharing(kmers, [0, 1, 2], [True, True, False], [2])[2]


def test_partner_columns_cover_the_twin_versus_fragment_split() -> None:
    assert {"frac_of_partner", "length_ratio"} <= set(pks.PARTNER_COLUMNS)
    assert {"frac_unique", "frac_sibling", "frac_other_group"} <= set(pks.BASELINE_COLUMNS)


def test_group_of_prefers_ictv_genus_and_falls_back(tmp_path) -> None:
    import zipfile

    ns = "xmlns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'"
    values = {
        "P": "Coltivirus",
        "R": "Coltivirus eyachense",
        "U": "Eyach virus",
        "X": "NC_003698.1",
    }
    strings = f"<sst {ns}>" + "".join(f"<si><t>{v}</t></si>" for v in values.values()) + "</sst>"
    cells = "".join(f'<c r="{col}2" t="s"><v>{i}</v></c>' for i, col in enumerate(values))
    sheet = f"<worksheet {ns}><row>{cells}</row></worksheet>"
    path = tmp_path / "vmr.xlsx"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("xl/sharedStrings.xml", strings)
        z.writestr("xl/worksheets/sheet2.xml", sheet)
    vmr = pks.load_vmr(path)
    by_name = {"species": "Eyach virus segment 9, complete genome", "genus": "Old"}
    assert pks.group_of(by_name, vmr, "NC_000000") == ("Coltivirus", "ictv_genus")  # by name
    assert pks.group_of({"species": "x"}, vmr, "NC_003698") == (
        "Coltivirus",
        "ictv_genus",
    )  # by accession
    assert pks.group_of({"species": "x", "genus": "G"}, vmr, "Z") == ("G", "genus")
    assert pks.group_of({"species": "x"}, vmr, "Z") == ("x", "species")
    assert pks.group_of({"species": "Eyach virus", "sibling_group": "S"}, vmr, "Z") == (
        "S",
        "sibling_group",
    )
    assert pks.group_of({"species": "Eyach virus"}, None, "Z") == ("Eyach virus", "species")
