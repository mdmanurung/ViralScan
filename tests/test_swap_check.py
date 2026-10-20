"""M-03 tier 1: scripts/swap_check.py barcode-overlap test, hand-computed hypergeometric values."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "swap_check", Path(__file__).parent.parent / "scripts" / "swap_check.py"
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def _run(root, name, called, rows):
    """A finished-run stand-in: ten observed barcodes b0..b9, the called ones, per-cell viral rows."""
    out = root / name
    (out / "kb-python" / "counts_unfiltered").mkdir(parents=True)
    (out / "results").mkdir()
    (out / "kb-python" / "counts_unfiltered" / "cells_x_genes.barcodes.txt").write_text(
        "".join(f"b{i}-1\n" for i in range(10))
    )
    (out / "results" / "called_cells.tsv").write_text(
        "barcode\n" + "".join(f"{b}-1\n" for b in sorted(called))
    )
    (out / "results" / "per_cell_viral.tsv").write_text(
        "barcode\tvirus_name\tviral_molecules_total_est\tmolecules_total_est\tviral_fraction"
        "\tis_called_cell\n"
        + "".join(f"{b}-1\t{v}\t1\t100\t0.01\t{b in called}\n" for b, v in rows)
    )
    return out


def _manifest(tmp_path, entries):
    path = tmp_path / "pool.tsv"
    path.write_text(
        "library_id\trun_dir\tpool_id\twhitelist\n" + "".join("\t".join(e) + "\n" for e in entries)
    )
    return path


def test_hypergeometric_columns_by_hand():
    # M=10 observed non-called barcodes, K=2 of them donor barcodes, n=3 recipient viral barcodes,
    # both of the donor's inside: P(X >= 2) = C(2,2) C(8,1) / C(10,3) = 8/120.
    row = mod.overlap_row(n_s=3, n_donor=5, universe_size=10, n_u_in_d=2, overlap=2)
    assert row["overlap_expected"] == pytest.approx(3 * 2 / 10)
    assert row["overlap_p"] == pytest.approx(8 / 120, abs=1e-6)
    assert row["swap_assessment"] == "not_above_chance"  # 0.0667 >= 0.05
    # M=20, K=4, n=4, all four overlap: P = 1 / C(20,4) = 1/4845
    row = mod.overlap_row(4, 4, 20, 4, 4)
    assert row["overlap_p"] == pytest.approx(1 / 4845, rel=1e-4)
    assert row["swap_assessment"] == "overlap_above_chance"
    # no overlap: P(X >= 0) = 1; no donor barcode observed unclaimed: no power
    assert mod.overlap_row(3, 5, 10, 2, 0)["overlap_p"] == 1.0
    assert mod.overlap_row(3, 0, 10, 0, 0)["swap_assessment"] == "not_assessed"


def test_both_directions_of_a_pool(tmp_path):
    """Recipient R: called b0,b1; viral b0 (called), b2,b3,b4 (not). Donor D: called b2; viral b2,b3,b9.

    R<-D: U_R = b2..b9 (8), D's virus barcodes {b2,b3,b9} are all in U_R (K=3), S_R = {b2,b3,b4},
    overlap {b2,b3} = 2: expected 3*3/8 = 1.125, p = (C(3,2) C(5,1) + C(3,3)) / C(8,3) = 16/56.
    D<-R: U_D = b0,b1,b3..b9 (9), R's virus barcodes {b0,b2,b3,b4} -> K = |{b0,b3,b4}| = 3,
    S_D = {b3,b9} (b2 is called in D), overlap {b3} = 1: expected 2*3/9, p = 1 - C(6,2)/C(9,2) = 21/36.
    """
    r = _run(tmp_path, "R", {"b0", "b1"}, [("b0", "V"), ("b2", "V"), ("b3", "V"), ("b4", "V")])
    d = _run(tmp_path, "D", {"b2"}, [("b2", "V"), ("b3", "V"), ("b9", "V")])
    manifest = _manifest(tmp_path, [("R", r.name, "pool1", "wl"), ("D", d.name, "pool1", "wl")])
    table = mod.swap_rows(mod.read_manifest(manifest))
    rows = {(x["recipient_library"], x["donor_library"]): x for x in table}
    assert set(rows) == {("R", "D"), ("D", "R")}
    rd = rows[("R", "D")]
    assert rd["n_noncalled_virus_barcodes"] == 3 and rd["n_noncalled_in_donor"] == 3
    assert rd["overlap_with_donor"] == 2
    assert rd["overlap_expected"] == pytest.approx(1.125)
    assert rd["overlap_p"] == pytest.approx(16 / 56, abs=1e-6)
    dr = rows[("D", "R")]
    assert dr["n_noncalled_virus_barcodes"] == 2 and dr["n_noncalled_in_donor"] == 3
    assert dr["overlap_with_donor"] == 1
    assert dr["overlap_expected"] == pytest.approx(2 * 3 / 9, abs=1e-6)
    assert dr["overlap_p"] == pytest.approx(21 / 36, abs=1e-6)


def test_different_pool_or_whitelist_is_never_compared(tmp_path):
    r = _run(tmp_path, "R", {"b0"}, [("b2", "V")])
    d = _run(tmp_path, "D", {"b0"}, [("b2", "V")])
    for other in (("D", d.name, "pool2", "wl"), ("D", d.name, "pool1", "other_wl")):
        manifest = _manifest(tmp_path, [("R", r.name, "pool1", "wl"), other])
        assert mod.swap_rows(mod.read_manifest(manifest)) == []


def test_main_writes_the_table_and_the_header_is_stable(tmp_path):
    r = _run(tmp_path, "R", {"b0"}, [("b2", "V")])
    d = _run(tmp_path, "D", {"b0"}, [("b2", "V")])
    manifest = _manifest(tmp_path, [("R", r.name, "p", "w"), ("D", d.name, "p", "w")])
    out = tmp_path / "out" / "swap.tsv"
    assert mod.main([str(manifest), "--out", str(out)]) == 0
    with open(out, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        assert reader.fieldnames == mod.COLUMNS
        assert len(list(reader)) == 2


@pytest.mark.parametrize(
    "entries, message",
    [
        ([("L", "a", "p", "w"), ("L", "b", "p", "w")], "duplicate library_id"),
        ([("L1", "a", "p", "w"), ("L2", "a", "p", "w")], "share one run_dir"),
        ([("L1", "a", "", "w")], "empty field"),
    ],
)
def test_manifest_errors_are_loud(tmp_path, entries, message):
    with pytest.raises(ValueError, match=message):
        mod.read_manifest(_manifest(tmp_path, entries))


def test_manifest_needs_its_columns(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("library_id\trun_dir\nL\ta\n")
    with pytest.raises(ValueError, match="missing column"):
        mod.read_manifest(path)


def test_a_run_dir_that_is_not_finished_is_an_error(tmp_path):
    (tmp_path / "R").mkdir()
    manifest = _manifest(tmp_path, [("R", "R", "p", "w"), ("D", "R2", "p", "w")])
    with pytest.raises(FileNotFoundError, match="not a finished ViralScan run"):
        mod.swap_rows(mod.read_manifest(manifest))


def test_tier2_planted_overlap_null_changes_and_different_target():
    rec = {
        "keys": {("C0", "U"): {"V"}, **{(f"C{i}", "H"): set() for i in range(1, 100)}},
        "barcodes": [f"C{i}" for i in range(100)],
    }
    don = {"keys": {("C0", "U"): {"V"}}}
    result = mod.molecule_overlap(rec, don, "V")
    assert result["observed_shared_keys"] == 1
    assert 0 < result["expected_shared_keys"] < 1
    assert 0 < result["p_value"] < 0.05
    assert result == mod.molecule_overlap(rec, don, "V")
    don["keys"][("C0", "U")] = {"other"}
    result = mod.molecule_overlap(rec, don, "V")
    assert result["observed_shared_keys"] == 0


def test_tier2_unrelated_pools_never_read_bus(monkeypatch):
    monkeypatch.setattr(mod, "load_molecule_keys", lambda _: pytest.fail("unadmitted run read"))
    manifest = [
        {"library_id": name, "pool_id": pool, "whitelist": "w", "run_dir": name}
        for name, pool in [("A", "one"), ("B", "two")]
    ]
    assert mod.tier2_rows(manifest)[0] == []
    manifest[1].update(pool_id="one", whitelist="other")
    assert mod.tier2_rows(manifest)[0] == []


def test_tier2_missing_bus_and_degenerate_null(tmp_path):
    assert mod.load_molecule_keys(tmp_path)["reason_codes"] == "bus_missing"
    run = {"keys": {("A", "U"): {"V"}}, "barcodes": ["A"]}
    assert mod.molecule_overlap(run, run, "V")["status"] == "not_assessed"
    assert mod.molecule_overlap(run, run, "V")["p_value"] == "not_assessed"


def test_native_tier2_bus_preference_fallback_and_target_resolution(tmp_path):
    import shutil
    import struct

    from viralscan.runconfig import RunConfig
    from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

    if not shutil.which("bustools"):
        pytest.skip("native bustools required")
    kb = tmp_path / "kb-python"
    counts = kb / "counts_unfiltered"
    counts.mkdir(parents=True)
    (counts / "cells_x_genes.genes.txt").write_text("V\nH\n")
    (kb / "transcripts.txt").write_text("tv\nth\n")
    (kb / "matrix.ec").write_text("0\t0\n1\t1\n")
    t2g = tmp_path / "t2g.txt"
    t2g.write_text("tv\tV\nth\tH\n")
    RunConfig(output=str(tmp_path), transcripts=str(t2g)).to_yaml(tmp_path / "config.yaml")
    VirusIdentityTable(
        genes=(
            GeneIdentity(
                gene_id="V",
                genome_accession="V.1",
                status="catalogued",
                viral=True,
                virus_name="Virus",
                virus_key="v",
            ),
        )
    ).write_tsv(tmp_path / "results" / "virus_identity.tsv")
    data = b"BUS\x00" + struct.pack("<IIII", 1, 16, 12, 0)
    data += struct.pack("<QQiIII", 0, 0, 0, 1, 0, 0) + struct.pack("<QQiIII", 1, 1, 1, 1, 0, 0)
    resolved = kb / "output.resolved.sorted.bus"
    corrected = kb / "output.corrected.bus"
    resolved.write_bytes(data)
    corrected.write_bytes(b"invalid, must not read when resolved exists")
    run = mod.load_molecule_keys(tmp_path)
    assert run["bus"] == str(resolved)
    assert run["keys"][("A" * 16, "A" * 12)] == {"Virus"}
    assert run["keys"][("A" * 15 + "C", "A" * 11 + "C")] == set()
    corrected.write_bytes(data)
    resolved.unlink()
    assert mod.load_molecule_keys(tmp_path)["bus"] == str(corrected)


def test_tier2_receipt_does_not_overwrite_tsv_named_json(tmp_path):
    a = _run(tmp_path, "A", {"b0"}, [("b2", "V")])
    b = _run(tmp_path, "B", {"b0"}, [("b2", "V")])
    manifest = _manifest(tmp_path, [("A", a.name, "p", "w"), ("B", b.name, "p", "w")])
    out = tmp_path / "fingerprint.json"
    assert mod.main([str(manifest), "--tier2", "--out", str(out)]) == 0
    assert out.read_text().startswith("pool_id\t")
    assert (tmp_path / "fingerprint.json.json").is_file()
