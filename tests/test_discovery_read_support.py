"""DISC-13: the route (ii) ``source=minimap2`` read-support producer, with a stubbed aligner."""

import gzip
import json
from pathlib import Path

import pytest

from viralscan import discovery

# ordinal, read_id, corrected_cb, umi
READS = [
    (0, "r0", "C1", "U1"),  # T10 resolved called: unique candidate win
    (1, "r1", "C1", "U1"),  # same molecule: candidate win
    (2, "r2", "C1", "U2"),  # tie with host
    (3, "r3", "C2", "U3"),  # host win
    (4, "r4", "C2", "U4"),  # unmapped
    (5, "r5", "C3", "U5"),  # T10 group but C3 is not called: never extracted
    (6, "r6", "C1", "U6"),  # T20 resolved called
    (7, "r7", "C2", "U7"),  # ambiguous group: ignored
    (8, "r8", "C2", "U8"),  # T10 candidate win in a second cell
]
GROUPS = [
    ("C1", "U1", 10, "resolved"),
    ("C1", "U2", 10, "resolved"),
    ("C2", "U3", 10, "resolved"),
    ("C2", "U4", 10, "resolved"),
    ("C3", "U5", 10, "resolved"),
    ("C1", "U6", 20, "resolved"),
    ("C2", "U7", 10, "ambiguous"),
    ("C2", "U8", 10, "resolved"),
]
# read -> {reference: AS}
HITS = {
    0: {"VIRUS|CAND_v10": 90, "HOST|chr1": 40},
    1: {"VIRUS|CAND_v10": 80},
    2: {"VIRUS|CAND_v10": 70, "HOST|chr1": 70},
    3: {"VIRUS|CAND_v10": 50, "HOST|chr1": 60},
    8: {"VIRUS|CAND_v10": 60, "VIRUS|PANEL_p": 30},
    6: {"VIRUS|CAND_v20": 99},
}


def _setup(tmp_path: Path, space: str = "X"):
    disc, sample = tmp_path / "disc", tmp_path / "sample"
    (sample / "results").mkdir(parents=True)
    disc.mkdir()
    (disc / discovery.CANDIDATES_TSV).write_text("taxid\n10\n20\n30\n")
    with gzip.open(disc / discovery.READS_TSV, "wt") as h:
        h.write("ordinal\tread_id\n")
        h.writelines(f"{o}\t{rid}\n" for o, rid, _, _ in READS)
    with gzip.open(disc / discovery.BARCODES_TSV, "wt") as h:
        h.write("ordinal\tcorrected_cb\tumi\tcb_status\n")
        h.writelines(f"{o}\t{cb}\t{u}\tonlist\n" for o, _, cb, u in READS)
    with gzip.open(disc / discovery.UMI_GROUPS_TSV, "wt") as h:
        h.write("\t".join(discovery.UMI_GROUP_COLUMNS) + "\n")
        for cb, umi, tax, cls in GROUPS:
            h.write(f"{cb}\t{umi}\t1\t1\t0\t{tax}\t{tax}\t{cls}\t1\t\t0\n")
    (sample / "results" / "called_cells.tsv").write_text("barcode\nC1-1\nC2-1\n")
    (sample / "results" / "cell_calling_input.json").write_text(json.dumps({"count_layer": space}))
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    r1.write_text("".join(f"@{rid}/1\nAAAA\n+\nIIII\n" for _, rid, _, _ in READS))
    r2.write_text("".join(f"@{rid}/2\nACGTACGT\n+\nIIIIIIII\n" for _, rid, _, _ in READS))
    fastas = {}
    for name in ("host", "c10", "c20", "panel"):
        fastas[name] = tmp_path / f"{name}.fa"
        fastas[name].write_text(
            f">{'chr1' if name == 'host' else 'v' + name[1:] if name != 'panel' else 'p'}\nACGT\n"
        )
    return disc, sample, r1, r2, fastas


def _stubs(seen: list):
    def aligner(reads_fa, ref_fa, bam):
        seen.append((Path(reads_fa).read_text(), Path(ref_fa).read_text()))
        return reads_fa

    def hit_reader(reads_fa):
        out = {}
        for line in Path(reads_fa).read_text().splitlines():
            if line.startswith(">"):
                name = line[1:]
                number = int(name.split("|")[0].rsplit("_", 1)[1])
                if number in HITS:
                    out[name] = HITS[number]
        return out

    return aligner, hit_reader


def _run(tmp_path, cands=None, **kw):
    disc, sample, r1, r2, f = _setup(tmp_path)
    seen: list = []
    aligner, hit_reader = _stubs(seen)
    cands = cands or {10: f["c10"]}
    out = discovery.build_read_support_evidence(
        disc, sample, r1, r2, f["host"], cands, aligner=aligner, hit_reader=hit_reader, **kw
    )
    return disc, out, seen, f


def test_counts_and_roundtrip(tmp_path):
    disc, out, seen, _ = _run(tmp_path)
    assert out == disc / discovery.READ_SUPPORT_EVIDENCE_TSV
    row = discovery.read_evidence_table(out)[(10, "read_support")]
    # reads 0,1,2,3,4,8 extracted (5 uncalled, 7 ambiguous, 6 other taxon); 0,1,8 best;
    # tie (2), host win (3) and unmapped (4) stay adjudicated only.
    assert (row["n_adjudicated"], row["n_taxon_best"]) == (6, 3)
    assert (row["n_shared_cells"], row["n_shared_umis"]) == (2, 2)  # C1/U1 and C2/U8
    assert row["source"] == "minimap2" and row["source"] not in discovery.CORROBORATION_SOURCES
    fasta = seen[0][0]
    assert ">C1_U1_0|r0\n" in fasta and "r5" not in fasta and "r7" not in fasta
    assert ">HOST|chr1" in seen[0][1] and ">VIRUS|CAND_v10" in seen[0][1]


def test_drives_read_supported_observation(tmp_path):
    disc, out, _, _ = _run(tmp_path, panel_fasta=None)
    # 3/6 best is below the 80% rule: negative, not not_assessed
    cand = {"taxid": "10", "name": "n", "rank": "S", "host_lineage": "0", "viral_only_db": "0"}
    ev = discovery.read_evidence_table(out)
    assert discovery.assign_state(cand, ev)["read_support"] == "negative"
    ev[(10, "read_support")].update(n_taxon_best=6, n_shared_umis=3)
    assert discovery.assign_state(cand, ev)["state"] == "read_supported_observation"
    assert discovery.assign_state({**cand, "taxid": "30"}, ev)["read_support"] == "not_assessed"


def test_panel_competitor_and_no_row_without_reads(tmp_path):
    disc, sample, r1, r2, f = _setup(tmp_path)
    seen: list = []
    aligner, hit_reader = _stubs(seen)
    out = discovery.build_read_support_evidence(
        disc,
        sample,
        r1,
        r2,
        f["host"],
        {10: f["c10"], 30: f["c10"]},
        panel_fasta=f["panel"],
        aligner=aligner,
        hit_reader=hit_reader,
    )
    assert ">VIRUS|PANEL_p" in seen[0][1]
    assert set(discovery.read_evidence_table(out)) == {(10, "read_support")}  # 30: not_assessed


def test_failures(tmp_path):
    with pytest.raises(discovery.DiscoveryError, match="max_candidates"):
        _run(tmp_path / "a", {1: Path("x"), 2: Path("x"), 3: Path("x"), 4: Path("x")})
    with pytest.raises(discovery.DiscoveryError, match="max_candidates=1"):
        disc, sample, r1, r2, f = _setup(tmp_path / "b")
        discovery.build_read_support_evidence(
            disc, sample, r1, r2, f["host"], {10: f["c10"], 20: f["c20"]}, max_candidates=1
        )
    disc, sample, r1, r2, f = _setup(tmp_path / "c")
    (disc / discovery.BARCODES_TSV).unlink()
    with pytest.raises(discovery.DiscoveryError, match="is missing"):
        discovery.build_read_support_evidence(disc, sample, r1, r2, f["host"], {10: f["c10"]})
    disc, sample, r1, r2, f = _setup(tmp_path / "d")
    with pytest.raises(discovery.DiscoveryError, match="not in candidates"):
        discovery.build_read_support_evidence(disc, sample, r1, r2, f["host"], {99: f["c10"]})
    disc, sample, r1, r2, f = _setup(tmp_path / "e", space="GeneFull.raw")
    with pytest.raises(discovery.DiscoveryError, match="not kb"):
        discovery.build_read_support_evidence(disc, sample, r1, r2, f["host"], {10: f["c10"]})
    disc, sample, r1, r2, f = _setup(tmp_path / "f")
    r2.write_text(r2.read_text().replace("@r0/2", "@zz/2"))
    with pytest.raises(ValueError, match="mismatch|differs"):
        discovery.build_read_support_evidence(disc, sample, r1, r2, f["host"], {10: f["c10"]})
