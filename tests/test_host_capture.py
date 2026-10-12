"""DLIST-02: opt-in ``--host-capture`` split of each virus's reads by host compatibility."""

from __future__ import annotations

import struct
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest
import yaml

from viralscan import evidence
from viralscan.menu import _write_run_config, create_help
from viralscan.run_safety import build_run_manifest
from viralscan.runconfig import RunConfig
from viralscan.scripts import host_capture as hc
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable

# Transcripts: H1, H2 are host; V1, V2 are viral. Classes: 0 = {V1}, 1 = {V2}, 2 = {H1},
# 3 = {H1, V2} (a read compatible with host and virus), 4 = {H2}.
TRANSCRIPTS = ["H1", "H2", "V1", "V2"]
EC_LINES = "0\t2\n1\t3\n2\t0\n3\t0,3\n4\t1\n"


def _write_bus(path: Path, ecs: list[int]) -> None:
    """A BUS file with one record per class in *ecs*; the record's flag is its read number."""
    with open(path, "wb") as fh:
        fh.write(b"BUS\x00" + struct.pack("<IIII", 1, 16, 12, 0))
        for read, ec in enumerate(ecs):
            fh.write(struct.pack("<QQiIII", read, read, ec, 1, read, 0))


def test_complement_capture_drops_reads_whose_class_holds_a_host_transcript(tmp_path):
    """bustools semantics the split rests on: a class is captured when it holds ANY listed
    transcript, so the complement keeps only classes with none of them."""
    bustools = evidence.tool_path("bustools")
    if bustools is None:
        pytest.skip("bustools not on PATH")
    (tmp_path / "transcripts.txt").write_text("\n".join(TRANSCRIPTS) + "\n")
    (tmp_path / "matrix.ec").write_text(EC_LINES)
    _write_bus(tmp_path / "in.bus", [0, 1, 2, 3, 4])
    viral = evidence._capture_flagged_text(
        bustools, tmp_path, tmp_path / "in.bus", ["V1", "V2"], tmp_path / "viral", 1
    )
    only = evidence._capture_flagged_text(
        bustools,
        tmp_path,
        tmp_path / "viral.bus",
        ["H1", "H2"],
        tmp_path / "only",
        1,
        complement=True,
    )

    def flags(path):
        return [int(line.split("\t")[4]) for line in path.read_text().splitlines()]

    assert flags(viral) == [0, 1, 3]  # the mixed class is a viral read
    assert flags(only) == [0, 1]  # ... and the complement removes it


def test_replay_chain_captures_viral_then_complements_the_host(tmp_path, monkeypatch):
    commands: list[list[str]] = []
    monkeypatch.setattr(evidence, "have_tools", lambda tools: [])
    monkeypatch.setattr(evidence, "tool_path", lambda name: name)
    monkeypatch.setattr(evidence, "_run", lambda cmd, **kw: commands.append(list(cmd)) or b"")
    texts = evidence.replay_host_capture_bus(
        index="index.idx",
        technology="10xv3",
        r1_path="R1.fastq.gz",
        r2_path="R2.fastq.gz",
        viral_transcripts=["V2", "V1"],
        host_transcripts=["H1", "H2"],
        workdir=str(tmp_path),
        threads=1,
    )
    assert next(c for c in commands if c[0] == "kallisto").count("-n") == 1
    captures = [c for c in commands if c[:2] == ["bustools", "capture"]]
    assert len(captures) == 2
    first, second = captures
    assert "-x" not in first and first[-1].endswith("output.bus")
    assert "-x" in second and second[-1] == str(tmp_path / "viral_touching.bus")
    assert (tmp_path / "viral_touching_transcripts.txt").read_text() == "V1\nV2\n"
    assert (tmp_path / "virus_only_transcripts.txt").read_text() == "H1\nH2\n"
    assert [t.name for t in texts] == [
        "viral_touching.by_flag.bus.txt",
        "virus_only.by_flag.bus.txt",
    ]


def test_target_replay_still_names_its_files_as_before(tmp_path, monkeypatch):
    commands: list[list[str]] = []
    monkeypatch.setattr(evidence, "have_tools", lambda tools: [])
    monkeypatch.setattr(evidence, "tool_path", lambda name: name)
    monkeypatch.setattr(evidence, "_run", lambda cmd, **kw: commands.append(list(cmd)) or b"")
    out = evidence.replay_exact_target_bus(
        index="i",
        technology="10xv3",
        r1_path="a",
        r2_path="b",
        target_transcripts=["t2", "t1"],
        workdir=str(tmp_path),
        threads=1,
    )
    assert out == tmp_path / "target.by_flag.bus.txt"
    assert (tmp_path / "target_transcripts.txt").read_text() == "t1\nt2\n"


def _rows(ecs_by_read: dict[int, int]) -> list[str]:
    return [f"AAAA\tCCCC\t{ec}\t1\t{read}\n" for read, ec in ecs_by_read.items()]


def test_read_viruses_and_summarise(tmp_path):
    ec_file = tmp_path / "matrix.ec"
    ec_file.write_text(EC_LINES + "5\t2,3\n")  # class 5 holds V1 and V2, no host
    ec_tx = hc.read_ec_transcripts(ec_file)
    tx_virus = {"V1": "Virus one", "V2": "Virus two"}
    viral = hc.read_viruses(_rows({0: 0, 1: 1, 2: 3, 3: 3, 4: 5}), ec_tx, TRANSCRIPTS, tx_virus)
    only = hc.read_viruses(_rows({0: 0, 1: 1, 4: 5}), ec_tx, TRANSCRIPTS, tx_virus)
    assert viral[4] == {"Virus one", "Virus two"}
    rows = {r["virus_name"]: r for r in hc.summarise(viral, only)}
    assert rows["Virus one"] == {
        "virus_name": "Virus one",
        "reads_total": 2,
        "reads_virus_only": 2,
        "reads_also_host": 0,
        "host_capture_status": "measured",
    }
    # Virus two: reads 1 and 4 are virus-only, reads 2 and 3 (host + V2) are also-host.
    assert (rows["Virus two"]["reads_virus_only"], rows["Virus two"]["reads_also_host"]) == (2, 2)
    assert rows["Virus two"]["reads_total"] == 4


def _identity() -> VirusIdentityTable:
    def gene(gid, acc, viral, name):
        return GeneIdentity(gid, acc, "catalogued" if viral else "host", viral, name, name)

    return VirusIdentityTable(
        (
            gene("gV1", "NC_1.1", True, "Virus one"),
            gene("gV2", "NC_2.1", True, "Virus two"),
            gene("gH1", "", False, ""),
            gene("gH2", "", False, ""),
        )
    )


def _ctx(tmp_path: Path, transcripts: list[str], t2g: dict[str, str]) -> SimpleNamespace:
    out = tmp_path / "out"
    kb = out / "kb-python"
    kb.mkdir(parents=True)
    (kb / "transcripts.txt").write_text("\n".join(transcripts) + "\n")
    t2g_path = tmp_path / "t2g.txt"
    t2g_path.write_text("".join(f"{tx}\t{g}\n" for tx, g in t2g.items()))
    config = SimpleNamespace(
        index="index.idx",
        transcripts=str(t2g_path),
        output=str(out),
        cores=1,
        technology="10xv3",
        strand=None,
        whitelist=None,
        kb_r1="R1.fastq.gz",
        kb_r2="R2.fastq.gz",
    )
    return SimpleNamespace(
        config=config, outputs=SimpleNamespace(transcripts_txt=kb / "transcripts.txt", root=kb)
    )


def test_run_writes_the_split_from_the_replays_own_ec_map(tmp_path, monkeypatch):
    ctx = _ctx(
        tmp_path,
        TRANSCRIPTS,
        {"H1": "gH1", "H2": "gH2", "V1": "gV1", "V2": "gV2"},
    )

    def fake_replay(**kw):
        work = Path(kw["workdir"])
        (work / "lineage_bus").mkdir(parents=True)
        (work / "lineage_bus" / "transcripts.txt").write_text("\n".join(TRANSCRIPTS) + "\n")
        (work / "lineage_bus" / "matrix.ec").write_text(EC_LINES)
        viral, only = work / "v.txt", work / "o.txt"
        viral.write_text("".join(_rows({0: 0, 1: 1, 2: 3})))
        only.write_text("".join(_rows({0: 0, 1: 1})))
        assert sorted(kw["host_transcripts"]) == ["H1", "H2"]
        assert sorted(kw["viral_transcripts"]) == ["V1", "V2"]
        return viral, only

    monkeypatch.setattr(hc, "replay_host_capture_bus", fake_replay)
    tsv, done = tmp_path / "hc.tsv", tmp_path / "hc.done"
    hc.run(ctx, _identity(), str(tsv), str(done))
    lines = [ln.split("\t") for ln in tsv.read_text().splitlines()]
    assert lines[0] == hc.COLUMNS
    assert {ln[0]: ln[1:4] for ln in lines[1:]} == {
        "Virus one": ["1", "1", "0"],
        "Virus two": ["2", "1", "1"],
    }
    assert done.read_text().startswith("measured")


def test_virus_only_index_is_not_applicable_never_zero(tmp_path, monkeypatch):
    ctx = _ctx(tmp_path, ["V1", "V2"], {"V1": "gV1", "V2": "gV2"})

    def no_replay(**kw):  # a replay on an index with no host would call every read virus-only
        raise AssertionError("must not replay without host transcripts")

    monkeypatch.setattr(hc, "replay_host_capture_bus", no_replay)
    tsv, done = tmp_path / "hc.tsv", tmp_path / "hc.done"
    hc.run(ctx, _identity(), str(tsv), str(done))
    rows = [ln.split("\t") for ln in tsv.read_text().splitlines()[1:]]
    assert [r[0] for r in rows] == ["Virus one", "Virus two"]
    assert all(r[1:4] == ["not_applicable"] * 3 for r in rows)
    assert "0" not in {c for r in rows for c in r[1:4]}
    assert done.read_text().startswith("not_applicable")


# ── flag, config, fingerprint, DAG ───────────────────────────────────────────────────────────
def _cfg(tmp_path: Path, **extra) -> dict:
    return {
        "output": str(tmp_path / "out") + "/",
        "index": str(tmp_path / "ref" / "index.idx"),
        "transcripts": "t",
        "sample1": "s1",
        "sample2": "s2",
        "technology": "10xv3",
        "visual": "True",
        "multimapping": "True",
        "gtf": "None",
        "fasta": "None",
        "f1": "None",
        "reference": "False",
        "umap": "False",
        "whitelist": "None",
        "emptydrops_seed": 100,
        **extra,
    }


def test_config_is_typed_and_default_off(tmp_path):
    assert RunConfig.from_snakemake_config(_cfg(tmp_path)).host_capture is False
    rc = RunConfig.from_snakemake_config(_cfg(tmp_path, host_capture="true"))
    assert rc.host_capture is True
    path = _write_run_config(rc)
    assert yaml.safe_load(path.read_text())["host_capture"] is True
    assert RunConfig.from_yaml(path).host_capture is True


def test_flag_parses_and_default_run_keeps_its_fingerprint():
    with patch("sys.argv", ["viralscan"]):
        args = create_help()
    assert args.host_capture is None
    assert "host_capture" not in build_run_manifest(args)["options"]
    with patch("sys.argv", ["viralscan", "--host-capture"]):
        assert create_help().host_capture is True
