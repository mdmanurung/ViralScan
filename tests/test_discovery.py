"""DISC-02/03/06: `viralscan discover` refuses unsafe inputs, never touches the source run, and
resolves the host-unmapped pair fail-closed."""

import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from viralscan import discovery
from viralscan.discovery import DiscoveryError


def _run(tmp_path: Path, samples=("S1",), complete: bool = True) -> Path:
    root = tmp_path / "run"
    for sample in samples:
        (root / sample / "results").mkdir(parents=True)
        (root / sample / "config.yaml").write_text("output: x\n")
        (root / sample / "results" / "viral_summary.tsv").write_text("virus\tdetected\nEBV\tTrue\n")
    (root / "run_manifest.json").write_text("{}\n")
    if complete:
        (root / "run_complete.json").write_text("{}\n")
    return root


def _fastq_pair(r1: Path, r2: Path, n: int) -> None:
    r1.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        for i in range(n):
            f1.write(f"@frag{i} 1:N\n{'A' * 28}\n+\n{'I' * 28}\n")
            f2.write(f"@frag{i} 2:N\n{'C' * 90}\n+\n{'I' * 90}\n")


def _host_filtered(
    hf: Path, n: int, audit_n: int | None = None, lineage_n: int | None = None
) -> None:
    """The host-filter layout: retained pair, audit and per-fragment lineage."""
    _fastq_pair(hf / "R1.fastq.gz", hf / "R2.fastq.gz", n)
    (hf / "host_filter_audit.tsv").write_text(
        "category\tfragments\tinterpretation\n"
        f"input\t{n + 5}\tx\nretained_host_unmapped\t{n if audit_n is None else audit_n}\tx\n"
        "star_param_set\tpinned\tx\nstar_param:outFilterMultimapNmax\t20\tx\n"
    )
    kept = n if lineage_n is None else lineage_n
    with gzip.open(hf / "fragment_lineage.tsv.gz", "wt") as handle:
        handle.write("read_id\tfilter_decision\treason\n")
        for i in range(kept + 5):
            decision = "retained" if i < kept else "removed"
            handle.write(f"frag{i}\t{decision}\thost\n")


def _args(**kw) -> SimpleNamespace:
    base = dict(host_subtract=None, host_index=None, host_filter_star_params="pinned", threads=1)
    return SimpleNamespace(**{**base, **kw})


def _db(tmp_path: Path, empty: str | None = None) -> Path:
    db = tmp_path / "k2db"
    db.mkdir()
    for name in discovery.KRAKEN2_DB_FILES:
        (db / name).write_bytes(b"" if name == empty else b"x")
    return db


def _inventory(root: Path) -> dict[str, str]:
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def test_resolves_single_sample(tmp_path):
    root = _run(tmp_path)
    assert discovery.resolve_sample_dir(root, None) == root / "S1"


@pytest.mark.parametrize(
    "samples,sample,complete,match",
    [
        (("S1",), None, False, "not a completed"),
        (("S1", "S2"), None, True, "several samples"),
        (("S1",), "S9", True, "not in run"),
    ],
)
def test_refuses_bad_run(tmp_path, samples, sample, complete, match):
    root = _run(tmp_path, samples, complete)
    with pytest.raises(DiscoveryError, match=match):
        discovery.resolve_sample_dir(root, sample)


def test_refuses_output_inside_run_or_non_empty(tmp_path):
    root = _run(tmp_path)
    for inside in (root, root / "S1" / "discovery"):
        with pytest.raises(DiscoveryError, match="outside --run-dir"):
            discovery.check_output_dir(inside, root)
    busy = tmp_path / "busy"
    busy.mkdir()
    (busy / "old.tsv").write_text("x")
    with pytest.raises(DiscoveryError, match="empty or absent"):
        discovery.check_output_dir(busy, root)
    # A leftover staging directory alone does not count as results.
    stale = tmp_path / "stale"
    (stale / discovery.PARTIAL).mkdir(parents=True)
    discovery.check_output_dir(stale, root)


def test_refuses_incomplete_db(tmp_path):
    with pytest.raises(DiscoveryError, match="does not exist"):
        discovery.validate_kraken2_db(tmp_path / "absent")
    with pytest.raises(DiscoveryError, match="hash.k2d"):
        discovery.validate_kraken2_db(_db(tmp_path, empty="hash.k2d"))


def test_refuses_missing_tools():
    with patch("viralscan.discovery.shutil.which", return_value=None):
        with pytest.raises(DiscoveryError, match="kraken2, k2"):
            discovery.check_tools()


def test_stage_and_publish(tmp_path):
    out = tmp_path / "out"
    partial = discovery.stage(out)
    (partial / "leftover.txt").write_text("from an interrupted run")
    partial = discovery.stage(out)  # restart discards the leftover
    (partial / "qc.json").write_text("{}")
    discovery.publish(out)
    assert sorted(p.name for p in out.iterdir()) == ["qc.json"]


def test_cli_leaves_source_run_untouched(tmp_path):
    """G2: the source inventory is byte-identical and the skeleton publishes nothing."""
    from viralscan.menu import main

    root = _run(tmp_path)
    _host_filtered(root / "S1" / "host_filtered", 3)
    before = _inventory(root)
    out = tmp_path / "disc"
    argv = [
        "viralscan",
        "discover",
        "--run-dir",
        str(root),
        "-o",
        str(out),
        "--kraken2-db",
        str(_db(tmp_path)),
    ]
    with (
        patch("sys.argv", argv),
        patch("viralscan.discovery.shutil.which", return_value="/bin/true"),
        pytest.raises(SystemExit) as exit_info,
    ):
        main()
    assert exit_info.value.code == 1
    assert _inventory(root) == before
    assert [p.name for p in out.iterdir()] == [discovery.PARTIAL]  # staged, never published
    staged = json.loads((out / discovery.PARTIAL / "input.json").read_text())
    assert (staged["source"], staged["n_fragments"]) == ("primary_host_filtered", 3)
    assert staged["stream"]["n_input"] == 3
    assert (out / discovery.PARTIAL / discovery.KRAKEN_INPUT).is_file()


def test_run_parser_gains_no_discovery_options():
    """G1: discovery adds nothing to the primary run's options (and so to its fingerprint)."""
    from viralscan.menu import build_parser

    with patch("sys.argv", ["viralscan"]):
        args = build_parser().parse_args([])
    assert not {"kraken2_db", "sample", "threads", "host_subtract"} & set(vars(args))


def test_reuses_primary_host_filtered(tmp_path):
    hf = tmp_path / "S1" / "host_filtered"
    _host_filtered(hf, 4)
    resolved = discovery.read_host_filtered(hf)
    assert resolved["n_fragments"] == 4
    assert resolved["star_param_set"] == "pinned"
    assert resolved["star_args"] == {"outFilterMultimapNmax": "20"}
    assert resolved["primary_verified"] is False
    assert resolved["sha256"]["r1"] == hashlib.sha256((hf / "R1.fastq.gz").read_bytes()).hexdigest()


def test_no_host_filter_output_is_none(tmp_path):
    assert discovery.read_host_filtered(tmp_path / "host_filtered") is None


@pytest.mark.parametrize("kw,match", [({"audit_n": 3}, "disagree"), ({"lineage_n": 5}, "disagree")])
def test_refuses_inconsistent_host_filter_audit(tmp_path, kw, match):
    hf = tmp_path / "host_filtered"
    _host_filtered(hf, 4, **kw)
    with pytest.raises(DiscoveryError, match=match):
        discovery.read_host_filtered(hf)


def test_refuses_partial_or_truncated_host_filter_output(tmp_path):
    hf = tmp_path / "host_filtered"
    _host_filtered(hf, 4)
    (hf / "fragment_lineage.tsv.gz").unlink()
    with pytest.raises(DiscoveryError, match="Incomplete"):
        discovery.read_host_filtered(hf)
    _host_filtered(hf, 4)
    _fastq_pair(hf / "R1.fastq.gz", tmp_path / "unused_R2.fastq.gz", 3)  # R1 one record short
    with pytest.raises(DiscoveryError, match="unreadable"):
        discovery.read_host_filtered(hf)


def test_input_resolution_refusals(tmp_path):
    root = _run(tmp_path)
    sample = root / "S1"
    with pytest.raises(DiscoveryError, match="--host-subtract starsolo --host-index"):
        discovery.resolve_input(root, sample, tmp_path / "p", _args())
    with pytest.raises(DiscoveryError, match="requires --host-index"):
        discovery.resolve_input(root, sample, tmp_path / "p", _args(host_subtract="starsolo"))
    _host_filtered(sample / "host_filtered", 2)
    with pytest.raises(DiscoveryError, match="drop --host-subtract"):
        discovery.resolve_input(root, sample, tmp_path / "p", _args(host_subtract="starsolo"))


def _star_run(tmp_path: Path, tamper: bool = False, delete: bool = False):
    """A run without host filtering whose manifest hashes its sample FASTQs."""
    root = _run(tmp_path)
    r1, r2 = tmp_path / "raw" / "R1.fastq.gz", tmp_path / "raw" / "R2.fastq.gz"
    _fastq_pair(r1, r2, 6)
    manifest = {
        "input_paths": {"sample1:0": str(r1.resolve()), "sample2:0": str(r2.resolve())},
        "input_fingerprints": {
            "sample1:0": hashlib.sha256(r1.read_bytes()).hexdigest(),
            "sample2:0": hashlib.sha256(r2.read_bytes()).hexdigest(),
        },
    }
    (root / "run_manifest.json").write_text(json.dumps(manifest))
    if tamper:
        _fastq_pair(r1, r2, 5)
    if delete:
        r1.unlink()
    config = SimpleNamespace(sample1=str(r1), sample2=str(r2), technology="10xv3", whitelist=None)
    host_index = tmp_path / "star_index"
    host_index.mkdir()
    return root, config, host_index


def _fake_starsolo(
    r1, r2, host_index, technology, whitelist, out_dir, fr1, fr2, threads, param_set
):
    _host_filtered(Path(out_dir), 2)


def test_discovery_star_pass_writes_only_under_output(tmp_path):
    root, config, host_index = _star_run(tmp_path)
    before = _inventory(root)
    args = _args(host_subtract="starsolo", host_index=str(host_index))
    with (
        patch("viralscan.runconfig.RunConfig.from_yaml", return_value=config),
        patch("viralscan.scripts.host_filter.check_host_filter_tools"),
        patch("viralscan.scripts.host_filter._starsolo_filter", side_effect=_fake_starsolo) as star,
    ):
        resolved = discovery.resolve_input(root, root / "S1", tmp_path / "out" / ".partial", args)
    assert star.call_args.args[:5] == (
        config.sample1,
        config.sample2,
        str(host_index),
        "10xv3",
        None,
    )
    assert star.call_args.args[-1] == "pinned"
    assert resolved["source"] == "discovery_star_pass"
    assert resolved["n_fragments"] == 2
    assert Path(resolved["r2"]).parent == tmp_path / "out" / ".partial" / "host_subtraction"
    assert _inventory(root) == before


@pytest.mark.parametrize("kw,match", [({"tamper": True}, "changed"), ({"delete": True}, "absent")])
def test_discovery_star_pass_checks_primary_fastqs(tmp_path, kw, match):
    root, config, host_index = _star_run(tmp_path, **kw)
    args = _args(host_subtract="starsolo", host_index=str(host_index))
    with (
        patch("viralscan.runconfig.RunConfig.from_yaml", return_value=config),
        patch("viralscan.scripts.host_filter._starsolo_filter") as star,
        pytest.raises(DiscoveryError, match=match),
    ):
        discovery.resolve_input(root, root / "S1", tmp_path / "out" / ".partial", args)
    star.assert_not_called()


def _pair(tmp_path: Path, records, r1_name="R1.fastq.gz", r2_name="R2.fastq.gz"):
    """Write (read_id, r1_seq, r2_seq) records; returns the two paths."""
    r1, r2 = tmp_path / r1_name, tmp_path / r2_name
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        for rid, s1, s2 in records:
            f1.write(f"@{rid} 1:N\n{s1}\n+\n{'I' * len(s1)}\n")
            f2.write(f"@{rid} 2:N\n{s2}\n+\n{'I' * len(s2)}\n")
    return str(r1), str(r2)


def _rows(path: Path) -> list[dict]:
    import csv

    with gzip.open(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_stream_reads_writes_ordinal_fasta_and_identity(tmp_path):
    cb, umi = "AAAACCCCGGGGTTTT", "ACGTACGTACGT"
    r1, r2 = _pair(
        tmp_path,
        [
            ("fragA/1", cb + umi + "TTTTT", "CCCCCCCCCC"),
            ("fragB", "AAAACCCCGGGGTTTN" + "ACGTACGTACGN", "GGGGGGGGGG"),
            ("fragC", "SHORT", "TTTTTTTTTT"),
        ],
    )
    counts = discovery.stream_reads(r1, r2, 16, 12, tmp_path)
    assert counts == {
        "n_input": 3,
        "n_bad_geometry": 1,
        "n_cb_has_n": 1,
        "n_umi_has_n": 1,
        "read_filter:not_run": 3,
    }
    with gzip.open(tmp_path / discovery.KRAKEN_INPUT, "rt") as handle:
        assert handle.read() == ">0\nCCCCCCCCCC\n>1\nGGGGGGGGGG\n>2\nTTTTTTTTTT\n"
    rows = _rows(tmp_path / discovery.READS_TSV)
    assert [r["read_id"] for r in rows] == ["fragA", "fragB", "fragC"]
    assert (rows[0]["raw_cb"], rows[0]["raw_umi"], rows[0]["geometry_ok"]) == (cb, umi, "1")
    assert (rows[1]["cb_has_n"], rows[1]["umi_has_n"]) == ("1", "1")
    assert (rows[2]["raw_cb"], rows[2]["geometry_ok"]) == ("", "0")


def test_stream_reads_empty_pair(tmp_path):
    r1, r2 = _pair(tmp_path, [])
    assert discovery.stream_reads(r1, r2, 16, 12, tmp_path)["n_input"] == 0


@pytest.mark.parametrize(
    "records2,match",
    [
        ([("x0", "C"), ("x9", "C")], "mate mismatch"),
        ([("x0", "C")], "Truncated or unpaired"),
    ],
)
def test_stream_reads_refuses_broken_pairs(tmp_path, records2, match):
    r1, _ = _pair(tmp_path, [("x0", "A" * 28, "C"), ("x1", "A" * 28, "C")])
    _, r2 = _pair(tmp_path, [(rid, "A" * 28, s) for rid, s in records2], r1_name="unused.fq.gz")
    with pytest.raises(DiscoveryError, match=match):
        discovery.stream_reads(r1, r2, 16, 12, tmp_path)


def test_stream_reads_refuses_duplicate_read_ids(tmp_path):
    r1, r2 = _pair(
        tmp_path, [("dup", "A" * 28, "C"), ("ok", "A" * 28, "C"), ("dup", "A" * 28, "G")]
    )
    with pytest.raises(DiscoveryError, match="Duplicate read IDs"):
        discovery.stream_reads(r1, r2, 16, 12, tmp_path)


def _lineage(path: Path, rows) -> Path:
    with gzip.open(path, "wt") as handle:
        handle.write("read_id\tfilter_decision\treason\n")
        for rid, decision, reason in rows:
            handle.write(f"{rid}\t{decision}\t{reason}\n")
    return path


def test_stream_reads_joins_read_filter_lineage_in_lockstep(tmp_path):
    r1, r2 = _pair(tmp_path, [("a", "A" * 28, "C"), ("b", "A" * 28, "C")])
    lineage = _lineage(
        tmp_path / "lin.tsv.gz", [("a", "retained", "complex_body"), ("b", "removed", "r2_reagent")]
    )
    counts = discovery.stream_reads(r1, r2, 16, 12, tmp_path, lineage)
    assert [r["primary_read_filter"] for r in _rows(tmp_path / discovery.READS_TSV)] == [
        "retained",
        "removed:r2_reagent",
    ]
    assert counts["read_filter:retained"] == counts["read_filter:removed"] == 1


@pytest.mark.parametrize(
    "lineage_rows,match",
    [
        ([("b", "retained", "x"), ("a", "retained", "x")], "out of step"),
        ([("a", "retained", "x")], "out of step"),
        ([("a", "retained", "x"), ("b", "retained", "x"), ("c", "retained", "x")], "more rows"),
    ],
)
def test_stream_reads_refuses_misaligned_lineage(tmp_path, lineage_rows, match):
    r1, r2 = _pair(tmp_path, [("a", "A" * 28, "C"), ("b", "A" * 28, "C")])
    lineage = _lineage(tmp_path / "lin.tsv.gz", lineage_rows)
    with pytest.raises(DiscoveryError, match=match):
        discovery.stream_reads(r1, r2, 16, 12, tmp_path, lineage)


@pytest.mark.parametrize(
    "technology,expected",
    [
        ("10xv3", (16, 12)),
        ("10xv2", (16, 10)),
        ("dropseq", (12, 8)),
        ("0,0,16:0,16,28:1,0,0", (16, 12)),
    ],
)
def test_barcode_geometry_supported(technology, expected):
    assert discovery.barcode_geometry(technology) == expected


@pytest.mark.parametrize("technology", ["10xv1", "0,0,16:1,0,12:1,0,0", "nonsense"])
def test_barcode_geometry_refused(technology):
    with pytest.raises(DiscoveryError):
        discovery.barcode_geometry(technology)


def test_read_filter_lineage_choice(tmp_path):
    star = {"source": "discovery_star_pass"}
    primary = {"source": "primary_host_filtered"}
    assert discovery.read_filter_lineage_for(tmp_path, primary, "off") == (None, "not_run")
    assert discovery.read_filter_lineage_for(tmp_path, star, "artefact") == (None, "not_joined")
    with pytest.raises(DiscoveryError, match="read-filter artefact"):
        discovery.read_filter_lineage_for(tmp_path, primary, "artefact")
    path = _lineage(tmp_path / "unused.tsv.gz", [])
    (tmp_path / "read_filtered").mkdir()
    path.rename(tmp_path / "read_filtered" / "fragment_lineage.tsv.gz")
    assert discovery.read_filter_lineage_for(tmp_path, primary, "artefact")[0] is not None
