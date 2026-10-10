"""Hand-computable competition contracts; no clinical data or remote resources."""

import csv
import gzip
import hashlib
import random
from pathlib import Path

import pytest

from viralscan import evidence as ev
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable


def inputs(tmp_path):
    seqs = {"target": ("t", "ACGTTGCA"), "host": ("h", "CCCCGGGG"), "competitor": ("r", "GATTACAG")}
    rows = []
    for role, (name, seq) in seqs.items():
        (tmp_path / f"{role}.fa").write_text(f">{name}\n{seq}\n")
        rows.append(
            dict(
                reference_id=name,
                source_role=role,
                **{
                    "class": {"target": "target", "host": "host", "competitor": "related_virus"}[
                        role
                    ]
                },
                virus_key={"target": "taxid:1", "host": "", "competitor": "taxid:2"}[role],
                reporting_group={"target": "taxid:1", "host": "", "competitor": "taxid:2"}[role],
                accession=name,
                display_name=name,
                sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                source="synthetic",
                relationship=role,
            )
        )
    path = tmp_path / "manifest.tsv"
    write_manifest(path, rows)
    identity = VirusIdentityTable(
        (GeneIdentity("g", "t", "catalogued", True, "taxid:1", "Target", sibling_group="siblings"),)
    )
    kw = dict(
        target_fasta=str(tmp_path / "target.fa"),
        host_fasta=str(tmp_path / "host.fa"),
        competitor_fasta=str(tmp_path / "competitor.fa"),
        identity=identity,
        target_genes=["g"],
    )
    return path, rows, kw


def write_manifest(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_manifest_and_encoded_fasta(tmp_path):
    path, _, kw = inputs(tmp_path)
    rows = ev.validate_competitor_manifest(path, **kw)
    out = tmp_path / "merged.fa"
    ev.write_competitor_fasta(
        kw["host_fasta"], kw["target_fasta"], kw["competitor_fasta"], rows, str(out)
    )
    assert set(line for line in out.read_text().splitlines() if line.startswith(">")) == {
        ">TARGET|t",
        ">RELATED|r",
        ">HOST|h",
    }


@pytest.mark.parametrize(
    "change",
    [
        "duplicate",
        "missing",
        "extra",
        "hash",
        "class",
        "role",
        "target_key",
        "same_group",
        "identity",
        "sequence_class",
    ],
)
def test_bad_manifest_rejected(tmp_path, change):
    path, rows, kw = inputs(tmp_path)
    if change == "duplicate":
        rows.append(dict(rows[0]))
    elif change == "missing":
        rows.pop()
    elif change == "extra":
        rows.append({**rows[2], "reference_id": "extra"})
    elif change == "hash":
        rows[0]["sequence_sha256"] = "0" * 64
    elif change == "class":
        rows[2]["class"] = "virus"
    elif change == "role":
        rows[2]["source_role"] = "host"
    elif change == "target_key":
        rows[0]["virus_key"] = "taxid:other"
    elif change == "same_group":
        rows[2]["class"] = "same_reporting_group"
    elif change == "identity":
        kw["identity"] = None
    elif change == "sequence_class":
        (tmp_path / "competitor.fa").write_text(">r\nACGTTGCA\n")
        rows[2]["sequence_sha256"] = rows[0]["sequence_sha256"]
    write_manifest(path, rows)
    with pytest.raises(ValueError):
        ev.validate_competitor_manifest(path, **kw)


def test_gzip_case_wrapping_normalization(tmp_path):
    path, rows, kw = inputs(tmp_path)
    src = tmp_path / "target.fa.gz"
    with gzip.open(src, "wt") as handle:
        handle.write(">t desc\nac gt\ntgca\n")
    kw["target_fasta"] = str(src)
    assert (
        ev.validate_competitor_manifest(path, **kw)["t"]["sequence_sha256"]
        == rows[0]["sequence_sha256"]
    )


def hit(query, subject, score=100, identity=99, coverage=100):
    return f"{query}\t{subject}\t{identity}\t50\t{coverage}\t1e-20\t{score}"


@pytest.mark.parametrize(
    "classes,expected",
    [
        (["TARGET"], "target_specific"),
        (["SAMEGROUP"], "within_group_support"),
        (["RELATED"], "related_preferred"),
        (["HOST"], "host_preferred"),
        (["DECOY"], "decoy_preferred"),
        (["TARGET", "SAMEGROUP"], "within_group_ambiguous"),
        (["TARGET", "RELATED"], "related_virus_ambiguous"),
        (["TARGET", "HOST"], "host_competitive"),
        (["TARGET", "DECOY"], "multi_class_ambiguous"),
        (["HOST", "RELATED"], "multi_class_ambiguous"),
    ],
)
def test_exhaustive_winner_classes(classes, expected):
    rows = ev.parse_competitor_blast_output(
        "\n".join(hit("q", f"{cls}|{i}") for i, cls in enumerate(classes)), ["q"]
    )
    assert rows[0]["diagnostic_class"] == expected
    assert rows[0]["n_top_tied_hits"] == len(classes)


def test_ties_distinct_subjects_order_independent_and_no_hits():
    lines = [
        hit("q", "TARGET|a", 101),
        hit("q", "TARGET|b", 101),
        hit("q", "TARGET|a", 101, identity=95),
        hit("q", "RELATED|r", 100),
        hit("q", "HOST|h", 90),
    ]
    a = ev.parse_competitor_blast_output("\n".join(lines), ["q", "absent"], tie_delta=1)
    random.Random(42).shuffle(lines)
    assert ev.parse_competitor_blast_output("\n".join(lines), ["q", "absent"], tie_delta=1) == a
    # Two target subjects still constitute one class; the only other class is related.
    assert a[0]["diagnostic_class"] == "related_virus_ambiguous"
    assert a[0]["n_top_tied_hits"] == 3
    assert a[0]["target_winner_hits"] == "TARGET|a;TARGET|b"
    assert a[0]["target_minus_host_bitscore"] == 11
    assert a[1]["diagnostic_class"] == "no_hits" and a[1]["target_bitscore"] is None
    assert a[1]["target_minus_host_bitscore"] is None
    assert a[1]["no_target_support"] is True


def test_incomplete_search_precedes_apparent_specificity():
    row = ev.parse_competitor_blast_output(hit("q", "TARGET|a"), ["q"], search_complete=False)[0]
    assert row["diagnostic_class"] == "incomplete_search"


@pytest.mark.parametrize("delta", [-1, float("inf"), float("nan")])
def test_invalid_tie_delta(delta):
    with pytest.raises(ValueError):
        ev.parse_competitor_blast_output("", [], tie_delta=delta)


def test_sampled_read_and_molecule_denominators_are_separate():
    lineage = [
        dict(read_number=i, read_id=f"id{i}", cb="cell", ub=umi)
        for i, umi in enumerate(["u1", "u1", "u2", "u3"])
    ]
    ids = ["cell_u1_0|id0", "cell_u1_1|id1", "cell_u2_2|id2"]
    rows = ev.parse_competitor_blast_output(
        hit(ids[0], "TARGET|a") + "\n" + hit(ids[1], "HOST|h"), ids
    )
    summary = ev.competitor_summary(rows, lineage, target_id="t", manifest_sha256="a" * 64)
    assert summary["n_reads_sampled_blast"] == 3
    assert summary["n_molecules"] == 3 and summary["n_sampled_molecules"] == 2
    assert summary["n_unsampled_reads"] == 1 and summary["n_unsampled_molecules"] == 1
    assert summary["fraction_target_specific"] == pytest.approx(1 / 3)
    assert summary["n_no_hit_reads"] == 1
    assert summary["n_unresolved_sampled_reads"] == 0


def test_qc_keeps_all_explicit_classes():
    refs = ["TARGET|t", "SAMEGROUP|s", "RELATED|r", "HOST|h", "DECOY|d"]
    sam = "\n".join(
        f"c_u_{i}|r\t0\t{ref}\t1\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0" for i, ref in enumerate(refs)
    )
    header = "\n".join(f"@SQ\tSN:{ref}\tLN:8" for ref in refs)
    expected = {"target", "same_reporting_group", "related_virus", "host", "decoy"}
    assert {r["reference_class"] for r in ev._alignment_qc_from_text(header, sam, "")} == expected
    assert {r["reference_class"] for r in ev._per_cell_qc_from_text(sam)} == expected


def test_absent_blast_qc_are_not_assessed():
    rows = ev.interpretation_flags([], [], [])
    assert {
        row["status"]
        for row in rows
        if row["flag"] in {"host_homology", "low_complexity", "coverage_hotspot"}
    } == {"not_assessed"}


def test_molecule_candidates_separate_secondary_evidence_and_unresolved_lineage():
    q = "cell_u1_0|id0"
    records = [
        f"{q}\t0\tTARGET|t\t1\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0",
        f"{q}\t256\tRELATED|r\t1\t20\t4M\t*\t0\t0\tACGT\t*\tNM:i:1",
    ]
    lineage = [dict(cb="cell", ub="u1"), dict(cb="cell", ub="u2")]
    a = ev.molecule_competition_rows("\n".join(records), lineage)
    b = ev.molecule_competition_rows("\n".join(reversed(records)), lineage)
    assert a == b and len(a) == 2
    assert a[0]["selected_class"] == "target" and a[0]["n_candidate_alignments"] == 2
    assert a[0]["candidate_classes"] == "related_virus;target"
    assert a[1]["status"] == "unresolved" and a[1]["selected_reference_id"] is None


def test_live_blast_failure_keeps_every_sampled_query_and_raw_failure(tmp_path, monkeypatch):
    fasta = tmp_path / "reads.fa"
    fasta.write_text(">q1\nACGT\n>q2\nGCTA\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">TARGET|t\nACGT\n")

    def fail(*args, **kwargs):
        raise RuntimeError("toy blast error")

    monkeypatch.setattr(ev, "_run", fail)
    rows = ev.competitive_blast_identity(
        str(fasta), str(ref), str(tmp_path / "blast"), multi_class=True
    )
    assert len(rows) == 2 and {r["diagnostic_class"] for r in rows} == {"failed"}
    assert all(r["target_bitscore"] is None for r in rows)
    assert "toy blast error" in (tmp_path / "blast/blast_failure.txt").read_text()


def test_complete_subject_cap_and_all_queries(tmp_path, monkeypatch):
    fasta = tmp_path / "reads.fa"
    fasta.write_text(">q\nACGT\n>nohit\nNNNN\n")
    ref = tmp_path / "ref.fa"
    ref.write_text("".join(f">TARGET|t{i}\nACGT\n" for i in range(25)))
    commands = []

    def fake_run(cmd, **kwargs):
        commands.append(cmd)
        if cmd[0] == "blastn":
            return ("\n".join(hit("q", f"TARGET|t{i}") for i in range(25)) + "\n").encode()
        return b""

    monkeypatch.setattr(ev, "_run", fake_run)
    rows = ev.competitive_blast_identity(
        str(fasta), str(ref), str(tmp_path / "blast"), multi_class=True
    )
    cmd = next(c for c in commands if c[0] == "blastn")
    assert cmd[cmd.index("-max_target_seqs") + 1] == "25"
    assert next(r for r in rows if r["read"] == "q")["n_top_tied_hits"] == 25
    assert next(r for r in rows if r["read"] == "nohit")["diagnostic_class"] == "no_hits"


def stub_run(tmp_path, monkeypatch, *, zero=True):
    """Stub replay only: references/identity/schema/output preflight remain real."""
    from types import SimpleNamespace

    from viralscan.scripts import evidence_run as runner

    path, rows, kw = inputs(tmp_path)
    run = tmp_path / "run"
    kb = run / "kb-python"
    kb.mkdir(parents=True)
    (run / "config.yaml").write_text("technology: 10xv3\n")
    for file in ["genes.txt", "transcripts.txt", "output.resolved.sorted.bus.txt", "matrix.ec"]:
        (kb / file).write_text("g\n" if file == "genes.txt" else "")
    (kb / "counts_unfiltered").mkdir()
    (kb / "counts_unfiltered/cells_x_genes.genes.txt").write_text("g\n")
    for file in ["R1.fq", "R2.fq", "index.idx"]:
        (tmp_path / file).write_text("")
    cfg = SimpleNamespace(
        technology="10xv3",
        index=str(tmp_path / "index.idx"),
        transcripts=str(kb / "transcripts.txt"),
        sample1=str(tmp_path / "R1.fq"),
        sample2=str(tmp_path / "R2.fq"),
        kb_r1=None,
        kb_r2=None,
        whitelist=None,
        strand=None,
        multimap_method="equal",
    )
    monkeypatch.setattr(runner.RunConfig, "from_yaml", lambda p: cfg)
    monkeypatch.setattr(runner, "load_transcripts", lambda *a: (["tx"], {"tx": "g"}))
    monkeypatch.setattr(runner, "load_run_identity", lambda *a: kw["identity"])
    monkeypatch.setattr(runner, "_replay_ec_map", lambda *a: {})
    replay = []

    def fake_replay(**args):
        replay.append(args)
        flagged = tmp_path / "flagged.txt"
        flagged.write_text("")
        return flagged

    monkeypatch.setattr(runner, "replay_exact_target_bus", fake_replay)
    monkeypatch.setattr(runner, "parse_flagged_target_bus", lambda *a: {})

    def fake_extract(*args):
        Path(args[3]).write_text("")
        with gzip.open(args[4], "wt") as handle:
            handle.write("read_number\tread_id\tcb\tub\n")
        return SimpleNamespace(viral_reads=0, total_reads=0)

    monkeypatch.setattr(runner, "extract_exact_reads_by_number", fake_extract)
    args = SimpleNamespace(
        run_dir=str(run),
        output=str(tmp_path / "out"),
        virus="g",
        cores=1,
        verbose=False,
        quiet=True,
        viral_fasta=kw["target_fasta"],
        host_fasta=kw["host_fasta"],
        competitor_fasta=kw["competitor_fasta"],
        competitor_manifest=str(path),
        blast=True,
        blast_tie_delta=0,
        sampling_seed=42,
    )
    return runner, args, replay, path, rows


def test_zero_read_path_validates_new_inputs_before_replay(tmp_path, monkeypatch):
    runner, args, replay, path, rows = stub_run(tmp_path, monkeypatch)
    rows[0]["sequence_sha256"] = "0" * 64
    write_manifest(path, rows)
    with pytest.raises(SystemExit):
        runner.run_evidence(args)
    assert replay == [] and not Path(args.output).exists()


def test_zero_read_competition_outputs_have_headers_status_and_schema(tmp_path, monkeypatch):
    import json

    from viralscan.validation import require_schema_valid

    runner, args, replay, _, _ = stub_run(tmp_path, monkeypatch)
    runner.run_evidence(args)
    out = Path(args.output)
    assert len(replay) == 1
    summary = list(csv.DictReader((out / "competitor_summary.tsv").open(), delimiter="\t"))[0]
    assert summary["n_reads_sampled_blast"] == "0" and summary["blast_status"] == "not_assessed"
    assert summary["fraction_target_specific"] == ""
    manifest = json.loads((out / "evidence_manifest.json").read_text())
    assert manifest["competition"]["extraction_status"] == "empty"
    assert (
        manifest["competition"]["manifest_sha256"]
        == hashlib.sha256(Path(args.competitor_manifest).read_bytes()).hexdigest()
    )
    assert set(manifest["tool_binaries"]) == {"kallisto", "bustools"}
    require_schema_valid(manifest, "evidence_manifest.schema.json")
    assert "selected_class" in (out / "molecule_competition.tsv").read_text()
    assert "target_bitscore" in (out / "blast_identity.tsv").read_text()


def test_new_mode_refuses_stale_output_before_replay(tmp_path, monkeypatch):
    runner, args, replay, _, _ = stub_run(tmp_path, monkeypatch)
    out = Path(args.output)
    out.mkdir()
    (out / "previous.tsv").write_text("keep")
    with pytest.raises(SystemExit):
        runner.run_evidence(args)
    assert replay == [] and (out / "previous.tsv").read_text() == "keep"


def test_zero_read_manifest_records_no_molecule_verdict(tmp_path, monkeypatch):
    import json

    runner, args, _, _, _ = stub_run(tmp_path, monkeypatch)
    runner.run_evidence(args)
    manifest = json.loads((Path(args.output) / "evidence_manifest.json").read_text())
    assert manifest["competition"]["molecule_verdict"] == "not_produced_in_competitor_mode"


def _parse_evidence(*extra):
    from unittest.mock import patch

    from viralscan.menu import create_help

    argv = ["viralscan", "evidence", "--run-dir", "r/", "-o", "o/", "--virus", "EBV", *extra]
    with patch("sys.argv", argv):
        return create_help()


def test_cli_flags_reach_args_with_defaults():
    args = _parse_evidence()
    assert args.competitor_fasta is None and args.competitor_manifest is None
    assert args.blast_tie_delta == 0.0
    args = _parse_evidence(
        "--competitor-fasta", "c.fa", "--competitor-manifest", "c.tsv", "--blast-tie-delta", "2.5"
    )
    assert (args.competitor_fasta, args.competitor_manifest) == ("c.fa", "c.tsv")
    assert args.blast_tie_delta == 2.5


@pytest.mark.parametrize(
    "extra",
    [
        ["--competitor-fasta", "c.fa"],
        ["--competitor-manifest", "c.tsv"],
        ["--blast-tie-delta", "-1"],
        ["--blast-tie-delta", "nan"],
        ["--competitor-fasta", "c.fa", "--competitor-manifest", "c.tsv"],  # no target/host FASTA
    ],
)
def test_cli_rejects_bad_competitor_flags_before_touching_run(extra, tmp_path):
    from viralscan.scripts.evidence_run import run_evidence

    args = _parse_evidence(*extra)
    args.run_dir, args.output = str(tmp_path / "missing_run"), str(tmp_path / "out")
    with pytest.raises(SystemExit):
        run_evidence(args)
    assert not (tmp_path / "out").exists()
