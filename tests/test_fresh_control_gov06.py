"""GOV-06 attempt-3: highmem tier, explicit-GTF v3 rows, parity, evidence, comparison."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from scripts.compare_legacy_v2_v3 import compare_fresh_vs_archive
from scripts.freeze_fresh_control_packet import freeze_packet
from scripts.prepare_fresh_controls import (
    EXPECTED_IDS,
    FreshControlPreparationError,
    build_panel_gtf,
    gtf_t2g_parity,
    prepare_tasks,
)
from scripts.run_fresh_control import build_command
from scripts.run_fresh_control_evidence import largest_non_target, sample_root

GTF_A = 'chr\tx\tgene\t1\t9\t.\t+\t.\tgene_id "VIR_A1"; transcript_id "t";\n'
GTF_B = 'chr\tx\tgene\t1\t9\t.\t+\t.\tgene_id "VIR_B1"; transcript_id "t";'  # no trailing \n


def _inputs(tmp_path: Path) -> dict[str, Path]:
    raw = tmp_path / "raw.tsv"
    fields = [
        "sample_id",
        "read1_path",
        "read2_path",
        "read1_storage_bytes",
        "read2_storage_bytes",
        "read1_storage_sha256",
        "read2_storage_sha256",
    ]
    with raw.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for sid in sorted(EXPECTED_IDS):
            for mate in (1, 2):
                (tmp_path / f"{sid}_{mate}.fq").write_bytes(b"x" * 10)
            writer.writerow(
                {
                    "sample_id": sid,
                    "read1_path": tmp_path / f"{sid}_1.fq",
                    "read2_path": tmp_path / f"{sid}_2.fq",
                    "read1_storage_bytes": 10,
                    "read2_storage_bytes": 10,
                    "read1_storage_sha256": "1" * 64,
                    "read2_storage_sha256": "2" * 64,
                }
            )
    paths = {"raw_manifest": raw}
    for name in ("v2_viralscan", "v3_viralscan", "index", "t2g", "whitelist"):
        paths[name] = tmp_path / name
        paths[name].write_text(name)
    return paths


def _panel(tmp_path: Path) -> tuple[Path, str]:
    data = tmp_path / "v2data"
    data.mkdir()
    (data / "b.gtf").write_text(GTF_B)
    (data / "a.gtf").write_text(GTF_A)
    gtf = tmp_path / "ref" / "panel.gtf"
    return gtf, build_panel_gtf(data, gtf)


def test_panel_gtf_is_sorted_deterministic_and_newline_terminated(tmp_path: Path) -> None:
    gtf, sha = _panel(tmp_path)
    assert gtf.read_text() == GTF_A + GTF_B + "\n"
    assert sha == hashlib.sha256(gtf.read_bytes()).hexdigest()


def test_parity_passes_on_subset_and_fails_on_gtf_only_gene(tmp_path: Path) -> None:
    gtf, _ = _panel(tmp_path)
    t2g = tmp_path / "t2g.txt"
    t2g.write_text("t1\tVIR_A1\nt2\tVIR_B1\nt3\tENSG1\nt4\tVIR_C1\n")
    ok = gtf_t2g_parity(gtf, t2g)
    assert ok["passed"] and ok["gtf_genes"] == 2 and ok["intersection"] == 2
    assert ok["t2g_genes"] == 4 and ok["t2g_non_ensembl_not_in_gtf"] == 1
    t2g.write_text("t1\tVIR_A1\n")
    bad = gtf_t2g_parity(gtf, t2g)
    assert not bad["passed"] and bad["gtf_only"] == 1 and bad["gtf_only_genes"] == ["VIR_B1"]


def test_highmem_tier_and_cache_free_explicit_gtf_rows(tmp_path: Path) -> None:
    paths = _inputs(tmp_path)
    gtf, sha = _panel(tmp_path)
    packet = tmp_path / "packet"
    tasks = prepare_tasks(
        **paths,
        output_root=tmp_path / "out",
        task_manifest=packet / "tasks.tsv",
        cores=8,
        attempt_id="attempt3",
        gtf=gtf,
        gtf_sha256=sha,
        highmem_task_ids=("v2__SRR6825024",),
    )
    by_id = {t["task_id"]: t for t in tasks}
    assert by_id["v2__SRR6825024"]["tier"] == "highmem"
    assert by_id["v3__SRR6825024"]["tier"] == "small"
    assert by_id["v3__SRR6825024"]["gtf_sha256"] == sha
    assert by_id["v2__SRR6825024"]["gtf_path"] == "-"
    assert by_id["v3__SRR6825024"]["viralscan_cache_path"] == "-"
    highmem = (packet / "tasks.highmem.tsv").read_text().splitlines()
    assert len(highmem) == 2 and highmem[1].startswith("v2__SRR6825024\t")


def test_v3_rows_without_cache_or_gtf_are_refused(tmp_path: Path) -> None:
    with pytest.raises(FreshControlPreparationError, match="cache or an explicit GTF"):
        prepare_tasks(
            **_inputs(tmp_path),
            output_root=tmp_path / "out",
            task_manifest=tmp_path / "p/tasks.tsv",
            cores=8,
            attempt_id="a",
        )


def test_v2_only_packet_needs_no_cache_and_unknown_highmem_is_refused(tmp_path: Path) -> None:
    paths = _inputs(tmp_path)
    tasks = prepare_tasks(
        **paths,
        output_root=tmp_path / "out",
        task_manifest=tmp_path / "p/tasks.tsv",
        cores=8,
        attempt_id="a",
        stacks=("v2",),
        sample_ids=("SRR6825024",),
        highmem_task_ids=("v2__SRR6825024",),
    )
    assert [t["tier"] for t in tasks] == ["highmem"]
    with pytest.raises(FreshControlPreparationError, match="unknown highmem"):
        prepare_tasks(
            **paths,
            output_root=tmp_path / "out2",
            task_manifest=tmp_path / "p2/tasks.tsv",
            cores=8,
            attempt_id="a",
            stacks=("v2",),
            highmem_task_ids=("v2__nope",),
        )


def test_build_command_passes_gtf_only_to_v3(tmp_path: Path) -> None:
    kw = {
        "viralscan": Path("vs"),
        "output": Path("o"),
        "read1": Path("r1"),
        "read2": Path("r2"),
        "index": Path("i"),
        "t2g": Path("t"),
        "whitelist": Path("w"),
        "cores": 2,
        "gtf": Path("p.gtf"),
    }
    assert build_command(stack="v3", **kw)[-2:] == ["-gtf", "p.gtf"]
    assert "-gtf" not in build_command(stack="v2", **kw)


def test_freeze_hashes_optional_highmem_gtf_and_evidence(tmp_path: Path) -> None:
    packet = tmp_path / "packet"
    (packet / "reference").mkdir(parents=True)
    for name in (
        "run_task.sh",
        "tasks.tsv",
        "tasks.small.tsv",
        "tasks.large.tsv",
        "tasks.highmem.tsv",
        "control_inputs.raw.tsv",
        "reference/v2_panel.gtf",
        "reference/gtf_t2g_parity.json",
    ):
        (packet / name).write_text(name)
    runner = tmp_path / "r.py"
    runner.write_text("r")
    evidence = tmp_path / "e.py"
    evidence.write_text("e")
    manifest = freeze_packet(packet_root=packet, runner_source=runner, evidence_source=evidence)
    names = {line.split("  ", 1)[1] for line in manifest.read_text().splitlines()}
    assert {
        "tasks.highmem.tsv",
        "reference/v2_panel.gtf",
        "reference/gtf_t2g_parity.json",
        "source/run_fresh_control_evidence.py",
    } <= names


def test_largest_non_target_skips_target_and_zero_rows(tmp_path: Path) -> None:
    summary = tmp_path / "s.tsv"
    summary.write_text(
        "virus_name\tviral_molecules_total_est\nEBV\t900\nHHV6B\t5\nTTV\t7\nZero\t0\n"
    )
    assert largest_non_target(summary, "ebv") == "TTV"
    summary.write_text("virus_name\tviral_molecules_total_est\nEBV\t900\n")
    assert largest_non_target(summary, "EBV") is None


def _tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.write_text(
        "\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows)
    )


def test_compare_fresh_vs_archive_matches_mismatches_and_retains_failures(tmp_path: Path) -> None:
    archive = tmp_path / "archive"
    archive.mkdir()
    rid = "SRR12682296__SRR12682296"
    _tsv(
        archive / "legacy_reproduction.tsv",
        ["run_id", "record_type", "identifier", "summary_value"],
        [
            [rid, "gene", "G1", 5.0],
            [rid, "virus", "EBV", 5.0],
        ],
    )
    _tsv(archive / "run_metrics.tsv", ["run_id", "legacy_total_viral_load"], [[rid, 5.0]])
    _tsv(
        archive / "virus_metrics.tsv",
        ["run_id", "virus_name", "v3_host_conservative"],
        [[rid, "EBV", 10.0]],
    )
    packet = tmp_path / "packet"
    (packet / "status").mkdir(parents=True)
    (packet / "revalidation").mkdir()
    v2root = tmp_path / "v2out"
    v2root.mkdir()
    (v2root / "summary.txt").write_text(
        "Gene ID; Gene Count\nG1;5.0\n\nEBV has a viral load of: 6.0 UMIs.\nTotal amount of viral load found: 5.0\n"
    )
    (packet / "status" / "v2__SRR12682296.json").write_text(
        json.dumps({"status": "failed", "exit_code": 65})
    )
    (packet / "revalidation" / "v2__SRR12682296.json").write_text(
        json.dumps({"status": "success", "v2_output_root": str(v2root)})
    )
    v3out = tmp_path / "v3out"
    (v3out / "SRR12682296" / "results").mkdir(parents=True)  # CLI nests by sample
    (v3out / "SRR12682296" / "config.yaml").write_text("x: 1\n")
    _tsv(
        v3out / "SRR12682296" / "results" / "viral_summary.tsv",
        ["virus_name", "viral_molecules_total_est"],
        [["EBV", 10.0], ["TTV", 3.0]],
    )
    (packet / "status" / "v3__SRR12682296.json").write_text(
        json.dumps({"status": "success", "output": str(v3out)})
    )
    report = compare_fresh_vs_archive([packet], archive)
    assert report["rows_compared"] == 2
    assert report["rows_not_compared"] == 8  # every missing row is kept, never dropped
    verdicts = {(r["stack"], r["identifier"]): r["verdict"] for r in report["records"]}
    assert verdicts[("v2", "G1")] == "match"
    assert verdicts[("v2", "EBV")] == "mismatch"
    assert verdicts[("v3", "EBV")] == "match"
    assert verdicts[("v3", "TTV")] == "mismatch"  # absent archived side counts as 0
    sources = {r["task_id"]: r["record_source"] for r in report["run_status"]}
    assert sources["v2__SRR12682296"] == "packet/revalidation"


def test_sample_root_resolves_nested_cli_layout(tmp_path: Path) -> None:
    (tmp_path / "S" ).mkdir()
    (tmp_path / "S" / "config.yaml").write_text("x: 1\n")
    assert sample_root(tmp_path) == tmp_path / "S"


def test_later_packet_overrides_and_malformed_row_does_not_sink_report(tmp_path: Path) -> None:
    archive = tmp_path / "archive"
    archive.mkdir()
    _tsv(archive / "legacy_reproduction.tsv", ["run_id", "record_type", "identifier", "summary_value"], [])
    _tsv(archive / "run_metrics.tsv", ["run_id", "legacy_total_viral_load"], [])
    _tsv(archive / "virus_metrics.tsv", ["run_id", "virus_name", "v3_host_conservative"], [])
    old, new = tmp_path / "attempt2", tmp_path / "attempt3"
    for root in (old, new):
        (root / "status").mkdir(parents=True)
    bad = tmp_path / "bad"
    bad.mkdir()
    (bad / "summary.txt").write_text("Gene ID; Gene Count\nnot a gene row\n")
    (old / "status" / "v2__SRR6825024.json").write_text(json.dumps({"status": "failed"}))
    (new / "status" / "v2__SRR6825024.json").write_text(
        json.dumps({"status": "success", "v2_output_root": str(bad)})
    )
    report = compare_fresh_vs_archive([old, new], archive)
    entry = next(r for r in report["run_status"] if r["task_id"] == "v2__SRR6825024")
    assert entry["record_source"] == "attempt3/status"
    assert entry["status"] == "not_compared" and "malformed" in entry["reason"]


def _run_spool_row(tmp_path: Path, row: dict[str, str], tier: str = "highmem") -> list[str]:
    """Execute run_fresh_control_task.sh from a simulated Slurm spool against one row."""

    import os
    import subprocess

    from scripts.prepare_fresh_controls import TASK_FIELDS

    repo = Path(__file__).resolve().parents[1]
    runner = repo / "scripts/run_fresh_control_task.sh"
    packet = tmp_path / "legacy/packet"
    (packet / "source").mkdir(parents=True)
    argv_file = tmp_path / "argv"
    fake_python = tmp_path / "legacy/env_full/bin/python"
    fake_python.parent.mkdir(parents=True)
    fake_python.write_text(f'#!/bin/sh\nprintf "%s\\n" "$@" > "{argv_file}"\n')
    fake_python.chmod(0o755)
    (packet / "source/run_fresh_control.py").write_text("x")
    (packet / "run_task.sh").write_bytes(runner.read_bytes())
    (packet / "control_inputs.raw.tsv").write_text("x")
    for name in ("tasks.tsv", "tasks.small.tsv", "tasks.large.tsv", "tasks.highmem.tsv"):
        (packet / name).write_text(
            "\t".join(TASK_FIELDS) + "\n" + "\t".join(row[f] for f in TASK_FIELDS) + "\n"
        )
    names = sorted(p.relative_to(packet).as_posix() for p in packet.rglob("*") if p.is_file())
    (packet / "packet.sha256").write_text(
        "".join(f"{hashlib.sha256((packet / n).read_bytes()).hexdigest()}  {n}\n" for n in names)
    )
    spool = tmp_path / "spool/slurm_script"
    spool.parent.mkdir()
    spool.write_bytes(runner.read_bytes())
    spool.chmod(0o755)
    env = os.environ | {
        "FRESH_PACKET_ROOT": str(packet),
        "FRESH_TASK_MANIFEST": str(packet / f"tasks.{tier}.tsv"),
        "SLURM_ARRAY_TASK_ID": "0",
    }
    done = subprocess.run([str(spool)], capture_output=True, text=True, env=env, check=False)
    assert done.returncode == 0, done.stderr
    return argv_file.read_text().splitlines()


def _row(**over: str) -> dict[str, str]:
    from scripts.prepare_fresh_controls import TASK_FIELDS

    row = {f: f"v_{f}" for f in TASK_FIELDS} | {
        "viralscan_cache_path": "-",
        "viralscan_cache_manifest_sha256": "-",
        "gtf_path": "-",
        "gtf_sha256": "-",
        "cores": "8",
    }
    return row | over


def test_highmem_row_with_empty_sentinels_runs_without_cache_or_gtf(tmp_path: Path) -> None:
    argv = _run_spool_row(tmp_path, _row(task_id="v2__SRR6825024", tier="highmem"))
    assert "--viralscan-cache" not in argv and "--gtf" not in argv
    assert argv[argv.index("--attempt-id") + 1] == "v_attempt_id"


def test_gtf_row_forwards_explicit_gtf_and_digest(tmp_path: Path) -> None:
    argv = _run_spool_row(
        tmp_path, _row(task_id="v3__X", gtf_path="/p/v2_panel.gtf", gtf_sha256="a" * 64), "small"
    )
    assert argv[argv.index("--gtf") + 1] == "/p/v2_panel.gtf"
    assert argv[argv.index("--gtf-sha256") + 1] == "a" * 64


def test_extra_gtf_appended_and_parity_residual_acceptance_is_fail_closed(tmp_path: Path) -> None:
    from scripts.prepare_fresh_controls import main, residual_sha256

    paths = _inputs(tmp_path)
    data = tmp_path / "v2data"
    data.mkdir()
    (data / "a.gtf").write_text(GTF_A)
    extra = tmp_path / "extra.gtf"
    extra.write_text(GTF_B)
    # t2g has VIR_B1 (via the extra GTF) but lacks VIR_A1 -> a gtf-only residual remains
    paths["t2g"].write_text("t1\tVIR_B1\nt2\tVIR_C1\nt3\tENSG1\n")
    argv = [
        "--raw-manifest", str(paths["raw_manifest"]),
        "--output-root", str(tmp_path / "out"),
        "--v2-viralscan", str(paths["v2_viralscan"]),
        "--v3-viralscan", str(paths["v3_viralscan"]),
        "--index", str(paths["index"]),
        "--t2g", str(paths["t2g"]),
        "--whitelist", str(paths["whitelist"]),
        "--attempt-id", "a",
        "--v2-data-dir", str(data),
        "--extra-gtf", str(extra),
    ]
    packet = tmp_path / "p1"
    with pytest.raises(FreshControlPreparationError, match="--accept-parity-residual"):
        main([*argv, "--task-manifest", str(packet / "tasks.tsv")])
    ref = packet / "reference"
    assert (ref / "v2_panel.gtf").read_text() == GTF_A + GTF_B + "\n"
    parity = json.loads((ref / "gtf_t2g_parity.json").read_text())
    assert parity["gtf_only_genes"] == ["VIR_A1"] and parity["intersection"] == 1
    assert parity["residual_sha256"] == residual_sha256(parity)
    prov = json.loads((ref / "panel_provenance.json").read_text())
    assert prov["extra_gtfs"][0]["sha256"] == hashlib.sha256(extra.read_bytes()).hexdigest()
    with pytest.raises(FreshControlPreparationError):  # wrong digest stays fail-closed
        main([*argv, "--task-manifest", str(tmp_path / "p2" / "tasks.tsv"),
              "--accept-parity-residual", "0" * 64])
    assert main([*argv, "--task-manifest", str(tmp_path / "p3" / "tasks.tsv"),
                 "--accept-parity-residual", parity["residual_sha256"]]) == 0
    accepted = json.loads((tmp_path / "p3/reference/gtf_t2g_parity.json").read_text())
    assert accepted["accepted_residual"] is True
