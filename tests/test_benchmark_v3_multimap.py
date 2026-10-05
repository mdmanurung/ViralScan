"""Focused tests for the retained-BUS v3 benchmark harness."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scripts import benchmark_v3_multimap

pytestmark = pytest.mark.research


def _tree_snapshot(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def test_code_fingerprints_follow_loaded_wheel_modules_when_script_is_frozen(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    frozen_script = tmp_path / "frozen" / "benchmark_v3_multimap.py"
    frozen_script.parent.mkdir()
    frozen_script.write_bytes(b"frozen benchmark entry point\n")
    wheel_root = tmp_path / "site-packages" / "viralscan"
    module_files = {
        "multimapping_module": wheel_root / "multimapping.py",
        "multimap_worker_module": wheel_root / "scripts" / "multimap.py",
        "virus_grouping_module": wheel_root / "virus_grouping.py",
        "virus_constants_module": wheel_root / "constants.py",
    }
    for name, path in module_files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"wheel artifact for {name}\n".encode())
        monkeypatch.setattr(
            benchmark_v3_multimap,
            name,
            SimpleNamespace(__file__=str(path)),
            raising=False,
        )
    monkeypatch.setattr(benchmark_v3_multimap, "__file__", str(frozen_script))

    fingerprints = benchmark_v3_multimap._code_fingerprints()

    assert (
        fingerprints["benchmark_v3_multimap"]["sha256"]
        == hashlib.sha256(frozen_script.read_bytes()).hexdigest()
    )
    expected_keys = {
        "multimapping": "multimapping_module",
        "multimap_worker": "multimap_worker_module",
        "virus_grouping": "virus_grouping_module",
        "virus_constants": "virus_constants_module",
    }
    for fingerprint_name, module_name in expected_keys.items():
        path = module_files[module_name]
        assert (
            fingerprints[fingerprint_name]["sha256"]
            == hashlib.sha256(path.read_bytes()).hexdigest()
        )


@dataclass(frozen=True)
class FakeAudit:
    input_molecules: int = 1
    resolved_molecules: int = 1
    unique_molecules: int = 1
    ambiguous_molecules: int = 0
    unresolved_molecules: int = 0
    ignored_read_multiplicity: int = 0

    def validate(self, _mass: float) -> None:
        return None


def test_cli_requires_whitelist_and_accepts_an_explicit_sample_id(tmp_path: Path) -> None:
    args = benchmark_v3_multimap.parse_args(
        [
            "--sample-id",
            "skin-run-17",
            "--sample-dir",
            str(tmp_path / "source"),
            "--whitelist",
            str(tmp_path / "3M-february-2018.txt.gz"),
            "--t2g",
            str(tmp_path / "t2g.txt"),
            "--output",
            str(tmp_path / "result.json"),
        ]
    )

    assert args.sample_id == "skin-run-17"
    assert args.whitelist == tmp_path / "3M-february-2018.txt.gz"

    with pytest.raises(SystemExit):
        benchmark_v3_multimap.parse_args(
            [
                "--sample-id",
                "skin-run-17",
                "--sample-dir",
                str(tmp_path / "source"),
                "--t2g",
                str(tmp_path / "t2g.txt"),
                "--output",
                str(tmp_path / "result.json"),
            ]
        )


def test_main_applies_the_declared_whitelist_to_the_raw_bus(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sample_dir = tmp_path / "source"
    kb = sample_dir / "kb-python"
    counts = kb / "counts_unfiltered"
    log = sample_dir / "log"
    counts.mkdir(parents=True)
    log.mkdir()
    required_files = [
        kb / "output.bus",
        kb / "matrix.ec",
        kb / "transcripts.txt",
        counts / "cells_x_genes.barcodes.txt",
        counts / "cells_x_genes.genes.txt",
        counts / "cells_x_genes.genes.names.txt",
        counts / "adata.h5ad",
    ]
    for path in required_files:
        path.write_text("", encoding="utf-8")
    (log / "analysis.txt").write_text("gene\n", encoding="utf-8")
    (kb / "run_info.json").write_text('{"n_pseudoaligned": 7}\n', encoding="utf-8")
    whitelist = tmp_path / "whitelist.txt"
    whitelist.write_text("AAAA\n", encoding="utf-8")
    t2g = tmp_path / "t2g.txt"
    t2g.write_text("tx\tgene\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_prepare(
        raw_bus,
        resolved_bus,
        resolved_text,
        *,
        whitelist,
        threads,
        corrected_bus,
    ) -> None:
        Path(resolved_bus).parent.mkdir(parents=True, exist_ok=True)
        Path(resolved_bus).write_bytes(b"resolved")
        Path(resolved_text).write_text("AAAA\tUMI\t0\t1\n", encoding="utf-8")
        Path(corrected_bus).write_bytes(b"corrected")
        captured.update(
            raw_bus=raw_bus,
            resolved_bus=resolved_bus,
            resolved_text=resolved_text,
            corrected_bus=corrected_bus,
            whitelist=whitelist,
            threads=threads,
        )

    monkeypatch.setattr(benchmark_v3_multimap, "prepare_resolved_bus", fake_prepare)
    source_adata = ad.AnnData(
        X=sparse.csr_matrix([[1.0]]),
        obs=pd.DataFrame(index=["AAAA-1"]),
        var=pd.DataFrame(index=["gene"]),
    )
    read_h5ad = ad.read_h5ad
    monkeypatch.setattr(
        benchmark_v3_multimap.ad,
        "read_h5ad",
        lambda path, **kwargs: (
            source_adata.copy()
            if Path(path) == counts / "adata.h5ad"
            else read_h5ad(path, **kwargs)
        ),
    )
    monkeypatch.setattr(benchmark_v3_multimap, "load_barcodes", lambda _path: ({"AAAA": 0}, 1))
    monkeypatch.setattr(benchmark_v3_multimap, "load_genes", lambda *_args: (["gene"], ["Gene"]))
    monkeypatch.setattr(
        benchmark_v3_multimap,
        "load_transcripts",
        lambda *_args: (["tx"], {"tx": "gene"}),
    )
    monkeypatch.setattr(benchmark_v3_multimap, "read_ec", lambda *_args: {0: [0]})
    audit = FakeAudit()
    layers = SimpleNamespace(
        corrected=sparse.csr_matrix([[0.0]]),
        unique=sparse.csr_matrix([[1.0]]),
        equal=sparse.csr_matrix([[0.0]]),
        host_conservative=sparse.csr_matrix([[0.0]]),
        audit=audit,
        method_diagnostics={"method": "host-conservative"},
    )

    def fake_build(*_args, **_kwargs):
        captured["allocator_calls"] = int(captured.get("allocator_calls", 0)) + 1
        return layers

    monkeypatch.setattr(benchmark_v3_multimap, "build_multimap_layers", fake_build)

    output = tmp_path / "out" / "result.json"
    scratch_root = tmp_path / "scratch"
    scratch_root.mkdir()
    monkeypatch.setenv("TMPDIR", str(scratch_root))
    source_before = _tree_snapshot(sample_dir)
    assert (
        benchmark_v3_multimap.main(
            [
                "--sample-id",
                "skin-run-17",
                "--sample-dir",
                str(sample_dir),
                "--whitelist",
                str(whitelist),
                "--t2g",
                str(t2g),
                "--output",
                str(output),
                "--threads",
                "3",
            ]
        )
        == 0
    )

    assert captured["raw_bus"] == kb / "output.bus"
    assert captured["resolved_bus"] == output.parent / "output.resolved.sorted.bus"
    assert Path(captured["resolved_text"]).parent.parent == scratch_root
    assert Path(captured["corrected_bus"]).parent == Path(captured["resolved_text"]).parent
    for generated in (
        captured["resolved_bus"],
        captured["resolved_text"],
        captured["corrected_bus"],
    ):
        assert not Path(generated).resolve().is_relative_to(sample_dir.resolve())
    assert not Path(captured["resolved_text"]).parent.exists()
    assert captured["whitelist"] == str(whitelist)
    assert captured["threads"] == 3
    assert captured["allocator_calls"] == 1
    result = json.loads(output.read_text())
    assert result["sample"] == "skin-run-17"
    assert result["primary_endpoint"] == "v3-host-conservative"
    assert set(result["comparison_endpoints"]) == {
        "v3-unique",
        "v3-equal",
        "v3-host-conservative",
    }
    exported = read_h5ad(output.parent / "adata_v3.h5ad")
    assert exported.layers["counts_unique"].toarray().tolist() == [[1.0]]
    assert exported.layers["counts_ambiguous_allocated"].toarray().tolist() == [[0.0]]
    assert exported.X.toarray().tolist() == [[1.0]]
    assert pd.read_csv(output.parent / "per_virus.tsv", sep="\t").to_dict("records") == [
        {
            "virus_name": "gene",
            "n_features": 1,
            "v3_unique": 1.0,
            "v3_equal": 1.0,
            "v3_host_conservative": 1.0,
        }
    ]
    assert pd.read_csv(output.parent / "per_cell.tsv", sep="\t").to_dict("records") == [
        {
            "barcode": "AAAA",
            "v3_unique_viral": 1.0,
            "v3_equal_viral": 1.0,
            "v3_host_conservative_viral": 1.0,
        }
    ]
    assert json.loads((output.parent / "count_audit.json").read_text())["input_molecules"] == 1
    hashes = json.loads((output.parent / "hashes.json").read_text())
    assert len(hashes["outputs"]["adata_v3.h5ad"]["sha256"]) == 64
    assert len(hashes["outputs"]["output.resolved.sorted.bus"]["sha256"]) == 64
    assert hashes["inputs_before"] == hashes["inputs_after"]
    assert set(hashes["code"]) >= {
        "benchmark_v3_multimap",
        "multimapping",
        "multimap_worker",
        "virus_grouping",
        "virus_constants",
    }
    assert len(hashes["virus_grouping_sha256"]) == 64
    status = json.loads((output.parent / "status.json").read_text())
    assert status["status"] == "success"
    assert status["exit_code"] == 0
    assert status["stage"] == "complete"
    assert status["scratch_retained"] is False
    assert status["scratch_bytes_before_cleanup"] > 0
    assert status["total_wall_seconds"] >= result["allocation_seconds"]
    assert len(status["attempt_id"]) >= 16
    assert len(status["scientific_parameter_hash"]) == 64
    assert not (output.parent / "failure.json").exists()
    assert _tree_snapshot(sample_dir) == source_before


def test_endpoint_metrics_report_all_frozen_v3_counting_arms() -> None:
    layers = SimpleNamespace(
        unique=np.array([[2.0, 0.0]]),
        equal=np.array([[1.5, 1.5]]),
        host_conservative=np.array([[3.0, 0.0]]),
    )

    assert benchmark_v3_multimap.endpoint_mass_metrics(layers) == {
        "v3-unique": {
            "allocation_method": "unique-only",
            "unique_molecule_mass": 2.0,
            "allocated_ambiguous_mass": 0.0,
            "selected_matrix_mass": 2.0,
        },
        "v3-equal": {
            "allocation_method": "equal",
            "unique_molecule_mass": 2.0,
            "allocated_ambiguous_mass": 3.0,
            "selected_matrix_mass": 5.0,
        },
        "v3-host-conservative": {
            "allocation_method": "host-conservative",
            "unique_molecule_mass": 2.0,
            "allocated_ambiguous_mass": 3.0,
            "selected_matrix_mass": 5.0,
        },
    }


def test_main_records_machine_readable_failure_status(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fail(*_args, **_kwargs):
        raise ValueError("counts_unique contains non-finite values")

    monkeypatch.setattr(benchmark_v3_multimap, "_run", fail)
    output = tmp_path / "failed" / "result.json"

    with pytest.raises(ValueError, match="non-finite"):
        benchmark_v3_multimap.main(
            [
                "--sample-id",
                "bad-row",
                "--sample-dir",
                str(tmp_path / "source"),
                "--whitelist",
                str(tmp_path / "whitelist.txt"),
                "--t2g",
                str(tmp_path / "t2g.txt"),
                "--output",
                str(output),
            ]
        )

    failure = json.loads((output.parent / "failure.json").read_text())
    assert failure["schema_version"] == "3.0.0"
    assert failure["sample_id"] == "bad-row"
    assert failure["status"] == "failed"
    assert failure["error_type"] == "ValueError"
    assert failure["message"] == "counts_unique contains non-finite values"
    assert failure["stage"] == "initialization"
    assert failure["exit_code"] == 1
    assert failure["command"][1:3] == ["--sample-id", "bad-row"]
    assert len(failure["attempt_id"]) >= 16
    assert len(failure["scientific_parameter_hash"]) == 64
    assert Path(failure["stderr_path"]).is_file()
    assert failure["total_wall_seconds"] >= 0
    assert json.loads((output.parent / "status.json").read_text()) == failure


@pytest.mark.parametrize(
    ("matrix", "message"),
    [
        (sparse.csr_matrix([[np.nan]]), "non-finite"),
        (sparse.csr_matrix([[-1.0]]), "negative"),
        (sparse.csr_matrix([[1.0, 0.0]]), "shape"),
    ],
)
def test_matrix_validation_fails_closed(matrix, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        benchmark_v3_multimap._validate_matrix("counts_unique", matrix, (1, 1))


def test_tmpdir_preflight_requires_three_times_raw_bus_space(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    raw_bus = tmp_path / "output.bus"
    raw_bus.write_bytes(b"0123456789")
    monkeypatch.setattr(
        benchmark_v3_multimap.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(total=100, used=71, free=29),
    )

    with pytest.raises(OSError, match="requires at least 30 bytes.*29 bytes"):
        benchmark_v3_multimap._require_tmpdir_capacity(raw_bus, tmp_path)


def test_serialized_h5ad_validation_checks_all_layers_and_x_equality(tmp_path: Path) -> None:
    path = tmp_path / "corrupt.h5ad"
    adata = ad.AnnData(
        X=sparse.csr_matrix([[9.0, 0.0]]),
        obs=pd.DataFrame(index=["BC"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    adata.layers["counts_unique"] = sparse.csr_matrix([[1.0, 0.0]])
    adata.layers["counts_ambiguous_allocated"] = sparse.csr_matrix([[0.0, 1.0]])
    adata.layers["counts_multimap_equal"] = sparse.csr_matrix([[0.5, 0.5]])
    adata.layers["counts_multimap_host_conservative"] = sparse.csr_matrix([[1.0, 0.0]])
    adata.write_h5ad(path)

    with pytest.raises(ValueError, match="X does not equal"):
        benchmark_v3_multimap._validate_written_h5ad(path, (1, 2))


def test_failed_bus_preparation_retains_scratch_boundary_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sample_dir = tmp_path / "source"
    kb = sample_dir / "kb-python"
    counts = kb / "counts_unfiltered"
    log = sample_dir / "log"
    counts.mkdir(parents=True)
    log.mkdir()
    for path in (
        kb / "output.bus",
        kb / "matrix.ec",
        kb / "transcripts.txt",
        counts / "cells_x_genes.barcodes.txt",
        counts / "cells_x_genes.genes.txt",
        counts / "cells_x_genes.genes.names.txt",
        counts / "adata.h5ad",
        kb / "run_info.json",
        log / "analysis.txt",
    ):
        path.write_bytes(b"x")
    whitelist = tmp_path / "whitelist.txt"
    t2g = tmp_path / "t2g.txt"
    whitelist.write_text("AAAA\n")
    t2g.write_text("tx\tgene\n")
    scratch_root = tmp_path / "scratch"
    scratch_root.mkdir()
    monkeypatch.setenv("TMPDIR", str(scratch_root))

    def fail_prepare(
        _raw_bus,
        _resolved_bus,
        resolved_text,
        *,
        whitelist,
        threads,
        corrected_bus,
    ) -> None:
        assert whitelist
        assert threads == 8
        Path(corrected_bus).write_bytes(b"corrected")
        Path(resolved_text).write_text("partial\n")
        raise RuntimeError("bustools sort failed")

    monkeypatch.setattr(benchmark_v3_multimap, "prepare_resolved_bus", fail_prepare)
    output = tmp_path / "out" / "result.json"
    with pytest.raises(RuntimeError, match="bustools sort failed"):
        benchmark_v3_multimap.main(
            [
                "--sample-id",
                "failed-row",
                "--sample-dir",
                str(sample_dir),
                "--whitelist",
                str(whitelist),
                "--t2g",
                str(t2g),
                "--output",
                str(output),
            ]
        )

    failure = json.loads((output.parent / "failure.json").read_text())
    scratch_path = Path(failure["scratch_path"])
    assert failure["stage"] == "prepare_resolved_bus"
    assert failure["scratch_retained"] is True
    assert scratch_path.is_dir()
    assert (scratch_path / "output.corrected.bus").read_bytes() == b"corrected"
    assert (scratch_path / "output.resolved.sorted.bus.txt").read_text() == "partial\n"
