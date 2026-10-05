import json
from pathlib import Path

import pytest

from scripts.revalidate_fresh_control import (
    REQUIRED_FAILURE_FIELDS,
    RevalidationError,
    revalidate,
)

V2_REQUIRED_ARTIFACTS = (
    "config.yaml",
    "summary.txt",
    "kb-python/run_info.json",
    "kb-python/output.bus",
    "kb-python/counts_unfiltered/adata.h5ad",
    "kb-python/counts_unfiltered/adata_multimap.h5ad",
)


def _write_v2_tree(root: Path) -> None:
    for relative in V2_REQUIRED_ARTIFACTS:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(relative, encoding="utf-8")


def _write_original_status(
    tmp_path: Path,
    *,
    output: Path,
    workflow_exit_code: int = 0,
    stack: str = "v2",
) -> Path:
    status = tmp_path / "status" / "v2__SRR12682296.json"
    status.parent.mkdir(parents=True, exist_ok=True)
    status.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "sample_id": "SRR12682296",
                "stack": stack,
                "status": "failed",
                "exit_code": 65,
                "workflow_exit_code": workflow_exit_code,
                "artifact_validation_errors": list(V2_REQUIRED_ARTIFACTS),
                "command": ["viralscan", "-i", "reads"],
                "command_sha256": "b" * 64,
                "output": str(output),
                "stdout": str(tmp_path / "logs/out.log"),
                "stderr": str(tmp_path / "logs/err.log"),
            }
        ),
        encoding="utf-8",
    )
    return status


def test_revalidation_clears_the_phantom_artifact_errors(tmp_path: Path) -> None:
    """The four attempt-2 rows had complete nested trees despite exit 65."""

    output = tmp_path / "runs/v2/SRR12682296"
    _write_v2_tree(output / "SRR12682296")
    status = _write_original_status(tmp_path, output=output)
    record = tmp_path / "revalidation/v2__SRR12682296.json"

    payload = revalidate(status_path=status, record_path=record, attempt_id="attempt2-revalidated")

    assert payload["status"] == "success"
    assert payload["exit_code"] == 0
    assert payload["artifact_validation_errors"] == []
    assert payload["original_exit_code"] == 65
    assert payload["original_artifact_validation_errors"] == list(V2_REQUIRED_ARTIFACTS)
    assert payload["v2_output_root"] == str((output / "SRR12682296").resolve())
    assert record.is_file()


def test_revalidation_records_every_required_field(tmp_path: Path) -> None:
    output = tmp_path / "runs/v2/SRR12682296"
    _write_v2_tree(output / "SRR12682296")
    status = _write_original_status(tmp_path, output=output)
    record = tmp_path / "revalidation/v2__SRR12682296.json"

    payload = revalidate(status_path=status, record_path=record, attempt_id="attempt2-revalidated")

    for field in REQUIRED_FAILURE_FIELDS:
        assert field in payload, field
    assert payload["stage"] == "artifact_validation"
    assert payload["attempt_id"] == "attempt2-revalidated"
    assert payload["scientific_parameter_hash"] == "b" * 64


def test_revalidation_never_overwrites_the_original_status(tmp_path: Path) -> None:
    output = tmp_path / "runs/v2/SRR12682296"
    _write_v2_tree(output / "SRR12682296")
    status = _write_original_status(tmp_path, output=output)
    before = status.read_text(encoding="utf-8")
    record = tmp_path / "revalidation/v2__SRR12682296.json"

    payload = revalidate(status_path=status, record_path=record, attempt_id="attempt2-revalidated")

    assert status.read_text(encoding="utf-8") == before
    assert payload["supersedes"] == str(status.resolve())


def test_revalidation_refuses_to_replace_an_existing_record(tmp_path: Path) -> None:
    output = tmp_path / "runs/v2/SRR12682296"
    _write_v2_tree(output / "SRR12682296")
    status = _write_original_status(tmp_path, output=output)
    record = tmp_path / "revalidation/v2__SRR12682296.json"
    record.parent.mkdir(parents=True)
    record.write_text("prior record\n", encoding="utf-8")

    with pytest.raises(RevalidationError, match="already exists"):
        revalidate(status_path=status, record_path=record, attempt_id="attempt2-revalidated")

    assert record.read_text(encoding="utf-8") == "prior record\n"


def test_revalidation_refuses_a_genuinely_failed_workflow(tmp_path: Path) -> None:
    """The out-of-memory row must not be laundered into a success."""

    output = tmp_path / "runs/v2/SRR6825024"
    status = _write_original_status(tmp_path, output=output, workflow_exit_code=1)
    record = tmp_path / "revalidation/v2__SRR6825024.json"

    with pytest.raises(RevalidationError, match="workflow did not succeed"):
        revalidate(status_path=status, record_path=record, attempt_id="attempt2-revalidated")

    assert not record.exists()


def test_revalidation_still_fails_on_a_genuinely_incomplete_tree(tmp_path: Path) -> None:
    output = tmp_path / "runs/v2/SRR12682296"
    (output / "SRR12682296").mkdir(parents=True)
    (output / "SRR12682296" / "config.yaml").write_text("config", encoding="utf-8")
    status = _write_original_status(tmp_path, output=output)
    record = tmp_path / "revalidation/v2__SRR12682296.json"

    payload = revalidate(status_path=status, record_path=record, attempt_id="attempt2-revalidated")

    assert payload["status"] == "failed"
    assert payload["exit_code"] == 65
    assert payload["artifact_validation_errors"]
