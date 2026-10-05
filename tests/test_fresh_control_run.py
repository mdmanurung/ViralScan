from pathlib import Path

import pytest

from scripts.run_fresh_control import (
    FreshControlError,
    build_command,
    ensure_fresh_output,
    resolve_v2_output_root,
    run_control,
    sha256_file,
    verify_frozen_fastq,
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
    """Create a complete, non-empty legacy 2.2.0 result tree at ``root``."""

    for relative in V2_REQUIRED_ARTIFACTS:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(relative, encoding="utf-8")


def test_fresh_control_commands_share_frozen_inputs_and_fix_stack_policies(
    tmp_path: Path,
) -> None:
    inputs = {
        "read1": tmp_path / "sample_1.fastq.gz",
        "read2": tmp_path / "sample_2.fastq.gz",
        "index": tmp_path / "original.idx",
        "t2g": tmp_path / "original.t2g",
        "whitelist": tmp_path / "10xv2.txt",
    }
    v2 = build_command(
        stack="v2",
        viralscan=tmp_path / "v2/bin/viralscan",
        output=tmp_path / "runs/v2/sample",
        cores=8,
        **inputs,
    )
    v3 = build_command(
        stack="v3",
        viralscan=tmp_path / "v3/bin/viralscan",
        output=tmp_path / "runs/v3/sample",
        cores=8,
        **inputs,
    )

    for value in inputs.values():
        assert str(value) in v2
        assert str(value) in v3
    assert v2[v2.index("-x") + 1] == "10xv2"
    assert v2[v2.index("-mm") + 1] == "True"
    assert v3[v3.index("--multimap-method") + 1] == "host-conservative"
    assert v3[v3.index("--cell-calling") + 1] == "auto"
    assert "--multimapping" in v3


def test_fresh_control_refuses_an_existing_output_tree(tmp_path: Path) -> None:
    output = tmp_path / "existing"
    output.mkdir()

    with pytest.raises(FreshControlError, match="already exists"):
        ensure_fresh_output(output)


def test_fresh_control_refuses_fastq_storage_size_drift(tmp_path: Path) -> None:
    fastq = tmp_path / "sample.fastq.gz"
    fastq.write_bytes(b"drifted")

    with pytest.raises(FreshControlError, match="stored byte count drifted"):
        verify_frozen_fastq(
            fastq,
            expected_bytes=99,
            expected_sha256="a" * 64,
        )


def test_fresh_control_refuses_a_same_size_content_swap(tmp_path: Path) -> None:
    """A stale restore or repointed symlink preserves the byte count."""

    audited = tmp_path / "audited.fastq.gz"
    audited.write_bytes(b"the frozen input")
    audited_sha256 = sha256_file(audited)

    swapped = tmp_path / "swapped.fastq.gz"
    swapped.write_bytes(b"a totally other!")
    assert swapped.stat().st_size == audited.stat().st_size

    with pytest.raises(FreshControlError, match="content drifted"):
        verify_frozen_fastq(
            swapped,
            expected_bytes=audited.stat().st_size,
            expected_sha256=audited_sha256,
        )


def test_fresh_control_accepts_an_unchanged_frozen_input(tmp_path: Path) -> None:
    fastq = tmp_path / "sample.fastq.gz"
    fastq.write_bytes(b"the frozen input")

    result = verify_frozen_fastq(
        fastq,
        expected_bytes=fastq.stat().st_size,
        expected_sha256=sha256_file(fastq),
    )

    assert result["runtime_sha256_verified"] is True
    assert result["audited_storage_sha256"] == sha256_file(fastq)


def test_fresh_control_retains_a_failed_stack_status(tmp_path: Path) -> None:
    viralscan = tmp_path / "bin/viralscan"
    viralscan.parent.mkdir()
    viralscan.write_text("#!/bin/sh\nexit 3\n", encoding="utf-8")
    viralscan.chmod(0o755)
    inputs = {}
    for name in ("read1", "read2", "index", "t2g", "whitelist"):
        inputs[name] = tmp_path / name
        inputs[name].write_text(name, encoding="utf-8")
    status = tmp_path / "attempt/status.json"

    payload = run_control(
        sample_id="SRR12682296",
        stack="v2",
        viralscan=viralscan,
        output=tmp_path / "output",
        cores=8,
        status_path=status,
        stdout_path=tmp_path / "attempt/stdout.log",
        stderr_path=tmp_path / "attempt/stderr.log",
        attempt_id="attempt-test",
        **inputs,
    )

    assert payload["status"] == "failed"
    assert payload["exit_code"] == 3
    assert status.is_file()


@pytest.mark.parametrize("artifact_name", ("status", "stdout", "stderr"))
def test_fresh_control_refuses_to_overwrite_attempt_evidence(
    tmp_path: Path,
    artifact_name: str,
) -> None:
    viralscan = tmp_path / "bin/viralscan"
    viralscan.parent.mkdir()
    viralscan.write_text("#!/bin/sh\nexit 3\n", encoding="utf-8")
    viralscan.chmod(0o755)
    inputs = {}
    for name in ("read1", "read2", "index", "t2g", "whitelist"):
        inputs[name] = tmp_path / name
        inputs[name].write_text(name, encoding="utf-8")
    attempt_paths = {
        "status": tmp_path / "attempt/status.json",
        "stdout": tmp_path / "attempt/stdout.log",
        "stderr": tmp_path / "attempt/stderr.log",
    }
    attempt_paths[artifact_name].parent.mkdir(parents=True, exist_ok=True)
    attempt_paths[artifact_name].write_text("prior evidence\n", encoding="utf-8")

    with pytest.raises(FreshControlError, match="attempt evidence already exists"):
        run_control(
            sample_id="SRR12682296",
            stack="v2",
            viralscan=viralscan,
            output=tmp_path / "output",
            cores=8,
            status_path=attempt_paths["status"],
            stdout_path=attempt_paths["stdout"],
            stderr_path=attempt_paths["stderr"],
            attempt_id="attempt-test",
            **inputs,
        )

    assert attempt_paths[artifact_name].read_text(encoding="utf-8") == "prior evidence\n"


def test_fresh_v3_status_fails_when_validate_run_fails(tmp_path: Path) -> None:
    viralscan = tmp_path / "bin/viralscan"
    viralscan.parent.mkdir()
    viralscan.write_text(
        '#!/bin/sh\nif [ "$1" = "validate-run" ]; then exit 4; fi\nexit 0\n',
        encoding="utf-8",
    )
    viralscan.chmod(0o755)
    inputs = {}
    for name in ("read1", "read2", "index", "t2g", "whitelist"):
        inputs[name] = tmp_path / name
        inputs[name].write_text(name, encoding="utf-8")

    payload = run_control(
        sample_id="SRR12682296",
        stack="v3",
        viralscan=viralscan,
        output=tmp_path / "output",
        cores=8,
        status_path=tmp_path / "attempt/status.json",
        stdout_path=tmp_path / "attempt/stdout.log",
        stderr_path=tmp_path / "attempt/stderr.log",
        attempt_id="attempt-test",
        **inputs,
    )

    assert payload["workflow_exit_code"] == 0
    assert payload["validation_exit_code"] == 4
    assert payload["status"] == "failed"


def test_fresh_v2_status_fails_when_required_outputs_are_missing(
    tmp_path: Path,
) -> None:
    viralscan = tmp_path / "bin/viralscan"
    viralscan.parent.mkdir()
    viralscan.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    viralscan.chmod(0o755)
    inputs = {}
    for name in ("read1", "read2", "index", "t2g", "whitelist"):
        inputs[name] = tmp_path / name
        inputs[name].write_text(name, encoding="utf-8")

    payload = run_control(
        sample_id="SRR12682296",
        stack="v2",
        viralscan=viralscan,
        output=tmp_path / "output",
        cores=8,
        status_path=tmp_path / "attempt/status.json",
        stdout_path=tmp_path / "attempt/stdout.log",
        stderr_path=tmp_path / "attempt/stderr.log",
        attempt_id="attempt-test",
        **inputs,
    )

    assert payload["workflow_exit_code"] == 0
    assert payload["status"] == "failed"
    assert payload["artifact_validation_errors"]
    assert payload["stage"] == "artifact_validation"


def test_resolve_v2_output_root_finds_the_nested_sample_directory(tmp_path: Path) -> None:
    """The legacy 2.2.0 CLI nests its tree one level under the sample identifier."""

    output = tmp_path / "output"
    _write_v2_tree(output / "SRR12682296")

    root, errors = resolve_v2_output_root(output)

    assert errors == []
    assert root == output / "SRR12682296"


def test_resolve_v2_output_root_accepts_a_flat_tree(tmp_path: Path) -> None:
    output = tmp_path / "output"
    _write_v2_tree(output)

    root, errors = resolve_v2_output_root(output)

    assert errors == []
    assert root == output


def test_resolve_v2_output_root_refuses_an_ambiguous_layout(tmp_path: Path) -> None:
    """Two candidate roots must fail rather than silently pick one."""

    output = tmp_path / "output"
    _write_v2_tree(output / "SRR12682296")
    _write_v2_tree(output / "SRR12682297")

    root, errors = resolve_v2_output_root(output)

    assert root == output
    assert len(errors) == 1
    assert "ambiguous" in errors[0]


def test_resolve_v2_output_root_refuses_a_missing_tree(tmp_path: Path) -> None:
    root, errors = resolve_v2_output_root(tmp_path / "absent")

    assert root == tmp_path / "absent"
    assert len(errors) == 1
    assert "missing v2 output directory" in errors[0]


def test_fresh_v2_status_succeeds_on_a_nested_output_tree(tmp_path: Path) -> None:
    """Regression: attempt 2 recorded exit 65 on four rows that had in fact succeeded."""

    viralscan = tmp_path / "bin/viralscan"
    viralscan.parent.mkdir()
    output = tmp_path / "output"
    nested = output / "SRR12682296"
    writes = "\n".join(
        f'mkdir -p "$(dirname "{nested / relative}")" && echo x > "{nested / relative}"'
        for relative in V2_REQUIRED_ARTIFACTS
    )
    viralscan.write_text(f"#!/bin/sh\n{writes}\nexit 0\n", encoding="utf-8")
    viralscan.chmod(0o755)
    inputs = {}
    for name in ("read1", "read2", "index", "t2g", "whitelist"):
        inputs[name] = tmp_path / name
        inputs[name].write_text(name, encoding="utf-8")

    payload = run_control(
        sample_id="SRR12682296",
        stack="v2",
        viralscan=viralscan,
        output=output,
        cores=8,
        status_path=tmp_path / "attempt/status.json",
        stdout_path=tmp_path / "attempt/stdout.log",
        stderr_path=tmp_path / "attempt/stderr.log",
        attempt_id="attempt-test",
        **inputs,
    )

    assert payload["workflow_exit_code"] == 0
    assert payload["artifact_validation_errors"] == []
    assert payload["exit_code"] == 0
    assert payload["status"] == "success"
    assert payload["stage"] == "complete"
    assert payload["v2_output_root"] == str((output / "SRR12682296").resolve())


def test_fresh_v3_validate_run_receives_the_unnested_output_path(tmp_path: Path) -> None:
    """v3 writes its manifest at the requested output dir; pin that assumption."""

    viralscan = tmp_path / "bin/viralscan"
    viralscan.parent.mkdir()
    recorded = tmp_path / "validate_argv"
    viralscan.write_text(
        f'#!/bin/sh\nif [ "$1" = "validate-run" ]; then echo "$2" > "{recorded}"; fi\nexit 0\n',
        encoding="utf-8",
    )
    viralscan.chmod(0o755)
    inputs = {}
    for name in ("read1", "read2", "index", "t2g", "whitelist"):
        inputs[name] = tmp_path / name
        inputs[name].write_text(name, encoding="utf-8")
    output = tmp_path / "output"

    payload = run_control(
        sample_id="SRR12682296",
        stack="v3",
        viralscan=viralscan,
        output=output,
        cores=8,
        status_path=tmp_path / "attempt/status.json",
        stdout_path=tmp_path / "attempt/stdout.log",
        stderr_path=tmp_path / "attempt/stderr.log",
        attempt_id="attempt-test",
        **inputs,
    )

    assert payload["status"] == "success"
    assert recorded.read_text(encoding="utf-8").strip() == str(output)


def test_fresh_status_records_every_required_failure_field(tmp_path: Path) -> None:
    """protocol.yaml failure_policy.required_failure_fields must all be present."""

    viralscan = tmp_path / "bin/viralscan"
    viralscan.parent.mkdir()
    viralscan.write_text("#!/bin/sh\nexit 3\n", encoding="utf-8")
    viralscan.chmod(0o755)
    inputs = {}
    for name in ("read1", "read2", "index", "t2g", "whitelist"):
        inputs[name] = tmp_path / name
        inputs[name].write_text(name, encoding="utf-8")

    payload = run_control(
        sample_id="SRR12682296",
        stack="v2",
        viralscan=viralscan,
        output=tmp_path / "output",
        cores=8,
        status_path=tmp_path / "attempt/status.json",
        stdout_path=tmp_path / "attempt/stdout.log",
        stderr_path=tmp_path / "attempt/stderr.log",
        attempt_id="attempt3",
        **inputs,
    )

    for field in (
        "run_id",
        "stage",
        "command",
        "exit_code",
        "stderr_path",
        "attempt_id",
        "scientific_parameter_hash",
    ):
        assert field in payload, field
    assert payload["run_id"] == "SRR12682296"
    assert payload["attempt_id"] == "attempt3"
    assert payload["stage"] == "workflow"
    assert payload["scientific_parameter_hash"] == payload["command_sha256"]
