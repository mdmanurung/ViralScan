import csv
import json
from pathlib import Path

from scripts.freeze_control_inputs import freeze_control_inputs

EBV_IDS = ("SRR12682296", "SRR12682297", "SRR12682298")
HIV_IDS = ("SRR6825024", "SRR6825025")


def _write_pair_audit(
    path: Path,
    *,
    sample_id: str,
    stored_bytes: tuple[int, int],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "sample_id": sample_id,
                "chemistry": "10xv2",
                "required_r1_bases": 26,
                "status": "valid",
                "records": 7,
                "read1": {
                    "bytes": 101,
                    "stored_bytes": stored_bytes[0],
                    "compression": "gzip",
                    "min_sequence_bases": 26,
                    "max_sequence_bases": 26,
                    "md5": "1" * 32,
                    "sha256": "1" * 64,
                },
                "read2": {
                    "bytes": 202,
                    "stored_bytes": stored_bytes[1],
                    "compression": "gzip",
                    "min_sequence_bases": 90,
                    "max_sequence_bases": 90,
                    "md5": "2" * 32,
                    "sha256": "2" * 64,
                },
            }
        ),
        encoding="utf-8",
    )


def test_freeze_control_inputs_emits_one_sanitized_five_sample_manifest(
    tmp_path: Path,
) -> None:
    ena_root = tmp_path / "ena"
    (ena_root / "status").mkdir(parents=True)
    (ena_root / "files").mkdir()
    with (ena_root / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sample_id", "mate", "url", "md5", "bytes"),
            delimiter="\t",
        )
        writer.writeheader()
        for sample_id in EBV_IDS:
            for mate, size in ((1, 11), (2, 22)):
                writer.writerow(
                    {
                        "sample_id": sample_id,
                        "mate": mate,
                        "url": f"https://example.org/{sample_id}_{mate}.fastq.gz",
                        "md5": str(mate) * 32,
                        "bytes": size,
                    }
                )
                (ena_root / "files" / f"{sample_id}_{mate}.fastq.gz").write_bytes(b"x" * size)
                (ena_root / "status" / f"{sample_id}_{mate}.tsv").write_text(
                    "sample_id\tmate\tstatus\tbytes\tmd5\tsha256\n"
                    f"{sample_id}\t{mate}\tvalid\t{size}\t"
                    f"{str(mate) * 32}\t{str(mate) * 64}\n",
                    encoding="utf-8",
                )

    ebv_audits = tmp_path / "ebv-audits"
    for sample_id in EBV_IDS:
        _write_pair_audit(
            ebv_audits / f"{sample_id}.json",
            sample_id=sample_id,
            stored_bytes=(11, 22),
        )

    retained_manifest = tmp_path / "retained.tsv"
    retained_audits = tmp_path / "retained-audits"
    with retained_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sample_id", "chemistry", "read1", "read2"),
            delimiter="\t",
        )
        writer.writeheader()
        for sample_id in HIV_IDS:
            read1 = tmp_path / f"{sample_id}_1.fastq"
            read2 = tmp_path / f"{sample_id}_2.fastq"
            read1.write_bytes(b"x" * 33)
            read2.write_bytes(b"x" * 44)
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "chemistry": "10xv2",
                    "read1": read1,
                    "read2": read2,
                }
            )
            _write_pair_audit(
                retained_audits / f"{sample_id}.json",
                sample_id=sample_id,
                stored_bytes=(33, 44),
            )

    raw_output = tmp_path / "ignored" / "control_inputs.raw.tsv"
    tracked_output = tmp_path / "analysis" / "control_inputs.tsv"
    rows = freeze_control_inputs(
        ena_root=ena_root,
        ebv_audit_root=ebv_audits,
        retained_manifest=retained_manifest,
        retained_audit_root=retained_audits,
        raw_output=raw_output,
        tracked_output=tracked_output,
    )

    assert [row["sample_id"] for row in rows] == sorted((*EBV_IDS, *HIV_IDS))
    assert len(tracked_output.read_text(encoding="utf-8").splitlines()) == 6
    assert str(tmp_path) in raw_output.read_text(encoding="utf-8")
    assert str(tmp_path) not in tracked_output.read_text(encoding="utf-8")
