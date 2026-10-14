import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_fastq_pair import FastqAuditError, _open_fastq, _record  # noqa: E402
from subset_fastq_pair import subset  # noqa: E402


def inputs(tmp_path):
    paths = [tmp_path / "source_R1.fastq", tmp_path / "source_R2.fastq"]
    for mate, path in enumerate(paths, 1):
        path.write_text("".join(f"@read{i}/{mate}\nACGT\n+\nIIII\n" for i in range(3)))
    return paths


def test_prefix_is_paired_and_full_records_only(tmp_path):
    read1, read2 = inputs(tmp_path)
    output = tmp_path / "subset"
    receipt = subset(read1, read2, output, 2)
    assert receipt["audit"]["records"] == 2
    with _open_fastq(output / "R1.fastq.gz") as handle:
        assert _record(handle, read1, 1)[0] == b"@read0/1\n"
        assert _record(handle, read1, 2)[0] == b"@read1/1\n"
        assert _record(handle, read1, 3) is None
    with pytest.raises(ValueError, match="new directory"):
        subset(read1, read2, output, 1)


@pytest.mark.parametrize("failure", ["mismatch", "truncated", "short"])
def test_bad_pair_leaves_no_published_subset(tmp_path, failure):
    read1, read2 = inputs(tmp_path)
    if failure == "mismatch":
        read2.write_text(read2.read_text().replace("read1/2", "other/2"))
    elif failure == "truncated":
        read2.write_text(read2.read_text()[:-6])
    output = tmp_path / "subset"
    with pytest.raises(FastqAuditError):
        subset(read1, read2, output, 4 if failure == "short" else 3)
    assert not output.exists()
