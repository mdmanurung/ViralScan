"""Raw output and failure parity against the frozen read filter."""

import gzip
import importlib.util
from pathlib import Path

import pytest

from viralscan import anello_align as aa
from viralscan.scripts import read_filter as candidate

spec = importlib.util.spec_from_file_location(
    "original_filter", Path(__file__).parent / "_frozen_oracles/frozen_read_filter_d056986.py"
)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)
CB = "ACGTACGTACGTACGTAAACCCGGGTTT"
BODY = "ACGT" * 6


def _run(fn, root, r1, r2, technology, mode):
    root.mkdir(exist_ok=True)
    p1, p2 = root.parent / "input1.fastq", root.parent / "input2.fastq"
    p1.write_text(r1)
    p2.write_text(r2)
    try:
        counts = fn(str(p1), str(p2), root, technology, mode)
        error = None
    except Exception as exc:
        counts, error = None, (type(exc), str(exc))
    return counts, error, {p.name: p.read_bytes() for p in root.iterdir()}


def _records(seqs, mate):
    return "".join(f"@read{i}/{mate}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(seqs))


@pytest.mark.parametrize("technology", ["10xv2", "10xv3", "10xv4", "dropseq"])
@pytest.mark.parametrize("mode", ["artefact", "tso-trim"])
@pytest.mark.parametrize("failure", [None, "mismatch", "truncated"])
def test_raw_files_counts_and_late_failure_match(tmp_path, monkeypatch, technology, mode, failure):
    monkeypatch.setattr(gzip.time, "time", lambda: 1234567890)
    seqs = [BODY, "ACgtACGT" * 3, "G" * 90, aa.TSO + "T" * 65, BODY * 25] * 4
    r1 = _records([CB + ("G" * i) for i in range(len(seqs))], 1)
    r2 = _records(seqs, 2)
    if failure == "mismatch":
        r2 = r2.replace(f"@read{len(seqs) - 1}/2", "@different/2")
    elif failure == "truncated":
        r2 = r2[:-3].rsplit("\n", 1)[0]
    # Identical absolute output paths also preserve gzip filename headers.
    out = tmp_path / "output"
    expected = _run(original.filter_pairs, out, r1, r2, technology, mode)
    for p in out.iterdir():
        p.unlink()
    assert expected == _run(original.filter_pairs, out, r1, r2, technology, mode)
    for p in out.iterdir():
        p.unlink()
    assert expected == _run(candidate.filter_pairs, out, r1, r2, technology, mode)


def test_cache_uses_raw_case_prefix_bounded_keys_and_invocation_lifetime(tmp_path, monkeypatch):
    calls = []
    real = candidate.classify

    def classify(r1, r2, length):
        calls.append((r1, r2, length))
        return real(r1, r2, length)

    monkeypatch.setattr(candidate, "classify", classify)
    assert real(CB, BODY, 28) != real(CB, "ACgtACGT" * 3, 28)
    seqs = [BODY, BODY, "ACgtACGT" * 3, BODY * 25, BODY * 25]
    # A suffix outside CB+UMI does not affect classification.
    r1 = _records([CB + "A" * i for i in range(5)], 1)
    r2 = _records(seqs, 2)
    _run(candidate.filter_pairs, tmp_path / "out", r1, r2, "10xv3", "artefact")
    assert len(calls) == 4
    assert all(a == CB for a, _, _ in calls)
    calls.clear()
    _run(candidate.filter_pairs, tmp_path / "again", r1, r2, "10xv3", "artefact")
    assert len(calls) == 4


def test_entry_bound_evicts_old_keys(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(candidate, "classify", lambda a, b, n: calls.append(b) or "")
    seqs = ["ACGT" + str(i) for i in range(257)]
    seqs += [seqs[-1], seqs[0]]
    _run(
        candidate.filter_pairs,
        tmp_path / "out",
        _records([CB] * len(seqs), 1),
        _records(seqs, 2),
        "10xv3",
        "artefact",
    )
    assert len(calls) == 258
