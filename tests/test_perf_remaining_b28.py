"""Single immutable GenBank parse retains grouped GTF bytes and errors."""

import importlib.util
import warnings
from pathlib import Path

import pytest

from viralscan.scripts import ncbi_fetch as candidate

spec = importlib.util.spec_from_file_location(
    "original_genbank", Path(__file__).parent / "_frozen_oracles/frozen_genbank_gtf_d056986.py"
)
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)


def record(features):
    return (
        "LOCUS       NC_TEST 5000 bp DNA circular\nVERSION     NC_TEST.1\nFEATURES             Location/Qualifiers\n"
        + features
        + "ORIGIN\n        1 acgt\n//\n"
    )


def capture(fn, text):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result, error = fn(text, "NC_TEST.1").encode(), None
        except Exception as exc:
            result, error = None, (type(exc), str(exc))
    return result, error, [(w.category, str(w.message)) for w in caught]


@pytest.mark.parametrize(
    "location",
    [
        "1..300",
        "complement(join(4000..5000,1..30))",
        "join(<1..20,40..>70)",
        "bogus",
        "join(1..30,bogus)",
        "1^2",
        "",
    ],
)
@pytest.mark.parametrize("symbols", [True, False])
def test_group_order_tokens_intervals_warnings_errors(location, symbols):
    parts = []
    for key in ["misc_RNA", "CDS", "source", "misc_RNA", "CDS"]:
        parts.append(f"     {key:<16}{location}\n")
        if symbols:
            parts.append('                     /gene="same/name"\n')
    text = record("".join(parts))
    expected = capture(original._genbank_to_gtf, text)
    assert expected == capture(original._genbank_to_gtf, text)
    assert expected == capture(candidate._genbank_to_gtf, text)


def test_no_features_and_single_traversal(monkeypatch):
    assert capture(original._genbank_to_gtf, record("")) == capture(
        candidate._genbank_to_gtf, record("")
    )
    real, calls = candidate.iter_features, []

    def observed(text):
        calls.append(text)
        return real(text)

    monkeypatch.setattr(candidate, "iter_features", observed)
    candidate._genbank_to_gtf(record("     CDS             1..300\n"), "NC_TEST.1")
    assert len(calls) == 1
