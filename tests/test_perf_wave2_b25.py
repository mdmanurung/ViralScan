"""B25 strict oracle: streamed reference-input concatenation == frozen whole-read baseline.

The frozen pre-B25 module is loaded by path (read-only, a verbatim copy of
``src/viralscan/menu.py`` in the wave-2 frozen baseline, commit ``faf5667``) and the
two implementations are compared byte for byte over the frozen edge cases the card
names: CR/LF/CRLF runs, trailing CR/LF runs of every length (including runs that span
a chunk boundary), empty and blank-only inputs, multibyte UTF-8, invalid UTF-8 bytes,
and multi-file combinations. Failure and partial-output behaviour is compared too.

The last test pins the memory contract itself: the streamed form must stay bounded
while the frozen whole-read form grows with the input, so the guard fails on the
pre-change code rather than merely describing it.
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
import tracemalloc
from pathlib import Path

import pytest

from viralscan import menu as cand

BASELINE = Path(__file__).resolve().parent / "_frozen_oracles" / "frozen_menu_faf5667.py"
assert BASELINE.is_file(), BASELINE

CHUNK = 1 << 20


def _load_baseline():
    spec = importlib.util.spec_from_file_location("viralscan_b25_frozen_menu", BASELINE)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


base = _load_baseline()
assert Path(cand.__file__).resolve() != BASELINE


EDGE_CASES = {
    "lf": b"acgt\nacgt\n",
    "crlf": b"acgt\r\nacgt\r\ncr-lf-tails\n",
    "cr-only": b"acgt\racgt\r",
    "trailing-run": b"acgt\n\r\n\r",
    "trailing-run-mixed": b"acgt\r\r\n\n\r\n",
    "only-newlines": b"\n\n\r\n",
    "only-crlf-and-cr": b"\r\n\r\n\r",
    "empty": b"",
    "no-final-newline": b"acgt",
    "blank-lines-mid": b"acgt\n\n\nacgt\n",
    "leading-blank": b"\nacgt\n",
    "lowercase": b"acgt\ntgca\n",
    "descriptions": b">chr1 human chromosome 1\nacgt\n",
    "utf8-multibyte": ">chr\u00e9 \u00fc\nacgt\u00e9\n".encode(),
    "invalid-utf8": b">chr1\nac\xffgt\n",
    "gtf-shaped": b'chr1\tsrc\tgene\t1\t9\t.\t+\t.\tgene_id "g1";\n',
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def edge_files(tmp_path: Path) -> dict[str, Path]:
    files = {}
    for name, data in EDGE_CASES.items():
        p = tmp_path / f"{name}.fa"
        p.write_bytes(data)
        files[name] = p
    # trailing run that spans the chunk boundary: 1 MiB of content then a long run
    big = tmp_path / "chunk-boundary-run.fa"
    big.write_bytes(b"acgt" * (1 << 18) + b"\r\n" * (1 << 19))
    files["chunk-boundary-run"] = big
    return files


COMBOS = [
    ["lf"],
    ["crlf"],
    ["cr-only"],
    ["trailing-run"],
    ["trailing-run-mixed"],
    ["only-newlines"],
    ["only-crlf-and-cr"],
    ["empty"],
    ["no-final-newline"],
    ["blank-lines-mid"],
    ["leading-blank"],
    ["utf8-multibyte"],
    ["invalid-utf8"],
    ["gtf-shaped"],
    ["chunk-boundary-run"],
    ["lf", "crlf"],
    ["empty", "lf"],
    ["empty", "empty"],
    ["only-newlines", "only-newlines"],
    ["trailing-run", "lf"],
    ["chunk-boundary-run", "lf"],
    ["lf", "crlf", "cr-only", "trailing-run", "empty", "utf8-multibyte", "gtf-shaped"],
    list(EDGE_CASES) + ["chunk-boundary-run"],
]


@pytest.mark.parametrize("combo", COMBOS, ids=lambda c: "+".join(c)[:60])
def test_byte_identical_to_frozen_baseline(edge_files, tmp_path, combo):
    paths = [str(edge_files[name]) for name in combo]
    got = cand._materialize_reference_input(paths, tmp_path / "cand" / "input.fa")
    want = base._materialize_reference_input(paths, tmp_path / "base" / "input.fa")
    assert Path(got).read_bytes() == Path(want).read_bytes()
    assert _sha(Path(got)) == _sha(Path(want))


def test_single_file_agreement(edge_files, tmp_path):
    """A one-file list still goes through the same contract."""
    for name in EDGE_CASES:
        paths = [str(edge_files[name])]
        got = cand._materialize_reference_input(paths, tmp_path / f"c-{name}.fa")
        want = base._materialize_reference_input(paths, tmp_path / f"b-{name}.fa")
        assert Path(got).read_bytes() == Path(want).read_bytes(), name


def test_trailing_run_lengths(edge_files, tmp_path):
    """Every trailing-run length 0..8, alone and after content, matches rstrip."""
    for content in (b"", b"acgt", b"acgt\n"):
        for run_len in range(9):
            for run in (b"\n" * run_len, b"\r" * run_len, (b"\r\n" * run_len)):
                src = tmp_path / "run.fa"
                src.write_bytes(content + run)
                paths = [str(src)]
                got = cand._materialize_reference_input(paths, tmp_path / "c.fa")
                want = base._materialize_reference_input(paths, tmp_path / "b.fa")
                assert Path(got).read_bytes() == Path(want).read_bytes(), (content, run)


def test_chunk_boundary_split_of_trailing_run(tmp_path, monkeypatch):
    """The run must be discarded even when chunks split it in every position."""
    content = b"acgt" * 1024
    run = b"\r\n" * 64
    src = tmp_path / "split.fa"
    src.write_bytes(content + run)
    want = base._materialize_reference_input([str(src)], tmp_path / "base.fa")
    for chunk in (3, 7, 64, 257, 4096):
        monkeypatch.setattr(cand, "_REFERENCE_COPY_CHUNK", chunk)
        got = cand._materialize_reference_input([str(src)], tmp_path / f"cand{chunk}.fa")
        assert Path(got).read_bytes() == Path(want).read_bytes(), chunk
    assert Path(want).read_bytes() == content + b"\n"


def test_internal_newline_at_chunk_boundary_is_kept(tmp_path, monkeypatch):
    """Only the *file's* trailing run may be dropped: an internal CR/LF that lands
    on a chunk boundary must survive, which a per-chunk rstrip would lose."""
    src = tmp_path / "internal.fa"
    for chunk in (3, 8, 64):
        data = b"".join(b"B" * chunk + b"\n" for _ in range(12))
        src.write_bytes(data)
        monkeypatch.setattr(cand, "_REFERENCE_COPY_CHUNK", chunk)
        want = base._materialize_reference_input([str(src)], tmp_path / "b.fa")
        got = cand._materialize_reference_input([str(src)], tmp_path / f"c{chunk}.fa")
        assert Path(got).read_bytes() == Path(want).read_bytes()
        # 12 newlines total, the last of which is the trailing run -> 11 kept + 1 added
        assert Path(want).read_bytes().count(b"\n") == 12


def test_missing_second_input_keeps_partial_output_and_error(tmp_path):
    """A later missing input must fail identically, leaving the same partial file."""
    first = tmp_path / "first.fa"
    first.write_bytes(b"acgt\nacgt\n")
    missing = tmp_path / "nope.fa"

    cand_out = tmp_path / "cand.fa"
    base_out = tmp_path / "base.fa"
    with pytest.raises(FileNotFoundError) as cand_exc:
        cand._materialize_reference_input([str(first), str(missing)], cand_out)
    with pytest.raises(FileNotFoundError) as base_exc:
        base._materialize_reference_input([str(first), str(missing)], base_out)
    assert str(cand_exc.value) == str(base_exc.value)
    assert type(cand_exc.value) is type(base_exc.value)
    assert cand_out.read_bytes() == base_out.read_bytes() == b"acgt\nacgt\n"


def test_first_input_missing_creates_empty_destination(tmp_path):
    missing = tmp_path / "nope.fa"
    cand_out = tmp_path / "cand.fa"
    base_out = tmp_path / "base.fa"
    with pytest.raises(FileNotFoundError) as cand_exc:
        cand._materialize_reference_input([str(missing)], cand_out)
    with pytest.raises(FileNotFoundError) as base_exc:
        base._materialize_reference_input([str(missing)], base_out)
    assert str(cand_exc.value) == str(base_exc.value)
    assert cand_out.exists() and base_out.exists()
    assert cand_out.read_bytes() == base_out.read_bytes() == b""


def test_directory_input_error_matches(tmp_path):
    directory = tmp_path / "a-directory"
    directory.mkdir()
    with pytest.raises(IsADirectoryError) as cand_exc:
        cand._materialize_reference_input([str(directory)], tmp_path / "c.fa")
    with pytest.raises(IsADirectoryError) as base_exc:
        base._materialize_reference_input([str(directory)], tmp_path / "b.fa")
    assert str(cand_exc.value) == str(base_exc.value)


def test_nested_destination_parent_is_created(tmp_path):
    dest = tmp_path / "deep" / "nested" / "input.fa"
    cand._materialize_reference_input([], dest)
    assert dest.is_file() and dest.read_bytes() == b""
    assert dest.parent.is_dir()


def test_memory_stays_bounded_and_baseline_grows(tmp_path):
    """The point of the change: bounded candidate peak, unbounded frozen baseline."""
    src = tmp_path / "big.fa"
    with open(src, "wb") as f:
        f.write(b"ACGT" * 512 + b"\n" * 0)
        remaining = 16 * (1 << 20)
        block = b"ACGT" * 512 + b"\n"
        while remaining > 0:
            f.write(block)
            remaining -= len(block)
        f.write(b"\r\n\r\n")

    def measure(fn):
        import gc

        gc.collect()
        tracemalloc.start()
        try:
            fn()
            return tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()

    size = src.stat().st_size
    cand_peak = measure(lambda: cand._materialize_reference_input([str(src)], tmp_path / "c.fa"))
    base_peak = measure(lambda: base._materialize_reference_input([str(src)], tmp_path / "b.fa"))

    assert Path(tmp_path / "c.fa").read_bytes() == Path(tmp_path / "b.fa").read_bytes()
    assert size > 16 * CHUNK
    # bounded: a few chunks regardless of input size
    assert cand_peak < 5 * CHUNK, (cand_peak, size)
    # the pre-change code holds whole-file copies, so the same guard fails there
    assert base_peak > 2 * size, (base_peak, size)
    assert cand_peak < base_peak / 4
