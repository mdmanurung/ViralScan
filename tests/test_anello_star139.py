"""STAR exit 139 (Solo crash with no whitelist and no mapped reads) is recovered, nothing else is."""

import subprocess

import pytest

from viralscan.scripts import anello_align

FINAL_ZERO = "Uniquely mapped reads number |\t0\nNumber of reads mapped to multiple loci |\t0\n"
FINAL_MAPPED = "Uniquely mapped reads number |\t3\nNumber of reads mapped to multiple loci |\t0\n"


def _cmd(whitelist):
    return [
        "STAR",
        "--runThreadN",
        "1",
        "--soloCBwhitelist",
        whitelist,
        "--outFileNamePrefix",
        "x/",
    ]


def _fake_run(monkeypatch, work, codes, final):
    calls = []

    def run(cmd, check=True):
        calls.append(list(cmd))
        code = codes[len(calls) - 1]
        if code:
            raise subprocess.CalledProcessError(code, cmd)
        (work / "Log.final.out").write_text(final)

    monkeypatch.setattr(anello_align.subprocess, "run", run)
    return calls


def test_zero_mapped_crash_is_retried_with_dummy_whitelist(tmp_path, monkeypatch):
    calls = _fake_run(monkeypatch, tmp_path, [139, 0], FINAL_ZERO)
    anello_align._run_star(_cmd("None"), tmp_path, 16)
    assert len(calls) == 2
    assert calls[1][calls[1].index("--soloCBwhitelist") + 1] == str(
        tmp_path / "dummy_whitelist.txt"
    )
    # the dummy must be as long as --soloCBlen, else STAR rejects the whitelist
    assert (tmp_path / "dummy_whitelist.txt").read_text() == "A" * 16 + "\n"
    assert (tmp_path / "star139_fallback.txt").is_file()


def test_retry_that_maps_reads_is_refused(tmp_path, monkeypatch):
    _fake_run(monkeypatch, tmp_path, [139, 0], FINAL_MAPPED)
    with pytest.raises(RuntimeError, match="refusing to drop barcodes"):
        anello_align._run_star(_cmd("None"), tmp_path, 12)


@pytest.mark.parametrize("code,whitelist", [(139, "wl.txt"), (1, "None")])
def test_other_failures_are_not_swallowed(tmp_path, monkeypatch, code, whitelist):
    calls = _fake_run(monkeypatch, tmp_path, [code], FINAL_ZERO)
    with pytest.raises(subprocess.CalledProcessError):
        anello_align._run_star(_cmd(whitelist), tmp_path, 12)
    assert len(calls) == 1


def test_main_passes_technology_barcode_length(tmp_path, monkeypatch):
    # explicit 16-base geometry: no on-list, so --soloCBlen 16 and the dummy must be 16 bases
    seen = {}
    monkeypatch.setattr(
        anello_align, "_run_star", lambda cmd, work, cb_len: seen.update(cb_len=cb_len) or 1 / 0
    )
    cfg = anello_align.RunConfig.__new__(anello_align.RunConfig)
    for k, v in dict(
        output=str(tmp_path),
        anello_index="idx",
        whitelist=None,
        technology="0,0,16:0,16,28:1,0,0",
        kb_r1="r1",
        kb_r2="r2",
    ).items():
        object.__setattr__(cfg, k, v)
    with pytest.raises(ZeroDivisionError):
        anello_align.main(cfg, "id.tsv", 1, "done")
    assert seen["cb_len"] == 16
