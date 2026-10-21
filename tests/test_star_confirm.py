"""STAR-04: ``viralscan run --star-confirm`` flags, config, rule script and summary columns."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
import yaml

from viralscan import star_hostvirus
from viralscan.cli.parsers import build_parser
from viralscan.menu import _write_run_config, errorhandler
from viralscan.runconfig import RunConfig
from viralscan.scripts import detection, star_confirm

from .test_auto_evidence import _cfg

STAR_COLS = list(star_hostvirus.STAR_SUMMARY_COLUMNS[1:])


def _args(tmp_path: Path, **kwargs) -> argparse.Namespace:
    (tmp_path / "R1.fastq.gz").touch()
    (tmp_path / "R2.fastq.gz").touch()
    defaults = {
        "ncbi_accession": None,
        "reference": False,
        "gtf": None,
        "fasta": None,
        "index": str(tmp_path / "ref" / "index.idx"),
        "transcripts": "/fake/t2g.txt",
        "f1": None,
        "host_filter": "starsolo",
        "host_index": "/fake/host",
        "sample1": str(tmp_path / "R1.fastq.gz"),
        "sample2": str(tmp_path / "R2.fastq.gz"),
        "star_confirm": True,
        "star_hostvirus_index": "/fake/star_idx",
    }
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def _check(args, problems=()):
    with (
        patch("os.path.exists", return_value=True),
        patch.object(star_hostvirus, "validate_index_manifest", return_value=list(problems)),
    ):
        errorhandler(args)


def test_parser_defaults_and_flags():
    ns = build_parser().parse_args(["run", "--r1", "a", "--r2", "b", "-o", "o"])
    assert ns.star_confirm is False and ns.star_hostvirus_index is None
    ns = build_parser().parse_args(
        ["run", "--r1", "a", "--r2", "b", "-o", "o", "--star-confirm"] + ["--star-hostvirus-index", "/x"]
    )
    assert ns.star_confirm is True and ns.star_hostvirus_index == "/x"


def test_valid_pairing_passes(tmp_path):
    _check(_args(tmp_path))


@pytest.mark.parametrize(
    "bad",
    [
        {"star_hostvirus_index": None},
        {"host_filter": None, "host_index": None},
        {"star_confirm": False},  # index without the flag
    ],
)
def test_pairing_errors(tmp_path, bad):
    with pytest.raises(SystemExit) as exc:
        _check(_args(tmp_path, **bad))
    assert exc.value.code != 0


def test_index_without_manifest_dies(tmp_path, caplog):
    with pytest.raises(SystemExit):
        _check(_args(tmp_path), problems=["star_hostvirus_index.json missing"])
    assert "star_hostvirus_index.json" in caplog.text


def test_runconfig_round_trip(tmp_path):
    off = RunConfig.from_snakemake_config(_cfg(tmp_path))
    assert off.star_confirm is False and off.star_hostvirus_index is None
    rc = RunConfig.from_snakemake_config(
        _cfg(tmp_path, star_confirm="true", star_hostvirus_index="/idx")
    )
    path = _write_run_config(rc)
    assert yaml.safe_load(path.read_text())["star_confirm"] is True
    back = RunConfig.from_yaml(path)
    assert back.star_confirm is True and back.star_hostvirus_index == "/idx"
    with pytest.raises(ValueError, match="star_hostvirus_index"):
        RunConfig.from_snakemake_config(_cfg(tmp_path, star_confirm="true"))


def _identity(path: Path) -> Path:
    path.write_text("genome_accession\tvirus_name\nNC_1.1\tVirA\nNC_2.1\tVirA\nNC_3.1\tVirB\n")
    return path


def _summary(path: Path) -> Path:
    star_hostvirus.write_summary(
        [
            dict(accession="NC_1.1", star_reads_primary=3, star_reads_any=4, star_molecules=2,
                 star_host_competing=1, star_other_virus=0),
            dict(accession="NC_2.1", star_reads_primary=1, star_reads_any=1, star_molecules=1,
                 star_host_competing=0, star_other_virus=1),
        ],
        path,
    )  # fmt: skip
    return path


def test_by_virus_sums_accessions(tmp_path):
    rows = star_confirm.by_virus(
        _summary(tmp_path / "s.tsv"), star_confirm._acc_to_virus(_identity(tmp_path / "i.tsv"))
    )
    assert rows == {
        "VirA": {
            "star_reads_primary": 4,
            "star_reads_any": 5,
            "star_molecules": 3,
            "star_host_competing": 1,
            "star_other_virus": 1,
        }
    }


def test_script_main_writes_tsv_and_renames_bam(tmp_path, monkeypatch):
    out = tmp_path / "out"
    (out / "results").mkdir(parents=True)
    (out / "log").mkdir()
    ident = _identity(out / "results" / "virus_identity.tsv")
    seen = {}

    def fake_align(r1, r2, index, work, technology, **kw):
        seen.update(kw, index=index, technology=technology)
        work.mkdir(parents=True, exist_ok=True)
        _summary(work / "star_support_by_virus.tsv")
        (work / "Aligned.viral.bam").write_text("b")
        (work / "Aligned.viral.bam.bai").write_text("i")
        return work / "Aligned.viral.bam"

    monkeypatch.setattr(star_confirm, "align_host_unmapped", fake_align)
    cfg = RunConfig.from_snakemake_config(
        _cfg(tmp_path, star_confirm="true", star_hostvirus_index="/idx", strand="reverse")
        | {"output": f"{out}/"}
    )
    star_confirm.main(cfg, str(ident), 4, str(out / "log" / "star_confirm.done"))
    assert seen["strand"] == "reverse" and seen["prune"] is True and seen["threads"] == 4
    assert (out / "star_hostvirus" / "viral.bam").is_file()
    assert (out / "star_hostvirus" / "viral.bam.bai").is_file()
    assert (out / "log" / "star_confirm.done").is_file()
    with open(out / "results" / "star_support_by_virus.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows == [
        {"virus_name": "VirA", **{c: str(v) for c, v in zip(STAR_COLS, (4, 5, 3, 1, 1))}}
    ]


def _stats():
    return {
        "viral_molecules_total_est": 5.0,
        "infected_cells": 2,
        "total_cells": 100,
        "pct_infected": 2.0,
        "viral_molecules_per_10k_est": 1.0,
        "n_called_cells": 50,
        "infected_called": 1,
        "pct_infected_called": 2.0,
        "infected_comparable": 1,
        "n_comparable_cells": 40,
        "pct_infected_comparable": 2.5,
        "accession_breadth": 1.0,
        "host_viral_ambig_fraction": 0.0,
    }


class _Cfg:
    def __init__(self, star):
        self.star_confirm = star


def _summary_table(run: Path) -> pd.DataFrame:
    return pd.read_csv(
        run / "results" / "viral_summary.tsv", sep="\t", dtype=str, keep_default_na=False
    )


def test_detection_merge_adds_columns_only_when_run(tmp_path):
    stats = {"VirA": _stats(), "VirB": _stats()}
    detection.write_tsv_outputs(stats, pd.DataFrame(), str(tmp_path / "off"))
    off = _summary_table(tmp_path / "off")
    assert not set(STAR_COLS) & set(off.columns)

    (tmp_path / "on" / "results").mkdir(parents=True)
    (tmp_path / "on" / "results" / "star_support_by_virus.tsv").write_text(
        "virus_name\t" + "\t".join(STAR_COLS) + "\nVirA\t4\t5\t3\t1\t1\n"
    )
    star = detection.star_evidence(_Cfg(True), str(tmp_path / "on"))
    detection.write_tsv_outputs(stats, pd.DataFrame(), str(tmp_path / "on"), star=star)
    on = _summary_table(tmp_path / "on").set_index("virus_name")
    assert list(on.columns)[: len(off.columns) - 1] == list(off.columns)[1:]  # old order kept
    assert list(on.columns)[-len(STAR_COLS) :] == STAR_COLS
    assert on.loc["VirA", STAR_COLS].tolist() == ["4", "5", "3", "1", "1"]
    assert on.loc["VirB", STAR_COLS].tolist() == ["0"] * 5  # ran, saw nothing


def test_star_evidence_none_when_off_and_missing_file_raises(tmp_path):
    assert detection.star_evidence(_Cfg(False), str(tmp_path)) is None
    with pytest.raises(FileNotFoundError, match="star_confirm"):
        detection.star_evidence(_Cfg(True), str(tmp_path))
