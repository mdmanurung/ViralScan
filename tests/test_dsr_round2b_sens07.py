"""DSR-16 round 2b / SENS-07: the planted-sensitivity array script's env block, config template and dry run."""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "slurm_dsr_round2b_sens07_array.sh"
TEMPLATE = REPO / "config" / "dsr_round2b_sens07.template.tsv"
M = "/exports/archive/hg-funcgenom-research/mdmanurung"


def _env_block(env: dict[str, str]) -> dict[str, str]:
    block = re.search(r"# >>> dsr-env.*?\n(.*?)# <<< dsr-env", SCRIPT.read_text(), re.S).group(1)
    names = ["ROUND2", "R", "CODE", "B", "CONFIG", "EXACT"]
    probe = block + "".join(f'echo "{n}=${n}"\n' for n in names)
    out = subprocess.run(
        ["bash", "-c", probe], env=env, capture_output=True, text=True, check=True
    ).stdout
    return dict(line.split("=", 1) for line in out.splitlines())


def test_defaults_are_round2b_on_the_q49_pin():
    got = _env_block({"PATH": "/usr/bin:/bin"})
    assert got["ROUND2"] == f"{M}/viralscan_work/dsr_round2"
    assert got["R"] == f"{M}/viralscan_work/dsr_round2b"
    assert got["CODE"] == f"{M}/vs_pinned/35de249"
    assert got["B"] == f"{M}/viral_ref_cat42d/build"
    assert got["CONFIG"] == f"{M}/viralscan_work/dsr_round2b/config.tsv"
    assert got["EXACT"] == f"{M}/viral_ref_cat42d/build/viral.fa"


def test_overrides_take_effect():
    got = _env_block(
        {"PATH": "/usr/bin:/bin", "ROUND_DIR": "/r", "CODE_DIR": "/c", "REF_BUILD": "/b/build"}
    )
    assert (got["R"], got["CODE"], got["B"]) == ("/r", "/c", "/b/build")
    assert got["CONFIG"] == "/r/config.tsv"
    assert got["EXACT"] == "/b/build/viral.fa"


def test_bash_syntax_and_shellcheck():
    subprocess.run(["bash", "-n", str(SCRIPT)], check=True)
    if shutil.which("shellcheck"):
        subprocess.run(["shellcheck", "-S", "error", str(SCRIPT)], check=True)


def test_config_template_shape():
    with open(TEMPLATE, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert list(rows[0]) == ["dataset", "sample", "run", "expected_family", "n_background", "note"]
    assert len(rows) == 17 and len({r["sample"] for r in rows}) == 17
    assert {r["run"] for r in rows} <= {"yes", "no"}
    assert all(r["note"] != "-" for r in rows if r["run"] == "no")


@pytest.fixture
def world(tmp_path):
    """A fake round 2 (two manifest rows + a skipped one), code dir, FASTA and identity table."""
    r2, code, b = tmp_path / "round2", tmp_path / "code", tmp_path / "build"
    for d in (r2 / "chem", code / "config", b):
        d.mkdir(parents=True)
    (code / "config" / "track_b_denylist.txt").write_text("GSE189536\n")
    head = "dataset\tsample\tgrp\tr1\tr2\twhitelist\tpairs\thint\n"
    rows = [
        ("hsv1", "S1", "/d/S1_1.fq.gz", "/d/S1_2.fq.gz", "-"),
        ("covid", "S2", "/d/S2_R1.fq.gz", "/d/S2_R2.fq.gz", "/wl.txt"),
        ("hhv6b", "S3", "/d/S3_1.fq.gz", "/d/S3_2.fq.gz", "-"),
        ("trackb", "S4", "/d/GSE189536_1.fq.gz", "/d/S4_2.fq.gz", "-"),
    ]
    (r2 / "manifest.tsv").write_text(
        head + "".join(f"{d}\t{s}\tS\t{a}\t{c}\t{w}\t1\t-\n" for d, s, a, c, w in rows)
    )
    for _, s, *_ in rows:
        (r2 / "chem" / f"{s}.json").write_text(
            '{"chemistry": {"chemistry": "10xv3"}, "strand_pilot": {"choice": "forward"}}'
        )
    (b / "viral.fa").write_text(">AB1.1 a\nACGT\n>AB2.1 b\nACGT\n>AB3.1 c\nACGT\n>AB4.1 d\nACGT\n")
    ident = ["gene_id\tgenome_accession\tstatus\tviral\tfamily"]
    for acc, fam in [
        ("AB1.1", "Herpesviridae"),
        ("AB2.1", "Papillomaviridae"),
        ("AB3.1", "Polyomaviridae"),
        ("AB4.1", "Coronaviridae"),
    ]:
        ident.append(f"g_{acc}\t{acc}\tcatalogued\ttrue\t{fam}")
    (r2 / "identity.tsv").write_text("\n".join(ident) + "\n")
    config = (
        "dataset\tsample\trun\texpected_family\tn_background\tnote\n"
        "hsv1\tS1\tyes\t-\t-\t-\n"
        "covid\tS2\tyes\tCoronaviridae\t1000\t-\n"
        "hhv6b\tS3\tno\t-\t-\tR1 too long\n"
        "trackb\tS4\tyes\t-\t-\t-\n"
    )
    (tmp_path / "config.tsv").write_text(config)
    return tmp_path


def _dry(world, rows, **extra):
    env = {
        "PATH": f"{Path(sys.executable).parent}:/usr/bin:/bin",
        "DRY_RUN": "1",
        "ROUND2_DIR": str(world / "round2"),
        "ROUND_DIR": str(world / "round2b"),
        "CODE_DIR": str(world / "code"),
        "REF_BUILD": str(world / "build"),
        "CONFIG_TSV": str(world / "config.tsv"),
        "FAMILY_IDENTITY": str(world / "round2" / "identity.tsv"),
        "PLANT_FAMILIES": "Herpesviridae Papillomaviridae Polyomaviridae",
        **extra,
    }
    return subprocess.run(
        ["bash", str(SCRIPT), *map(str, rows)], env=env, capture_output=True, text=True
    )


def test_dry_run_prints_commands_and_creates_nothing(world):
    done = _dry(world, [0, 1, 2])
    out = done.stdout
    assert done.returncode == 0, done.stderr
    assert not (world / "round2b").exists()
    # row 0: default families, exact plants only, Q48 parameters
    assert "families: 1 Herpesviridae 1 Papillomaviridae 1 Polyomaviridae" in out
    assert "--per-family 2 --seed 20261012 --n-background 50000000" in out
    assert "--heldout-fasta" not in out
    # row 1: its own expected family is added, per-row n_background, whitelist, same arm flags as combined_off
    assert "1 Coronaviridae" in out and "--n-background 1000" in out
    assert "-w /wl.txt" in out and "--cell-calling emptydrops --anellovirus-gene-ids --yes" in out
    assert "--read-filter" not in out and "--host-filter" not in out
    # planted + twin each go through viralscan and the ledger; then score
    assert out.count("viralscan.menu") == 4 and "--arm combined_off_planted" in out
    assert "--arm combined_off_unplanted" in out and "sens07_planted_sensitivity.py score" in out
    # row 2: config run=no
    assert "row 2 hhv6b/S3 skipped (config run=no): R1 too long" in out


def test_parameters_override_defaults(world):
    out = _dry(
        world,
        [0],
        PER_FAMILY="3",
        SEED="7",
        N_BACKGROUND_DEFAULT="99",
        PLANT_EXTRA="--levels 10 100",
    ).stdout
    assert "--per-family 3 --seed 7 --n-background 99 --levels 10 100" in out


def test_unknown_family_and_denylisted_row_are_refused(world):
    bad = _dry(world, [0], PLANT_FAMILIES="Nonexistentviridae")
    assert bad.returncode != 0 and "Nonexistentviridae" in bad.stderr
    deny = _dry(world, [3])
    assert deny.returncode == 4 and "deny-listed" in deny.stderr


def test_script_is_executable():
    assert os.access(SCRIPT, os.X_OK)
