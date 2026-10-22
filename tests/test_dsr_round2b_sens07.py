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
M = "/site/root"


def _env_block(env: dict[str, str]) -> dict[str, str]:
    block = re.search(r"# >>> dsr-env.*?\n(.*?)# <<< dsr-env", SCRIPT.read_text(), re.S).group(1)
    names = ["ROUND2", "R", "CODE", "B", "CONFIG", "EXACT"]
    probe = block + "".join(f'echo "{n}=${n}"\n' for n in names)
    out = subprocess.run(
        ["bash", "-c", probe], env=env, capture_output=True, text=True, check=True
    ).stdout
    return dict(line.split("=", 1) for line in out.splitlines())


def test_defaults_are_round2b_on_the_q49_pin():
    got = _env_block({"PATH": "/usr/bin:/bin", "VIRALSCAN_SITE_ROOT": M})
    assert got["ROUND2"] == f"{M}/viralscan_work/dsr_round2"
    assert got["R"] == f"{M}/viralscan_work/dsr_round2b"
    assert got["CODE"] == f"{M}/vs_pinned/35de249"
    assert got["B"] == f"{M}/viral_ref_cat42d/build"
    assert got["CONFIG"] == f"{M}/viralscan_work/dsr_round2b/config.tsv"
    assert got["EXACT"] == f"{M}/viral_ref_cat42d/build/viral.fa"


def test_overrides_take_effect():
    got = _env_block(
        {
            "PATH": "/usr/bin:/bin",
            "VIRALSCAN_SITE_ROOT": M,
            "ROUND_DIR": "/r",
            "CODE_DIR": "/c",
            "REF_BUILD": "/b/build",
        }
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
        "VIRALSCAN_SITE_ROOT": M,
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
    # plants come from the index's sense transcripts, in the library's orientation (forward here)
    assert f"--transcripts-fasta {world / 'build' / 'cdna.fa'} --library-strand forward" in out
    # row 1: its own expected family is added, per-row n_background, whitelist, same arm flags as combined_off
    assert "1 Coronaviridae" in out and "--n-background 1000" in out
    assert "-w /wl.txt" in out and "--cell-calling emptydrops --anellovirus-gene-ids --yes" in out
    assert "--read-filter" not in out and "--host-filter" not in out
    # planted + twin each go through viralscan and the ledger; then score
    assert out.count("viralscan.menu") == 4 and "--arm combined_off_planted" in out
    assert "--arm combined_off_unplanted" in out and "sens07_planted_sensitivity.py score" in out
    # row 2: config run=no
    assert "row 2 hhv6b/S3 skipped (config run=no): R1 too long" in out


def test_reverse_stranded_row_plants_antisense_reads(world):
    (world / "round2" / "chem" / "S1.json").write_text(
        '{"chemistry": {"chemistry": "10xv3"}, "strand_pilot": {"choice": "reverse"}}'
    )
    out = _dry(world, [0]).stdout
    assert "--library-strand reverse" in out and "--strand reverse" in out


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


# ── round 2c: per-family seeds (Q64) and trimmed HHV-6B (Q63), the same script ────────────────────
TEMPLATE_2C = REPO / "config" / "dsr_round2c_sens07.template.tsv"
GUARD_2C = REPO / "docs" / "plans" / "dsr_round2c_sens07_GUARD.md"


def _rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_round2c_template_is_round2b_plus_trimmed_hhv6b():
    old, new = _rows(TEMPLATE), _rows(TEMPLATE_2C)
    assert list(new[0]) == [
        "dataset", "sample", "run", "expected_family", "n_background", "trim_r1", "note",
    ]  # fmt: skip
    assert [r["sample"] for r in new] == [r["sample"] for r in old]
    assert {r["sample"] for r in new if r["trim_r1"] == "yes"} == {"SRR20710641"}
    assert {r["trim_r1"] for r in new} <= {"yes", "no"}
    for a, b in zip(old, new, strict=True):
        assert all(a[k] == b[k] for k in ("dataset", "expected_family", "n_background"))
        assert a["run"] == b["run"] or b["sample"] == "SRR20710641"
    by = {r["sample"]: r for r in new}
    assert by["SRR20710641"]["run"] == "yes" and "Q63" in by["SRR20710641"]["note"]
    assert by["GSM4682311"]["run"] == "no" and "Q63" in by["GSM4682311"]["note"]
    assert all(r["note"] != "-" for r in new if r["run"] == "no")


def test_round2c_guard_names_the_changes_and_environment():
    text = GUARD_2C.read_text()
    for needle in (
        "SELECTION_SEED_MODE=per-family",
        "--selection-seed-mode per-family",
        "--trim-r1-to-geometry",
        "SRR20710641",
        "GSM4682311",
        "dsr_round2c",
    ):
        assert needle in text, needle


def _config_2c(world: Path, trim_for_s3: str = "yes") -> Path:
    path = world / "config2c.tsv"
    path.write_text(
        "dataset\tsample\trun\texpected_family\tn_background\ttrim_r1\tnote\n"
        "hsv1\tS1\tyes\t-\t-\tno\t-\n"
        "covid\tS2\tyes\tCoronaviridae\t1000\tno\t-\n"
        f"hhv6b\tS3\tyes\t-\t-\t{trim_for_s3}\ttrim R1\n"
        "trackb\tS4\tyes\t-\t-\tno\t-\n"
    )
    return path


def test_round2c_dry_run_adds_the_flags_only_where_asked(world):
    done = _dry(
        world,
        [0, 2],
        CONFIG_TSV=str(_config_2c(world)),
        SELECTION_SEED_MODE="per-family",
        PLANT_EXTRA="--levels 10 100",
    )
    assert done.returncode == 0, done.stderr
    plants = [ln for ln in done.stdout.splitlines() if "sens07_planted_sensitivity.py plant" in ln]
    assert len(plants) == 2 and not (world / "round2b").exists()
    assert all("--selection-seed-mode per-family" in ln for ln in plants)
    assert "--trim-r1-to-geometry" not in plants[0] and "--trim-r1-to-geometry" in plants[1]
    # after PLANT_EXTRA, so round 2b's `--n-background N --levels ...` order is untouched
    assert "--n-background 50000000 --levels 10 100 --selection-seed-mode per-family" in plants[0]
    assert "row 2 hhv6b/S3" in done.stdout and "selection=per-family trim_r1=yes" in done.stdout


def test_round2b_command_line_is_unchanged_without_the_round2c_settings(world):
    out = _dry(world, [0, 1]).stdout  # the 6-column 2b config, no SELECTION_SEED_MODE
    assert "--selection-seed-mode" not in out and "--trim-r1-to-geometry" not in out
    assert "selection=global trim_r1=no" in out


def test_round2c_refuses_bad_settings_and_a_pin_without_the_flags(world):
    assert _dry(world, [0], SELECTION_SEED_MODE="shuffle").returncode == 2
    bad = _dry(world, [2], CONFIG_TSV=str(_config_2c(world, "maybe")))
    assert bad.returncode == 3 and "trim_r1 must be yes or no" in bad.stderr
    (world / "code" / "scripts").mkdir()
    (world / "code" / "scripts" / "sens07_planted_sensitivity.py").write_text("# 2b planter\n")
    old = _dry(world, [0], SELECTION_SEED_MODE="per-family")
    assert old.returncode == 2 and "has no --selection-seed-mode" in old.stderr
    assert _dry(world, [0]).returncode == 0  # a round-2b pin is fine without the flags


def test_round2c_settings_refuse_the_round2b_directory(world):
    cases = (
        (0, {"SELECTION_SEED_MODE": "per-family"}),
        (2, {"CONFIG_TSV": str(_config_2c(world))}),  # hhv6b row: trim_r1=yes
    )
    for row, env in cases:
        for round_dir in (f"{M}/viralscan_work/dsr_round2b", ""):  # explicit, or unset -> default
            done = _dry(world, [row], ROUND_DIR=round_dir, **env)
            assert done.returncode == 2 and "set ROUND_DIR to dsr_round2c" in done.stderr, env
    # round 2b's own settings may use that directory
    assert _dry(world, [0], ROUND_DIR=f"{M}/viralscan_work/dsr_round2b").returncode == 0
