"""HR-05: evaluation settings survive parser/config/Snakemake boundaries."""

import argparse

import pytest

from viralscan.menu import _build_run_config, _hostresponse_evaluation_kwargs, build_parser
from viralscan.run_safety import build_run_manifest
from viralscan.runconfig import RunConfig


def test_group_options_roundtrip_into_worker_config(tmp_path):
    args = build_parser().parse_args(
        [
            "-s1",
            "R1.fq",
            "-s2",
            "R2.fq",
            "-x",
            "10xv3",
            "--hostresponse-cv",
            "group",
            "--hostresponse-groups",
            "donor",
            "--hostresponse-cv-folds",
            "3",
            "--hostresponse-cell-type-column",
            "cell_type",
            "--hostresponse-cell-types",
            "CD4 T",
            "B cells",
            "--hostresponse-permutations",
            "7",
            "--hostresponse-permutation-unit",
            "group",
        ]
    )
    cfg = _build_run_config(args, "out/", "index.idx", "t2g.txt", None, "R1.fq", "R2.fq")
    restored = RunConfig.from_snakemake_config(cfg.to_dict())
    assert restored.hostresponse_cv == "group"
    assert restored.hostresponse_groups == "donor"
    assert restored.hostresponse_cv_folds == 3
    assert restored.hostresponse_cell_types == ["CD4 T", "B cells"]
    assert restored.hostresponse_permutations == 7
    assert restored.hostresponse_permutation_unit == "group"


def test_rerun_inherits_saved_group_settings_and_respects_explicit_zero():
    cfg = RunConfig(
        hostresponse_cv="group", hostresponse_groups="donor", hostresponse_permutations=9
    )
    args = build_parser().parse_args(
        ["hostresponse", "-o", "out", "--host-h5ad", "host.h5ad", "--permutations", "0"]
    )
    opts = _hostresponse_evaluation_kwargs(args, cfg)
    assert opts["cv_mode"] == "group" and opts["groups_column"] == "donor"
    assert opts["permutations"] == 0


def test_unset_evaluation_settings_preserve_old_fingerprint():
    old = argparse.Namespace(multimap_method="equal", cell_calling="none")
    new = argparse.Namespace(**vars(old), hostresponse_cv=None, hostresponse_cell_types=None)
    assert build_run_manifest(old)["run_fingerprint"] == build_run_manifest(new)["run_fingerprint"]


@pytest.mark.parametrize(
    "settings",
    [
        {"hostresponse_cv": "group"},
        {"hostresponse_cv_folds": 0},
        {"hostresponse_permutations": -1},
        {"hostresponse_min_groups": 0},
        {"hostresponse_cell_types": ["B", "B"]},
    ],
)
def test_config_rejects_invalid_evaluation_contract(settings):
    with pytest.raises(ValueError, match="[Hh]ostresponse|[Hh]ost-response"):
        RunConfig.from_snakemake_config(settings)


def test_ncbi_credentials_are_not_recorded_or_fingerprinted():
    """AUD-09 / C-02: the key and email stay out of run_manifest.json and the fingerprint."""
    base = argparse.Namespace(multimap_method="equal", cell_calling="none")
    with_creds = argparse.Namespace(
        **vars(base), ncbi_api_key="SECRETKEY123", ncbi_email="me@example.org"
    )
    manifest = build_run_manifest(with_creds)
    assert "SECRETKEY123" not in str(manifest) and "me@example.org" not in str(manifest)
    assert manifest["run_fingerprint"] == build_run_manifest(base)["run_fingerprint"]


# ── HR-02: the null specification is enforced before kb count, not inside the rule ──

_GROUP = {"hostresponse_cv": "group", "hostresponse_groups": "donor"}
_BAD_NULLS = [
    ({"hostresponse_permutations": 5, "hostresponse_permutation_unit": "group"}, "group CV"),
    ({**_GROUP, "hostresponse_permutations": 5}, "declared explicitly"),
    (
        {**_GROUP, "hostresponse_permutations": 5, "hostresponse_permutation_unit": "cell_within_block"},
        "requires a permutation block",
    ),
    (
        {
            **_GROUP,
            "hostresponse_permutations": 5,
            "hostresponse_permutation_unit": "group",
            "hostresponse_permutation_block": "batch",
        },
        "only applies",
    ),
]  # fmt: skip


@pytest.mark.parametrize(("settings", "reason"), _BAD_NULLS)
def test_config_rejects_an_invalid_null_specification(settings, reason):
    with pytest.raises(ValueError, match=reason):
        RunConfig.from_snakemake_config(settings)


@pytest.mark.parametrize(("settings", "reason"), _BAD_NULLS)
def test_errorhandler_rejects_an_invalid_null_specification(settings, reason, caplog):
    from unittest.mock import patch

    from viralscan.menu import errorhandler

    args = argparse.Namespace(
        ncbi_accession=None,
        reference=False,
        gtf=None,
        fasta=None,
        index="/fake/index.idx",
        transcripts="/fake/t2g.txt",
        f1=None,
        host_filter=None,
        host_index=None,
        sample1="/fake/R1.fastq.gz",
        sample2="/fake/R2.fastq.gz",
        **settings,
    )
    with patch("os.path.exists", return_value=True), pytest.raises(SystemExit) as exc:
        errorhandler(args)
    assert exc.value.code != 0
    assert reason in caplog.text


def test_a_declared_group_null_passes_the_config():
    args = build_parser().parse_args(
        ["-s1", "R1.fq", "-s2", "R2.fq", "-x", "10xv3"]
        + ["--hostresponse-cv", "group", "--hostresponse-groups", "donor"]
        + ["--hostresponse-permutations", "5", "--hostresponse-permutation-unit", "group"]
    )
    cfg = _build_run_config(args, "out/", "index.idx", "t2g.txt", None, "R1.fq", "R2.fq")
    assert cfg.hostresponse_permutations == 5


# ── HR-05: citation year and help text match what the code does ──


def test_help_cites_luebbert_2025_and_states_the_differential_family(capsys):
    from viralscan.scripts import hostresponse as hr

    texts = []
    for argv in (["hostresponse", "--help"], ["--help"]):
        with pytest.raises(SystemExit):
            build_parser().parse_args(argv)
        texts.append(capsys.readouterr().out)
    assert all("Luebbert et al. 2026" not in t for t in texts)
    assert any("Luebbert et al. 2025" in t for t in texts)
    with pytest.raises(SystemExit):
        hr._build_parser().parse_args(["--help"])
    standalone = " ".join(capsys.readouterr().out.split())
    assert "over ALL features" not in standalone and "not donor-adjusted" in standalone


# ── HR-04: the in-fold panel baseline needs group CV, checked before kb count ──


def test_panel_in_fold_without_group_cv_fails_before_kb_count():
    with pytest.raises(ValueError, match="panel_in_fold requires group CV"):
        RunConfig.from_snakemake_config({"hostresponse_panel_in_fold": True})
    args = build_parser().parse_args(
        ["-s1", "R1.fq", "-s2", "R2.fq", "-x", "10xv3", "--hostresponse-panel-in-fold"]
    )
    with pytest.raises(ValueError, match="panel_in_fold requires group CV"):
        _build_run_config(args, "out/", "index.idx", "t2g.txt", None, "R1.fq", "R2.fq")
    ok = build_parser().parse_args(
        ["-s1", "R1.fq", "-s2", "R2.fq", "-x", "10xv3", "--hostresponse-panel-in-fold"]
        + ["--hostresponse-cv", "group", "--hostresponse-groups", "donor"]
    )
    cfg = _build_run_config(ok, "out/", "index.idx", "t2g.txt", None, "R1.fq", "R2.fq")
    assert cfg.hostresponse_panel_in_fold is True


# ── HR-03: one --cell-types syntax on every entry point; strata need their column early ──


@pytest.mark.parametrize(
    "tokens", [["T_cells,B_cells"], ["T_cells", "B_cells"], ["T_cells,", "B_cells"], [" T_cells , B_cells"]]
)  # fmt: skip
def test_cell_types_syntax_is_the_same_on_every_entry_point(tokens):
    from viralscan.scripts import hostresponse as hr

    want = ["T_cells", "B_cells"]
    parser = build_parser()
    sub = parser.parse_args(["hostresponse", "-o", "out", "--host-h5ad", "h", "--cell-types", *tokens])
    main = parser.parse_args(
        ["-s1", "R1.fq", "-s2", "R2.fq", "-x", "10xv3", "--hostresponse-cell-types", *tokens]
    )
    flags = ["--virus-h5ad", "v", "--host-h5ad", "h", "--viral-accessions", "a", "--output", "o"]
    alone = hr._build_parser().parse_args([*flags, "--cell-types", *tokens])
    assert sub.cell_types == main.hostresponse_cell_types == alone.cell_types == want


def test_cell_types_without_a_column_fail_before_kb_count():
    with pytest.raises(ValueError, match="cell_type_column"):
        RunConfig.from_snakemake_config({"hostresponse_cell_types": ["B"]})
    args = build_parser().parse_args(
        ["-s1", "R1.fq", "-s2", "R2.fq", "-x", "10xv3", "--hostresponse-cell-types", "B"]
    )
    with pytest.raises(ValueError, match="cell_type_column"):
        _build_run_config(args, "out/", "index.idx", "t2g.txt", None, "R1.fq", "R2.fq")
