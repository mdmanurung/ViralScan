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
