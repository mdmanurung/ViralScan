"""REVIEW-FIX-B: an explicit 0 is kept, out-of-range 0 exits 2, huge -c is clamped."""

import os
from unittest.mock import patch

import pytest
from test_errorhandler import _always_exists, _args

from viralscan.menu import errorhandler
from viralscan.runconfig import RunConfig

BASE = {
    "output": "/tmp/o/", "index": "/tmp", "transcripts": "t", "sample1": "s1", "sample2": "s2",
    "technology": "10xv3", "visual": "True", "multimapping": "True", "gtf": "None",
    "fasta": "None", "f1": "None", "reference": "False", "umap": "False", "whitelist": "None",
}  # fmt: skip


def test_valid_zeros_round_trip():
    c = RunConfig.from_snakemake_config(
        {**BASE, "emptydrops_lower": 0, "hostresponse_stab_min_prob": 0.0}
    )
    assert c.emptydrops_lower == 0
    assert c.hostresponse_stab_min_prob == 0.0


@pytest.mark.parametrize(
    "key",
    [
        "emptydrops_fdr",
        "emptydrops_niters",
        "hostresponse_n_stab_iter",
        "hostresponse_top_n_genes",
        "knee_min_umi",
    ],
)
def test_config_rejects_zero_where_out_of_range(key):
    with pytest.raises(ValueError, match=key):
        RunConfig.from_snakemake_config({**BASE, key: 0})


@pytest.mark.parametrize("attr", ["emptydrops_fdr", "emptydrops_niters", "knee_min_umi"])
def test_cli_rejects_zero_with_exit_2(attr):
    with patch("os.path.exists", side_effect=_always_exists), pytest.raises(SystemExit) as e:
        errorhandler(_args(**{attr: 0}))
    assert e.value.code == 2


def test_huge_cores_clamped_with_warning(caplog):
    args = _args(cores=99999999999)
    with patch("os.path.exists", side_effect=_always_exists):
        errorhandler(args)
    assert args.cores == (os.cpu_count() or 1)
    assert "exceeds" in caplog.text
