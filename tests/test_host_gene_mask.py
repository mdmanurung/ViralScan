"""MECH-F: one reading of the gene-role contract for every consumer."""

import numpy as np
import pandas as pd

from viralscan.utils import host_gene_mask


def test_gene_role_wins_over_is_viral():
    var = pd.DataFrame({"gene_role": ["host", "latent", "unknown"], "is_viral": [True, True, True]})
    assert host_gene_mask(var).tolist() == [True, False, False]


def test_is_viral_and_legacy_viral_columns_are_complemented():
    assert host_gene_mask(pd.DataFrame({"is_viral": [False, True]})).tolist() == [True, False]
    assert host_gene_mask(pd.DataFrame({"viral": [True, False]})).tolist() == [False, True]


def test_no_role_columns_returns_none():
    assert host_gene_mask(pd.DataFrame({"gene_name": ["a"]})) is None
    assert isinstance(host_gene_mask(pd.DataFrame({"is_viral": [True]})), np.ndarray)
