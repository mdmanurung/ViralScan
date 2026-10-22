"""AUD-06 / AUD-07 / B-16: one infected rule, a called-cell detection gate, a per-cell label cut-off.

Audit findings B-03, B-04, B-10, B-16, B-23 (analysis/audit_2026-10/findings.tsv). User decisions
of 2026-10-11, recorded in PLAN WP4L rows AUD-06 and AUD-07.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from viralscan.defaults import DEFAULTS  # noqa: E402
from viralscan.enrichment import cell_type_enrichment  # noqa: E402
from viralscan.runconfig import RunConfig  # noqa: E402
from viralscan.scripts import detection as D  # noqa: E402
from viralscan.scripts.umap import label_infected_cells  # noqa: E402
from viralscan.utils import infected_masks, virus_molecules_per_cell  # noqa: E402


def _adata(X, genes, cells=None):
    a = ad.AnnData(np.asarray(X, dtype=float))
    a.var_names = genes
    a.obs_names = cells or [f"c{i}" for i in range(a.n_obs)]
    return a


def _stats(a, called=None, genes=("v1",), **kw):
    stats, per_cell = D.compute_stats(
        a, {}, {"V": list(genes)}, list(genes), called_mask=called, **kw
    )
    return stats["V"], per_cell


# ── the shared rule ──────────────────────────────────────────────────────────


class TestInfectedRule:
    def test_repro_r_det1_tiny_fraction_is_subthreshold_not_infected(self):
        """B-03: X=[[1e-9,0,500],[0,0,500]] used to give infected_cells=1."""
        a = _adata([[1e-9, 0, 500], [0, 0, 500]], ["v1", "v2", "h"])
        s, per_cell = _stats(a)
        assert s["infected_cells"] == 0
        assert s["infected_subthreshold_cells"] == 1
        assert per_cell.empty
        assert s["detected"] is False

    def test_threshold_is_inclusive_and_sums_over_the_virus_genes(self):
        a = _adata([[0.5, 0.5, 9], [0.5, 0, 9], [1, 0, 9]], ["v1", "v2", "h"])
        s, per_cell = _stats(a, genes=("v1", "v2"))
        assert s["infected_cells"] == 2  # c0 (0.5+0.5) and c2 (1.0)
        assert s["infected_subthreshold_cells"] == 1  # c1
        assert list(per_cell["barcode"]) == ["c0", "c2"]

    def test_float32_rounding_just_under_one_molecule_is_infected(self):
        """EM splits of one molecule can sum to 0.99999994 in float32: not "sub-threshold"."""
        weight = np.array([np.nextafter(np.float32(1), np.float32(0))], dtype=np.float32)
        assert float(weight[0]) < 1.0
        infected, sub = infected_masks(weight, 1.0)
        assert infected.tolist() == [True] and sub.tolist() == [False]
        # ... but a genuinely smaller weight still is
        assert infected_masks(np.array([0.99]), 1.0)[1].tolist() == [True]

    def test_custom_cutoff(self):
        a = _adata([[2, 0, 9], [3, 0, 9]], ["v1", "v2", "h"])
        s, _ = _stats(a, infected_min_molecules=3)
        assert (s["infected_cells"], s["infected_subthreshold_cells"]) == (1, 1)

    def test_enrichment_uses_the_same_rule(self, tmp_path):
        a = _adata([[0.4, 5], [0.4, 5], [2.0, 5], [0, 5]], ["v1", "h"], list("ABCD"))
        csv = tmp_path / "ct.csv"
        pd.DataFrame({"barcode": list("ABCD"), "cell_type": ["T", "T", "B", "B"]}).to_csv(
            csv, index=False
        )
        out = cell_type_enrichment(a, {"V": ["v1"]}, RunConfig(cell_types=str(csv)))
        n_inf = dict(zip(out["cell_type"], out["n_infected"]))
        assert n_inf == {"T": 0, "B": 1}  # the 0.4-weight cells are not infected
        loose = cell_type_enrichment(
            a, {"V": ["v1"]}, RunConfig(cell_types=str(csv), infected_min_molecules=0.1)
        )
        assert dict(zip(loose["cell_type"], loose["n_infected"])) == {"T": 2, "B": 1}

    def test_umap_labels_use_the_same_rule_and_the_detection_matrix(self):
        a = _adata([[0.4, 0, 5], [0.4, 0.4, 5], [1.0, 0, 5], [0, 0, 5]], ["v1", "v2", "h"])
        labels, counts = label_infected_cells(a, a.X, ["v1", "v2"], None, 1.0)
        # gene->virus uses the legacy prefix rules here; v1 and v2 are separate "viruses",
        # so c1 (0.4 in v1, 0.4 in v2) is positive for neither.
        assert labels[0] == "No Virus" and labels[1] == "No Virus"
        assert labels[2] != "No Virus" and labels[3] == "No Virus"
        assert counts.tolist() == pytest.approx([0.4, 0.8, 1.0, 0.0])

    def test_umap_sums_the_genes_of_one_virus(self):
        a = _adata([[0.6, 0.6, 5]], ["v1", "v2", "h"])
        table = MagicMock()
        table.by_gene.return_value = {
            g: MagicMock(virus_name="VirusX", viral=True) for g in ("v1", "v2")
        }
        labels, _ = label_infected_cells(a, a.X, ["v1", "v2"], table, 1.0)
        assert labels == ["VirusX"]  # 1.2 molecules over the virus, neither gene reaches 1 alone

    def test_helper_sums_sparse_and_dense(self):
        import scipy.sparse as sp

        a = _adata([[0.5, 0.5, 1]], ["v1", "v2", "h"])
        a.X = sp.csr_matrix(a.X)
        assert virus_molecules_per_cell(a, a.X, ["v1", "v2"]).tolist() == [1.0]


# ── the called-cell gate ─────────────────────────────────────────────────────


class TestCalledCellGate:
    def _droplet_fixture(self):
        # c0, c1 are real cells with no virus; c2 is an empty droplet holding one viral UMI.
        a = _adata([[0, 600], [0, 800], [1, 3]], ["v1", "h"])
        return a, np.array([True, True, False])

    def test_one_umi_in_an_empty_droplet_is_a_candidate_not_detected(self):
        a, called = self._droplet_fixture()
        s, _ = _stats(a, called)
        assert s["detected"] is False
        assert s["infected_called"] == 0
        assert s["infected_cells"] == 1  # the all-barcode screen still sees it
        assert s["background_umis"] == 1
        assert s["viral_molecules_total_est"] == 1

    def test_a_called_cell_with_the_umi_passes(self):
        a = _adata([[1, 600], [0, 800], [0, 3]], ["v1", "h"])
        s, _ = _stats(a, np.array([True, True, False]))
        assert s["detected"] is True and s["background_umis"] == 0

    def test_cell_calling_none_means_every_barcode_is_called(self):
        """No called_mask: all barcodes are cells, so integer-count behaviour is unchanged."""
        a, _ = self._droplet_fixture()
        s, _ = _stats(a, None)
        assert s["detected"] is True and s["background_umis"] == 0
        assert s["infected_called"] == s["infected_cells"] == 1

    def test_min_cells(self):
        a = _adata([[1, 600], [0, 800], [4, 3]], ["v1", "h"])
        called = np.array([True, True, True])
        assert _stats(a, called, detection_min_cells=2)[0]["detected"] is True
        assert _stats(a, called, detection_min_cells=3)[0]["detected"] is False

    def test_detection_threshold_counts_called_cell_molecules_only(self):
        a = _adata([[1, 600], [0, 800], [9, 3]], ["v1", "h"])
        s, _ = _stats(a, np.array([True, True, False]), detection_threshold=2)
        assert s["detected"] is False  # 1 molecule in called cells, 9 in the droplet
        assert s["background_umis"] == 9

    def test_candidate_is_kept_in_viral_summary_with_new_columns(self, tmp_path):
        a, called = self._droplet_fixture()
        stats, per_cell = D.compute_stats(a, {}, {"V": ["v1"]}, ["v1"], called_mask=called)
        D.write_tsv_outputs(stats, per_cell, str(tmp_path))
        row = pd.read_csv(tmp_path / "results" / "viral_summary.tsv", sep="\t").iloc[0]
        assert row["virus_name"] == "V"
        assert bool(row["detected"]) is False
        assert row["background_umis"] == 1
        assert row["infected_subthreshold_cells"] == 0
        assert row["infected_called"] == 0

    def test_summary_line_flags_a_candidate_with_no_infected_called_cell(self):
        a, called = self._droplet_fixture()
        s, _ = _stats(a, called)
        line = D._virus_summary_line("V", s)
        assert "CANDIDATE ONLY" in line and "no infected called cell" in line
        assert "background" in line
        assert ";" not in line  # compare_legacy_v2_v3._GENE_TOTAL_RE

    def test_summary_line_of_a_detected_virus_has_no_flag(self):
        a = _adata([[3, 600]], ["v1", "h"])
        s, _ = _stats(a)
        assert "CANDIDATE" not in D._virus_summary_line("V", s)

    def test_gate_note_states_the_rule_and_the_cell_calling_none_caveat(self):
        cfg = RunConfig(detection_min_cells=2, detection_threshold=3, infected_min_molecules=1.5)
        note = D._gate_note(cfg, "none", ["W"])
        assert "at least 2 called cell(s)" in note and "1.5" in note and "at least 3" in note
        assert "Cell calling is off" in note
        assert "Candidate only, not detected" in note and "W" in note
        assert ";" not in note
        assert "Cell calling is off" not in D._gate_note(cfg, "emptydrops", [])

    def test_barcode_suffix_mismatch_does_not_silently_zero_the_comparable_cells(self):
        """B-10: host sets are suffix-stripped, obs_names are raw."""
        a = _adata([[5, 500], [0, 500], [2, 600]], ["v1", "h"], ["AAA-1", "CCC-1", "GGG-1"])
        hc = {"called": {"AAA", "CCC", "GGG"}, "comparable": {"AAA", "CCC", "GGG"}}
        s, _ = _stats(a, np.ones(3, bool), host_cells=hc)
        assert s["n_comparable_cells"] == 3 and s["infected_comparable"] == 2


# ── plumbing ─────────────────────────────────────────────────────────────────


class TestOptions:
    def test_cli_defaults_are_unset_and_resolve_at_use(self):
        from viralscan.cli.parsers import build_parser

        args = build_parser().parse_args(["run", "-o", "o", "--r1", "a", "--r2", "b"])
        assert args.infected_min_molecules is None
        assert args.detection_min_cells is None
        assert args.hostresponse_min_umi is None
        cfg = RunConfig.from_snakemake_config(_cfg_in())
        assert cfg.infected_min_molecules == DEFAULTS["infected_min_molecules"] == 1.0
        assert cfg.detection_min_cells == DEFAULTS["detection_min_cells"] == 1
        assert cfg.hostresponse_min_umi == DEFAULTS["hostresponse_min_umi"] == 10

    def test_values_are_recorded_in_config_yaml(self, tmp_path):
        cfg = RunConfig.from_snakemake_config(
            _cfg_in(infected_min_molecules=2.5, detection_min_cells=3, hostresponse_min_umi=20)
        )
        cfg.to_yaml(tmp_path / "config.yaml")
        back = RunConfig.from_yaml(str(tmp_path / "config.yaml"))
        assert (
            back.infected_min_molecules,
            back.detection_min_cells,
            back.hostresponse_min_umi,
        ) == (
            2.5,
            3,
            20,
        )

    @pytest.mark.parametrize(
        "key, value",
        [
            ("infected_min_molecules", 0),
            ("infected_min_molecules", -1),
            ("infected_min_molecules", float("nan")),
            ("detection_min_cells", 0),
            ("hostresponse_min_umi", 0),
        ],
    )
    def test_ranges_are_validated(self, key, value):
        with pytest.raises(ValueError, match=key):
            RunConfig.from_snakemake_config(_cfg_in(**{key: value}))


def _cfg_in(**extra):
    base = {
        "output": "o/", "index": "i", "transcripts": "t", "sample1": "a", "sample2": "b",
        "gtf": None, "fasta": None, "visual": True, "f1": None, "reference": False,
        "umap": False, "technology": "10xv3", "whitelist": None, "multimapping": True,
    }  # fmt: skip
    base.update(extra)
    return base


# ── B-16: the host-response per-cell label ───────────────────────────────────


def _sample_dir(tmp_path: Path, **cfg_kw) -> Path:
    d = tmp_path / "sample"
    (d / "log").mkdir(parents=True)
    RunConfig(output=str(d) + "/", multimapping=True, **cfg_kw).to_yaml(d / "config.yaml")
    (d / "log" / "analysis.txt").write_text("NC_001806.2\n")
    return d


def _run_subcommand(sample_dir: Path, **overrides) -> dict:
    from viralscan.menu import _run_hostresponse_subcommand

    args = MagicMock()
    for k, v in dict(
        run_dir=str(sample_dir), output=None, host_h5ad="h.h5ad", n_seeds=None, n_stab_iter=None, use_hvg=True,
        stab_min_prob=None, top_n_genes=None, detection_threshold=None, hostresponse_min_umi=None,
        enrichment=False, enrichment_db=None, verbose=False, quiet=False,
    ).items():  # fmt: skip
        setattr(args, k, v)
    for k, v in overrides.items():
        setattr(args, k, v)
    kb = MagicMock()
    kb.current_adata.return_value = sample_dir / "virus.h5ad"
    (sample_dir / "virus.h5ad").touch()
    with (
        patch("viralscan.kb_outputs.KbCountOutputs.from_config_output", return_value=kb),
        patch("viralscan.scripts.hostresponse.run_hostresponse") as run,
    ):
        _run_hostresponse_subcommand(args)
    return run.call_args.kwargs


class TestHostresponseMinUmi:
    def test_label_threshold_defaults_to_10_and_ignores_the_sample_gate(self, tmp_path):
        d = _sample_dir(tmp_path, detection_threshold=1)
        assert _run_subcommand(d)["detection_threshold"] == 10

    def test_config_value_is_used(self, tmp_path):
        d = _sample_dir(tmp_path, detection_threshold=1, hostresponse_min_umi=25)
        assert _run_subcommand(d)["detection_threshold"] == 25

    def test_flag_overrides_config(self, tmp_path):
        d = _sample_dir(tmp_path, hostresponse_min_umi=25)
        assert _run_subcommand(d, hostresponse_min_umi=3)["detection_threshold"] == 3

    def test_the_removed_detection_threshold_alias_is_ignored(self, tmp_path):
        # CLI-UNIFY-01: the parser refuses `hostresponse --detection-threshold`; a stray
        # namespace attribute must not set the label either.
        d = _sample_dir(tmp_path)
        assert _run_subcommand(d, detection_threshold=4)["detection_threshold"] == 10
