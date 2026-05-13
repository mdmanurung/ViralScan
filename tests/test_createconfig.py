"""Tests for config YAML creation logic (createconfig.py logic).

createconfig.py runs inside Snakemake, so we unit-test the config-building
logic directly — simulating the dict that Snakemake would provide as
``snakemake.config`` — without invoking Snakemake itself.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any

import pytest
import yaml

from viralscan.defaults import DEFAULTS


# ---------------------------------------------------------------------------
# Helpers mirroring the config-building logic in createconfig.py
# ---------------------------------------------------------------------------


def _build_cfg(cfg_in: dict[str, Any]) -> dict[str, Any]:
    """Reproduce the dict-building logic from createconfig.py."""
    detection_threshold = int(cfg_in.get("detection_threshold", 1))
    if detection_threshold < 1:
        raise ValueError(
            f"detection_threshold must be >= 1, got {detection_threshold}. "
            "A threshold of 0 or below would flag every viral accession as detected."
        )
    multimap_pseudocount = float(
        cfg_in.get("multimap_pseudocount", DEFAULTS["multimap_pseudocount"])
    )
    if multimap_pseudocount <= 0:
        raise ValueError(f"multimap_pseudocount must be > 0, got {multimap_pseudocount}.")
    return {
        **DEFAULTS,
        "output": cfg_in["output"],
        "index": cfg_in["index"],
        "transcripts": cfg_in["transcripts"],
        "sample1": cfg_in["sample1"],
        "sample2": cfg_in["sample2"],
        "overwrite": "yes",
        "gtf": cfg_in["gtf"] or None,
        "fasta": cfg_in["fasta"] or None,
        "visual": bool(cfg_in["visual"]),
        "f1": cfg_in["f1"] or None,
        "reference": bool(cfg_in["reference"]),
        "umap": bool(cfg_in["umap"]),
        "technology": cfg_in["technology"],
        "whitelist": cfg_in["whitelist"] or None,
        "multimapping": bool(cfg_in["multimapping"]),
        "se_threshold": int(cfg_in.get("se_threshold", DEFAULTS["se_threshold"])),
        "detection_threshold": detection_threshold,
        "min_counts": int(cfg_in.get("min_counts", DEFAULTS["min_counts"])),
        "min_genes": int(cfg_in.get("min_genes", DEFAULTS["min_genes"])),
        "hvg_min_mean": float(cfg_in.get("hvg_min_mean", DEFAULTS["hvg_min_mean"])),
        "hvg_max_mean": float(cfg_in.get("hvg_max_mean", DEFAULTS["hvg_max_mean"])),
        "hvg_min_disp": float(cfg_in.get("hvg_min_disp", DEFAULTS["hvg_min_disp"])),
        "umap_n_neighbors": int(cfg_in.get("umap_n_neighbors", DEFAULTS["umap_n_neighbors"])),
        "multimap_method": cfg_in.get("multimap_method", DEFAULTS["multimap_method"]),
        "multimap_pseudocount": multimap_pseudocount,
        "multimap_primary_call": cfg_in.get(
            "multimap_primary_call", DEFAULTS["multimap_primary_call"]
        ),
        "cell_types": cfg_in.get("cell_types") or None,
    }


def _minimal_cfg_in(**overrides) -> dict[str, Any]:
    """Return a minimal valid cfg_in dict with optional overrides."""
    base: dict[str, Any] = {
        "output": "/out/",
        "index": "/ref/index.idx",
        "transcripts": "/ref/t2g.txt",
        "sample1": "R1.fastq.gz",
        "sample2": "R2.fastq.gz",
        "gtf": None,
        "fasta": None,
        "visual": True,
        "f1": None,
        "reference": False,
        "umap": False,
        "technology": "10xv3",
        "whitelist": None,
        "multimapping": True,
        "multimap_method": DEFAULTS["multimap_method"],
        "multimap_pseudocount": DEFAULTS["multimap_pseudocount"],
        "multimap_primary_call": DEFAULTS["multimap_primary_call"],
    }
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# Boolean round-trip (§1.6 regression)
# ---------------------------------------------------------------------------


class TestBooleanFields:
    """Ensure boolean config values are stored as proper Python bools, not strings."""

    @pytest.mark.parametrize(
        "field,truthy_input",
        [
            ("visual", True),
            ("visual", 1),
            ("visual", "True"),  # string "True" is truthy — must become Python True
            ("reference", True),
            ("umap", True),
            ("multimapping", True),
        ],
    )
    def test_truthy_becomes_true(self, field, truthy_input) -> None:
        cfg = _build_cfg(_minimal_cfg_in(**{field: truthy_input}))
        assert cfg[field] is True
        assert isinstance(cfg[field], bool)

    @pytest.mark.parametrize(
        "field,falsy_input",
        [
            ("visual", False),
            ("visual", 0),
            ("reference", False),
            ("umap", False),
            ("multimapping", False),
        ],
    )
    def test_falsy_becomes_false(self, field, falsy_input) -> None:
        cfg = _build_cfg(_minimal_cfg_in(**{field: falsy_input}))
        assert cfg[field] is False
        assert isinstance(cfg[field], bool)


# ---------------------------------------------------------------------------
# None / empty-string normalisation
# ---------------------------------------------------------------------------


class TestNoneNormalisation:
    """Empty strings and falsy values in optional fields must become None."""

    @pytest.mark.parametrize("field", ["gtf", "fasta", "f1", "whitelist"])
    def test_none_stays_none(self, field) -> None:
        cfg = _build_cfg(_minimal_cfg_in(**{field: None}))
        assert cfg[field] is None

    @pytest.mark.parametrize("field", ["gtf", "fasta", "f1", "whitelist"])
    def test_empty_string_becomes_none(self, field) -> None:
        cfg = _build_cfg(_minimal_cfg_in(**{field: ""}))
        assert cfg[field] is None

    def test_non_empty_gtf_preserved(self) -> None:
        cfg = _build_cfg(_minimal_cfg_in(gtf="/path/to/virus.gtf"))
        assert cfg["gtf"] == "/path/to/virus.gtf"

    def test_cell_types_none_by_default(self) -> None:
        cfg = _build_cfg(_minimal_cfg_in())
        assert cfg["cell_types"] is None


# ---------------------------------------------------------------------------
# Integer threshold fields
# ---------------------------------------------------------------------------


class TestIntegerThresholds:
    def test_defaults_applied_when_absent(self) -> None:
        cfg = _build_cfg(_minimal_cfg_in())
        assert cfg["se_threshold"] == 10
        assert cfg["detection_threshold"] == 1
        assert cfg["min_counts"] == 1000
        assert cfg["min_genes"] == 200
        assert cfg["hvg_min_mean"] == DEFAULTS["hvg_min_mean"]
        assert cfg["hvg_max_mean"] == DEFAULTS["hvg_max_mean"]
        assert cfg["hvg_min_disp"] == DEFAULTS["hvg_min_disp"]
        assert cfg["umap_n_neighbors"] == DEFAULTS["umap_n_neighbors"]
        assert cfg["multimap_method"] == DEFAULTS["multimap_method"]
        assert cfg["multimap_pseudocount"] == DEFAULTS["multimap_pseudocount"]
        assert cfg["multimap_primary_call"] == DEFAULTS["multimap_primary_call"]

    def test_custom_thresholds_stored(self) -> None:
        cfg = _build_cfg(
            _minimal_cfg_in(
                se_threshold=50,
                detection_threshold=3,
                min_counts=500,
                min_genes=100,
                hvg_min_mean=0.2,
                hvg_max_mean=4.0,
                hvg_min_disp=1.5,
                umap_n_neighbors=25,
                multimap_method="unique-weighted",
                multimap_pseudocount=0.25,
                multimap_primary_call="unique-only",
            )
        )
        assert cfg["se_threshold"] == 50
        assert cfg["detection_threshold"] == 3
        assert cfg["min_counts"] == 500
        assert cfg["min_genes"] == 100
        assert cfg["hvg_min_mean"] == 0.2
        assert cfg["hvg_max_mean"] == 4.0
        assert cfg["hvg_min_disp"] == 1.5
        assert cfg["umap_n_neighbors"] == 25
        assert cfg["multimap_method"] == "unique-weighted"
        assert cfg["multimap_pseudocount"] == 0.25
        assert cfg["multimap_primary_call"] == "unique-only"

    @pytest.mark.parametrize(
        "field",
        [
            "se_threshold",
            "detection_threshold",
            "min_counts",
            "min_genes",
            "umap_n_neighbors",
        ],
    )
    def test_threshold_is_int(self, field) -> None:
        cfg = _build_cfg(_minimal_cfg_in())
        assert isinstance(cfg[field], int)

    @pytest.mark.parametrize("field", ["hvg_min_mean", "hvg_max_mean", "hvg_min_disp"])
    def test_hvg_threshold_is_float(self, field) -> None:
        cfg = _build_cfg(_minimal_cfg_in())
        assert isinstance(cfg[field], float)

    def test_multimap_pseudocount_is_float(self) -> None:
        cfg = _build_cfg(_minimal_cfg_in())
        assert isinstance(cfg["multimap_pseudocount"], float)


# ---------------------------------------------------------------------------
# YAML round-trip
# ---------------------------------------------------------------------------


class TestYamlRoundTrip:
    """Config dict serialised to YAML and read back must preserve types."""

    def _roundtrip(self, cfg_in: dict[str, Any]) -> dict[str, Any]:
        cfg = _build_cfg(cfg_in)
        buf = io.StringIO()
        yaml.dump(cfg, buf)
        buf.seek(0)
        return yaml.safe_load(buf)

    def test_booleans_survive_roundtrip(self) -> None:
        rt = self._roundtrip(_minimal_cfg_in(visual=True, multimapping=False))
        assert rt["visual"] is True
        assert rt["multimapping"] is False

    def test_none_survives_roundtrip(self) -> None:
        rt = self._roundtrip(_minimal_cfg_in(gtf=None))
        assert rt["gtf"] is None

    def test_strings_survive_roundtrip(self) -> None:
        rt = self._roundtrip(_minimal_cfg_in(output="/my/out/"))
        assert rt["output"] == "/my/out/"

    def test_integers_survive_roundtrip(self) -> None:
        rt = self._roundtrip(_minimal_cfg_in(detection_threshold=5))
        assert rt["detection_threshold"] == 5
        assert isinstance(rt["detection_threshold"], int)

    def test_all_required_keys_present(self) -> None:
        rt = self._roundtrip(_minimal_cfg_in())
        required = {
            "output",
            "index",
            "transcripts",
            "sample1",
            "sample2",
            "overwrite",
            "gtf",
            "fasta",
            "visual",
            "f1",
            "reference",
            "umap",
            "technology",
            "whitelist",
            "multimapping",
            "se_threshold",
            "detection_threshold",
            "min_counts",
            "min_genes",
            "hvg_min_mean",
            "hvg_max_mean",
            "hvg_min_disp",
            "umap_n_neighbors",
            "multimap_method",
            "multimap_pseudocount",
            "multimap_primary_call",
        }
        assert required.issubset(rt.keys())


# ---------------------------------------------------------------------------
# load_config utility
# ---------------------------------------------------------------------------


class TestLoadConfig:
    """viralscan.utils.load_config must read the YAML back correctly."""

    def test_load_config_returns_dict(self, tmp_path: Path) -> None:
        from viralscan.utils import load_config

        cfg = _build_cfg(_minimal_cfg_in())
        p = tmp_path / "config.yaml"
        with open(p, "w") as f:
            yaml.dump(cfg, f)
        loaded = load_config(str(p))
        assert loaded["output"] == "/out/"
        assert loaded["visual"] is True
        assert loaded["detection_threshold"] == 1

    def test_load_config_preserves_booleans(self, tmp_path: Path) -> None:
        from viralscan.utils import load_config

        cfg = _build_cfg(_minimal_cfg_in(visual=False, multimapping=True))
        p = tmp_path / "config.yaml"
        with open(p, "w") as f:
            yaml.dump(cfg, f)
        loaded = load_config(str(p))
        assert loaded["visual"] is False
        assert loaded["multimapping"] is True

    def test_load_config_missing_file_raises(self, tmp_path: Path) -> None:
        from viralscan.utils import load_config

        with pytest.raises(Exception):
            load_config(str(tmp_path / "nonexistent.yaml"))


# ---------------------------------------------------------------------------
# Audit §3.3 — detection_threshold must be >= 1
# ---------------------------------------------------------------------------


class TestDetectionThresholdValidation:
    """Audit §3.3: detection_threshold=0 makes every virus appear detected.

    Regression for: audits/2026-05-08-full-pipeline.md §3.3
    """

    def test_threshold_zero_raises_value_error(self) -> None:
        """
        GIVEN: detection_threshold=0
        WHEN:  _build_cfg validates the config
        THEN:  ValueError is raised before the config is returned
        """
        with pytest.raises(ValueError, match="detection_threshold"):
            _build_cfg(_minimal_cfg_in(detection_threshold=0))

    def test_threshold_negative_raises_value_error(self) -> None:
        """
        GIVEN: detection_threshold=-1
        WHEN:  _build_cfg validates the config
        THEN:  ValueError is raised
        """
        with pytest.raises(ValueError, match="detection_threshold"):
            _build_cfg(_minimal_cfg_in(detection_threshold=-1))

    def test_threshold_one_is_accepted(self) -> None:
        """threshold=1 (the default) must be accepted without error."""
        cfg = _build_cfg(_minimal_cfg_in(detection_threshold=1))
        assert cfg["detection_threshold"] == 1

    def test_threshold_large_is_accepted(self) -> None:
        """Any threshold >= 1 must be accepted."""
        cfg = _build_cfg(_minimal_cfg_in(detection_threshold=1000))
        assert cfg["detection_threshold"] == 1000


class TestMultimapConfigValidation:
    def test_pseudocount_zero_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="multimap_pseudocount"):
            _build_cfg(_minimal_cfg_in(multimap_pseudocount=0))

    def test_pseudocount_negative_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="multimap_pseudocount"):
            _build_cfg(_minimal_cfg_in(multimap_pseudocount=-0.5))

    def test_pseudocount_positive_is_accepted(self) -> None:
        cfg = _build_cfg(_minimal_cfg_in(multimap_pseudocount=0.1))
        assert cfg["multimap_pseudocount"] == 0.1
