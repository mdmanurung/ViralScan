"""Where ViralScan-only artifacts live under a sample's output directory.

evonk's ViralScan 2.2.0 wrote ``config.yaml``, ``summary.txt``, ``log/``,
``kb-python/`` and ``plots/``. Everything v3 added lives under ``results/`` (and
``read_filtered/``). This module names those paths once, like
:mod:`viralscan.kb_outputs` does for ``kb-python/``, and performs no I/O.

``LAYOUT_VERSION`` is stamped in ``run_manifest.json``; bump it when a published
path moves. Version 2 adds the reserved ``posthoc/`` subfolder, where post-hoc verbs
(``evidence``, ``rerun-multimap``, ``rerun-programs``, ``hostresponse``) write by default.
Readers accept every version in ``SUPPORTED_LAYOUT_VERSIONS``; version 1 trees stay valid.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union

LAYOUT_VERSION = 2
SUPPORTED_LAYOUT_VERSIONS = (1, 2)
POSTHOC_DIRNAME = "posthoc"


def posthoc_dir(run_dir: Union[str, Path], *parts: str) -> Path:
    """``<run_dir>/posthoc/<parts>``: the only place post-hoc verbs write by default."""
    return Path(run_dir).joinpath(POSTHOC_DIRNAME, *parts)


def is_posthoc(path: Path, run_dir: Path) -> bool:
    """True when ``path`` is inside ``run_dir``'s reserved ``posthoc/`` folder."""
    return path.relative_to(run_dir).parts[:1] == (POSTHOC_DIRNAME,)


@dataclass(frozen=True)
class V3Outputs:
    """Named v3 paths beneath one sample's output directory."""

    output: Path

    @classmethod
    def from_config_output(cls, output: Union[str, Path]) -> V3Outputs:
        """Build from a ``config["output"]`` value (which may carry a trailing sep)."""
        return cls(Path(output))

    @property
    def results_dir(self) -> Path:
        return self.output / "results"

    @property
    def posthoc_dir(self) -> Path:
        return posthoc_dir(self.output)

    @property
    def read_filtered_dir(self) -> Path:
        return self.output / "read_filtered"

    @property
    def count_audit(self) -> Path:
        return self.results_dir / "count_audit.tsv"
