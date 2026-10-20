"""Where ViralScan-only artifacts live under a sample's output directory.

evonk's ViralScan 2.2.0 wrote ``config.yaml``, ``summary.txt``, ``log/``,
``kb-python/`` and ``plots/``. Everything v3 added lives under ``results/`` (and
``read_filtered/``). This module names those paths once, like
:mod:`viralscan.kb_outputs` does for ``kb-python/``, and performs no I/O.

``LAYOUT_VERSION`` is stamped in ``run_manifest.json``; bump it when a published
path moves.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union

LAYOUT_VERSION = 1


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
    def read_filtered_dir(self) -> Path:
        return self.output / "read_filtered"

    @property
    def count_audit(self) -> Path:
        return self.results_dir / "count_audit.tsv"
