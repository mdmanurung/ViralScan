"""The on-disk layout that ``kb count`` writes under a Run's output directory.

This module owns the kb-python file convention so it lives in exactly one place
instead of being hand-built as f-strings in ``multimap.py``, ``detection.py``,
``umap.py`` and the ``Snakefile``. It answers "what / where" only — it performs
no I/O (consistent with the passive Run Context). See ``CONTEXT.md``
("Kb Count Outputs").

After ``kb_count`` runs, the Snakefile moves kb-python's products into
``<output>/kb-python/``; ``counts_unfiltered/`` holds the count matrices.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union


@dataclass(frozen=True)
class KbCountOutputs:
    """Named paths within one Run's ``kb-python`` output tree.

    ``output`` is the per-sample output directory (the same value the rules read
    as ``config["output"]``); ``kb-python/`` lives directly beneath it.
    """

    output: Path

    @classmethod
    def from_config_output(cls, output: Union[str, Path]) -> KbCountOutputs:
        """Build from a ``config["output"]`` value (which may carry a trailing sep)."""
        return cls(Path(output))

    # ── directories ───────────────────────────────────────────────────────
    @property
    def root(self) -> Path:
        return self.output / "kb-python"

    @property
    def counts_dir(self) -> Path:
        return self.root / "counts_unfiltered"

    # ── files under kb-python/ ────────────────────────────────────────────
    @property
    def bus(self) -> Path:
        return self.root / "output.bus"

    @property
    def bus_txt(self) -> Path:
        return self.root / "output.bus.txt"

    @property
    def ec(self) -> Path:
        return self.root / "matrix.ec"

    @property
    def transcripts_txt(self) -> Path:
        """kb-python's ``transcripts.txt`` — distinct from the t2g in config."""
        return self.root / "transcripts.txt"

    # ── files under kb-python/counts_unfiltered/ ──────────────────────────
    @property
    def adata(self) -> Path:
        return self.counts_dir / "adata.h5ad"

    @property
    def adata_multimap(self) -> Path:
        return self.counts_dir / "adata_multimap.h5ad"

    @property
    def barcodes(self) -> Path:
        return self.counts_dir / "cells_x_genes.barcodes.txt"

    @property
    def genes(self) -> Path:
        return self.counts_dir / "cells_x_genes.genes.txt"

    @property
    def gene_names(self) -> Path:
        return self.counts_dir / "cells_x_genes.genes.names.txt"

    # ── resolution ────────────────────────────────────────────────────────
    def current_adata(self, *, multimapping: bool) -> Path:
        """The adata a Run's detection / umap step should read.

        Resolved from the ``multimapping`` flag, **not** file existence: if
        multimapper correction ran it wrote :attr:`adata_multimap`, so that file
        must be present. A missing file then surfaces as a read error rather than
        a silent downgrade to the uncorrected matrix.
        """
        return self.adata_multimap if multimapping else self.adata
