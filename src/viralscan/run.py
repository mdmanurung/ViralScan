"""Read-only access to a finished ViralScan run: ``open_run(path)``.

A thin facade over the existing readers and path names. It resolves a run root (or one
sample directory) and returns pandas and AnnData objects. It never writes into the run
directory and never changes an output. Stable names: ``open_run`` and ``Run``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import anndata
    import pandas

__all__ = ["Run", "open_run"]

SUMMARY_TSV = Path("results") / "viral_summary.tsv"
PER_CELL_TSV = Path("results") / "per_cell_viral.tsv"


class Run:
    """One sample of a finished run. Build it with :func:`open_run`.

    ``root`` is the run root (holds ``run_manifest.json``), ``sample_dir`` is the sample
    directory, and ``sample`` is its name.
    """

    def __init__(self, root: Path, sample_dir: Path) -> None:
        self.root = Path(root)
        self.sample_dir = Path(sample_dir)
        self.sample = self.sample_dir.name

    def __repr__(self) -> str:
        return f"Run(root={str(self.root)!r}, sample={self.sample!r})"

    def manifest(self) -> dict[str, Any]:
        """The run's ``run_manifest.json`` as a dict (options, versions, fingerprint)."""
        from viralscan.run_safety import RUN_MANIFEST

        manifest: dict[str, Any] = json.loads(
            (self.root / RUN_MANIFEST).read_text(encoding="utf-8")
        )
        return manifest

    def summary(self) -> pandas.DataFrame:
        """``results/viral_summary.tsv``: one row per candidate virus, as a DataFrame."""
        return self._read_tsv(SUMMARY_TSV)

    def cells(self, called_only: bool = False) -> pandas.DataFrame:
        """``results/per_cell_viral.tsv``: one row per infected cell and virus.

        ``called_only=True`` keeps rows whose ``is_called_cell`` is true, the
        denominator the headline numbers use.
        """
        table = self._read_tsv(PER_CELL_TSV)
        if called_only:
            flag = table["is_called_cell"].astype(str).str.lower() == "true"
            table = table[flag].reset_index(drop=True)
        return table

    def adata(self) -> anndata.AnnData:
        """The selected-method count matrix, ``adata_multimap.h5ad``.

        Falls back to the uncorrected ``adata.h5ad`` when no multimapping step wrote one.
        """
        import anndata

        from viralscan.kb_outputs import KbCountOutputs

        kb = KbCountOutputs(self.sample_dir)
        path = kb.adata_multimap if kb.adata_multimap.is_file() else kb.adata
        if not path.is_file():
            raise FileNotFoundError(f"No count matrix under {kb.counts_dir}")
        return anndata.read_h5ad(path)

    def _read_tsv(self, relative: Path) -> pandas.DataFrame:
        import pandas as pd

        path = self.sample_dir / relative
        if not path.is_file():
            raise FileNotFoundError(f"{path} is missing; the run did not write it.")
        return pd.read_csv(path, sep="\t")


def open_run(path: str | Path, sample: str | None = None) -> Run:
    """Open a finished run read-only.

    ``path`` is the run root (the ``-o`` directory) or one sample directory. ``sample``
    names the sample when the root holds several. Raises ``ValueError`` when the path is
    not a finished run or the sample is ambiguous.
    """
    from viralscan.discovery import DiscoveryError, resolve_sample_dir

    path = Path(path)
    if (path / "config.yaml").is_file() and sample is None:
        root = path.parent  # a sample directory
        sample = path.name
    else:
        root = path
    try:
        sample_dir = resolve_sample_dir(root, sample)
    except DiscoveryError as err:
        raise ValueError(str(err)) from err
    return Run(root, sample_dir)
