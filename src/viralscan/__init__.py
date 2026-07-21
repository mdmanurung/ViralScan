"""ViralScan — quantify viral load from single-cell RNA-seq.

Bulk RNA-seq is not supported: the pipeline requires a cell-barcode/UMI geometry
(``evidence.cb_umi_geometry`` has no BULK entry, and the host-filter / evidence
paths raise ``ValueError`` without one). Use a droplet or plate single-cell
chemistry (10x v1/v2/v3, 10x 5', Drop-seq, or an explicit ``bc:umi:seq`` string).

The single source of truth for the package version. ``pyproject.toml`` reads this
value at build time via ``[tool.setuptools.dynamic] version = {attr =
"viralscan.__version__"}``, so bumping the release version means editing only this line.
"""

__version__ = "2.7.0"
