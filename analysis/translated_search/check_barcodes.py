#!/usr/bin/env python3
"""Compare the barcodes of an arm-P count matrix with a frozen called-cell list.

Reports overlap before and after normalising (upper case, trailing ``-1``-style suffix
removed), suffix/length/case differences, and the overlap if one side were reverse-complemented.
Exits 2 when the share of frozen cells present in the observed barcodes is below
``--abort-below`` (default 0.9; 0 disables), so a format mismatch cannot silently shrink a
cell-restricted summary.

``OBSERVED`` and ``FROZEN`` are ``.h5ad`` files (obs names) or text files (first
tab-separated column, one barcode per line, optionally gzip-compressed).
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

_SUFFIX = re.compile(r"-\d+$")
_COMP = str.maketrans("ACGTN", "TGCAN")


def normalize_barcode(barcode: str) -> str:
    return _SUFFIX.sub("", barcode.strip().upper())


def read_barcodes(path: Path) -> list[str]:
    if str(path).endswith(".h5ad"):
        import anndata as ad

        adata = ad.read_h5ad(path, backed="r")
        try:
            return [str(b) for b in adata.obs_names]
        finally:
            adata.file.close()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return [
            line.split("\t")[0].strip()
            for line in handle
            if line.strip() and not line.startswith("#")
        ]


def _shape(barcodes: list[str]) -> dict:
    return {
        "n": len(barcodes),
        "with_suffix": sum(bool(_SUFFIX.search(b)) for b in barcodes),
        "lowercase": sum(b != b.upper() for b in barcodes),
        "lengths": dict(Counter(len(normalize_barcode(b)) for b in barcodes).most_common(5)),
    }


def compare(observed: list[str], frozen: list[str]) -> dict:
    obs, fro = {normalize_barcode(b) for b in observed}, {normalize_barcode(b) for b in frozen}
    raw = len(set(observed) & set(frozen))
    both = len(obs & fro)
    rc = len({b.translate(_COMP)[::-1] for b in obs} & fro)
    return {
        "observed": _shape(observed),
        "frozen": _shape(frozen),
        "raw_overlap": raw,
        "normalized_overlap": both,
        "frozen_found_fraction": both / len(fro) if fro else 0.0,
        "observed_in_frozen_fraction": both / len(obs) if obs else 0.0,
        "revcomp_overlap": rc,
        "format_differs": raw != both,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("observed", type=Path)
    parser.add_argument("frozen", type=Path)
    parser.add_argument(
        "--abort-below", type=float, default=0.9, help="minimum frozen_found_fraction (0 disables)"
    )
    args = parser.parse_args()
    report = compare(read_barcodes(args.observed), read_barcodes(args.frozen))
    print(json.dumps(report, indent=1))
    if report["frozen_found_fraction"] < args.abort_below:
        print(
            f"ABORT: only {report['frozen_found_fraction']:.1%} of frozen cells found "
            f"(threshold {args.abort_below:.0%})",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
