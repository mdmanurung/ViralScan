"""Barcode-whitelist match-rate preflight.

Chemistry/whitelist mismatch is a silent failure mode. If the wrong whitelist
(or the wrong ``--technology``) is used, ``bustools correct`` discards nearly
every read and kb produces an all-empty-droplet matrix with **no error** —
finding F-005: a GEM-X 5' library mislabeled as ``10xv3`` lost 96.5% of its
reads and initially produced an invalidated result after a full compute job.

This module samples the first reads of R1, extracts the cell barcode using the
technology's CB geometry, and reports the fraction that match the whitelist, so
a chemistry mismatch is caught in seconds *before* a wasted run.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from typing import IO

from viralscan.evidence import cb_umi_geometry

DEFAULT_SAMPLE = 100_000
DEFAULT_MIN_MATCH_RATE = 0.5


def _open_text(path: str) -> IO[str]:
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p)


def load_whitelist(path: str) -> set[str]:
    """Read a barcode whitelist (one barcode per line; optionally gzipped)."""
    wl: set[str] = set()
    with _open_text(path) as fh:
        for line in fh:
            bc = line.strip()
            if bc:
                wl.add(bc)
    return wl


def whitelist_match_rate(
    r1_fastq: str, whitelist: str | set[str], cb_len: int, n_sample: int = DEFAULT_SAMPLE
) -> tuple[float, int]:
    """Fraction of sampled R1 barcodes (bases ``0:cb_len``) present in the whitelist.

    ``whitelist`` may be a path (str) or a preloaded ``set``. Reads at most
    ``n_sample`` records. Returns ``(match_rate, n_sampled)``; the rate is 0.0
    when no full-length barcode could be sampled.
    """
    wl = whitelist if isinstance(whitelist, set) else load_whitelist(whitelist)
    if not wl:
        raise ValueError("whitelist is empty")
    n = 0
    matched = 0
    with _open_text(r1_fastq) as fh:
        for i, line in enumerate(fh):
            # In a 4-line FASTQ record the sequence is the 2nd line (index % 4 == 1).
            if i % 4 != 1:
                continue
            bc = line.strip()[:cb_len]
            if len(bc) < cb_len:
                continue
            n += 1
            if bc in wl:
                matched += 1
            if n >= n_sample:
                break
    return (matched / n if n else 0.0), n


@dataclass
class WhitelistCheck:
    match_rate: float
    n_sampled: int
    cb_len: int
    ok: bool
    message: str


def check_whitelist(
    r1_fastq: str,
    whitelist: str | set[str],
    technology: str,
    min_match_rate: float = DEFAULT_MIN_MATCH_RATE,
    n_sample: int = DEFAULT_SAMPLE,
) -> WhitelistCheck:
    """Run the barcode match-rate check for one R1 file against ``whitelist``.

    ``cb_len`` is resolved from ``technology`` via :func:`cb_umi_geometry`, so a
    wrong ``--technology`` (which changes the barcode length/offset) also shows up
    as a low match rate.
    """
    cb_len, _ = cb_umi_geometry(technology)
    rate, n = whitelist_match_rate(r1_fastq, whitelist, cb_len, n_sample=n_sample)
    ok = rate >= min_match_rate
    if ok:
        msg = (
            f"{rate:.1%} of {n} sampled R1 barcodes match the whitelist "
            f"(technology={technology}, CB={cb_len}bp) — OK."
        )
    else:
        msg = (
            f"Only {rate:.1%} of {n} sampled R1 barcodes match the whitelist "
            f"(technology={technology}, CB={cb_len}bp). This usually means "
            f"--technology / --whitelist do not match the library chemistry; "
            f"bustools will discard most reads and produce an all-empty matrix "
            f"(finding F-005). Verify the chemistry (3' vs 5', v2 vs v3, GEM-X)."
        )
    return WhitelistCheck(match_rate=rate, n_sampled=n, cb_len=cb_len, ok=ok, message=msg)
