"""Shared STAR helpers: gunzipped whitelist and the exit-139 Solo fallback (ANDET-09-STAR139).

Used by ``scripts/anello_align`` and ``star_hostvirus``. Python standard library only.
"""

import gzip
import logging
import re
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)

_MAPPED_LINE = re.compile(
    r"(?:Uniquely mapped reads number|Number of reads mapped to multiple loci)\s*\|\s*(\d+)"
)


def plain_whitelist(whitelist: str | None, work_dir: Path) -> str | None:
    """Return a plain-text copy of a gzipped *whitelist*; STAR cannot read .gz."""
    if not whitelist or not whitelist.endswith(".gz"):
        return whitelist
    plain = work_dir / "whitelist.txt"
    with gzip.open(whitelist, "rb") as source, open(plain, "wb") as target:
        shutil.copyfileobj(source, target)
    return str(plain)


def run_star(cmd: list[str], work: Path, cb_len: int = 12) -> None:
    """Run STAR; recover from its Solo crash on a no-whitelist run that mapped nothing.

    STAR 2.7.11b segfaults (exit 139) in the Solo step with ``--soloCBwhitelist None``
    when no read maps. The retry gives it a one-barcode dummy whitelist (``cb_len`` As,
    which must equal ``--soloCBlen``) and is accepted only if ``Log.final.out`` reports
    zero mapped reads, so no barcode is dropped from a library that has alignments. A
    retry is recorded in ``star139_fallback.txt``. *work* is the ``--outFileNamePrefix``
    directory.
    """
    try:
        subprocess.run(cmd, check=True)
        return
    except subprocess.CalledProcessError as exc:
        no_whitelist = (
            "--soloCBwhitelist" in cmd and cmd[cmd.index("--soloCBwhitelist") + 1] == "None"
        )
        if exc.returncode != 139 or not no_whitelist:
            raise
    for leftover in ("_STARtmp", "_STARgenome", "_STARpass1", "Solo.out"):
        shutil.rmtree(work / leftover, ignore_errors=True)
    dummy = work / "dummy_whitelist.txt"
    dummy.write_text("A" * cb_len + "\n", encoding="utf-8")
    retry = list(cmd)
    retry[retry.index("--soloCBwhitelist") + 1] = str(dummy)
    subprocess.run(retry, check=True)
    mapped = sum(
        int(n) for n in _MAPPED_LINE.findall((work / "Log.final.out").read_text(encoding="utf-8"))
    )
    if mapped:
        raise RuntimeError(
            f"STAR exited 139 with no whitelist and the dummy-whitelist retry mapped {mapped} "
            "reads; refusing to drop barcodes"
        )
    (work / "star139_fallback.txt").write_text(
        "STAR exited 139 with --soloCBwhitelist None; retried with a dummy whitelist; "
        "zero reads mapped, so the result is empty\n",
        encoding="utf-8",
    )
    log.warning("STAR crashed with no whitelist and no mapped reads; empty result")
