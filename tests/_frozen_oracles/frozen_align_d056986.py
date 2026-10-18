import subprocess
from pathlib import Path
from typing import Optional

from viralscan.evidence import MINIMAP2_ALIGN_OPTIONS


def _run(cmd: list[str], *, stdin: Optional[bytes] = None, capture: bool = False) -> bytes:
    """Run *cmd* (list form, no shell), surfacing the tool's stderr on failure.

    Returns stdout bytes when *capture* is True. On a non-zero exit raises
    ``RuntimeError`` that includes the tool name and its stderr — never silently
    swallows the error (the wrappers used to ``stderr=DEVNULL`` which hid the
    real cause behind a bare CalledProcessError).
    """
    proc = subprocess.run(
        cmd,
        input=stdin,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE,
    )
    if proc.returncode != 0:
        err = (proc.stderr or b"").decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"{cmd[0]} failed (exit {proc.returncode}): {err or 'no stderr output'}")
    return proc.stdout if capture else b""


def align_reads_to_viral(reads_fasta: str, viral_fasta: str, out_bam: str, threads: int = 4) -> str:
    """minimap2 short-read align *reads_fasta* to *viral_fasta* -> sorted, indexed BAM.

    Secondary alignments within 0.8 of the primary score (at most 5) are kept, so the host-versus-virus
    comparison sees both sides; see :data:`MINIMAP2_ALIGN_OPTIONS`.
    """
    out_bam = str(out_bam)
    Path(out_bam).parent.mkdir(parents=True, exist_ok=True)
    sam = _run(
        ["minimap2", *MINIMAP2_ALIGN_OPTIONS, "-t", str(threads), viral_fasta, reads_fasta],
        capture=True,
    )
    _run(["samtools", "sort", "-@", str(threads), "-o", out_bam, "-"], stdin=sam)
    _run(["samtools", "index", out_bam])
    return out_bam
