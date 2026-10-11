"""_write_evidence_manifest hashes references (kallisto index, host FASTA) with Path.read_bytes():
peak memory ~ file size. run_safety.sha256_file streams in 1 MiB blocks."""
import tempfile
import tracemalloc
from pathlib import Path

from viralscan.scripts import evidence_run

tmp = Path(tempfile.mkdtemp())
big = tmp / "index.idx"
with open(big, "wb") as h:
    h.truncate(300 * 1024 * 1024)  # 300 MiB stand-in for a multi-GB kallisto index
out = tmp / "ev"
out.mkdir()
tracemalloc.start()
evidence_run._write_evidence_manifest(
    out, target_label="X", target_genes=["g"], method="equal", run_fingerprint=None,
    references={"index": str(big)},
)
peak = tracemalloc.get_traced_memory()[1]
print(f"file 300 MiB, tracemalloc peak {peak / 2**20:.0f} MiB")
