"""Pull the reads behind a viral call out of a finished run, for BLAST, alignment or any other tool.

The ViralScan counterpart of kallisto's "extract reads that aligned to specific virus IDs"
(https://kallisto.readthedocs.io/en/latest/translated/notebooks/virus_detection_sc.html).
Extraction is `viralscan evidence` without a reference: the exact (CB, UMI) molecules assigned to
the target are replayed through kallisto, so the reads are the ones the call was counted from.

    from viralscan.reads import extract_virus_reads, align_reads, blast_reads

    got = extract_virus_reads("out/S1/S1", ["Human papillomavirus 16"], "out/S1/reads")
    seqs = got["Human papillomavirus 16"].sequences()      # {read name: sequence}
    align_reads(got[...].fasta, "viral.fa", "hpv16.bam")   # minimap2 -> sorted, indexed BAM
    blast_reads(got[...].fasta, "viral.fa", "out/S1/blast")  # local BLAST, best hit per read
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from viralscan.evidence import align_reads_to_viral as align_reads
from viralscan.evidence import blast_identity as blast_reads

__all__ = ["VirusReads", "align_reads", "blast_reads", "extract_virus_reads", "read_fasta"]


def read_fasta(path: str | Path) -> dict[str, str]:
    """{record name: sequence} from a FASTA file (read names are `<CB>_<UMI>_<n>`)."""
    seqs: dict[str, list[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None and line:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


@dataclass(frozen=True)
class VirusReads:
    virus: str
    out_dir: Path

    @property
    def fasta(self) -> Path:
        return self.out_dir / "viral_reads.fasta"

    @property
    def lineage_tsv(self) -> Path:
        return self.out_dir / "read_lineage.tsv.gz"

    def sequences(self) -> dict[str, str]:
        return read_fasta(self.fasta)

    def lineage(self) -> list[dict[str, str]]:
        """One row per read: which cell/UMI and equivalence class it came from."""
        with gzip.open(self.lineage_tsv, "rt") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))


def extract_virus_reads(
    run_dir: str | Path, viruses: Iterable[str] | str, out_dir: str | Path, *, cores: int = 4
) -> dict[str, VirusReads]:
    """Extract the reads assigned to each virus into `<out_dir>/<virus slug>/viral_reads.fasta`.

    `viruses` takes what `viralscan evidence --virus` takes: an accession or gene ID, a canonical
    label, or a detected call from viral_summary.tsv. Needs the FASTQs and kallisto index the run used.
    """
    from viralscan.scripts.evidence_run import run_evidence

    if isinstance(viruses, str):
        viruses = [viruses]
    out: dict[str, VirusReads] = {}
    for virus in viruses:
        target = Path(out_dir) / re.sub(r"[^A-Za-z0-9]+", "_", virus).strip("_")
        args = argparse.Namespace(
            run_dir=str(run_dir),
            output=str(target),
            virus=virus,
            viral_fasta=None,
            host_fasta=None,
            blast=False,
            read_start_profile=False,
            cell_tags=False,
            dedup="umi",
            bin_size=1,
            sampling_seed=42,
            cores=cores,
            verbose=False,
            quiet=True,
        )
        try:
            run_evidence(args)
        except SystemExit as exc:  # the CLI path dies with sys.exit; a library call should raise
            raise RuntimeError(
                f"read extraction failed for {virus!r} (exit {exc.code}); see the log above"
            ) from None
        out[virus] = VirusReads(virus, target)
    return out
