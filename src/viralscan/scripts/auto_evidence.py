"""Snakemake rule ``auto_evidence``: read-level host confirmation (PLAN ANDET-04).

Opt-in (``--auto-evidence``). After ``detection``, runs ``viralscan evidence`` with
``--viral-fasta`` and ``--host-fasta`` for every Anelloviridae genus listed in
``results/viral_summary.tsv``, into ``results/evidence/<genus>/``. The work is
:func:`viralscan.reads.extract_virus_reads` -> ``run_evidence``; this script only
picks the genera. A failed genus fails the rule (nothing is skipped silently).

``log/auto_evidence.done`` lists the genera processed, or ``none``, so an empty
result is distinguishable from a rule that never ran.
"""

import csv
import logging
from pathlib import Path
from typing import Any, Union

from viralscan.reads import extract_virus_reads
from viralscan.run_context import RunContext
from viralscan.virus_grouping import load_run_identity
from viralscan.virus_identity import VirusIdentityTable

log = logging.getLogger(__name__)

snakemake: Any


def anello_genera(summary_tsv: Union[str, Path], identity: VirusIdentityTable) -> list[str]:
    """Detected Anelloviridae genera in *summary_tsv*, in file order.

    ``alignment_only`` rows are skipped: the evidence replay only sees reads in
    kallisto equivalence classes, so those rows would extract nothing.
    """
    names = {g.virus_name for g in identity.genes if g.viral and g.family == "Anelloviridae"}
    with open(summary_tsv, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [
        r["virus_name"]
        for r in rows
        if r["virus_name"] in names and r.get("detection_source") != "alignment_only"
    ]


def run(ctx: RunContext, done_file: str) -> None:
    """Snakemake entry point for one Run."""
    config = ctx.config
    run_dir = Path(config.output)
    identity = load_run_identity(run_dir)
    if identity is None:
        raise RuntimeError(
            f"No results/virus_identity.tsv in {run_dir}; cannot tell which detected viruses "
            "are Anelloviridae."
        )
    genera = anello_genera(run_dir / "results" / "viral_summary.tsv", identity)
    log.info("auto_evidence: %d detected anellovirus genera: %s", len(genera), genera)
    # ponytail: one `evidence` run per genus, each copying the host genome into its own
    # competitive FASTA and re-reading the FASTQs; share one index across genera if slow.
    extract_virus_reads(
        run_dir,
        genera,
        run_dir / "results" / "evidence",
        cores=config.cores,
        viral_fasta=config.viral_fasta,
        host_fasta=config.host_fasta,
    )
    Path(done_file).write_text("\n".join(genera) + "\n" if genera else "none\n")


if "snakemake" in globals():
    run(
        RunContext.from_yaml(snakemake.params.configfile),  # noqa: F821
        snakemake.output[0],  # noqa: F821
    )
