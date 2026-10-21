"""Snakemake rule ``star_confirm``: STARsolo host+virus confirmation arm (PLAN STAR-04).

Aligns the reads kb counted (host-filtered, ``config.kb_r1/kb_r2``) to the ``viralscan
star-index`` index and writes ``results/star_support_by_virus.tsv``: the ``star_hostvirus``
summary per ``virus_name`` (accessions mapped through ``virus_identity.tsv``), which
``detection`` merges into ``viral_summary.tsv`` as labels only. The pruned VIRUS|-only BAM
stays in ``star_hostvirus/``.

``summarise_by_virus`` aggregates by virus itself: a read that places on several accessions of
one virus counts once for that virus, and ``star_other_virus`` needs a different virus.
"""

import csv
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

from viralscan.runconfig import RunConfig
from viralscan.star_hostvirus import STAR_SUMMARY_COLUMNS, align_host_unmapped

log = logging.getLogger(__name__)

STAR_TSV = "results/star_support_by_virus.tsv"
STAR_VIRUS_COLUMNS = ("virus_name", *STAR_SUMMARY_COLUMNS[1:])


def _acc_to_virus(identity_tsv: Path) -> dict[str, str]:
    """genome_accession -> virus_name for every row of virus_identity.tsv."""
    with open(identity_tsv, encoding="utf-8") as fh:
        return {
            r["genome_accession"]: r["virus_name"]
            for r in csv.DictReader(fh, delimiter="\t")
            if r.get("genome_accession")
        }


def read_summary(summary_tsv: Path) -> dict[str, dict[str, int]]:
    """The per-virus ``star_support_by_virus.tsv`` (its ``accession`` column holds virus_name)."""
    with open(summary_tsv, encoding="utf-8") as fh:
        return {
            row["accession"]: {c: int(row[c]) for c in STAR_VIRUS_COLUMNS[1:]}
            for row in csv.DictReader(fh, delimiter="\t")
        }


def main(config: RunConfig, identity_tsv: str, n_threads: int, done_path: str) -> None:
    out = Path(config.output)
    work = out / "star_hostvirus"
    index = config.star_hostvirus_index
    assert index, "star_confirm runs only with a star_hostvirus_index"
    virus_of = _acc_to_virus(Path(identity_tsv))
    pruned = align_host_unmapped(
        config.kb_r1,
        config.kb_r2,
        index,
        work,
        config.technology,
        whitelist=config.whitelist,
        threads=n_threads,
        strand=config.strand or "unstranded",  # kb's effective default is unstranded
        prune=True,
        virus_of=virus_of,
    )
    # The rule's declared outputs are viral.bam(.bai); phase 1 names the pruned file Aligned.viral.bam.
    pruned.replace(work / "viral.bam")
    pruned.with_name(pruned.name + ".bai").replace(work / "viral.bam.bai")
    rows = read_summary(work / "star_support_by_virus.tsv")
    with open(out / STAR_TSV, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(STAR_VIRUS_COLUMNS)
        for virus in sorted(rows):
            writer.writerow([virus, *(rows[virus][c] for c in STAR_VIRUS_COLUMNS[1:])])
    log.info("star_confirm: %d viruses with STAR alignment evidence", len(rows))
    Path(done_path).touch()


if TYPE_CHECKING:
    snakemake: Any  # injected by the Snakemake runner; declared for the type checker only

if "snakemake" in globals():
    main(
        config=RunConfig.from_yaml(snakemake.params.configfile),  # noqa: F821
        identity_tsv=str(snakemake.input.virus_identity),  # noqa: F821
        n_threads=snakemake.threads,  # noqa: F821
        done_path=str(snakemake.output.done),  # noqa: F821
    )
