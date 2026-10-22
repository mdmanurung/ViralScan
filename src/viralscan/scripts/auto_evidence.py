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
from viralscan.runconfig import RunConfig
from viralscan.virus_grouping import detected_summary_rows, load_run_identity
from viralscan.virus_identity import VirusIdentityTable

log = logging.getLogger(__name__)

snakemake: Any


REVIEW_FIELDS = (
    "virus_name",
    "detection_source",
    "detected",
    "evidence_source",
    "n_reads",
    "n_unique_cb_umi",
    "complex_body_fraction",
    "reagent_fraction",
    "status",
    "reason_codes",
)


def candidate_genera(
    summary_tsv: str | Path,
    identity: VirusIdentityTable,
    include_below_gate: bool = False,
    include_alignment_only: bool = False,
) -> list[str]:
    """Anelloviridae selection, with legacy defaults retained."""
    names = {g.virus_name for g in identity.genes if g.viral and g.family == "Anelloviridae"}
    with open(summary_tsv, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not include_below_gate:
        rows = detected_summary_rows(rows)
    return list(
        dict.fromkeys(
            r["virus_name"]
            for r in rows
            if r["virus_name"] in names
            and (include_alignment_only or r.get("detection_source") != "alignment_only")
        )
    )


def review_candidates(
    run_dir: str | Path,
    config: RunConfig,
    evidence_dir: str | Path,
    *,
    viral_fasta: str | Path | None = None,
    host_fasta: str | Path | None = None,
) -> list[dict[str, Any]]:
    """One diagnostic row per candidate, preserving source and primary outputs."""
    import os
    import tempfile

    from viralscan.anello_align import has_reagent, is_complex_body
    from viralscan.scripts.evidence_run import star_evidence

    run_dir, evidence_dir = Path(run_dir), Path(evidence_dir)
    identity = load_run_identity(run_dir)
    if identity is None:
        raise ValueError("Candidate review requires results/virus_identity.tsv")
    summary = run_dir / "results" / "viral_summary.tsv"
    names = set(candidate_genera(summary, identity, True, True))
    with summary.open(newline="") as handle:
        candidates = [r for r in csv.DictReader(handle, delimiter="\t") if r["virus_name"] in names]
    rows = []
    for candidate in candidates:
        name = candidate["virus_name"]
        row: dict[str, Any] = dict.fromkeys(REVIEW_FIELDS, "not_assessed")
        row.update(
            virus_name=name,
            detected=candidate.get("detected", ""),
            detection_source=candidate.get("detection_source", "kallisto"),
        )
        if row["detection_source"] == "alignment_only":
            targets = {
                g.genome_accession for g in identity.genes if g.viral and g.virus_name == name
            }
            row.update(
                star_evidence(targets, run_dir / "anello_align" / "Aligned.sortedByCoord.out.bam")
            )
        elif row["detection_source"] in {"kallisto", "kallisto+alignment", ""}:
            extracted = extract_virus_reads(
                run_dir,
                [name],
                evidence_dir,
                cores=config.cores,
                viral_fasta=viral_fasta,
                host_fasta=host_fasta,
            )[name]
            sequences = extracted.sequences()
            lineage = extracted.lineage()
            # Lineage keys are corrected replay keys; absent tags are not zero molecules.
            keys = {(r.get("cb"), r.get("ub")) for r in lineage}
            keys_valid = all(cb and umi for cb, umi in keys)
            row.update(
                evidence_source="kallisto_replay",
                status="reviewed",
                n_reads=len(sequences),
                n_unique_cb_umi=len(keys) if keys_valid else "not_assessed",
                complex_body_fraction=sum(is_complex_body(seq) for seq in sequences.values())
                / len(sequences)
                if sequences
                else "not_assessed",
                reagent_fraction=sum(has_reagent(seq) for seq in sequences.values())
                / len(sequences)
                if sequences
                else "not_assessed",
                reason_codes="" if keys_valid else "corrected_cb_umi_missing",
            )
        else:
            row.update(
                evidence_source="none",
                status="unsupported_source",
                reason_codes="unsupported_detection_source",
            )
        rows.append(row)
    output = evidence_dir / "evidence_review.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", dir=output.parent, delete=False, newline=""
    ) as handle:
        temporary = handle.name
        writer = csv.DictWriter(handle, fieldnames=REVIEW_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, output)
    return rows


def anello_genera(summary_tsv: Union[str, Path], identity: VirusIdentityTable) -> list[str]:
    """Detected Anelloviridae genera in *summary_tsv*, in file order.

    Rows with ``detected == False`` (candidates failing the called-cell gate) are skipped,
    as are ``alignment_only`` rows: the evidence replay only sees reads in
    kallisto equivalence classes, so those rows would extract nothing.
    """
    return candidate_genera(summary_tsv, identity)


def run(ctx: RunContext, done_file: str) -> None:
    """Snakemake entry point for one Run."""
    config = ctx.config
    run_dir = Path(config.output)
    if config.review_candidates:
        review_candidates(
            run_dir,
            config,
            run_dir / "results" / "evidence",
            viral_fasta=config.viral_fasta,
            host_fasta=config.host_fasta,
        )
        Path(done_file).write_text("candidate_review\n")
        return
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
        RunContext.from_yaml(snakemake.params.configfile).with_cores(snakemake.threads),  # noqa: F821
        snakemake.output[0],  # noqa: F821
    )
