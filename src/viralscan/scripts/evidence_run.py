"""Orchestrator for ``viralscan evidence``.

Operates on a *completed* ViralScan run directory: traces the reads that drove
the viral calls, extracts them, optionally re-aligns to a viral genome for IGV,
and scores evidence quality (coverage + BLAST identity).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import sys
from pathlib import Path
from typing import NoReturn

from viralscan.evidence import (
    align_reads_to_viral,
    alignment_qc_table,
    competitive_blast_identity,
    coverage_table,
    deduplicate_bam,
    extract_exact_reads_by_number,
    have_tools,
    interpretation_flags,
    parse_flagged_target_bus,
    per_cell_alignment_qc,
    plot_coverage_comparison,
    read_start_distribution,
    replay_exact_target_bus,
    resolve_viral_target,
    write_competitive_fasta,
    write_igv_session,
    write_tagged_bam,
)
from viralscan.kb_outputs import KbCountOutputs
from viralscan.runconfig import RunConfig
from viralscan.scripts.multimap import load_transcripts, read_ec
from viralscan.utils import configure_logging

log = logging.getLogger("viralscan")


def _die(msg: str) -> NoReturn:
    log.error(msg)
    sys.exit(1)


def _write_evidence_manifest(
    output: Path,
    *,
    target_label: str,
    target_genes: list[str],
    method: str,
    run_fingerprint: str | None,
    references: dict[str, str | None],
) -> Path:
    """Write hashes for every completed evidence artifact."""
    outputs: dict[str, str] = {}
    for path in sorted(output.iterdir()):
        if path.name == "evidence_manifest.json" or not path.is_file():
            continue
        outputs[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    reference_hashes = {
        name: hashlib.sha256(Path(path).read_bytes()).hexdigest() if path else None
        for name, path in references.items()
    }
    manifest = {
        "schema_version": "3.0.0",
        "target": {"label": target_label, "genes": sorted(target_genes)},
        "method": method,
        "run_fingerprint": run_fingerprint,
        "references": reference_hashes,
        "outputs": outputs,
    }
    path = output / "evidence_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def run_evidence(args: argparse.Namespace) -> None:
    configure_logging(
        verbose=bool(getattr(args, "verbose", False)),
        quiet=bool(getattr(args, "quiet", False)),
    )
    run_dir = Path(args.run_dir).resolve()
    cfg_path = run_dir / "config.yaml"
    if not cfg_path.exists():
        _die(f"No config.yaml in {run_dir}; point --run-dir at a completed ViralScan run.")
    config = RunConfig.from_yaml(cfg_path)
    kb = KbCountOutputs.from_config_output(os.path.join(str(run_dir), ""))

    bus_text = kb.resolved_bus_txt
    if not bus_text.exists():
        _die(
            f"Missing v3 resolved BUS text {bus_text}. Run or rerun the v3 multimap step; "
            "raw output.bus is not a valid molecule lineage source."
        )

    for required in (kb.genes, kb.transcripts_txt, kb.ec):
        if not Path(required).exists():
            _die(f"Missing kb-python output {required}; --run-dir is not a completed run.")
    with open(kb.genes) as fh:
        gene_ids = [line.strip() for line in fh]
    transcripts, t2g_map = load_transcripts(str(kb.transcripts_txt), config.transcripts)  # type: ignore[no-untyped-call]
    ec_map = read_ec(str(kb.ec), transcripts, t2g_map, gene_ids)  # type: ignore[no-untyped-call]

    analysis = run_dir / "log" / "analysis.txt"
    if not analysis.exists():
        _die(f"No log/analysis.txt in {run_dir}; was the run completed?")
    with open(analysis) as fh:
        viral_ids = {line.strip() for line in fh}
    summary_path = run_dir / "results" / "viral_summary.tsv"
    detected_names: list[str] = []
    if summary_path.exists():
        with summary_path.open(newline="") as handle:
            detected_names = [row["virus_name"] for row in csv.DictReader(handle, delimiter="\t")]
    try:
        target_label, target_genes = resolve_viral_target(
            getattr(args, "virus", "") or "", viral_ids, detected_virus_names=detected_names
        )
    except ValueError as exc:
        _die(str(exc))
    viral_idx = {i for i, g in enumerate(gene_ids) if g in set(target_genes)}
    log.info("Tracing exact target %s (%d genes) …", target_label, len(viral_idx))

    technology = config.technology
    if not technology:
        _die("config.yaml has no 'technology'; cannot resolve barcode geometry. Re-run the sample.")
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    run_manifest_path = run_dir / "run_manifest.json"
    run_fingerprint = None
    if run_manifest_path.exists():
        try:
            run_fingerprint = json.loads(run_manifest_path.read_text()).get("run_fingerprint")
        except (OSError, json.JSONDecodeError):
            run_fingerprint = None
    if not config.index:
        _die("config.yaml has no kallisto index; exact read-lineage replay is impossible.")
    target_gene_set = set(target_genes)
    target_transcripts = [tx for tx in transcripts if t2g_map.get(tx) in target_gene_set]
    if not target_transcripts:
        _die(f"No transcripts resolve to exact target {target_label!r}.")
    try:
        flagged_text = replay_exact_target_bus(
            index=config.index,
            technology=technology,
            r1_path=config.sample1,
            r2_path=config.sample2,
            ec_file=str(kb.ec),
            transcripts_file=str(kb.transcripts_txt),
            target_transcripts=target_transcripts,
            workdir=str(out),
            threads=int(args.cores),
        )
        with flagged_text.open() as handle:
            lineage_by_number = parse_flagged_target_bus(
                handle,
                ec_map,
                viral_idx,
                {i for i, gene in enumerate(gene_ids) if gene in viral_ids},
                config.multimap_method,
            )
    except (RuntimeError, ValueError) as exc:
        _die(str(exc))
    log.info("%d exact target read records found by replay", len(lineage_by_number))

    ev_fasta = out / "viral_reads.fasta"
    stats = extract_exact_reads_by_number(
        config.sample1,
        config.sample2,
        lineage_by_number,
        str(ev_fasta),
        str(out / "read_lineage.tsv.gz"),
    )
    log.info(
        "Extracted %d viral-assigned reads (of %d) -> %s",
        stats.viral_reads,
        stats.total_reads,
        ev_fasta,
    )
    if stats.viral_reads == 0:
        log.warning("No reads extracted; nothing to align or BLAST.")
        _write_evidence_manifest(
            out,
            target_label=target_label,
            target_genes=target_genes,
            method=config.multimap_method,
            run_fingerprint=run_fingerprint,
            references={"index": config.index, "transcripts": config.transcripts},
        )
        return

    if args.viral_fasta:
        if not getattr(args, "host_fasta", None):
            _die(
                "--viral-fasta requires --host-fasta in v3 so evidence is aligned and "
                "BLASTed competitively against the full host genome."
            )
        missing = have_tools(["minimap2", "samtools"])
        if missing:
            _die(f"--viral-fasta given but missing tools: {', '.join(missing)}")
        competitive_fasta = write_competitive_fasta(
            args.host_fasta, args.viral_fasta, str(out / "competitive_host_target.fasta")
        )
        bam = align_reads_to_viral(
            str(ev_fasta),
            competitive_fasta,
            str(out / "competitive_reads.raw.bam"),
            int(args.cores),
        )
        dedup_mode = getattr(args, "dedup", "umi")
        dedup_bam = deduplicate_bam(
            bam,
            str(out / f"competitive_reads.{dedup_mode}_dedup.bam"),
            dedup_mode,
            int(args.cores),
        )
        raw_cov = coverage_table(bam)
        dedup_cov = coverage_table(dedup_bam)
        raw_qc = alignment_qc_table(bam)
        dedup_qc = alignment_qc_table(dedup_bam)
        dedup_reads = {str(row["reference"]): int(row["reads"]) for row in dedup_qc}
        qc_rows: list[dict[str, object]] = []
        for layer, rows in (("raw", raw_qc), ("deduplicated", dedup_qc)):
            for row in rows:
                raw_reads = int(row["reads"])
                if layer == "raw" and raw_reads:
                    duplicate_fraction = 1 - dedup_reads.get(str(row["reference"]), 0) / raw_reads
                else:
                    duplicate_fraction = 0.0
                qc_rows.append(
                    {"count_layer": layer, "duplicate_fraction": duplicate_fraction, **row}
                )
        qc_path = out / "alignment_qc.tsv"
        with qc_path.open("w", newline="") as handle:
            fieldnames = list(qc_rows[0]) if qc_rows else ["count_layer", "reference"]
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(qc_rows)
        cell_qc_rows: list[dict[str, object]] = []
        for layer, layer_bam in (("raw", bam), ("deduplicated", dedup_bam)):
            cell_qc_rows.extend(
                {"count_layer": layer, **row} for row in per_cell_alignment_qc(layer_bam)
            )
        with (out / "per_cell_alignment_qc.tsv").open("w", newline="") as handle:
            fields = list(cell_qc_rows[0]) if cell_qc_rows else ["count_layer", "cell_barcode"]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(cell_qc_rows)
        coverage_fields = [
            "rname",
            "startpos",
            "endpos",
            "numreads",
            "covbases",
            "coverage",
            "meandepth",
            "meanbaseq",
            "meanmapq",
        ]
        for rows, cov_path in (
            (raw_cov, out / "coverage.raw.tsv"),
            (dedup_cov, out / "coverage.deduplicated.tsv"),
        ):
            with open(cov_path, "w", newline="") as fh:
                w = csv.DictWriter(
                    fh, fieldnames=list(rows[0].keys()) if rows else coverage_fields, delimiter="\t"
                )
                w.writeheader()
                w.writerows(rows)
        plot_coverage_comparison(bam, dedup_bam, str(out / "coverage.raw_vs_deduplicated.png"))
        if not raw_cov:
            log.warning(
                "0 reads aligned to %s — viral call has no competitive support.",
                competitive_fasta,
            )
        else:
            log.info("Aligned -> %s; deduplicated -> %s", bam, dedup_bam)
            igv_bams = [bam, dedup_bam]

            if getattr(args, "read_start_profile", False):
                profile = read_start_distribution(
                    dedup_bam,
                    dedup="none",
                    bin_size=int(getattr(args, "bin_size", 1)),
                )
                prof_path = out / "read_start_profile.tsv"
                with open(prof_path, "w", newline="") as fh:
                    w = csv.DictWriter(
                        fh,
                        fieldnames=["reference", "position", "n_read_starts", "n_reads"],
                        delimiter="\t",
                    )
                    w.writeheader()
                    w.writerows(profile)
                log.info(
                    "Read-start profile (dedup=%s) -> %s (%d positions)",
                    dedup_mode,
                    prof_path,
                    len(profile),
                )

            if getattr(args, "cell_tags", False):
                tagged = write_tagged_bam(
                    dedup_bam, str(out / "competitive_reads.deduplicated.tagged.bam")
                )
                igv_bams.append(tagged)
                log.info("Cell-tagged BAM -> %s (IGV: group by tag CB)", tagged)
            session = write_igv_session(
                competitive_fasta, igv_bams, str(out / "viralscan_evidence.igv.xml")
            )
            log.info("IGV session -> %s", session)
        for r in raw_cov[:10]:
            log.info(
                "  %s: reads=%s coverage=%s%% meandepth=%s",
                r.get("rname"),
                r.get("numreads"),
                r.get("coverage"),
                r.get("meandepth"),
            )

        blast_rows: list[dict[str, str]] = []
        if getattr(args, "blast", False):
            blast_missing = have_tools(["blastn", "makeblastdb"])
            if blast_missing:
                _die(f"--blast requires {', '.join(blast_missing)} on PATH (install blast+).")
            else:
                blast_rows = competitive_blast_identity(
                    str(ev_fasta),
                    competitive_fasta,
                    str(out / "blast"),
                    threads=int(args.cores),
                    seed=int(args.sampling_seed),
                    sampling_manifest=str(out / "blast_sampling.json"),
                )
                bpath = out / "blast_identity.tsv"
                with open(bpath, "w", newline="") as fh:
                    fieldnames = (
                        list(blast_rows[0])
                        if blast_rows
                        else [
                            "read",
                            "top_viral_hit",
                            "viral_identity",
                            "viral_query_coverage",
                            "viral_evalue",
                            "viral_bitscore",
                            "top_host_hit",
                            "host_identity",
                            "host_query_coverage",
                            "host_evalue",
                            "host_bitscore",
                            "viral_minus_host_bitscore",
                            "low_complexity",
                        ]
                    )
                    w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
                    w.writeheader()
                    w.writerows(blast_rows)
                log.info("Competitive BLAST: %d reads -> %s", len(blast_rows), bpath)
        flag_rows = interpretation_flags(qc_rows, blast_rows, lineage_by_number.values())
        with (out / "interpretation_flags.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(flag_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(flag_rows)
    elif getattr(args, "blast", False):
        _die("--blast requires --viral-fasta (to build the local BLAST database).")
    elif getattr(args, "read_start_profile", False):
        _die("--read-start-profile requires --viral-fasta.")
    elif getattr(args, "cell_tags", False):
        _die("--cell-tags requires --viral-fasta.")

    log.info("Evidence outputs written under %s", out)
    _write_evidence_manifest(
        out,
        target_label=target_label,
        target_genes=target_genes,
        method=config.multimap_method,
        run_fingerprint=run_fingerprint,
        references={
            "index": config.index,
            "transcripts": config.transcripts,
            "viral_fasta": args.viral_fasta,
            "host_fasta": getattr(args, "host_fasta", None),
        },
    )
