"""Orchestrator for ``viralscan evidence``.

Operates on a *completed* ViralScan run directory: traces the reads that drove
the viral calls, extracts them, optionally re-aligns to a viral genome for IGV,
and scores evidence quality (coverage + BLAST identity).
"""

from __future__ import annotations

import csv
import logging
import subprocess
import sys
from pathlib import Path

from viralscan.evidence import (
    align_reads_to_viral,
    blast_identity,
    cb_umi_geometry,
    coverage_table,
    extract_viral_reads,
    have_tools,
    viral_assigned_keys,
    viral_equivalence_classes,
)
from viralscan.kb_outputs import KbCountOutputs
from viralscan.scripts.multimap import load_transcripts, read_ec
from viralscan.utils import configure_logging, load_config

log = logging.getLogger("viralscan")


def _die(msg: str) -> None:
    log.error(msg)
    sys.exit(1)


def run_evidence(args) -> None:
    configure_logging(
        verbose=bool(getattr(args, "verbose", False)),
        quiet=bool(getattr(args, "quiet", False)),
    )
    run_dir = Path(args.run_dir).resolve()
    cfg_path = run_dir / "config.yaml"
    if not cfg_path.exists():
        _die(f"No config.yaml in {run_dir}; point --run-dir at a completed ViralScan run.")
    config = load_config(str(cfg_path))
    kb = KbCountOutputs.from_config_output(str(run_dir) + "/")

    if not kb.bus_txt.exists():
        if have_tools(["bustools"]):
            _die("output.bus.txt missing and 'bustools' not on PATH to generate it.")
        log.info("Converting BUS to text …")
        subprocess.run(["bustools", "text", "-o", str(kb.bus_txt), str(kb.bus)], check=True)

    for required in (kb.genes, kb.transcripts_txt, kb.ec):
        if not Path(required).exists():
            _die(f"Missing kb-python output {required}; --run-dir is not a completed run.")
    gene_ids = [line.strip() for line in open(kb.genes)]
    transcripts, t2g_map = load_transcripts(str(kb.transcripts_txt), config["transcripts"])
    ec_map = read_ec(str(kb.ec), transcripts, t2g_map, gene_ids)

    analysis = run_dir / "log" / "analysis.txt"
    if not analysis.exists():
        _die(f"No log/analysis.txt in {run_dir}; was the run completed?")
    viral_ids = {line.strip() for line in open(analysis)}
    viral_idx = {i for i, g in enumerate(gene_ids) if g in viral_ids}
    if getattr(args, "virus", None):
        needle = args.virus.lower()
        viral_idx = {i for i in viral_idx if needle in gene_ids[i].lower()}
        if not viral_idx:
            _die(f"No viral genes match --virus {args.virus!r} in this run.")
    log.info("Tracing reads for %d viral genes …", len(viral_idx))

    viral_ecs = viral_equivalence_classes(ec_map, viral_idx)
    with open(kb.bus_txt) as fh:
        keys = viral_assigned_keys(fh, viral_ecs)
    log.info("%d viral-assigned (barcode, UMI) pairs", len(keys))

    technology = config.get("technology")
    if not technology:
        _die("config.yaml has no 'technology'; cannot resolve barcode geometry. Re-run the sample.")
    cb_len, umi_len = cb_umi_geometry(technology)
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    ev_fasta = out / "viral_reads.fasta"
    stats = extract_viral_reads(
        config["sample1"], config["sample2"], keys, cb_len, umi_len, str(ev_fasta)
    )
    log.info(
        "Extracted %d viral-assigned reads (of %d) -> %s",
        stats.viral_reads,
        stats.total_reads,
        ev_fasta,
    )
    if stats.viral_reads == 0:
        log.warning("No reads extracted; nothing to align or BLAST.")
        return

    if args.viral_fasta:
        missing = have_tools(["minimap2", "samtools"])
        if missing:
            _die(f"--viral-fasta given but missing tools: {', '.join(missing)}")
        bam = align_reads_to_viral(
            str(ev_fasta), args.viral_fasta, str(out / "viral_reads.bam"), int(args.cores)
        )
        cov = coverage_table(bam)
        cov_path = out / "coverage.tsv"
        if not cov:
            log.warning("0 reads aligned to %s — viral call has no genome-level support.", args.viral_fasta)
        else:
            with open(cov_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(cov[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(cov)
            log.info("Aligned -> %s (+ .bai); coverage -> %s", bam, cov_path)
            log.info("Open in IGV: load %s as genome, then %s", args.viral_fasta, bam)
        for r in cov[:10]:
            log.info(
                "  %s: reads=%s coverage=%s%% meandepth=%s",
                r.get("rname"),
                r.get("numreads"),
                r.get("coverage"),
                r.get("meandepth"),
            )

        if getattr(args, "blast", False):
            blast_missing = have_tools(["blastn", "makeblastdb"])
            if blast_missing:
                _die(f"--blast requires {', '.join(blast_missing)} on PATH (install blast+).")
            else:
                rows = blast_identity(
                    str(ev_fasta), args.viral_fasta, str(out / "blast"), threads=int(args.cores)
                )
                bpath = out / "blast_identity.tsv"
                with open(bpath, "w", newline="") as fh:
                    w = csv.DictWriter(
                        fh, fieldnames=["read", "subject", "pident", "length"], delimiter="\t"
                    )
                    w.writeheader()
                    w.writerows(rows)
                if rows:
                    idents = sorted(float(r["pident"]) for r in rows)
                    med = idents[len(idents) // 2]
                    log.info(
                        "BLAST: %d reads, median identity %.1f%% -> %s", len(rows), med, bpath
                    )
    elif getattr(args, "blast", False):
        _die("--blast requires --viral-fasta (to build the local BLAST database).")

    log.info("Evidence outputs written under %s", out)
