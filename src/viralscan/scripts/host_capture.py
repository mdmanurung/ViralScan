"""Snakemake rule ``host_capture``: split each virus's reads by host compatibility (PLAN DLIST-02).

Opt-in (``--host-capture``). Replays the quantified FASTQs through ``kallisto bus -n`` on the
run's own index, captures every read whose equivalence class holds a viral transcript, then
applies ``bustools capture --complement`` with the host transcripts. Per virus:

* ``reads_virus_only``: reads whose class holds no host transcript;
* ``reads_also_host``: reads whose class holds the virus and a host transcript;
* ``reads_total``: their sum.

Counts are reads (not molecules): every read has its own BUS record. A read compatible with
two viruses counts for both. It needs an index that holds host transcripts (a combined
host + virus index). On a virus-only index the host list is empty, so every read would look
virus-only; the rule writes ``not_applicable`` for every virus instead (never 0).

Writes ``results/host_capture.tsv`` and ``log/host_capture.done``. Viruses with no read in
any viral class have no row.
"""

import csv
import logging
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, Union

from viralscan.evidence import replay_host_capture_bus
from viralscan.run_context import RunContext
from viralscan.scripts.evidence_run import _replay_whitelist, replay_fastqs
from viralscan.scripts.multimap import load_transcripts
from viralscan.virus_identity import VirusIdentityTable

log = logging.getLogger(__name__)

snakemake: Any

COLUMNS = [
    "virus_name",
    "reads_total",
    "reads_virus_only",
    "reads_also_host",
    "host_capture_status",
]
NOT_APPLICABLE = "not_applicable"
MEASURED = "measured"


def read_ec_transcripts(path: Union[str, Path]) -> dict[int, list[int]]:
    """``matrix.ec`` as ``{ec: [transcript index, ...]}`` (one ``ec<TAB>t1,t2`` row each)."""
    out: dict[int, list[int]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1]:
                out[int(parts[0])] = [int(t) for t in parts[1].split(",") if t.isdigit()]
    return out


def read_viruses(
    lines: Iterable[str],
    ec_transcripts: Mapping[int, list[int]],
    transcripts: list[str],
    transcript_virus: Mapping[str, str],
) -> dict[int, set[str]]:
    """``{read number: viruses its class touches}`` from ``bustools text -f`` rows.

    A read whose class holds no viral transcript is absent.
    """
    out: dict[int, set[str]] = {}
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        ec, flag = int(fields[2]), int(fields[4])
        names = {
            transcript_virus[transcripts[i]]
            for i in ec_transcripts.get(ec, [])
            if i < len(transcripts) and transcripts[i] in transcript_virus
        }
        if names:
            out[flag] = names
    return out


def summarise(
    viral_reads: Mapping[int, set[str]], virus_only_reads: Mapping[int, set[str]]
) -> list[dict[str, Any]]:
    """Per-virus read counts. A virus-only read is a virus-touching read the complement kept."""
    total: dict[str, int] = {}
    only: dict[str, int] = {}
    for flag, names in viral_reads.items():
        for name in names:
            total[name] = total.get(name, 0) + 1
            if flag in virus_only_reads:
                only[name] = only.get(name, 0) + 1
    return [
        {
            "virus_name": name,
            "reads_total": total[name],
            "reads_virus_only": only.get(name, 0),
            "reads_also_host": total[name] - only.get(name, 0),
            "host_capture_status": MEASURED,
        }
        for name in sorted(total)
    ]


def not_applicable_rows(identity: VirusIdentityTable) -> list[dict[str, Any]]:
    """One row per virus with every count ``not_applicable``: the index has no host to split on."""
    names = sorted({g.virus_name for g in identity.genes if g.viral and g.virus_name})
    return [
        {
            "virus_name": n,
            "reads_total": NOT_APPLICABLE,
            "reads_virus_only": NOT_APPLICABLE,
            "reads_also_host": NOT_APPLICABLE,
            "host_capture_status": "not_applicable_no_host_transcripts_in_index",
        }
        for n in names
    ]


def write_tsv(path: Union[str, Path], rows: list[dict[str, Any]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run(ctx: RunContext, identity: VirusIdentityTable, tsv_path: str, done_path: str) -> None:
    """Snakemake entry point for one Run."""
    config = ctx.config
    if not config.index:
        raise RuntimeError("config.yaml has no kallisto index; host_capture cannot replay.")
    transcripts, t2g_map = load_transcripts(str(ctx.outputs.transcripts_txt), config.transcripts)  # type: ignore[no-untyped-call]
    by_gene = identity.by_gene()
    transcript_virus = {
        tx: by_gene[gene].virus_name
        for tx, gene in t2g_map.items()
        if gene in by_gene and by_gene[gene].viral and by_gene[gene].virus_name
    }
    host = [tx for tx in transcripts if tx not in transcript_virus]
    viral = [tx for tx in transcripts if tx in transcript_virus]
    if not host or not viral:
        reason = "no host transcripts" if not host else "no viral transcripts"
        log.warning("host_capture: the index holds %s; writing not_applicable.", reason)
        write_tsv(tsv_path, not_applicable_rows(identity))
        Path(done_path).write_text(f"not_applicable: {reason} in the index\n")
        return
    r1, r2 = replay_fastqs(config)
    work = Path(config.output) / "host_capture"
    work.mkdir(parents=True, exist_ok=True)
    strand = getattr(config, "strand", None)
    viral_text, only_text = replay_host_capture_bus(
        index=config.index,
        technology=config.technology,
        r1_path=r1,
        r2_path=r2,
        viral_transcripts=viral,
        host_transcripts=host,
        workdir=str(work),
        threads=int(config.cores),
        whitelist=_replay_whitelist(config, ctx.outputs),
        strand=None if strand == "auto" else strand,
    )
    # kallisto numbers classes in discovery order: use the replay's own EC map and transcript list.
    replay = work / "lineage_bus"
    replay_transcripts = (replay / "transcripts.txt").read_text().split("\n")
    replay_transcripts = [t for t in replay_transcripts if t]
    ec_tx = read_ec_transcripts(replay / "matrix.ec")
    with open(viral_text, encoding="utf-8") as fh:
        viral_reads = read_viruses(fh, ec_tx, replay_transcripts, transcript_virus)
    with open(only_text, encoding="utf-8") as fh:
        only_reads = read_viruses(fh, ec_tx, replay_transcripts, transcript_virus)
    rows = summarise(viral_reads, only_reads)
    write_tsv(tsv_path, rows)
    Path(done_path).write_text(f"{MEASURED}: {len(rows)} viruses\n")


if "snakemake" in globals():
    _ctx = RunContext.from_yaml(snakemake.params.configfile).with_cores(snakemake.threads)  # noqa: F821
    _identity = VirusIdentityTable.read_tsv(snakemake.input.virus_identity)  # noqa: F821
    run(_ctx, _identity, snakemake.output.tsv, snakemake.output.done)  # noqa: F821
