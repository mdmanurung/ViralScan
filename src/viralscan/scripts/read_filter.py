"""
Opt-in read-artefact filter before ``kb count`` (PLAN DEF-01, R2.0).

Runs only with ``--read-filter artefact``. It reads the pair that ``kb count``
would otherwise see (the host-filtered pair when ``--host-filter`` is on) and
drops a fragment for the first of these reasons that applies:

``r1_tso``
    The 10x TSO sits where R1's barcode and UMI should be: a chimera by
    construction (``anello_align.has_r1_tso`` on the CB+UMI span only).
``r2_reagent``
    R2 carries reagent in a shape no genuine molecule has: ``TSO|poly-T`` or a
    forward TruSeq R1 chimera (``anello_align.has_reagent``).
``r2_no_complex_body``
    R2 has no templated body of >= 20 nt and >= 2.0 bits before its first
    homopolymer run (``anello_align.is_complex_body``). A genuine
    ``[viral body][poly-A]`` read keeps its body and is retained.

The classifiers are the F-019 measures that ``anello_align`` reports as labels.
The filter is off by default because a low-complexity viral body is removed
too (``.living/decisions.md``, 2026-10-04 and 2026-10-05). The reagent checks
use 10x oligos; on other chemistries they rarely fire, and the audit shows
their counts.

Output
    {output}read_filtered/R1.fastq.gz, R2.fastq.gz   retained pairs, input order
    {output}read_filtered/read_filter_audit.tsv      counts per reason + pinned constants
    {output}read_filtered/fragment_lineage.tsv.gz    one row per input fragment
"""

import csv
import gzip
from pathlib import Path

from viralscan import anello_align as aa
from viralscan.chemistry import cb_umi_geometry
from viralscan.evidence import _open_maybe_gzip
from viralscan.runconfig import RunConfig
from viralscan.scripts.host_filter import canonical_read_id
from viralscan.utils import setup_script_logging

log = setup_script_logging()

#: Fixed reason order: a pair that fails several checks gets the first one.
REASONS: tuple[str, ...] = ("r1_tso", "r2_reagent", "r2_no_complex_body")

#: Constants that decide the split, written to the audit.
PINNED: tuple[tuple[str, object], ...] = (
    ("homopolymer_run", aa.HOMOPOLYMER_RUN),
    ("min_body_len", aa.MIN_BODY_LEN),
    ("min_body_entropy", aa.MIN_BODY_ENTROPY),
    ("tso", aa.TSO),
    ("tso_max_mismatch", aa.TSO_MAX_MISMATCH),
    ("truseq_r1", aa.TRUSEQ_R1),
)

#: Speed over size: kb count reads the files once.
GZIP_LEVEL = 1


def classify(r1: str, r2: str, cb_umi_len: int) -> str:
    """Return the removal reason for a pair, or ``""`` to keep it."""
    if aa.has_r1_tso(r1[:cb_umi_len]):
        return "r1_tso"
    if aa.has_reagent(r2):
        return "r2_reagent"
    if not aa.is_complex_body(r2):
        return "r2_no_complex_body"
    return ""


def filter_pairs(r1_in: str, r2_in: str, out_dir: Path, technology: str) -> dict[str, int]:
    """Stream the pair, write retained records, the lineage and the audit."""
    # ponytail: one core, about 80 us per pair (1.15 core-h per 50 M pairs);
    # an ordered multiprocessing imap over chunks if wall time matters.
    cb_len, umi_len = cb_umi_geometry(technology)
    counts = dict.fromkeys(("input", "retained", *REASONS), 0)
    with (
        _open_maybe_gzip(r1_in) as fq1,
        _open_maybe_gzip(r2_in) as fq2,
        gzip.open(out_dir / "R1.fastq.gz", "wt", compresslevel=GZIP_LEVEL) as out1,
        gzip.open(out_dir / "R2.fastq.gz", "wt", compresslevel=GZIP_LEVEL) as out2,
        gzip.open(
            out_dir / "fragment_lineage.tsv.gz", "wt", newline="", compresslevel=GZIP_LEVEL
        ) as lineage,
    ):
        writer = csv.writer(lineage, delimiter="\t")
        writer.writerow(["read_id", "filter_decision", "reason"])
        while True:
            rec1 = [fq1.readline() for _ in range(4)]
            rec2 = [fq2.readline() for _ in range(4)]
            if not rec1[0] and not rec2[0]:
                break
            if not all(rec1) or not all(rec2):
                raise ValueError(f"Truncated or unpaired FASTQ record in {r1_in!r} / {r2_in!r}")
            read_id = canonical_read_id(rec1[0])
            if read_id != canonical_read_id(rec2[0]):
                raise ValueError(f"FASTQ mate mismatch: {rec1[0].strip()!r} vs {rec2[0].strip()!r}")
            counts["input"] += 1
            reason = classify(rec1[1].rstrip(), rec2[1].rstrip(), cb_len + umi_len)
            if reason:
                counts[reason] += 1
                writer.writerow([read_id, "removed", reason])
            else:
                counts["retained"] += 1
                out1.writelines(rec1)
                out2.writelines(rec2)
                writer.writerow([read_id, "retained", "complex_body"])
    _write_audit(out_dir, counts)
    return counts


def _write_audit(out_dir: Path, counts: dict[str, int]) -> None:
    with (out_dir / "read_filter_audit.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "fragments", "interpretation"])
        writer.writerow(["input", counts["input"], "paired fragments presented to the filter"])
        writer.writerow(["retained", counts["retained"], "reached kb count"])
        for reason in REASONS:
            writer.writerow([f"removed_{reason}", counts[reason], "first failing check"])
        pct = 100.0 * counts["retained"] / counts["input"] if counts["input"] else 0.0
        writer.writerow(["pct_retained", f"{pct:.2f}", "fraction of input reaching kb count"])
        for name, value in PINNED:
            writer.writerow([f"param:{name}", value, "pinned filter constant"])


# ── Entry point ───────────────────────────────────────────────────────────────
def main(config: RunConfig, r1_in: str, r2_in: str, done_path: str) -> None:
    out_dir = Path(config.output) / "read_filtered"
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("Read-artefact filter (%s): %s, %s", config.read_filter, r1_in, r2_in)
    counts = filter_pairs(r1_in, r2_in, out_dir, config.technology)
    log.info(
        "Read-artefact filter kept %d / %d pairs (r1_tso %d, r2_reagent %d, r2_no_complex_body %d).",
        counts["retained"],
        counts["input"],
        counts["r1_tso"],
        counts["r2_reagent"],
        counts["r2_no_complex_body"],
    )
    Path(done_path).touch()


# ── Snakemake wiring (only runs under snakemake) ─────────────────────────────
if "snakemake" in globals():
    main(
        config=RunConfig.from_yaml(snakemake.params.configfile),  # noqa: F821
        r1_in=str(snakemake.input.r1),  # noqa: F821
        r2_in=str(snakemake.input.r2),  # noqa: F821
        done_path=str(snakemake.output.done),  # noqa: F821
    )
