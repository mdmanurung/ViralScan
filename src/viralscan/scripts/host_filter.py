"""
Optional host-subtraction pre-step for ViralScan.

Runs before ``kb_count`` when the user supplies ``--host-filter`` and
``--host-index``.  The stable v3 workflow supports one mode:

starsolo
    STAR genome alignment (STARsolo barcode-aware mode).  Reads that do NOT
    align to the host genome are written as FASTQ by STAR via
    ``--outReadsUnmapped Fastx`` and collected as the filtered output.
    Requires a STAR genome directory built with ``STAR --runMode genomeGenerate``.

The former kallisto subtraction mode is intentionally unavailable in v3.
Kallisto BUS output does not retain exact source read identifiers; subtracting
all fragments that share a mapped CB–UMI can delete unrelated viral fragments.

Output
    {output}host_filtered/R1.fastq.gz   — barcode+UMI read (R1 in 10x convention)
    {output}host_filtered/R2.fastq.gz   — cDNA read (R2 in 10x convention)

These paths are pre-registered in config["kb_r1"] / config["kb_r2"] by
``createconfig.py``, so the downstream ``kb_count`` rule consumes them
transparently with no further changes.

CB/UMI geometry is resolved via ``viralscan.evidence.cb_umi_geometry``, which
handles 10x v1/v2/v3, Drop-seq, and explicit kallisto ``bc:umi:seq`` triplets,
and raises ``ValueError`` for unknown chemistries instead of silently
mis-slicing barcodes (fixes PLAN S1/S6).
"""

import csv
import gzip
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from viralscan.evidence import _open_maybe_gzip, cb_umi_geometry
from viralscan.runconfig import RunConfig
from viralscan.utils import setup_script_logging

log = setup_script_logging()


def required_host_filter_tools(aligner: str) -> tuple[str, ...]:
    """Return native executables required by a host-filter mode."""
    if aligner == "kallisto":
        raise ValueError(
            "--host-filter kallisto is not available in ViralScan v3: kallisto BUS "
            "output does not preserve exact read identifiers, so safe fragment-level "
            "subtraction cannot be guaranteed. Use --host-filter starsolo."
        )
    if aligner == "starsolo":
        return ("STAR",)
    raise ValueError(f"Unknown host_filter_aligner: {aligner!r}. Choose 'starsolo'.")


def check_host_filter_tools(aligner: str) -> None:
    """Fail before launching Snakemake work if native tools are unavailable."""
    missing = [tool for tool in required_host_filter_tools(aligner) if shutil.which(tool) is None]
    if missing:
        tools = ", ".join(missing)
        raise RuntimeError(
            f"--host-filter {aligner} requires {tools} on PATH. "
            "Install the ViralScan conda/container environment with STAR before running "
            "host filtering."
        )


# ── Helper: gzip-copy a plain-text file to a .gz destination ─────────────────
def _gzip_file(src: Path, dst: str) -> None:
    with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def canonical_read_id(header: str) -> str:
    """Return the fragment identifier used to compare paired FASTQ records."""
    token = header.strip().split(maxsplit=1)[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith(("/1", "/2")):
        token = token[:-2]
    return token


def iter_paired_fastq_ids(r1_path: str, r2_path: str):
    """Yield synchronized fragment IDs, failing on truncation or mate mismatch."""
    with _open_maybe_gzip(r1_path) as fq1, _open_maybe_gzip(r2_path) as fq2:
        record = 0
        while True:
            lines1 = [fq1.readline() for _ in range(4)]
            lines2 = [fq2.readline() for _ in range(4)]
            if not lines1[0] and not lines2[0]:
                return
            if not lines1[0] or not lines2[0]:
                raise ValueError(
                    f"Paired FASTQs have different record counts: {r1_path!r}, {r2_path!r}"
                )
            if not all(lines1) or not all(lines2):
                raise ValueError(f"Truncated FASTQ record in {r1_path!r} or {r2_path!r}")
            record += 1
            id1 = canonical_read_id(lines1[0])
            id2 = canonical_read_id(lines2[0])
            if id1 != id2:
                raise ValueError(f"FASTQ mate mismatch at record {record}: {id1!r} != {id2!r}")
            yield id1


def validate_paired_fastq_ids(r1_path: str, r2_path: str) -> int:
    """Validate pair synchronization and return the number of fragments."""
    return sum(1 for _ in iter_paired_fastq_ids(r1_path, r2_path))


def filter_fastq_pairs(*_args, **_kwargs):
    """Reject the pre-v3 CB–UMI-wide subtraction API."""
    raise RuntimeError(
        "CB–UMI-wide FASTQ subtraction was removed in ViralScan v3 because it can "
        "delete unrelated fragments. Use exact read identifiers via STAR host filtering."
    )


def _write_filter_audit(
    out_dir: Path,
    original_pairs: int,
    retained_pairs: int,
    filtered_r1: str,
    filtered_r2: str,
    input_r1: str,
    input_r2: str,
    param_set: str = "pinned",
) -> None:
    """Write aggregate filtering counts and per-fragment lineage (retained and removed)."""
    star_args = STAR_FILTER_PARAM_SETS[param_set]
    removed_pairs = original_pairs - retained_pairs
    with (out_dir / "host_filter_audit.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "fragments", "interpretation"])
        writer.writerow(["input", original_pairs, "paired fragments presented to STAR"])
        writer.writerow(["retained_host_unmapped", retained_pairs, "STAR unmapped mate pair"])
        writer.writerow(
            [
                "removed_host_aligned_or_ambiguous",
                removed_pairs,
                "not emitted by STAR --outReadsUnmapped; alignment subclass unavailable",
            ]
        )
        writer.writerow(
            [
                "pct_retained",
                f"{100.0 * retained_pairs / original_pairs:.2f}" if original_pairs else "0.00",
                "fraction of input reaching viral quant",
            ]
        )
        # Pin the parameters that decided the split, so a re-run that removes a
        # different fraction can be attributed to a parameter change rather than
        # to a STAR version difference.
        writer.writerow(["star_param_set", param_set, "--host-filter-star-params"])
        writer.writerow(
            [
                "cell_calling_matrix",
                str(out_dir / "star_tmp" / "Solo.out" / "GeneFull" / "raw"),
                "STARsolo host GeneFull raw counts for two-step EmptyDrops",
            ]
        )
        for i in range(0, len(star_args), 2):
            writer.writerow(
                [
                    f"star_param:{star_args[i].lstrip('-')}",
                    star_args[i + 1],
                    "STAR filter parameter",
                ]
            )

    with gzip.open(out_dir / "fragment_lineage.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["read_id", "filter_decision", "reason"])
        # STAR runs with --outSAMtype None, so removed = input IDs - retained IDs.
        # Only the retained side (the small one in practice) is held in memory;
        # input IDs are streamed in input order.
        # ponytail: memory is O(retained fragments); switch to a sorted merge if
        # retained ever approaches the input size.
        retained = set(iter_paired_fastq_ids(filtered_r1, filtered_r2))
        for read_id in iter_paired_fastq_ids(input_r1, input_r2):
            if read_id in retained:
                writer.writerow([read_id, "retained", "host_unmapped"])
            else:
                writer.writerow([read_id, "removed", "host_mapped"])


def lost_truth_counts(truth_tsv: str, lineage_path: str, viral_label: str = "viral") -> dict:
    """D15/D16-style loss of truth-viral material at the host-filter boundary.

    *truth_tsv* has a header with ``read_id`` and ``label`` columns, plus an
    optional ``molecule_id`` (canonical CB-UMI truth key). Returns:

    ``d15_truth_fragments`` / ``d15_removed_fragments`` / ``d15_loss_fraction``
        truth-viral fragments entering the boundary and those removed by it.
    ``d16_truth_molecules`` / ``d16_lost_molecules`` / ``d16_loss_fraction``
        truth molecules with no surviving fragment (``None`` without the
        ``molecule_id`` column).
    ``not_in_lineage``
        truth-viral reads the lineage never saw; non-zero means the row is not
        estimable (protocol: ``not-estimable-and-row-failed``), so the fractions
        are then ``None``.
    """
    with gzip.open(lineage_path, "rt", newline="") as handle:
        decision = {
            r["read_id"]: r["filter_decision"] for r in csv.DictReader(handle, delimiter="\t")
        }
    fragments = removed = missing = 0
    molecules: dict[str, bool] = {}
    has_molecule = False
    with open(truth_tsv, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        has_molecule = "molecule_id" in (reader.fieldnames or ())
        for row in reader:
            if row["label"] != viral_label:
                continue
            fragments += 1
            state = decision.get(row["read_id"])
            if state is None:
                missing += 1
                continue
            removed += state == "removed"
            if has_molecule:
                molecules[row["molecule_id"]] = molecules.get(row["molecule_id"], False) or (
                    state == "retained"
                )
    ok = missing == 0
    lost = sum(1 for survived in molecules.values() if not survived)
    return {
        "d15_truth_fragments": fragments,
        "d15_removed_fragments": removed,
        "d15_loss_fraction": removed / fragments if ok and fragments else None,
        "d16_truth_molecules": len(molecules) if has_molecule else None,
        "d16_lost_molecules": lost if has_molecule else None,
        "d16_loss_fraction": lost / len(molecules) if ok and molecules else None,
        "not_in_lineage": missing,
    }


#: Every STAR parameter the host filter depends on, pinned explicitly. This is
#: the ``pinned`` set of ``--host-filter-star-params`` and the CLI default.
#:
#: Before this existed the command set *no* alignment or filter options, so the
#: run inherited whatever the installed STAR defaulted to. Defaults below are
#: from ``STAR --help`` (2.7.4a and 2.7.11b agree); the pinned values differ
#: from them as follows:
#:
#:   ``--outFilterMismatchNmax`` 10 (STAR default) -> 4
#:       Tighter than the default, so a read with 5-10 mismatches is *not*
#:       called host and survives into viral quantification. Chosen together
#:       with ``--outFilterMismatchNoverReadLmax 0.05`` (default 1.0) to demand
#:       a near-identical host alignment before removing a read.
#:
#:   ``--outFilterMatchNminOverLread`` 0.66 (STAR default) -> 0.9
#:       The default accepts an alignment covering only two thirds of the read;
#:       0.9 requires near-full-length host alignment. Again stricter: partially
#:       host-aligned reads are kept for the viral index.
#:
#:   ``--outFilterMultimapNmax`` 10 (STAR default) -> 20
#:       A read mapping to more loci than this limit is not output as aligned,
#:       so it is reported as unmapped ("mapped to too many loci") and survives
#:       into viral quantification. Raising the limit to 20 lets reads with
#:       11-20 host loci be called host.
#:
#: The remaining pins (``--alignIntronMin 20`` vs default 21,
#: ``--alignSJoverhangMin 8`` vs 5, ``--alignSJDBoverhangMin 8`` vs 3,
#: ``--outSAMattributes None`` vs Standard) are small splice/output tweaks;
#: ``--outFilterType Normal`` and ``--outSAMprimaryFlag OneBestScore`` equal
#: the STAR defaults and are pinned only for reproducibility.
#:
#: ``--outSAMtype None`` is kept: it is what makes ``--outReadsUnmapped`` cheap.
#: The cost is that the alignment subclass of every removed fragment is lost, so
#: the audit can only report aggregate counts (see ``_write_filter_audit``).
STAR_FILTER_ARGS: tuple[str, ...] = (
    "--outFilterType",
    "Normal",
    "--outFilterMismatchNmax",
    "4",
    "--outFilterMatchNminOverLread",
    "0.9",
    "--outFilterMultimapNmax",
    "20",
    "--outFilterMismatchNoverReadLmax",
    "0.05",
    "--alignIntronMin",
    "20",
    "--alignSJoverhangMin",
    "8",
    "--alignSJDBoverhangMin",
    "8",
    "--outSAMattributes",
    "None",
    "--outSAMprimaryFlag",
    "OneBestScore",
)

#: The same flags at the STAR 2.7.11b ``--help`` defaults, still passed
#: explicitly so a different STAR install cannot change them. Against
#: ``pinned`` it calls more reads host on mismatches and partial coverage, but
#: fewer on multimapping: reads with 11-20 host loci go unmapped and survive.
STAR_DEFAULT_FILTER_ARGS: tuple[str, ...] = (
    "--outFilterType",
    "Normal",
    "--outFilterMismatchNmax",
    "10",
    "--outFilterMatchNminOverLread",
    "0.66",
    "--outFilterMultimapNmax",
    "10",
    "--outFilterMismatchNoverReadLmax",
    "1.0",
    "--alignIntronMin",
    "21",
    "--alignSJoverhangMin",
    "5",
    "--alignSJDBoverhangMin",
    "3",
    "--outSAMattributes",
    "Standard",
    "--outSAMprimaryFlag",
    "OneBestScore",
)

#: ``--host-filter-star-params`` choices (the protocol's star_alignment grid).
STAR_FILTER_PARAM_SETS: dict[str, tuple[str, ...]] = {
    "pinned": STAR_FILTER_ARGS,
    "star-default": STAR_DEFAULT_FILTER_ARGS,
}


def starsolo_barcode_args(technology: str, whitelist: Optional[str]) -> list[str]:
    """STARsolo barcode/UMI arguments for the host-filter pass (pure; SW-15).

    ``--soloBarcodeReadLength 0`` disables STAR's check that the barcode read is
    exactly CB+UMI long. 10x 5' libraries ship a 150 bp R1 (barcode, UMI, TSO,
    then cDNA), and without it STAR aborts with "barcode read length 150 not
    equal to expected 28". *whitelist* must already be plain text; ``None``
    passes STAR's literal ``None`` (accept every barcode).
    """
    cb_len, umi_len = cb_umi_geometry(technology)
    return [
        "--soloType",
        "CB_UMI_Simple",
        "--soloCBstart",
        "1",
        "--soloCBlen",
        str(cb_len),
        "--soloUMIstart",
        str(cb_len + 1),
        "--soloUMIlen",
        str(umi_len),
        "--soloBarcodeReadLength",
        "0",
        "--soloCBwhitelist",
        whitelist or "None",
    ]


def _plain_whitelist(whitelist: Optional[str], work_dir: Path) -> Optional[str]:
    """Return a plain-text copy of a gzipped *whitelist*; STAR cannot read .gz."""
    if not whitelist or not whitelist.endswith(".gz"):
        return whitelist
    plain = work_dir / "whitelist.txt"
    with gzip.open(whitelist, "rb") as source, open(plain, "wb") as target:
        shutil.copyfileobj(source, target)
    return str(plain)


# ── STARsolo mode ─────────────────────────────────────────────────────────────
def _starsolo_filter(
    r1: str,
    r2: str,
    host_index: str,
    technology: str,
    whitelist: Optional[str],
    out_dir: Path,
    filtered_r1: str,
    filtered_r2: str,
    n_threads: int,
    param_set: str = "pinned",
) -> None:
    """Run STARsolo with ``--outReadsUnmapped Fastx``.

    STAR writes unmapped mates to:
      Unmapped.out.mate1  — the *first* file given to --readFilesIn (= R2/cDNA in 10x)
      Unmapped.out.mate2  — the *second* file (= R1/barcode+UMI in 10x)

    In the 10x convention we pass R2 first and R1 second so that STAR's
    soloType CB_UMI_Simple knows which read carries the barcode.
    """
    original_pairs = validate_paired_fastq_ids(r1, r2)
    star_tmp = out_dir / "star_tmp"
    star_tmp.mkdir(exist_ok=True)

    read_files_cmd = "zcat" if r1.endswith(".gz") or r2.endswith(".gz") else "-"

    cmd = [
        "STAR",
        "--runThreadN",
        str(n_threads),
        "--genomeDir",
        host_index,
        # 10x convention: cDNA read (R2) first, barcode+UMI read (R1) second
        "--readFilesIn",
        r2,
        r1,
        "--readFilesCommand",
        read_files_cmd,
        *starsolo_barcode_args(technology, _plain_whitelist(whitelist, star_tmp)),
        "--soloFeatures",
        "GeneFull",
        "--outSAMtype",
        "None",
        "--outReadsUnmapped",
        "Fastx",
        "--outFileNamePrefix",
        str(star_tmp) + os.sep,
    ]
    cmd += list(STAR_FILTER_PARAM_SETS[param_set])

    log.info("Running STARsolo host filter...")
    subprocess.run(cmd, check=True)

    # Mate1 = cDNA (→ R2), Mate2 = barcode+UMI (→ R1)
    unmapped_cdna = star_tmp / "Unmapped.out.mate1"
    unmapped_bc = star_tmp / "Unmapped.out.mate2"

    log.info("Compressing unmapped reads → %s", out_dir)
    _gzip_file(unmapped_bc, filtered_r1)  # barcode+UMI → R1
    _gzip_file(unmapped_cdna, filtered_r2)  # cDNA       → R2

    retained_pairs = validate_paired_fastq_ids(filtered_r1, filtered_r2)
    if retained_pairs > original_pairs:
        raise RuntimeError(
            f"STAR returned more pairs than it received: {retained_pairs} > {original_pairs}"
        )
    _write_filter_audit(
        out_dir, original_pairs, retained_pairs, filtered_r1, filtered_r2, r1, r2, param_set
    )
    pct = 100.0 * retained_pairs / original_pairs if original_pairs else 0.0
    log.info(
        "STARsolo host filter complete: kept %d / %d read pairs (%.1f%% passed).",
        retained_pairs,
        original_pairs,
        pct,
    )


# ── Entry point ───────────────────────────────────────────────────────────────
def main(config: RunConfig, n_threads: int, done_path: str) -> None:
    output = config.output
    aligner = config.host_filter_aligner or "starsolo"
    host_index = config.host_index
    r1 = config.sample1
    r2 = config.sample2
    technology = config.technology
    whitelist = config.whitelist

    out_dir = Path(output) / "host_filtered"
    out_dir.mkdir(parents=True, exist_ok=True)
    filtered_r1 = str(out_dir / "R1.fastq.gz")
    filtered_r2 = str(out_dir / "R2.fastq.gz")

    param_set = config.host_filter_star_params
    log.info(
        "Host pre-subtraction: aligner=%s, host_index=%s, star_params=%s",
        aligner,
        host_index,
        param_set,
    )
    check_host_filter_tools(aligner)

    assert host_index is not None, "host_filter runs only when host_index is set"
    if aligner == "starsolo":
        _starsolo_filter(
            r1,
            r2,
            host_index,
            technology,
            whitelist,
            out_dir,
            filtered_r1,
            filtered_r2,
            n_threads,
            param_set,
        )
    else:
        required_host_filter_tools(aligner)

    Path(done_path).touch()


# ── Snakemake wiring (only runs under snakemake) ─────────────────────────────
if "snakemake" in globals():
    main(
        config=RunConfig.from_yaml(snakemake.params.configfile),  # noqa: F821
        n_threads=snakemake.threads,  # noqa: F821
        done_path=str(snakemake.output.done),  # noqa: F821
    )
