#!/usr/bin/env python3
"""
Build the ViralScan bundled-panel reference index — one-time setup.

Builds a combined HOST + VIRUS kallisto index for joint quantification.
Human cDNA is downloaded from Ensembl via fetch_host_cdna; viral FASTAs are
fetched from NCBI via ncbi_fetch._fetch_one.  Both host and viral transcripts
are quantified together, enabling coexpression analysis between host genes
(ENST*) and viral genes (ADENO_*, EPSTEIN_*, etc.).

Usage:
    NCBI_EMAIL=you@example.com PYTHONPATH=/path/to/ViralScan/src \\
        python scripts/build_bundled_panel_ref.py \\
        --out /exports/para-lipg-hpc/mdmanurung/viralscan_bulk_gse128078/ref

SLURM (activate test_viralscan conda env before sbatch):
    sbatch --job-name=build_panel_ref --cpus-per-task=2 --mem=16G --time=08:00:00 \\
        --wrap "NCBI_EMAIL=... PYTHONPATH=$PWD/src python $PWD/scripts/build_bundled_panel_ref.py --out $WORKDIR/ref"

Pre-conditions verified at startup:
  1. Bundled GTFs contain exon features (kb ref extracts cDNA from exon rows;
     missing exons yield an empty/partial index that wc -l alone won't catch).
  2. Every GTF seqname is a valid NCBI accession (flag non-matching names early).
  3. Every GTF seqname has a matching FASTA header after download (a suppressed
     or missing NCBI accession would otherwise let kb ref silently drop viruses).
  4. Every accession the virus catalogue claims reaches the assembled panel
     FASTA, or the miss carries a recorded decision in
     src/viralscan/data/index_exclusions.tsv.  Unexplained misses are written to
     catalogued_not_indexed.tsv and fail the build under --strict-reconciliation
     (PLAN CAT-31, finding F-015).  Runs before kb ref, so a panel that would
     ship a lower sensitivity bound fails in seconds instead of after ~64 GB and
     ~8 h of indexing.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


def _find_repo_root() -> Path:
    """Walk up from this script to find the repo root (contains src/viralscan/)."""
    here = Path(__file__).resolve().parent
    for candidate in [here.parent, here]:
        if (candidate / "src" / "viralscan").is_dir():
            return candidate
    sys.exit(
        "ERROR: cannot locate ViralScan repo root (src/viralscan/ not found). "
        "Run from inside the ViralScan repo."
    )


_ACC_RE = re.compile(r"^[A-Za-z]{1,3}_?\d+(\.\d+)?$")


def _extract_gtf_seqnames(gtf_files: list[Path]) -> tuple[set[str], set[str]]:
    """Return (all_seqnames, feature_types) from the bundled GTF corpus."""
    all_seqnames: set[str] = set()
    feature_types: set[str] = set()
    for gtf in gtf_files:
        for line in gtf.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            all_seqnames.add(cols[0])
            feature_types.add(cols[2])
    return all_seqnames, feature_types


def _dustmask_dir() -> str | None:
    for cand in ("dustmasker", "dustmasker_ng"):
        found = shutil.which(cand)
        if found:
            return found
    return None


def _dustmask_fasta_file(target: Path, level: int, windows: tuple[int, ...]) -> None:
    """Hardmask a FASTA in place via NCBI dustmasker, merging several windows.

    The CAT-17 gate requires a masked panel: the deployed reference was 99.99%
    unmasked yet carried 170 pure-homopolymer 31-mers, and a poly-A/poly-T 10x
    tail matches one of those exactly (F-014).  A kallisto D-list cannot fix this
    -- it filters host-homologous k-mers, not self-similarity inside a viral
    contig.  Masking must therefore happen before the gate can ever pass.

    One pass is not enough.  The gate scores 31-mers, so a homopolymer run that
    dustmasker's default window (30) tolerates can still supply 31 consecutive
    identical bases.  Masking therefore runs at each window in *windows* and the
    results are merged as a union of masked positions, which is what the gate's
    own failure message prescribes ("windows 64 and 30, merged, masked to N").
    """
    binary = _dustmask_dir()
    if binary is None:
        sys.exit(
            "ERROR: dustmasker not on PATH. The CAT-17 low-complexity gate requires a "
            "masked panel. Install blast (provides dustmasker) or pass "
            "--no-dustmask to build unmasked and accept the gate failure."
        )
    original = target.with_suffix(target.suffix + ".premask")
    target.replace(original)
    passes: list[Path] = []
    for w in windows:
        out_path = target.with_suffix(target.suffix + f".w{w}")
        subprocess.run(
            [
                binary,
                "-infmt",
                "fasta",
                "-in",
                str(original),
                "-outfmt",
                "fasta",
                "-level",
                str(level),
                "-window",
                str(w),
                "-out",
                str(out_path),
            ],
            check=True,
            capture_output=True,
        )
        passes.append(out_path)
    if len(passes) == 1:
        passes[0].replace(target)
    else:
        _merge_masked_fastas(original, passes, target)
        for p in passes:
            p.unlink()
    original.unlink()


def _mask_homopolymer_runs(target: Path, run_length: int) -> int:
    """Replace every run of >= *run_length* identical bases with N, in place.

    dustmasker is not sufficient for the class the CAT-17 gate actually fails on.
    Measured on this panel it masked 785 of 11,048,660 bases (0.01 %) and left all
    85 pure-homopolymer 31-mers in NC_001479.1 intact: DUST scores compositional
    complexity over a window and is not designed to strip a homopolymer tract, so
    the poly-A/poly-T 10x tail that manufactured 1.44 % of R2 reads (F-014) passes
    straight through it.

    A "pure homopolymer 31-mer" is by definition a run of >= 31 identical bases, but
    masking only at the k length is not enough (F-021, PLAN CAT-42, 2026-10-01). A
    29-T head on KP343824.1 still yields 31-mers of 29 T plus 2 flank bases, which
    host poly-A tails match: 10,037 molecules in an anellovirus-free cell line, with
    0 reads on any other viral 31-mer. Every Gamma/Betatorquevirus molecule of that
    run landed on a genome with a >= 20-nt run, so the default is 20: a k-mer then
    needs >= 11 specific flank bases. The cost is bounded and auditable: 38 runs in
    38 of 2,343 genomes, mostly poly-A tracts and genomic poly-A tails.

    Returns the number of bases masked.
    """
    records = _read_fasta(target)
    masked_total = 0
    out: list[tuple[str, str]] = []
    for name, seq in records:
        pieces: list[str] = []
        i = 0
        n = len(seq)
        while i < n:
            j = i + 1
            while j < n and seq[j] == seq[i]:
                j += 1
            length = j - i
            if length >= run_length and seq[i] != "N":
                pieces.append("N" * length)
                masked_total += length
            else:
                pieces.append(seq[i:j])
            i = j
        out.append((name, "".join(pieces)))
    if masked_total:
        with open(target, "w") as fh:
            for name, seq in out:
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i : i + 60] + "\n")
    return masked_total


def _anellovirus_accessions() -> set[str]:
    """Versioned and bare accessions of every Anelloviridae record the package knows."""
    from viralscan.anellovirus import load_accession_table
    from viralscan.virus_catalog import accession_metadata

    accs = {a for a, row in accession_metadata().items() if row.get("family") == "Anelloviridae"}
    accs.update(r["accession"].strip() for r in load_accession_table())
    return accs


def _lowcomplexity_window_mask(seq: str, window: int = 31) -> list[bool]:
    """Per-base mask for windows dominated by one or two bases (CAT-42).

    A window is masked whole when it has <= 2 distinct ACGT bases and a base count
    >= 20, or any base count >= 28.  Catches tracts a pure-run mask misses because
    interruptions break the run (KP343822.1 A-rich tract, KP343842.1 G/C tail).
    """
    n = len(seq)
    mask = [False] * n
    if n < window:
        return mask
    cum = {}
    for b in "ACGT":
        c = [0]
        for ch in seq:
            c.append(c[-1] + (ch == b))
        cum[b] = c
    for i in range(n - window + 1):
        counts = [cum[b][i + window] - cum[b][i] for b in "ACGT"]
        top = max(counts)
        if top >= 28 or (top >= 20 and sum(1 for c in counts if c) <= 2):
            for j in range(i, i + window):
                mask[j] = True
    return mask


def _mask_lowcomplexity_kmers(
    target: Path, accessions: set[str], window: int = 31
) -> dict[str, int]:
    """N-mask low-complexity windows in the *accessions* records of *target*, in place.

    Anelloviridae only: applied panel-wide this rule masks real herpes gene sequence
    (EBNA-2, HSV-1 s-genes).  Returns {accession: newly masked bases} for records
    with any.
    """
    out: list[tuple[str, str]] = []
    masked: dict[str, int] = {}
    for name, seq in _read_fasta(target):
        acc = name.split()[0]
        if acc in accessions or acc.split(".")[0] in accessions:
            m = _lowcomplexity_window_mask(seq, window)
            new = [i for i, f in enumerate(m) if f and seq[i] != "N"]
            if new:
                chars = list(seq)
                for i in new:
                    chars[i] = "N"
                seq = "".join(chars)
                masked[acc] = len(new)
        out.append((name, seq))
    if masked:
        with open(target, "w") as fh:
            for name, seq in out:
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i : i + 60] + "\n")
    return masked


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def _merge_masked_fastas(original: Path, masked_passes: list[Path], out_path: Path) -> None:
    """Write *original* with any position masked in ANY pass replaced by N."""
    base = dict(_read_fasta(original))
    merged: dict[str, str] = {}
    for pass_path in masked_passes:
        for name, seq in _read_fasta(pass_path):
            prior = merged.get(name, base.get(name, seq))
            if len(prior) != len(seq):
                continue
            merged[name] = "".join("N" if (a == "N" or b == "N") else a for a, b in zip(prior, seq))
    for name, seq in base.items():
        merged.setdefault(name, seq)
    with open(out_path, "w") as fh:
        for name, seq in merged.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def _fasta_seq_ids(fasta_paths: list[Path]) -> set[str]:
    """Collect all sequence IDs (text before first space on '>' lines) from FASTAs."""
    ids: set[str] = set()
    for fp in fasta_paths:
        for line in fp.read_text().splitlines():
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def _build_arg_parser() -> argparse.ArgumentParser:
    """Return the CLI parser, so the flag contract is testable without a build."""
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--gtf-manifest",
        type=Path,
        default=_find_repo_root() / "analysis" / "panel_expansion" / "gtf_corpus_manifest.tsv",
        help="Local GTF corpus hash manifest, checked before downloads or indexing.",
    )
    p.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output directory for panel.idx, panel.t2g, cdna.fa, combined.*",
    )
    p.add_argument(
        "--anello-star-only",
        action="store_true",
        help="Only (re)build the STAR anellovirus index anello_star/ from an existing "
        "<out>/viral.fa, skipping every other step (PLAN ANDET-09a).",
    )
    p.add_argument(
        "--star",
        default="STAR",
        help="STAR binary for the anello_star/ index (default: STAR on PATH).",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for STAR genomeGenerate (default: 4).",
    )
    p.add_argument(
        "--host-species",
        default="human",
        help="Ensembl host species for the combined host+viral reference (default: human). "
        "Host transcripts (ENST*) are quantified alongside viral ones for coexpression analysis.",
    )
    p.add_argument(
        "--ncbi-email",
        default=os.environ.get("NCBI_EMAIL"),
        help="Contact email for NCBI E-utilities (or set NCBI_EMAIL env var)",
    )
    p.add_argument(
        "--ncbi-api-key",
        default=os.environ.get("NCBI_API_KEY"),
        help="NCBI API key for higher rate limits (optional)",
    )
    p.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Override the NCBI fetch cache directory",
    )
    p.add_argument(
        "--max-low-complexity-fraction",
        type=float,
        default=0.05,
        metavar="FRAC",
        help="Per-record ceiling on the fraction of 31-mers that are low-complexity "
        "(homopolymer run > 11, fewer than 3 distinct bases, or a perfect tandem "
        "repeat of a unit <= 5 bp).  Default 0.05.  This is a backstop on the "
        "long_run/few_bases/tandem classes, NOT the primary control: a tandem repeat "
        "is legitimate sequence, and measured across this panel the worst record sits "
        "at 0.0205 (median 0.0035, p90 0.0100), so a 0.0 ceiling is unachievable by "
        "any masking and only guaranteed a red build.  The class that actually "
        "manufactures poly-A false reads is pure_homopolymer, gated absolutely by "
        "--max-pure-homopolymer-kmers (PLAN CAT-17, finding F-014).",
    )
    p.add_argument(
        "--max-pure-homopolymer-kmers",
        type=int,
        default=0,
        metavar="N",
        help="Per-record ceiling on 31-mers that are one base repeated.  This is the "
        "class a 10x poly-A/poly-T tail matches exactly: the deployed panel carried "
        "170 of them across 9 records and they manufactured >1 %% of R2 reads in the "
        "EBV LCL.  A D-list cannot fix this -- it filters host-homologous k-mers, not "
        "self-similarity inside a viral contig (PLAN CAT-17).",
    )
    p.add_argument(
        "--genome-dlist",
        type=Path,
        default=None,
        metavar="GENOME_FA",
        help="Path to a genome-level FASTA (e.g. GRCh38 primary assembly) to use as "
        "a kallisto D-list.  k-mers shared between the D-list and any viral sequence "
        "are masked in the index, preventing reads with host-genomic sequence homology "
        "(including intronic/intergenic regions) from being counted as viral.  "
        "This resolves the cDNA-only artefact described in finding F-005.  "
        "Building with a 3 GB genome D-list requires ~64 GB RAM and ~6 h.",
    )
    p.add_argument(
        "--catalogue",
        type=Path,
        default=None,
        metavar="CATALOGUE_TSV",
        help="Virus catalogue declaring the intended detection panel, columns including "
        "accession/family/species.  Default: the packaged "
        "src/viralscan/data/virus_catalog.tsv.",
    )
    p.add_argument(
        "--index-exclusions",
        type=Path,
        default=None,
        metavar="EXCLUSIONS_TSV",
        help="Reviewed allowlist of deliberate omissions, columns accession/reason/decided_by. "
        "A miss listed here is reported as 'intentional'; every other miss is "
        "'unexplained'.  Default: the packaged src/viralscan/data/index_exclusions.tsv.",
    )
    p.add_argument(
        "--no-dustmask",
        action="store_true",
        help="Skip the NCBI dustmasker pass over the viral records.  Off by default: the "
        "CAT-17 low-complexity gate requires a masked panel, and an unmasked build "
        "manufactures poly-A/poly-T false reads (finding F-014).  Passing this flag will "
        "normally make the gate fail, which is the intended signal.",
    )
    p.add_argument(
        "--dustmask-level",
        type=int,
        default=30,
        metavar="N",
        help="dustmasker -level (default: 30).",
    )
    p.add_argument(
        "--homopolymer-run-length",
        type=int,
        default=20,
        metavar="N",
        help="Mask runs of N or more identical bases to N, after dustmasker. Defaults "
        "to 20. A run of 31 (the k length) is a pure-homopolymer 31-mer, but a 29-nt "
        "run still let host poly-A tails through its 29+2 k-mers (F-021, CAT-42), so "
        "the default leaves at least 11 specific flank bases in any such k-mer. "
        "dustmasker alone does not do this (it masked 0.01 %% of this panel and left "
        "all 85 pure-homopolymer 31-mers in NC_001479.1).  See findings F-014, F-021.",
    )
    p.add_argument(
        "--lowcomplexity-kmer-mask",
        action="store_true",
        help="After the homopolymer mask, N-mask every 31-nt window of an Anelloviridae "
        "record that has <= 2 distinct bases and a base count >= 20, or any base count "
        ">= 28 (PLAN CAT-42: interrupted poly-A / G-C tracts in KP343822.1, KP343842.1). "
        "OFF by default so the run-length-only build stays reproducible. Anelloviridae "
        "only: panel-wide it masks real herpes genes. Writes lowcomplexity_mask.tsv.",
    )
    p.add_argument(
        "--dustmask-windows",
        type=int,
        nargs="+",
        default=[64, 30],
        metavar="N",
        help="dustmasker -window values, run in sequence and merged as a union of "
        "masked positions (default: 64 30).  A single window is not sufficient: the "
        "CAT-17 gate scores 31-mers, so a run dustmasker tolerates at -window 30 can "
        "still yield 31 identical bases.  See finding F-014.",
    )
    p.add_argument(
        "--strict-reconciliation",
        action="store_true",
        help="Fail the build when a catalogued accession is missing from the panel with no "
        "recorded decision.  OFF by default: the catalogue (2,249 accessions) is ahead of "
        "the index by design, so a default-on gate would make the build unusable rather "
        "than honest.  The miss report is written either way, and unexplained misses are "
        "printed prominently regardless.  Turn it on in CI, where a silently-dropped "
        "accession is a sensitivity regression (PLAN CAT-31, finding F-015).",
    )
    return p


def _run_reconciliation(panel_fasta: Path, args: argparse.Namespace) -> None:
    """Reconcile the assembled panel against the catalogue and enforce the verdict.

    Compared on the version-stripped, underscore-normalised base accession, because
    the same genome appears here as NC_007605.1, there as 'NC 007605.1', and in a
    GenBank-only catalogue row with no version at all.
    """
    from viralscan.scripts.build_reference import (
        RECONCILIATION_REPORT_NAME,
        format_reconciliation_summary,
        reconcile_reference_panel,
        reconciliation_failure,
    )

    report = args.out / RECONCILIATION_REPORT_NAME
    try:
        result = reconcile_reference_panel(
            panel_fasta,
            catalogue=args.catalogue,
            exclusions=args.index_exclusions,
            report=report,
        )
    except ValueError as exc:
        sys.exit(f"ERROR: catalogue<->index reconciliation could not run: {exc}")

    print(format_reconciliation_summary(result, report))
    failure = reconciliation_failure(result, strict=args.strict_reconciliation, report=report)
    if failure:
        sys.exit(failure)


def _build_anello_star(viral_fa: Path, out_dir: Path, star: str, threads: int) -> None:
    """Build the STAR anellovirus index anello_star/ (PLAN ANDET-09a).

    Same sequences as the kallisto index (the panel's own masked Anelloviridae
    records), so a disagreement between the two branches is an algorithm
    difference, not a reference difference.
    """
    import hashlib
    import json

    from viralscan.anello_align import genome_generate_cmd, write_star_reference

    print(f"Building STAR anellovirus index → {out_dir}")
    fasta, gtf, n = write_star_reference(viral_fa, _anellovirus_accessions(), out_dir)
    cmd = genome_generate_cmd(star, out_dir, fasta, gtf, threads)
    subprocess.run(cmd, check=True, cwd=out_dir)  # noqa: S603
    if not (out_dir / "SA").is_file():
        sys.exit(f"ERROR: STAR genomeGenerate wrote no SA file in {out_dir}")
    version = subprocess.run(
        [star, "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    manifest = {
        "source_fasta": str(viral_fa),
        "anello_fa_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
        "n_contigs": n,
        "star_version": version,
        "command": cmd,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"  {n} anellovirus contigs indexed (STAR {version})")


def _check_gtf_manifest(repo: Path, manifest: Path) -> None:
    """Fail before fetching when local GTFs differ from the reviewed corpus."""
    subprocess.run(
        [
            sys.executable,
            str(repo / "scripts" / "write_gtf_manifest.py"),
            "--data-dir",
            str(repo / "src" / "viralscan" / "data"),
            "--check",
            str(manifest),
        ],
        check=True,
    )


def main() -> None:
    repo = _find_repo_root()
    sys.path.insert(0, str(repo / "src"))

    from viralscan.scripts.build_reference import validate_reference_records

    args = _build_arg_parser().parse_args()

    if args.anello_star_only:
        viral_fa = args.out / "viral.fa"
        if not viral_fa.is_file():
            sys.exit(f"ERROR: --anello-star-only needs an existing {viral_fa}")
        _build_anello_star(viral_fa, args.out / "anello_star", args.star, args.threads)
        return

    if not args.ncbi_email:
        sys.exit("ERROR: NCBI requires an email. Pass --ncbi-email or set NCBI_EMAIL.")

    _check_gtf_manifest(repo, args.gtf_manifest)

    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    from viralscan.anellovirus import gtf_text_for as _anello_gtf_text  # noqa: E402
    from viralscan.anellovirus import load_accession_table as _load_anello_table  # noqa: E402
    from viralscan.gtf_normalise import normalise_gtf_file  # noqa: E402
    from viralscan.scripts.build_reference import (  # noqa: E402
        fetch_host_cdna,
        host_cdna_as_gtf,
    )
    from viralscan.scripts.ncbi_fetch import (  # noqa: E402
        DEFAULT_CACHE_DIR,
        NCBIFetchError,
        _fetch_one,
    )

    cache_dir: Path = args.cache_dir or DEFAULT_CACHE_DIR

    # ── 1. Download human (host) reference from Ensembl ──────────────────────
    print(f"Step 1/8  Downloading {args.host_species} cDNA + GTF from Ensembl …")
    host_cache = Path.home() / ".cache" / "viralscan" / "ensembl" / args.host_species
    host_fasta_gz, host_gtf_gz = fetch_host_cdna(args.host_species, out / "host", host_cache)
    print(f"  cDNA  : {host_fasta_gz}")
    print(f"  GTF   : {host_gtf_gz}")

    # ── 2. Discover bundled GTFs ──────────────────────────────────────────────
    print("Step 2/8  Scanning bundled viral GTFs …")
    gtf_dir = repo / "src" / "viralscan" / "data"
    gtf_files = sorted(gtf_dir.glob("*.gtf"))
    if not gtf_files:
        sys.exit(f"ERROR: no *.gtf files found in {gtf_dir}")
    print(f"  Found {len(gtf_files)} bundled GTFs in {gtf_dir}")

    # ── 3. Pre-checks ─────────────────────────────────────────────────────────
    print("Step 3/8  Pre-checks …")
    all_seqnames, feature_types = _extract_gtf_seqnames(gtf_files)
    if "exon" not in feature_types:
        sys.exit(
            f"ERROR: bundled GTFs contain no 'exon' features (found: {sorted(feature_types)}). "
            "kb ref extracts cDNA from exon rows — the index would be empty."
        )
    print(f"  'exon' present. All feature types: {sorted(feature_types)}")

    non_accession = sorted(s for s in all_seqnames if not _ACC_RE.match(s))
    if non_accession:
        sys.exit(
            f"ERROR: {len(non_accession)} GTF seqname(s) do not look like NCBI accessions:\n"
            + "\n".join(f"  {s}" for s in non_accession)
        )
    accessions = sorted(all_seqnames)
    print(f"  {len(accessions)} unique NCBI accessions in GTF seqnames")

    # ── 4. Download viral FASTAs from NCBI ───────────────────────────────────
    print(f"Step 4/8  Downloading {len(accessions)} viral FASTAs from NCBI …")
    fasta_paths: list[Path] = []
    errors: list[str] = []
    for i, acc in enumerate(accessions, 1):
        print(f"  [{i:3d}/{len(accessions)}] {acc} … ", end="", flush=True)
        try:
            fasta_path, _ = _fetch_one(acc, cache_dir, args.ncbi_email, args.ncbi_api_key)
            fasta_paths.append(fasta_path)
            print("OK")
        except NCBIFetchError as exc:
            errors.append(f"{acc}: {exc}")
            print(f"FAILED: {exc}")

    if errors:
        sys.exit(
            f"\nERROR: {len(errors)} accession(s) failed to download "
            "(fix and re-run — successful downloads are cached):\n"
            + "\n".join(f"  {e}" for e in errors)
        )

    # ── 4b. Fetch anellovirus panel (clareaulab accessions) ──────────────────
    print("Step 4b/8  Fetching anellovirus panel (clareaulab accessions) …")
    anello_rows = _load_anello_table()
    # clareaulab rows plus any other table row no bundled GTF already covers (the
    # CAT-05 canonical NC_038359.1 is source "viralscan-refseq; CAT-05 canonical").
    anello_accs = sorted(
        row["accession"].strip()
        for row in anello_rows
        if row.get("source", "").strip() == "clareaulab"
        or row["accession"].strip() not in all_seqnames
    )
    print(f"  {len(anello_accs)} anellovirus accessions to fetch")

    anello_fasta_texts: list[str] = []
    anello_gtf_texts: list[str] = []
    anello_errors: list[str] = []
    anello_fasta_by_acc: dict[str, str] = {}

    for i, acc in enumerate(anello_accs, 1):
        if i % 250 == 0:
            print(f"  Anellovirus fetch progress: {i} / {len(anello_accs)}", flush=True)
        try:
            fasta_path, _ = _fetch_one(acc, cache_dir, args.ncbi_email, args.ncbi_api_key)
            text = fasta_path.read_text()
            if text and not text.endswith("\n"):
                text += "\n"
            anello_fasta_texts.append(text)
            anello_fasta_by_acc[acc] = text
        except NCBIFetchError as exc:
            anello_errors.append(f"{acc}: {exc}")

    if anello_errors:
        fail_frac = len(anello_errors) / len(anello_accs) if anello_accs else 0.0
        print(
            f"  WARNING: {len(anello_errors)} / {len(anello_accs)} anellovirus accessions "
            f"failed to download ({fail_frac:.0%}). "
            "Skipping failures — re-run to retry (successful downloads are cached).",
            flush=True,
        )
        if fail_frac > 0.5:
            sys.exit(
                f"ERROR: >50% of anellovirus accessions failed "
                f"({len(anello_errors)}/{len(anello_accs)}). "
                "Check NCBI connectivity and re-run (cached downloads will be reused)."
            )

    # Real CDS structure from the packaged gene catalogue, not one placeholder
    # gene per genome. A whole-genome transcript shares sequence with every other
    # genome in the panel, so reads cross-map in proportion to conservation and
    # the most conserved genome absorbs the panel's entire viral signal.
    anello_gtf_texts.append(_anello_gtf_text(anello_accs, fasta_texts=anello_fasta_by_acc))

    n_anello_genes = anello_gtf_texts[0].count("\texon\t") if anello_gtf_texts else 0
    print(
        f"  {len(anello_fasta_texts)} anellovirus FASTAs fetched, "
        f"{n_anello_genes} exon rows written "
        f"({n_anello_genes / max(1, len(anello_accs)):.1f} genes/genome)",
        flush=True,
    )

    # ── 5. Seqname coverage check ─────────────────────────────────────────────
    print("Step 5/8  Seqname coverage check …")
    fasta_ids = _fasta_seq_ids(fasta_paths)
    missing = set(accessions) - fasta_ids
    if missing:
        sys.exit(
            f"ERROR: {len(missing)} GTF seqname(s) have no matching FASTA record "
            "(kb ref would silently drop them):\n" + "\n".join(f"  {s}" for s in sorted(missing))
        )
    print(f"  OK: all {len(accessions)} GTF seqnames have FASTA records")

    # ── 6. Concatenate: human (gzip) + viral FASTAs → combined.fa ────────────
    print("Step 6/8  Concatenating references …")
    combined_fa = out / "combined.fa"
    combined_gtf = out / "combined.gtf"

    # CAT-17: hardmask every viral record before it reaches the index. Done on a
    # scratch copy so the shared NCBI cache stays pristine, and before the
    # low-complexity gate below so the gate measures what will actually be
    # indexed rather than failing on sequences the builder had not yet fixed.
    # One dustmasker call over the assembled panel, not one per record: dustmasker
    # scores each sequence independently, so the result is identical and 2,343
    # process spawns become 1.
    viral_fa = out / "viral.fa"
    with open(viral_fa, "wb") as fh:
        for fp in fasta_paths:
            data = fp.read_bytes()
            fh.write(data)
            if not data.endswith(b"\n"):
                fh.write(b"\n")
        for text in anello_fasta_texts:
            fh.write(text.encode())

    if not args.no_dustmask:
        print(
            f"  dustmasking {len(fasta_paths) + len(anello_fasta_texts)} viral records "
            f"(level {args.dustmask_level}, windows {'+'.join(map(str, args.dustmask_windows))}) …",
            flush=True,
        )
        _dustmask_fasta_file(viral_fa, args.dustmask_level, tuple(args.dustmask_windows))
        print("  dustmask complete")
        n_masked = _mask_homopolymer_runs(viral_fa, args.homopolymer_run_length)
        print(
            f"  homopolymer runs >= {args.homopolymer_run_length} masked to N: {n_masked:,} base(s)"
        )
        if args.lowcomplexity_kmer_mask:
            lc = _mask_lowcomplexity_kmers(viral_fa, _anellovirus_accessions())
            with open(out / "lowcomplexity_mask.tsv", "w") as fh:
                fh.write("accession\tmasked_bases\n")
                for acc, n in sorted(lc.items()):
                    fh.write(f"{acc}\t{n}\n")
            print(
                f"  low-complexity 31-nt windows masked to N (Anelloviridae): "
                f"{sum(lc.values()):,} base(s) in {len(lc)} record(s)"
            )
    else:
        print("  WARNING: --no-dustmask; the CAT-17 gate will evaluate unmasked sequence")
    print(f"  viral.fa     → {viral_fa}")

    with open(combined_fa, "wb") as fh:
        # Human cDNA first (gzip-encoded from Ensembl)
        with gzip.open(host_fasta_gz, "rb") as gz:
            shutil.copyfileobj(gz, fh)
        # Dustmasked viral records
        with open(viral_fa, "rb") as viral_fh:
            shutil.copyfileobj(viral_fh, fh)
    print(f"  combined.fa  → {combined_fa}")

    # The Ensembl companion GTF (host_gtf_gz) is *chromosomal* (seqnames 1/2/X) and does
    # NOT match the cDNA FASTA headers (ENST…) — handing that pair to kb ref makes it hang
    # forever at "Splitting genome". Generate a cDNA-level host GTF from the FASTA instead.
    host_cdna_gtf = out / "host" / "host_cdna.gtf"
    host_cdna_gtf.parent.mkdir(parents=True, exist_ok=True)
    n_host_tx = host_cdna_as_gtf(host_fasta_gz, host_cdna_gtf)
    print(f"  host cDNA GTF → {host_cdna_gtf} ({n_host_tx} transcripts)")

    with open(combined_gtf, "wb") as fh:
        # Host cDNA-level GTF first (seqname = ENST, matches the cDNA FASTA)
        with open(host_cdna_gtf, "rb") as host_fh:
            shutil.copyfileobj(host_fh, fh)
        # Bundled viral GTFs, REF-13 normalised (gene-scoped transcripts with exons, introns
        # removed for single-protein spliced CDS); the shipped files are not rewritten.
        for gtf in gtf_files:
            fh.write(("\n".join(normalise_gtf_file(gtf, cache_dir)) + "\n").encode())
        # Anellovirus GTF (real NCBI CDS structure where the packaged gene
        # catalogue covers the accession; whole-genome placeholder otherwise)
        for gtf_text in anello_gtf_texts:
            encoded = gtf_text.encode()
            fh.write(encoded)
            if not encoded.endswith(b"\n"):
                fh.write(b"\n")
    print(f"  combined.gtf → {combined_gtf}")

    # ── 7. Catalogue↔index reconciliation (CAT-31) ────────────────────────────
    # The panel is assembled; the catalogue still names 2,249 genomes the project
    # claims to quantify against. A catalogued genome that never reached
    # `viral.fa` is undetectable, and the index content is decided entirely by
    # the bundled GTFs plus the anellovirus fetch — so nothing else in this build
    # would ever notice. Reconcile before kb ref, while failing is still cheap.
    print("Step 7/8  Catalogue<->index reconciliation …")
    _run_reconciliation(viral_fa, args)

    # ── 8. Run kb ref ─────────────────────────────────────────────────────────
    print("Step 8/8  Running kb ref …")
    kb_bin = shutil.which("kb")
    if kb_bin is None:
        sys.exit(
            "ERROR: 'kb' not found on PATH. "
            "Activate the test_viralscan conda env before running this script."
        )

    panel_idx = out / "panel.idx"
    panel_t2g = out / "panel.t2g"
    cdna_fa = out / "cdna.fa"

    cmd = [
        kb_bin,
        "ref",
        "-i",
        str(panel_idx),
        "-g",
        str(panel_t2g),
        "-f1",
        str(cdna_fa),
        "--overwrite",
    ]
    if args.genome_dlist:
        if not args.genome_dlist.exists():
            sys.exit(f"ERROR: --genome-dlist path does not exist: {args.genome_dlist}")
        cmd += ["--d-list", str(args.genome_dlist)]
        print(
            f"  genome D-list: {args.genome_dlist}\n"
            f"  k-mers shared with genome will be masked (F-005 fix)"
        )
    cmd += [str(combined_fa), str(combined_gtf)]
    print(f"  {' '.join(cmd)}")

    # Low-complexity k-mer gate (CAT-17). This builder never called dustmasker, so
    # nothing upstream guarantees the panel is free of homopolymer k-mers — and
    # those k-mers match the poly-A/poly-T tails that dominate 10x R2 reads. A
    # D-list cannot remove them: it filters host-homologous k-mers, not
    # self-similarity within a viral contig. Check before spending ~64 GB and
    # ~8 h on a kallisto build that would bake the artifact into the index.
    viral_fa = out / "viral.fa"
    if viral_fa.exists():
        try:
            validate_reference_records(
                viral_fa,
                max_low_complexity_fraction=args.max_low_complexity_fraction,
                max_pure_homopolymer_kmers=args.max_pure_homopolymer_kmers,
            )
        except ValueError as exc:
            sys.exit(
                f"ERROR: viral panel failed the low-complexity k-mer gate: {exc}\n"
                f"  Mask it first: dustmasker -window 64 -level 30 (and -window 30), "
                f"merged, masked to N. Override with --max-pure-homopolymer-kmers / "
                f"--max-low-complexity-fraction if this is deliberate."
            )
        print(
            f"  low-complexity k-mer gate passed (pure-homopolymer <= "
            f"{args.max_pure_homopolymer_kmers}, fraction <= {args.max_low_complexity_fraction})"
        )
    else:
        print(f"  WARNING: {viral_fa} not found; low-complexity k-mer gate skipped")

    subprocess.run(cmd, check=True)  # noqa: S603

    # ── Verify output ─────────────────────────────────────────────────────────
    for path in [panel_idx, panel_t2g, cdna_fa]:
        if not path.exists() or path.stat().st_size == 0:
            sys.exit(f"ERROR: expected output missing or empty: {path}")

    # PROG-13: HSV-1 / VZV repeat copies are 31-mer-identical, so merge them to one gene in
    # the t2g (the index is untouched). Idempotent; src/ is on sys.path since step 1.
    from viralscan.repeat_merges import apply_repeat_merges

    print(f"  repeat-copy gene merges (PROG-13): {apply_repeat_merges(panel_t2g)} t2g row(s)")

    t2g_lines = sum(1 for _ in panel_t2g.open())
    human_lines = sum(1 for ln in panel_t2g.open() if ln.startswith("ENST"))
    viral_lines = t2g_lines - human_lines
    print("\nIndex build complete:")
    print(f"  panel.idx : {panel_idx}  ({panel_idx.stat().st_size // (1024 * 1024)} MB)")
    print(
        f"  panel.t2g : {panel_t2g}  ({t2g_lines:,} total  |  {human_lines:,} human ENST*  |  {viral_lines:,} viral)"
    )
    print(f"  cdna.fa   : {cdna_fa}")

    epstein = sum(1 for ln in panel_t2g.open() if "EPSTEIN" in ln)
    print(f"  Spot check: {epstein} EPSTEIN (EBV) entries in panel.t2g", end="")
    if epstein == 0:
        print("  ← WARNING: expected >0; check EBV GTF/FASTA inclusion")
    else:
        print("  ✓")

    anello_t2g = sum(1 for ln in panel_t2g.open() if "_gene" in ln and not ln.startswith("ENST"))
    print(f"  Spot check: {anello_t2g} anellovirus-style (_geneN) entries in panel.t2g", end="")
    if anello_fasta_texts and anello_t2g == 0:
        print("  ← WARNING: expected >0; check Step 4b anellovirus fetch and GTF append")
    else:
        print("  ✓")

    _build_anello_star(out / "viral.fa", out / "anello_star", args.star, args.threads)

    print(
        "\nNext step: run format probe on one sample to confirm kb count -x BULK output layout.\n"
        "The output h5ad/mtx will contain both ENST* (host) and viral gene_ids —\n"
        "bulk_viral_summarize.py emits a viral summary TSV; the h5ad is the full\n"
        "host+viral matrix for coexpression analysis.\n"
        "See the plan notes in scripts/bulk_viral_summarize.py."
    )


if __name__ == "__main__":
    main()
