"""Reference-strategy benchmark contracts and validators.

This module keeps the 12-row ViralScan/STARsolo benchmark auditable without
requiring the heavy alignments to run inside unit tests.  CLI wrappers in
``scripts/`` use these helpers to freeze reference provenance, generate command
manifests, and validate final benchmark tables.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from viralscan.evidence import cb_umi_geometry

PANEL_ID = "serratus_plus_expanded_anellovirus"
# STAR executable for the benchmark harness. Defaults to ``STAR`` on ``PATH``;
# override with the ``VIRALSCAN_STAR_BIN`` environment variable to pin a specific
# build. (Previously hardcoded to an author-specific absolute path, which broke
# the benchmark on any other machine.)
STAR_BIN = os.environ.get("VIRALSCAN_STAR_BIN", "STAR")

DATASETS: tuple[dict[str, str], ...] = (
    {
        "dataset": "hhv6b",
        "sample_id": "hhv6_carT_ref",
        "srr": "SRR20710641",
        "technology": "10xv3",
        "target_virus": "HHV-6B",
        "target_regex": r"(?i)(hhv-?6b|human[_ -]herpesvirus[_ -]?6b|herpesvirus[_ -]?6b|hum[_ -]herp6b)",
        "off_target_regex": r"(?i)(hhv-?6a|human[_ -]herpesvirus[_ -]?6(?!b)|herpesvirus[_ -]?6(?!b)|hum[_ -]herp6(?!b))",
    },
    {
        "dataset": "ebv",
        "sample_id": "lcl_5lines",
        "srr": "SRR12682296",
        "technology": "10xv2",
        "target_virus": "EBV",
        "target_regex": r"(?i)(epstein|ebv|hhv-?4|epstein[_ -]hhv4)",
        "off_target_regex": r"(?i)(hhv-?[1235678]|herpesvirus[_ -]?[1235678])",
    },
    {
        "dataset": "hsv1",
        "sample_id": "hsv1_fibroblast",
        "srr": "SRR8315713",
        "technology": "DROPSEQ",
        "target_virus": "HSV-1",
        # HSV-1 == Human herpesvirus 1. STARsolo's Serratus gene_ids use the `HHV1gp…`
        # form (and the NC_001806 accession), which the hsv-1/herpesvirus-1 aliases do not
        # match — unlike EBV (`hhv-?4`) and HHV-6B (`hhv-?6b`), whose regexes already include
        # the HHV-N form. Without `hhv-?1` here, STARsolo HSV-1 is silently counted as 0.
        # `(?![0-9])` prevents `hhv1` from matching a hypothetical HHV-1x.
        "target_regex": r"(?i)(hsv-?1|hhv-?1(?![0-9])|human[_ -]herpesvirus[_ -]?1|herpesvirus[_ -]?1|hum[_ -]herp1|nc_?001806)",
        "off_target_regex": r"(?i)(hsv-?2|hhv-?2(?![0-9])|human[_ -]herpesvirus[_ -]?2|herpesvirus[_ -]?2|hum[_ -]herp2)",
    },
)

METHOD_ROWS: tuple[dict[str, str], ...] = (
    {"method": "starsolo", "reference_strategy": "combined"},
    {"method": "starsolo", "reference_strategy": "two_step"},
    {"method": "viralscan", "reference_strategy": "combined"},
    {"method": "viralscan", "reference_strategy": "two_step"},
)

REQUIRED_RESULT_COLUMNS: tuple[str, ...] = (
    "dataset",
    "srr",
    "technology",
    "method",
    "reference_strategy",
    "target_virus",
    "reference_hash_id",
    "command_id",
    "slurm_job_id",
    "count_layer",
    "barcode_universe",
    "denominator_fixed_barcodes",
    "denominator_method_called_cells",
    "denominator_shared_anchor_cells",
    "target_positive_cells_fixed",
    "target_umi_fixed",
    "related_off_target_positive_cells_fixed",
    "related_off_target_umi_fixed",
    "combined_vs_two_step_delta_target_umi",
    "status",
    "failure_reason",
)

REQUIRED_AUDIT_TARGETS: tuple[str, ...] = ("HHV-6B", "EBV", "HSV-1")

MANIFEST_PATH_FIELDS: frozenset[str] = frozenset(
    {
        "anellovirus_accession_table",
        "fastq_root",
        "genes_gtf",
        "genome_dir",
        "genome_fasta",
        "genome_gtf",
        "gtf",
        "human_cdna",
        "human_transcriptome_fasta",
        "kallisto_index",
        "t2g",
    }
)

REQUIRED_PROVENANCE_FIELDS: tuple[str, ...] = (
    "build_commands",
    "provenance",
    "viral_panel.anellovirus_expected_count",
    "viral_panel.anellovirus_fetched_count",
    "viral_panel.anellovirus_missing_count",
)

REQUIRED_BUILD_COMMANDS: tuple[str, ...] = (
    "starsolo_human_only",
    "starsolo_all_virus",
    "starsolo_combined",
    "kallisto_human_only",
    "kallisto_all_virus",
    "kallisto_combined",
)

STAR_INDEX_REQUIRED_FILES: tuple[str, ...] = (
    "Genome",
    "SA",
    "SAindex",
    "genomeParameters.txt",
)

_PLACEHOLDER_BUILD_COMMANDS: frozenset[str] = frozenset(
    {"todo", "tbd", "template", "placeholder", "set_me", "fill_me"}
)

FORBIDDEN_REFERENCE_PATTERNS: tuple[str, ...] = (
    "fasta_split",
    "split_gtf",
    "GRCh38_EBV",
    "starsolo_p22_6",
    "transcriptomev3_noHuman",
    "Epstein_Barr_virus_NC_007605",
    "Human_herpesvirus_1_NC_001806",
    "Human_herpesvirus_6_NC_001664",
)


class BenchmarkContractError(ValueError):
    """Raised when a manifest, command manifest, or result table violates the contract."""


@dataclass(frozen=True)
class PathAudit:
    key: str
    path: str
    exists: bool
    size_bytes: int | str
    mtime: str
    sha256: str
    feature_count: int | str
    expected_virus_presence: str


@dataclass(frozen=True)
class FastqAudit:
    srr: str
    mate: str
    path: str
    exists: bool
    size_bytes: int | str
    sha256: str
    records: int | str
    valid_fastq: bool
    source_url: str
    failure_reason: str


@dataclass
class ParsedMetrics:
    counts_by_barcode: dict[str, tuple[float, float]]
    anchor_barcodes: set[str]
    status: str
    failure_reason: str = ""


def _nested_get(mapping: dict[str, Any], dotted_key: str) -> Any:
    value: Any = mapping
    for part in dotted_key.split("."):
        if not isinstance(value, dict) or part not in value:
            return None
        value = value[part]
    return value


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open_fastq_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def _fastq_stats(path: Path) -> tuple[int, str]:
    records = 0
    try:
        with _open_fastq_text(path) as handle:
            while True:
                header = handle.readline()
                if not header:
                    return records, ""
                seq = handle.readline()
                plus = handle.readline()
                qual = handle.readline()
                if not (seq and plus and qual):
                    return records, "FASTQ ends mid-record"
                if not header.startswith("@"):
                    return records, "FASTQ header does not start with @"
                if not plus.startswith("+"):
                    return records, "FASTQ plus line does not start with +"
                if len(seq.rstrip("\n\r")) != len(qual.rstrip("\n\r")):
                    return records, "quality string length is not equal to sequence length"
                records += 1
    except OSError as exc:
        return records, f"gzip or file integrity error: {exc}"


def _fastq_source_url(manifest: dict[str, Any], dataset: dict[str, str], mate: str) -> str:
    fastqs = manifest.get("fastqs", {})
    srr = dataset["srr"]
    if isinstance(fastqs, dict):
        value = fastqs.get(srr) or fastqs.get(dataset["dataset"]) or {}
        if isinstance(value, dict):
            source_url = value.get("source_url")
            mate_value = value.get(mate) or value.get(f"{srr}_{mate}.fastq.gz") or {}
            if isinstance(mate_value, dict):
                return str(mate_value.get("source_url") or source_url or "")
            if isinstance(source_url, str):
                return source_url
    return ""


def fastq_paths_for_dataset(manifest: dict[str, Any], dataset: dict[str, str]) -> tuple[Path, Path]:
    fastqs = manifest.get("fastqs", {})
    srr = dataset["srr"]
    if isinstance(fastqs, dict):
        value = fastqs.get(srr) or fastqs.get(dataset["dataset"]) or {}
        if isinstance(value, dict):
            r1 = value.get("R1") or value.get("r1") or value.get("mate1")
            r2 = value.get("R2") or value.get("r2") or value.get("mate2")
            if isinstance(r1, dict):
                r1 = r1.get("path")
            if isinstance(r2, dict):
                r2 = r2.get("path")
            if r1 and r2:
                return Path(str(r1)), Path(str(r2))
    fastq_root = Path(
        manifest.get("fastq_root", "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data")
    )
    sample_fastq_dir = fastq_root / dataset["sample_id"] / srr
    return sample_fastq_dir / f"{srr}_1.fastq.gz", sample_fastq_dir / f"{srr}_2.fastq.gz"


def audit_fastqs(
    manifest: dict[str, Any],
    *,
    base_dir: Path | None = None,
    compute_sha256: bool = True,
) -> list[FastqAudit]:
    root = Path(base_dir or Path.cwd())
    expected: list[tuple[dict[str, str], str, Path]] = []
    for dataset in DATASETS:
        for mate, raw_path in zip(("R1", "R2"), fastq_paths_for_dataset(manifest, dataset)):
            expected.append((dataset, mate, _normalize_manifest_path(raw_path, base_dir=root)))
    fail_fast_missing = any(not path.exists() for _dataset, _mate, path in expected)
    rows: list[FastqAudit] = []
    pair_counts: dict[str, list[int | str]] = {}
    for dataset, mate, path in expected:
        exists = path.exists()
        records: int | str = ""
        reason = ""
        if not exists:
            reason = "missing FASTQ"
        elif fail_fast_missing:
            reason = ""
        elif path.is_file():
            records, reason = _fastq_stats(path)
            pair_counts.setdefault(dataset["srr"], []).append(records)
        rows.append(
            FastqAudit(
                srr=dataset["srr"],
                mate=mate,
                path=str(path),
                exists=exists,
                size_bytes=path.stat().st_size if exists else "",
                sha256=(
                    sha256_file(path)
                    if compute_sha256
                    and exists
                    and path.is_file()
                    and not reason
                    and not fail_fast_missing
                    else ""
                ),
                records=records,
                valid_fastq=exists and not reason,
                source_url=_fastq_source_url(manifest, dataset, mate),
                failure_reason=reason,
            )
        )
    for srr, counts in pair_counts.items():
        if len(counts) == 2 and counts[0] != counts[1]:
            for idx, row in enumerate(rows):
                if row.srr == srr:
                    reason = row.failure_reason or "paired FASTQ record counts differ"
                    rows[idx] = FastqAudit(
                        **{**row.__dict__, "valid_fastq": False, "failure_reason": reason}
                    )
    return rows


def write_fastq_audit(path: Path, rows: list[FastqAudit]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [field for field in FastqAudit.__dataclass_fields__]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row.__dict__)


def validate_fastq_audit_rows(
    rows: list[FastqAudit],
    *,
    require_hashes: bool = True,
    require_source_urls: bool = True,
) -> None:
    missing = [row.path for row in rows if not row.exists]
    if missing:
        raise BenchmarkContractError(f"FASTQ files are missing: {missing}")
    invalid = [row for row in rows if not row.valid_fastq]
    if invalid:
        reasons = {f"{row.srr}:{row.mate}": row.failure_reason for row in invalid}
        raise BenchmarkContractError(f"invalid FASTQ records: {reasons}")
    missing_hashes = [row.path for row in rows if not row.sha256]
    if require_hashes and missing_hashes:
        raise BenchmarkContractError(f"missing FASTQ SHA256 values: {missing_hashes}")
    if require_source_urls:
        missing_sources = [f"{row.srr}:{row.mate}" for row in rows if not row.source_url]
        if missing_sources:
            raise BenchmarkContractError(f"FASTQ source URLs are missing: {missing_sources}")


def iter_manifest_paths(manifest: dict[str, Any]) -> Iterable[tuple[str, Path]]:
    def walk(prefix: str, value: Any) -> Iterable[tuple[str, Path]]:
        if isinstance(value, dict):
            for key, nested in value.items():
                yield from walk(f"{prefix}.{key}" if prefix else key, nested)
        elif isinstance(value, list):
            for idx, nested in enumerate(value):
                yield from walk(f"{prefix}[{idx}]", nested)
        elif isinstance(value, str) and prefix.rsplit(".", 1)[-1] in MANIFEST_PATH_FIELDS:
            yield prefix, Path(value)

    yield from walk("", manifest)


def load_manifest(path: Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        manifest = json.load(handle)
    if not isinstance(manifest, dict):
        raise BenchmarkContractError(f"{path} must contain a JSON object")
    return manifest


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _read_text_sample(path: Path, max_bytes: int = 2_000_000) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="ignore") as handle:
        return handle.read(max_bytes)


def _read_text_edges(path: Path, max_bytes: int = 2_000_000) -> str:
    if path.suffix == ".gz":
        return _read_text_sample(path, max_bytes=max_bytes)
    size = path.stat().st_size
    with path.open("rb") as handle:
        head = handle.read(max_bytes)
        if size > max_bytes:
            handle.seek(max(size - max_bytes, 0))
            tail = handle.read(max_bytes)
        else:
            tail = b""
    return (head + b"\n" + tail).decode("utf-8", errors="ignore")


def _feature_count(
    path: Path, max_lines: int = 2_000_000, max_bytes: int = 100_000_000
) -> int | str:
    if not path.exists() or not path.is_file():
        return ""
    if path.suffix not in {".gtf", ".gff", ".gff3", ".tsv", ".txt"}:
        return ""
    if path.stat().st_size > max_bytes:
        return f"large_file>{max_bytes}"
    count = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="ignore") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                count += 1
                if count >= max_lines:
                    return f">={max_lines}"
    return count


def _presence(path: Path) -> str:
    if not path.exists() or not path.is_file():
        return "missing"
    sample = _read_text_edges(path)
    present = []
    for dataset in DATASETS:
        if re.search(dataset["target_regex"], sample):
            present.append(dataset["target_virus"])
    return ",".join(present) if present else "not_detected_in_sample"


def _should_scan_presence(key: str, path: Path) -> bool:
    if not path.is_file():
        return False
    if path.suffix not in {".fa", ".fasta", ".gtf", ".gff", ".gff3", ".tsv", ".txt"}:
        return False
    return any(token in key for token in ("all_virus", "combined", "viral_panel"))


def _missing_star_index_files(row: PathAudit) -> list[str]:
    if not row.key.endswith("genome_dir") or not Path(row.path).is_dir():
        return []
    genome_dir = Path(row.path)
    return [name for name in STAR_INDEX_REQUIRED_FILES if not (genome_dir / name).exists()]


def _normalize_manifest_path(value: Path, *, base_dir: Path) -> Path:
    return value if value.is_absolute() else (base_dir / value)


def _is_placeholder_value(value: Any) -> bool:
    if not isinstance(value, str):
        return True
    lowered = value.strip().lower()
    return not lowered or lowered in _PLACEHOLDER_BUILD_COMMANDS


def _validate_provenance(manifest: dict[str, Any]) -> None:
    build_commands = manifest.get("build_commands")
    if not isinstance(build_commands, dict) or not build_commands:
        raise BenchmarkContractError("build_commands must be a non-empty object")

    missing = [name for name in REQUIRED_BUILD_COMMANDS if name not in build_commands]
    if missing:
        raise BenchmarkContractError(f"build_commands missing required step names: {missing}")

    for name in REQUIRED_BUILD_COMMANDS:
        command = build_commands.get(name)
        if not isinstance(command, str) or _is_placeholder_value(command):
            raise BenchmarkContractError(f"build_commands[{name!r}] is missing or a placeholder")

    expected = _nested_get(manifest, "viral_panel.anellovirus_expected_count")
    fetched = _nested_get(manifest, "viral_panel.anellovirus_fetched_count")
    missing = _nested_get(manifest, "viral_panel.anellovirus_missing_count")
    for field, value in (
        ("viral_panel.anellovirus_expected_count", expected),
        ("viral_panel.anellovirus_fetched_count", fetched),
        ("viral_panel.anellovirus_missing_count", missing),
    ):
        if not isinstance(value, int):
            raise BenchmarkContractError(f"{field} must be an integer")
        if value < 0:
            raise BenchmarkContractError(f"{field} must be >= 0")

    if fetched + missing > expected:
        raise BenchmarkContractError(
            "viral_panel anellovirus counts are inconsistent "
            f"(expected={expected}, fetched={fetched}, missing={missing})"
        )


def audit_manifest(
    manifest: dict[str, Any],
    *,
    compute_sha256: bool = True,
    base_dir: Path | None = None,
) -> list[PathAudit]:
    root = Path(base_dir or Path.cwd())
    rows: list[PathAudit] = []
    for key, path in iter_manifest_paths(manifest):
        path = _normalize_manifest_path(path, base_dir=root)
        exists = path.exists()
        stat = path.stat() if exists else None
        rows.append(
            PathAudit(
                key=key,
                path=str(path),
                exists=exists,
                size_bytes=stat.st_size if stat else "",
                mtime=str(int(stat.st_mtime)) if stat else "",
                sha256=sha256_file(path) if exists and path.is_file() and compute_sha256 else "",
                feature_count=_feature_count(path),
                expected_virus_presence=_presence(path)
                if exists and _should_scan_presence(key, path)
                else "",
            )
        )
    return rows


def write_reference_audit(path: Path, rows: list[PathAudit]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [field for field in PathAudit.__dataclass_fields__]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row.__dict__)


def validate_manifest(
    manifest: dict[str, Any],
    *,
    expected_panel: str = PANEL_ID,
    fail_on_single_virus: bool = False,
    require_provenance: bool = False,
) -> None:
    panel = manifest.get("panel") or manifest.get("viral_panel", {}).get("id")
    if panel != expected_panel:
        raise BenchmarkContractError(f"expected panel {expected_panel!r}, found {panel!r}")

    human_source = manifest.get("human", {}).get("source_release")
    star_human = manifest.get("references", {}).get("starsolo", {}).get("human_source_release")
    vs_human = manifest.get("references", {}).get("viralscan", {}).get("human_source_release")
    releases = {x for x in (human_source, star_human, vs_human) if x}
    if len(releases) > 1:
        raise BenchmarkContractError(f"human source releases do not match: {sorted(releases)}")

    all_text = json.dumps(manifest, sort_keys=True)
    if fail_on_single_virus:
        hits = [p for p in FORBIDDEN_REFERENCE_PATTERNS if p in all_text]
        if hits:
            raise BenchmarkContractError(f"single-virus or historical references found: {hits}")

    star_refs = manifest.get("references", {}).get("starsolo", {})
    all_virus = star_refs.get("all_virus", {})
    if not (all_virus.get("genome_fasta") and all_virus.get("genome_gtf")):
        raise BenchmarkContractError("STAR-compatible all-virus genome FASTA/GTF is missing")

    if require_provenance:
        missing = [
            field
            for field in REQUIRED_PROVENANCE_FIELDS
            if _nested_get(manifest, field) in (None, "", {})
        ]
        if missing:
            raise BenchmarkContractError(f"missing required provenance fields: {missing}")
        _validate_provenance(manifest)


def validate_audit_rows(
    rows: list[PathAudit],
    *,
    require_hashes: bool = True,
    require_target_presence: bool = True,
) -> None:
    missing = [row for row in rows if not row.exists]
    if missing:
        raise BenchmarkContractError(f"{len(missing)} manifest paths are missing")

    if require_hashes:
        missing_hashes = [row.key for row in rows if Path(row.path).is_file() and not row.sha256]
        if missing_hashes:
            raise BenchmarkContractError(f"missing SHA256 for file artifacts: {missing_hashes}")

    incomplete_star_dirs = {
        row.key: _missing_star_index_files(row)
        for row in rows
        if row.key.endswith("genome_dir") and row.exists
    }
    incomplete_star_dirs = {key: files for key, files in incomplete_star_dirs.items() if files}
    if incomplete_star_dirs:
        raise BenchmarkContractError(
            f"incomplete STAR genomeGenerate directories: {incomplete_star_dirs}"
        )

    if require_target_presence:
        viral_rows = [
            row
            for row in rows
            if any(token in row.key for token in ("all_virus", "combined"))
            and Path(row.path).is_file()
            and Path(row.path).suffix in {".fa", ".fasta", ".gtf", ".gff", ".gff3", ".tsv", ".txt"}
        ]
        seen: set[str] = set()
        for row in viral_rows:
            seen.update(filter(None, row.expected_virus_presence.split(",")))
        missing_targets = sorted(set(REQUIRED_AUDIT_TARGETS) - seen)
        if missing_targets:
            raise BenchmarkContractError(
                f"expected viruses absent from audited viral files: {missing_targets}"
            )


def benchmark_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset in DATASETS:
        for method in METHOD_ROWS:
            row = {**dataset, **method}
            row["row_id"] = (
                f"{dataset['dataset']}__{method['method']}__{method['reference_strategy']}"
            )
            rows.append(row)
    return rows


def _thread_arg() -> str:
    return "__VS_THREADS__"


def starsolo_geometry_args(technology: str) -> list[str]:
    cb_len, umi_len = cb_umi_geometry(technology)
    return [
        "--soloCBstart",
        "1",
        "--soloCBlen",
        str(cb_len),
        "--soloUMIstart",
        str(cb_len + 1),
        "--soloUMIlen",
        str(umi_len),
    ]


def _starsolo_base_command(
    *,
    genome_dir: str,
    r1: str,
    r2: str,
    technology: str,
    out_prefix: str,
    out_reads_unmapped: str = "None",
    read_files_command: str | None = "zcat",
) -> list[str]:
    cmd = [
        STAR_BIN,
        "--genomeDir",
        genome_dir,
        "--readFilesIn",
        r2,
        r1,
    ]
    if read_files_command:
        cmd += ["--readFilesCommand", read_files_command]
    cmd += [
        "--soloType",
        "CB_UMI_Simple",
        *starsolo_geometry_args(technology),
        "--soloCBwhitelist",
        "None",
        "--soloBarcodeReadLength",
        "0",
        "--soloFeatures",
        "GeneFull",
        "--soloCellFilter",
        "CellRanger2.2",
        "--outReadsUnmapped",
        out_reads_unmapped,
        "--runThreadN",
        _thread_arg(),
        "--outFileNamePrefix",
        out_prefix,
    ]
    return cmd


def commands_for_row(
    row: dict[str, str], manifest: dict[str, Any], run_dir: Path
) -> list[list[str]]:
    refs = manifest.get("references", {})
    r1_path, r2_path = fastq_paths_for_dataset(manifest, row)
    r1 = str(r1_path)
    r2 = str(r2_path)
    out = str(run_dir / "runs" / row["row_id"])
    if row["method"] == "viralscan":
        vs = refs.get("viralscan", {})
        ref = vs["combined"] if row["reference_strategy"] == "combined" else vs["all_virus"]
        cmd: list[str] = [
            "python",
            "src/viralscan/menu.py",
            "-i",
            ref["kallisto_index"],
            "-t",
            ref["t2g"],
            "-gtf",
            ref["gtf"],
            "-o",
            out,
            "-s1",
            r1,
            "-s2",
            r2,
            "-x",
            row["technology"],
            "-c",
            _thread_arg(),
            "--yes",
        ]
        if row["reference_strategy"] == "two_step":
            cmd += ["--host-filter", "kallisto", "--host-index", vs["human_only"]["kallisto_index"]]
        return [cmd]

    star = refs.get("starsolo", {})
    if row["reference_strategy"] == "combined":
        return [
            _starsolo_base_command(
                genome_dir=star["combined"]["genome_dir"],
                r1=r1,
                r2=r2,
                technology=row["technology"],
                out_prefix=out + "/starsolo/",
            )
        ]

    host_prefix = out + "/starsolo_host/"
    virus_prefix = out + "/starsolo_virus/"
    host_cmd = _starsolo_base_command(
        genome_dir=star["human_only"]["genome_dir"],
        r1=r1,
        r2=r2,
        technology=row["technology"],
        out_prefix=host_prefix,
        out_reads_unmapped="Fastx",
    )
    # STAR writes unmapped mate1 for first readFilesIn (R2/cDNA) and mate2 for
    # second readFilesIn (R1/barcode+UMI).  Feed those to the viral pass in the
    # same STARsolo order: cDNA first, barcode read second.
    viral_cmd = _starsolo_base_command(
        genome_dir=star["all_virus"]["genome_dir"],
        r1=host_prefix + "Unmapped.out.mate2",
        r2=host_prefix + "Unmapped.out.mate1",
        technology=row["technology"],
        out_prefix=virus_prefix,
        read_files_command=None,
    )
    return [host_cmd, viral_cmd]


def command_for_row(row: dict[str, str], manifest: dict[str, Any], run_dir: Path) -> list[str]:
    return commands_for_row(row, manifest, run_dir)[0]


def write_commands(run_dir: Path, manifest: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    run_dir.mkdir(parents=True, exist_ok=True)
    commands_path = run_dir / "commands.jsonl"
    with open(commands_path, "w", encoding="utf-8") as handle:
        for row in benchmark_rows():
            commands = commands_for_row(row, manifest, run_dir)
            record = {**row, "commands": commands, "command": commands[0]}
            rows.append(record)
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    return rows


def write_slurm_array(run_dir: Path) -> Path:
    path = run_dir / "run_reference_strategy_array.sh"
    resolved_run_dir = run_dir.resolve()
    log_dir = run_dir.resolve() / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""#!/usr/bin/env bash
#SBATCH --job-name=vs_ref_strategy
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output={log_dir}/%x_%A_%a.out
#SBATCH --error={log_dir}/%x_%A_%a.err

set -euo pipefail

source /share/software/tools/miniconda/3.10/23.3.1/etc/profile.d/conda.sh
conda activate /exports/archive/hg-funcgenom-research/evonk/conda/envs/test_viralscan
KB_PYTHON_BIN_DIR=$(python - <<'PY'
from pathlib import Path
import kb_python

path = Path(kb_python.__file__).parent / "bins" / "linux" / "kallisto"
print(path if path.exists() else "")
PY
)
if [[ -n "$KB_PYTHON_BIN_DIR" ]]; then
    export PATH="$KB_PYTHON_BIN_DIR:$PATH"
fi
export PYTHONPATH=/exports/para-lipg-hpc/mdmanurung/ViralScan/src

missing_tools=()
for tool in python kb snakemake kallisto bustools; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        missing_tools+=("$tool")
    fi
done
if (( ${{#missing_tools[@]}} )); then
    printf 'ERROR: missing required benchmark runtime tools after conda activation: %s\\n' "${{missing_tools[*]}}" >&2
    printf 'PATH=%s\\n' "$PATH" >&2
    exit 127
fi

RUN_DIR={resolved_run_dir}
mkdir -p "$RUN_DIR/logs"
ROW=${{SLURM_ARRAY_TASK_ID:?submit with --array=0-11}}
export VS_THREADS=${{SLURM_CPUS_PER_TASK:-8}}
python - "$RUN_DIR/commands.jsonl" "$ROW" <<'PY'
import json
import os
import subprocess
import sys

path, idx = sys.argv[1], int(sys.argv[2])
os.chdir("/exports/para-lipg-hpc/mdmanurung/ViralScan")
with open(path) as handle:
    for i, line in enumerate(handle):
        if i == idx:
            record = json.loads(line)
            break
    else:
        raise SystemExit(f"row {{idx}} not found")

threads = os.environ.get("VS_THREADS", "8")
for step, cmd in enumerate(record["commands"], 1):
    materialized = [threads if arg == "__VS_THREADS__" else arg for arg in cmd]
    print(f"row={{idx}} step={{step}} command={{materialized}}", flush=True)
    subprocess.run(materialized, check=True)
PY
""",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def validate_commands(commands_path: Path, *, fail_on_single_virus: bool = True) -> None:
    text = commands_path.read_text(encoding="utf-8")
    if fail_on_single_virus:
        hits = [p for p in FORBIDDEN_REFERENCE_PATTERNS if p in text]
        if hits:
            raise BenchmarkContractError(f"forbidden references in command manifest: {hits}")
    rows = [json.loads(line) for line in text.splitlines() if line.strip()]
    expected = {(r["dataset"], r["method"], r["reference_strategy"]) for r in benchmark_rows()}
    found = {(r["dataset"], r["method"], r["reference_strategy"]) for r in rows}
    if found != expected:
        raise BenchmarkContractError(
            f"command rows mismatch; missing={sorted(expected - found)} extra={sorted(found - expected)}"
        )
    for row in rows:
        commands = row.get("commands") or [row.get("command")]
        if row["method"] == "starsolo" and row["reference_strategy"] == "two_step":
            if len(commands) != 2:
                raise BenchmarkContractError(
                    f"{row['row_id']} must have host and viral STAR commands"
                )
            if "Fastx" not in commands[0] or "Unmapped.out.mate" not in " ".join(commands[1]):
                raise BenchmarkContractError(
                    f"{row['row_id']} missing STAR host-unmapped or viral second pass"
                )
        for cmd in commands:
            if not isinstance(cmd, list) or not all(isinstance(arg, str) for arg in cmd):
                raise BenchmarkContractError(f"{row['row_id']} command must be an argv string list")
            if "${SLURM_CPUS_PER_TASK:-8}" in cmd:
                raise BenchmarkContractError(
                    f"{row['row_id']} contains unevaluated shell thread expression"
                )


def _read_lines(path: Path) -> list[str]:
    return [line.rstrip("\n") for line in path.read_text(encoding="utf-8").splitlines()]


def _match_feature(name: str, pattern: str) -> bool:
    return re.search(pattern, name or "") is not None


def parse_starsolo_metrics(
    raw_dir: Path,
    filtered_barcodes: Path,
    *,
    target_regex: str,
    off_target_regex: str,
) -> ParsedMetrics:
    features_path = raw_dir / "features.tsv"
    barcodes_path = raw_dir / "barcodes.tsv"
    matrix_path = raw_dir / "matrix.mtx"
    required = [features_path, barcodes_path, matrix_path, filtered_barcodes]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        return ParsedMetrics({}, set(), "incomplete", f"missing STARsolo outputs: {missing}")

    features = _read_lines(features_path)
    barcodes = _read_lines(barcodes_path)
    target_idx: set[int] = set()
    off_target_idx: set[int] = set()
    for idx, line in enumerate(features, 1):
        parts = line.split("\t")
        searchable = " ".join(parts)
        if _match_feature(searchable, target_regex):
            target_idx.add(idx)
        if _match_feature(searchable, off_target_regex):
            off_target_idx.add(idx)

    counts: dict[str, list[float]] = {barcode: [0.0, 0.0] for barcode in barcodes}
    with matrix_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            dims = line.split()
            if len(dims) == 3:
                break
        for line in handle:
            parts = line.split()
            if len(parts) != 3:
                continue
            feature_i = int(parts[0])
            barcode_i = int(parts[1])
            value = float(parts[2])
            if barcode_i < 1 or barcode_i > len(barcodes):
                continue
            barcode = barcodes[barcode_i - 1]
            if feature_i in target_idx:
                counts[barcode][0] += value
            if feature_i in off_target_idx:
                counts[barcode][1] += value

    anchor = set(_read_lines(filtered_barcodes))
    return ParsedMetrics(
        {barcode: (values[0], values[1]) for barcode, values in counts.items()},
        anchor,
        "complete",
    )


def _find_viralscan_results(row_dir: Path, filename: str) -> Path | None:
    direct = row_dir / "results" / filename
    if direct.exists():
        return direct
    matches = sorted(row_dir.glob(f"*/results/{filename}"))
    return matches[0] if matches else None


def parse_viralscan_metrics(
    row_dir: Path,
    *,
    target_regex: str,
    off_target_regex: str,
) -> ParsedMetrics:
    per_cell = _find_viralscan_results(row_dir, "per_cell_viral.tsv")
    summary = _find_viralscan_results(row_dir, "viral_summary.tsv")
    missing = [
        name
        for name, path in (("per_cell_viral.tsv", per_cell), ("viral_summary.tsv", summary))
        if path is None
    ]
    if missing:
        return ParsedMetrics({}, set(), "incomplete", f"missing ViralScan outputs: {missing}")
    assert per_cell is not None and summary is not None  # guaranteed by the `missing` check above

    counts: dict[str, list[float]] = {}
    with per_cell.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"barcode", "virus_name", "viral_umi"}
        if not required.issubset(reader.fieldnames or []):
            return ParsedMetrics({}, set(), "failed", "per_cell_viral.tsv missing required columns")
        for row in reader:
            barcode = row["barcode"]
            values = counts.setdefault(barcode, [0.0, 0.0])
            virus = row.get("virus_name", "")
            umi = float(row.get("viral_umi") or 0)
            if _match_feature(virus, target_regex):
                values[0] += umi
            if _match_feature(virus, off_target_regex):
                values[1] += umi

    anchor: set[str] = set(counts)
    with summary.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"virus_name", "total_cells"}.issubset(reader.fieldnames or []):
            return ParsedMetrics({}, set(), "failed", "viral_summary.tsv missing required columns")
        for row in reader:
            if row.get("total_cells"):
                # ViralScan writes the same denominator on all virus rows.
                break

    return ParsedMetrics(
        {barcode: (values[0], values[1]) for barcode, values in counts.items()},
        anchor,
        "complete",
    )


def _row_metrics(run_dir: Path, row: dict[str, str]) -> ParsedMetrics:
    row_dir = run_dir / "runs" / row["row_id"]
    if row["method"] == "viralscan":
        return parse_viralscan_metrics(
            row_dir,
            target_regex=row["target_regex"],
            off_target_regex=row["off_target_regex"],
        )

    if row["reference_strategy"] == "combined":
        solo_dir = row_dir / "starsolo" / "Solo.out" / "GeneFull"
        raw_dir = solo_dir / "raw"
        filtered = solo_dir / "filtered" / "barcodes.tsv"
    else:
        virus_dir = row_dir / "starsolo_virus" / "Solo.out" / "GeneFull"
        raw_dir = virus_dir / "raw"
        filtered = row_dir / "starsolo_host" / "Solo.out" / "GeneFull" / "filtered" / "barcodes.tsv"
    return parse_starsolo_metrics(
        raw_dir,
        filtered,
        target_regex=row["target_regex"],
        off_target_regex=row["off_target_regex"],
    )


def _latest_slurm_job_id(log_dir: Path, row_index: int) -> str:
    outs = sorted(log_dir.glob(f"vs_ref_strategy_*_{row_index}.out"))
    if not outs:
        return ""
    return outs[-1].name.replace("vs_ref_strategy_", "").removesuffix(".out")


def _read_run_status(run_dir: Path) -> dict[str, tuple[str, str]]:
    status_path = run_dir / "run_status.tsv"
    if not status_path.exists():
        return {}
    with status_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"dataset", "method", "reference_strategy", "status", "failure_reason"}
        if not required.issubset(reader.fieldnames or []):
            return {}
        return {
            f"{row['dataset']}__{row['method']}__{row['reference_strategy']}": (
                row.get("status", ""),
                row.get("failure_reason", ""),
            )
            for row in reader
        }


def parsed_benchmark_rows(run_dir: Path) -> list[dict[str, str]]:
    commands_path = run_dir / "commands.jsonl"
    audit_path = run_dir / "reference_audit.tsv"
    if not commands_path.exists():
        raise BenchmarkContractError(f"missing command manifest: {commands_path}")
    if not audit_path.exists():
        raise BenchmarkContractError(f"missing reference audit: {audit_path}")

    command_rows = [
        json.loads(line)
        for line in commands_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    metrics_by_id = {row["row_id"]: _row_metrics(run_dir, row) for row in command_rows}
    run_status = _read_run_status(run_dir)
    audit_hash = hashlib.sha256(audit_path.read_bytes()).hexdigest()

    fixed_universe_by_dataset: dict[str, set[str]] = {}
    shared_anchor_by_dataset: dict[str, set[str]] = {}
    for dataset in DATASETS:
        rows = [row for row in command_rows if row["dataset"] == dataset["dataset"]]
        combined_sets = [
            set(metrics_by_id[row["row_id"]].counts_by_barcode)
            for row in rows
            if row["reference_strategy"] == "combined"
            and metrics_by_id[row["row_id"]].status == "complete"
        ]
        fixed = set().union(*combined_sets) if combined_sets else set()
        fixed_universe_by_dataset[dataset["dataset"]] = fixed
        anchor_sets = [
            metrics_by_id[row["row_id"]].anchor_barcodes
            for row in rows
            if metrics_by_id[row["row_id"]].status == "complete"
        ]
        shared_anchor_by_dataset[dataset["dataset"]] = (
            set.intersection(*anchor_sets) if anchor_sets else set()
        )

    out_rows: list[dict[str, str]] = []
    target_by_dataset_strategy: dict[tuple[str, str, str], float] = {}
    for idx, row in enumerate(command_rows):
        metrics = metrics_by_id[row["row_id"]]
        status = metrics.status
        failure_reason = metrics.failure_reason if metrics.status != "complete" else ""
        status_evidence = run_status.get(row["row_id"])
        if metrics.status != "complete" and status_evidence:
            evidence_status, evidence_reason = status_evidence
            if evidence_status and evidence_status != "outputs_present_unparsed":
                status = evidence_status
                failure_reason = evidence_reason or failure_reason
        fixed = fixed_universe_by_dataset[row["dataset"]]
        shared = shared_anchor_by_dataset[row["dataset"]]
        target_umi = sum(metrics.counts_by_barcode.get(barcode, (0.0, 0.0))[0] for barcode in fixed)
        off_umi = sum(metrics.counts_by_barcode.get(barcode, (0.0, 0.0))[1] for barcode in fixed)
        target_pos = sum(
            1 for barcode in fixed if metrics.counts_by_barcode.get(barcode, (0.0, 0.0))[0] > 0
        )
        off_pos = sum(
            1 for barcode in fixed if metrics.counts_by_barcode.get(barcode, (0.0, 0.0))[1] > 0
        )
        target_by_dataset_strategy[(row["dataset"], row["method"], row["reference_strategy"])] = (
            target_umi
        )
        delta = ""
        pair_key = (
            row["dataset"],
            row["method"],
            "combined" if row["reference_strategy"] == "two_step" else "two_step",
        )
        if pair_key in target_by_dataset_strategy:
            if row["reference_strategy"] == "two_step":
                delta = str(target_umi - target_by_dataset_strategy[pair_key])
            else:
                delta = str(target_by_dataset_strategy[pair_key] - target_umi)

        out_rows.append(
            {
                "dataset": row["dataset"],
                "srr": row["srr"],
                "technology": row["technology"],
                "method": row["method"],
                "reference_strategy": row["reference_strategy"],
                "target_virus": row["target_virus"],
                "reference_hash_id": f"reference_audit_sha256:{audit_hash}",
                "command_id": row["row_id"],
                "slurm_job_id": _latest_slurm_job_id(run_dir / "logs", idx),
                "count_layer": "GeneFull.raw"
                if row["method"] == "starsolo"
                else "per_cell_viral.viral_umi",
                "barcode_universe": "combined_starsolo_raw_union_plus_viralscan_combined_barcodes",
                "denominator_fixed_barcodes": str(len(fixed)),
                "denominator_method_called_cells": str(len(metrics.anchor_barcodes)),
                "denominator_shared_anchor_cells": str(len(shared)),
                "target_positive_cells_fixed": str(target_pos),
                "target_umi_fixed": str(target_umi),
                "related_off_target_positive_cells_fixed": str(off_pos),
                "related_off_target_umi_fixed": str(off_umi),
                "combined_vs_two_step_delta_target_umi": delta or "0",
                "status": status,
                "failure_reason": failure_reason if status != "complete" else "",
            }
        )
    return out_rows


def write_parsed_benchmark_results(run_dir: Path, out_path: Path) -> None:
    rows = parsed_benchmark_rows(run_dir)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUIRED_RESULT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def validate_results(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_cols = [
            col for col in REQUIRED_RESULT_COLUMNS if col not in (reader.fieldnames or [])
        ]
        if missing_cols:
            raise BenchmarkContractError(f"missing required columns: {missing_cols}")
        rows = [dict(row) for row in reader]

    expected = {(r["dataset"], r["method"], r["reference_strategy"]) for r in benchmark_rows()}
    found = {(r["dataset"], r["method"], r["reference_strategy"]) for r in rows}
    if found != expected:
        raise BenchmarkContractError(
            f"benchmark rows mismatch; missing={sorted(expected - found)} extra={sorted(found - expected)}"
        )

    text = path.read_text(encoding="utf-8")
    hits = [p for p in FORBIDDEN_REFERENCE_PATTERNS if p in text]
    if hits:
        raise BenchmarkContractError(
            f"forbidden historical/single-virus references in results: {hits}"
        )

    for row in rows:
        if not row["count_layer"]:
            raise BenchmarkContractError(f"{row} has empty count_layer")
        if row["status"] == "complete":
            for col in ("reference_hash_id", "command_id", "slurm_job_id"):
                if not row[col]:
                    raise BenchmarkContractError(
                        f"{row['dataset']} {row['method']} {row['reference_strategy']} missing {col}"
                    )
            if row["failure_reason"]:
                raise BenchmarkContractError(
                    f"{row['dataset']} complete row has non-empty failure_reason"
                )
        for col in (
            "denominator_fixed_barcodes",
            "denominator_method_called_cells",
            "denominator_shared_anchor_cells",
            "target_positive_cells_fixed",
            "target_umi_fixed",
            "related_off_target_positive_cells_fixed",
            "related_off_target_umi_fixed",
            "combined_vs_two_step_delta_target_umi",
        ):
            if row[col] == "":
                raise BenchmarkContractError(
                    f"{row['dataset']} {row['method']} {row['reference_strategy']} missing {col}"
                )
            try:
                float(row[col])
            except ValueError as exc:
                raise BenchmarkContractError(f"{col} must be numeric, got {row[col]!r}") from exc
        if row["reference_strategy"] == "two_step" and row["barcode_universe"] in {
            "viral_only_called_cells",
            "",
        }:
            raise BenchmarkContractError(
                "two-step viral-only denominator must use host/combined/external barcode universe"
            )

    fixed_by_dataset: dict[str, str] = {}
    shared_by_dataset: dict[str, str] = {}
    for row in rows:
        dataset = row["dataset"]
        fixed = row["denominator_fixed_barcodes"]
        shared = row["denominator_shared_anchor_cells"]
        fixed_by_dataset.setdefault(dataset, fixed)
        shared_by_dataset.setdefault(dataset, shared)
        if fixed_by_dataset[dataset] != fixed:
            raise BenchmarkContractError(f"{dataset} has inconsistent fixed barcode denominator")
        if shared_by_dataset[dataset] != shared:
            raise BenchmarkContractError(f"{dataset} has inconsistent shared-anchor denominator")
    return rows


def write_run_status_from_results(results_path: Path, out_path: Path) -> None:
    rows = validate_results(results_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["dataset", "srr", "method", "reference_strategy", "status", "failure_reason"]
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def default_manifest() -> dict[str, Any]:
    """Return a conservative manifest template populated with known local defaults."""
    return {
        "panel": PANEL_ID,
        "created_for": "reference_strategy_benchmark",
        "human": {
            "source_release": "GRCh38-2024-A",
            "genome_fasta": "/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa",
            "genes_gtf": "/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf",
        },
        "viral_panel": {
            "id": PANEL_ID,
            "source": "Serratus plus ViralScan expanded anellovirus accession table",
            "anellovirus_expected_count": 2022,
            "anellovirus_accession_table": "src/viralscan/data/anellovirus_accessions.tsv",
        },
        "fastq_root": "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data",
        "fastqs": {
            "SRR20710641": {
                "source_url": "https://www.ebi.ac.uk/ena/browser/view/SRR20710641",
                "R1": "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data/hhv6_carT_ref/SRR20710641/SRR20710641_1.fastq.gz",
                "R2": "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data/hhv6_carT_ref/SRR20710641/SRR20710641_2.fastq.gz",
            },
            "SRR12682296": {
                "R1": {
                    "path": "/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/reference_strategy/SRR12682296/SRR12682296_1.fastq.gz",
                    "source_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/096/SRR12682296/SRR12682296_1.fastq.gz",
                    "md5": "dd1bfe5861d10c89e0f67b836d3ea262",
                    "bytes": "2273612250",
                },
                "R2": {
                    "path": "/exports/para-lipg-hpc/mdmanurung/ViralScan/benchmark_inputs/reference_strategy/SRR12682296/SRR12682296_2.fastq.gz",
                    "source_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/096/SRR12682296/SRR12682296_2.fastq.gz",
                    "md5": "54c5806e695152f7d03a01d8159f8261",
                    "bytes": "8672883135",
                },
            },
            "SRR8315713": {
                "source_url": "https://www.ebi.ac.uk/ena/browser/view/SRR8315713",
                "R1": "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data/hsv1_fibroblast/SRR8315713/SRR8315713_1.fastq.gz",
                "R2": "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/data/hsv1_fibroblast/SRR8315713/SRR8315713_2.fastq.gz",
            },
        },
        "references": {
            "starsolo": {
                "human_source_release": "GRCh38-2024-A",
                "human_only": {"genome_dir": "references/starsolo/human_GRCh38_2024A"},
                "all_virus": {
                    "genome_fasta": "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun/refs/merged/viral_serratus_plus_anellovirus.fa",
                    "genome_gtf": "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun/refs/merged/viral_serratus_plus_anellovirus.gtf",
                    "genome_dir": "references/starsolo/all_virus_serratus_plus_anellovirus",
                },
                "combined": {
                    "genome_dir": "references/starsolo/combined_GRCh38_2024A_serratus_plus_anellovirus"
                },
            },
            "viralscan": {
                "human_source_release": "GRCh38-2024-A",
                "human_only": {"kallisto_index": "references/kallisto/human_GRCh38_2024A.idx"},
                "all_virus": {
                    "kallisto_index": "references/kallisto/all_virus_serratus_plus_anellovirus.idx",
                    "t2g": "references/kallisto/all_virus_serratus_plus_anellovirus.t2g.tsv",
                    "gtf": "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun/refs/merged/viral_serratus_plus_anellovirus.gtf",
                },
                "combined": {
                    "kallisto_index": "references/kallisto/combined_GRCh38_2024A_serratus_plus_anellovirus.idx",
                    "t2g": "references/kallisto/combined_GRCh38_2024A_serratus_plus_anellovirus.t2g.tsv",
                    "gtf": "references/kallisto/combined_GRCh38_2024A_serratus_plus_anellovirus.gtf",
                },
            },
        },
        "provenance": {
            "note": "Template must be audited before running; missing relative reference paths are build targets, not completed artifacts.",
        },
    }
