#!/usr/bin/env python3
"""Package Serratus plus anellovirus references for STARsolo.

This creates STAR-compatible FASTA/GTF inputs only. Building the STAR genome
directories still requires a local STAR binary and should be done via the
commands recorded in ``reference_manifest.json``.
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path


DEFAULT_SERRATUS_FASTA = Path(
    "/exports/archive/hg-funcgenom-research/evonk/old/intern/fasta_viruses/Serratus/all_fastas.fasta"
)
DEFAULT_SERRATUS_FASTA_DIR = Path(
    "/exports/archive/hg-funcgenom-research/evonk/old/intern/fasta_viruses/Serratus/fasta_split"
)
DEFAULT_SERRATUS_GTF = Path(
    "/exports/archive/hg-funcgenom-research/evonk/old/intern/fasta_viruses/Serratus/all_gtf.gtf"
)
DEFAULT_ANELLO_FASTA = Path(
    "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun/refs/anellovirus/anellovirus.fa"
)
DEFAULT_ANELLO_GTF = Path(
    "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun/refs/anellovirus/anellovirus.gtf"
)
DEFAULT_HUMAN_FASTA = Path(
    "/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
)
DEFAULT_HUMAN_GTF = Path(
    "/exports/archive/hg-funcgenom-research/evonk/old/intern/human_reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
)
DEFAULT_VIRALSCAN_TRANSCRIPTOME = Path(
    "/exports/para-lipg-hpc/mdmanurung/viralscan_showcase/fullrun/refs/merged/transcriptome_plus_anellovirus.fa"
)


TARGET_PATTERNS = {
    "HHV-6": re.compile(
        r"(Human[_ -]herpesvirus[_ -]?6|Human[_ -]betaherpesvirus[_ -]?6|"
        r"HUM_HERP6[AB]?|hhv-?6[ab]?|NC_001664|HHV6)",
        re.I,
    ),
    "EBV": re.compile(r"(Epstein[_ -]Barr|EPSTEIN_HHV4|ebv|hhv-?4|NC_007605|HHV4)", re.I),
    "HSV-1": re.compile(
        r"(Human[_ -]herpesvirus[_ -]?1|HUM_HERP1|hsv-?1|NC_001806|HHV1)",
        re.I,
    ),
}


def fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def gtf_seqnames(path: Path) -> set[str]:
    seqnames: set[str] = set()
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            seqnames.add(line.split("\t", 1)[0])
    return seqnames


def copy_many(inputs: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as out:
        for source in inputs:
            with source.open("rb") as handle:
                shutil.copyfileobj(handle, out)
            out.write(b"\n")


def append_hhv6b_pseudocontigs(transcriptome: Path, fasta: Path, gtf: Path) -> int:
    """Append ViralScan/kallisto HHV-6B transcript targets as STAR contigs."""
    count = 0
    current_id: str | None = None
    current_chunks: list[str] = []

    def flush_record() -> None:
        nonlocal count, current_id, current_chunks
        if current_id is None or not current_id.startswith("HUM_HERP6B"):
            return
        sequence = "".join(current_chunks).replace(" ", "").replace("\t", "")
        if not sequence:
            return
        with fasta.open("a", encoding="utf-8") as fasta_out:
            fasta_out.write(f">{current_id} source:ViralScan_kallisto_HHV6B_pseudocontig\n")
            for start in range(0, len(sequence), 80):
                fasta_out.write(sequence[start : start + 80] + "\n")
        with gtf.open("a", encoding="utf-8") as gtf_out:
            attrs = f'gene_id "{current_id}"; transcript_id "{current_id}"; gene_name "{current_id}";'
            gtf_out.write(
                f"{current_id}\tViralScan_kallisto\tgene\t1\t{len(sequence)}\t.\t+\t.\t{attrs}\n"
            )
            gtf_out.write(
                f"{current_id}\tViralScan_kallisto\ttranscript\t1\t{len(sequence)}\t.\t+\t.\t{attrs}\n"
            )
            gtf_out.write(
                f"{current_id}\tViralScan_kallisto\texon\t1\t{len(sequence)}\t.\t+\t.\t{attrs} exon_number \"1\";\n"
            )
        count += 1

    with transcriptome.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                flush_record()
                current_id = line[1:].split()[0]
                current_chunks = []
            elif current_id is not None:
                current_chunks.append(line.strip())
    flush_record()
    return count


def fasta_inputs(single_fasta: Path | None, fasta_dir: Path | None) -> list[Path]:
    if fasta_dir is not None:
        paths = sorted(fasta_dir.glob("*.fasta")) + sorted(fasta_dir.glob("*.fa"))
        if not paths:
            raise SystemExit(f"No FASTA files found in {fasta_dir}")
        return paths
    if single_fasta is None:
        raise SystemExit("Either --serratus-fasta or --serratus-fasta-dir is required")
    return [single_fasta]


def require_existing(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise SystemExit("Missing source paths:\n" + "\n".join(missing))


def validate_pair(fasta: Path, gtf: Path) -> None:
    ids = fasta_ids(fasta)
    seqnames = gtf_seqnames(gtf)
    missing = sorted(seqnames - ids)
    if missing:
        preview = "\n".join(missing[:25])
        raise SystemExit(
            f"{gtf} contains {len(missing)} seqnames absent from {fasta}; first entries:\n{preview}"
        )


def filter_gtf_to_fasta(gtf: Path, fasta: Path, dropped_out: Path) -> None:
    ids = fasta_ids(fasta)
    dropped: set[str] = set()
    tmp = gtf.with_suffix(gtf.suffix + ".tmp")
    with gtf.open(encoding="utf-8", errors="ignore") as source, tmp.open(
        "w", encoding="utf-8"
    ) as out:
        for line in source:
            if not line.strip() or line.startswith("#"):
                out.write(line)
                continue
            seqname = line.split("\t", 1)[0]
            if seqname in ids:
                out.write(line)
            else:
                dropped.add(seqname)
    tmp.replace(gtf)
    dropped_out.parent.mkdir(parents=True, exist_ok=True)
    dropped_out.write_text("\n".join(sorted(dropped)) + ("\n" if dropped else ""), encoding="utf-8")


def target_presence(path: Path) -> dict[str, bool]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    return {name: bool(pattern.search(text)) for name, pattern in TARGET_PATTERNS.items()}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, default=Path("references/starsolo"))
    parser.add_argument("--serratus-fasta", type=Path, default=None)
    parser.add_argument("--serratus-fasta-dir", type=Path, default=DEFAULT_SERRATUS_FASTA_DIR)
    parser.add_argument("--serratus-gtf", type=Path, default=DEFAULT_SERRATUS_GTF)
    parser.add_argument("--anellovirus-fasta", type=Path, default=DEFAULT_ANELLO_FASTA)
    parser.add_argument("--anellovirus-gtf", type=Path, default=DEFAULT_ANELLO_GTF)
    parser.add_argument("--viralscan-transcriptome", type=Path, default=DEFAULT_VIRALSCAN_TRANSCRIPTOME)
    parser.add_argument("--human-fasta", type=Path, default=DEFAULT_HUMAN_FASTA)
    parser.add_argument("--human-gtf", type=Path, default=DEFAULT_HUMAN_GTF)
    args = parser.parse_args(argv)

    require_existing(
        [
            args.serratus_fasta_dir if args.serratus_fasta_dir is not None else args.serratus_fasta,
            args.serratus_gtf,
            args.anellovirus_fasta,
            args.anellovirus_gtf,
            args.viralscan_transcriptome,
            args.human_fasta,
            args.human_gtf,
        ]
    )

    all_virus_dir = args.out_root / "all_virus_serratus_plus_anellovirus"
    combined_dir = args.out_root / "combined_GRCh38_2024A_serratus_plus_anellovirus"
    viral_fasta = all_virus_dir / "viral_genome.fa"
    viral_gtf = all_virus_dir / "viral_genome.gtf"
    combined_fasta = combined_dir / "combined_genome.fa"
    combined_gtf = combined_dir / "combined.gtf"

    serratus_fastas = fasta_inputs(args.serratus_fasta, args.serratus_fasta_dir)
    copy_many([*serratus_fastas, args.anellovirus_fasta], viral_fasta)
    copy_many([args.serratus_gtf, args.anellovirus_gtf], viral_gtf)
    hhv6b_count = append_hhv6b_pseudocontigs(args.viralscan_transcriptome, viral_fasta, viral_gtf)
    if hhv6b_count == 0:
        raise SystemExit(f"No HUM_HERP6B records found in {args.viralscan_transcriptome}")
    filter_gtf_to_fasta(viral_gtf, viral_fasta, all_virus_dir / "dropped_gtf_seqnames.tsv")
    validate_pair(viral_fasta, viral_gtf)

    copy_many([args.human_fasta, viral_fasta], combined_fasta)
    copy_many([args.human_gtf, viral_gtf], combined_gtf)
    validate_pair(combined_fasta, combined_gtf)

    presence = target_presence(viral_gtf)
    missing_targets = [name for name, present in presence.items() if not present]
    if missing_targets:
        raise SystemExit(f"Packaged viral GTF is missing expected targets: {missing_targets}")

    print(f"Wrote {viral_fasta}")
    print(f"Wrote {viral_gtf}")
    print(f"Appended {hhv6b_count} ViralScan/kallisto HHV-6B pseudo-contigs")
    print(f"Wrote {combined_fasta}")
    print(f"Wrote {combined_gtf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
