#!/usr/bin/env python3
"""Write the curated reference manifest and reproducibility audit for a finished panel build.

Purpose: PLAN REF-07 and REF-04. `build_bundled_panel_ref.py` builds the panel but does not
         write `reference_manifest.json`. This script reads a finished build directory and writes
         the per-sequence manifest (accession version, taxonomy, snapshot, retrieval date,
         SHA-256, length, draft licence) for the viral panel, plus the digests of every build
         file. Run it on the frozen build and on a rebuild, then compare the two audits.
Run:     PYTHONPATH=src python scripts/make_panel_reference_manifest.py \\
             --build-dir <panel_builds/viral_ref_cat42d/build> --out-dir <fresh dir> \\
             [--genome-dlist <dlist.fa>]
         python scripts/compare_reference_reproducibility.py A/reference_reproducibility.json \\
             B/reference_reproducibility.json
Inputs (read only, never changed):
  <build-dir>/viral.fa               viral panel FASTA, one record per accession
  <build-dir>/viral_panel.gtf        viral annotation
  <build-dir>/combined.fa            host cDNA + viral FASTA given to kb ref
  <build-dir>/combined.gtf           matching GTF
  <build-dir>/cdna.fa                transcriptome FASTA kb ref wrote
  <build-dir>/panel.t2g              transcript-to-gene map
  <build-dir>/panel.idx              kallisto index (digest recorded, not claimed reproducible)
  <build-dir>/lowcomplexity_mask.tsv low-complexity mask table (optional)
  ~/.cache/viralscan/ncbi            cached GenBank/FASTA payloads and retrieval receipts
Outputs (in --out-dir, which must not exist):
  reference_manifest.json            viral rows with provenance; build-file digests
  reference_reproducibility.json     `content_files` digests + `content_sha256`, for the compare script
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from viralscan.run_safety import sha256_file  # noqa: E402
from viralscan.scripts import build_reference as br  # noqa: E402
from viralscan.scripts.ncbi_fetch import cached_accession_provenance  # noqa: E402

DEFAULT_CACHE = Path.home() / ".cache" / "viralscan" / "ncbi"

# Text outputs that must be byte-identical across rebuilds. panel.idx is binary and
# version-dependent, so it is recorded in build_receipt.index but kept out of the audit.
AUDITED_FILES = {
    "viral_fasta": "viral.fa",
    "viral_gtf": "viral_panel.gtf",
    "combined_fasta": "combined.fa",
    "combined_gtf": "combined.gtf",
    "cdna_fasta": "cdna.fa",
    "t2g": "panel.t2g",
    "lowcomplexity_mask": "lowcomplexity_mask.tsv",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--genome-dlist", type=Path, default=None)
    parser.add_argument("--host-species", default="Homo sapiens")
    args = parser.parse_args()

    # a run dir is never overwritten: mkdir fails if it exists
    args.out_dir.mkdir(parents=True)

    # viral rows: provenance comes from the NCBI cache, never from file timestamps
    viral_fasta = args.build_dir / "viral.fa"
    viral_ids = {name for name, _ in br._fasta_records(viral_fasta)}
    provenance = {
        name: cached_accession_provenance(name, args.cache_dir) for name in sorted(viral_ids)
    }
    manifest_path = br.write_reference_manifest(
        viral_fasta,
        args.out_dir / "reference_manifest.json",
        profile="curated",
        host_species=args.host_species,
        viral_identifiers=viral_ids,
        provenance=provenance,
        genome_dlist=args.genome_dlist,
        source_adapters=["panel_build_dir"],
    )

    # digests of the build files, then one content hash over rows + digests
    manifest = json.loads(manifest_path.read_text())
    files = {
        key: sha256_file(args.build_dir / name)
        for key, name in AUDITED_FILES.items()
        if (args.build_dir / name).is_file()
    }
    manifest["content_files"] = files
    index = args.build_dir / "panel.idx"
    manifest["build_receipt"]["index"] = {
        "path": str(index.resolve()),
        "sha256": sha256_file(index),
        "bytes": index.stat().st_size,
    }
    manifest["build_receipt"]["genome_dlist_applied"] = args.genome_dlist is not None
    manifest["content_sha256"] = br._reference_content_sha256(manifest)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    audit = {"content_files": files, "content_sha256": manifest["content_sha256"]}
    (args.out_dir / "reference_reproducibility.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n"
    )
    print(f"wrote {manifest_path} ({len(viral_ids)} viral rows, {len(files)} audited files)")


if __name__ == "__main__":
    main()
