#!/usr/bin/env python3
"""Run ``kb count --aa`` for one library of a ``make_configs.py`` config against a prebuilt index.

Same command as the Snakefile ``count`` rule (``--parity single``, technology and optional
on-list from the config), plus optional ``--num`` (read numbers in the BUS flag column, so
suspicious hits can be traced to reads later) and a receipt with versions, hashes, FASTQ paths,
timing and the pseudoalignment summary. The BUS files are kept.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml


def sha256_file(path: Path) -> str:  # standalone: compute nodes only see a copy of this script
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def command(cfg: dict, sample: str, index: Path, t2g: Path, out: Path, threads: int, num: bool):
    if sample not in cfg["samples"]:
        raise SystemExit(f"{sample!r} not in config samples")
    return [
        cfg["kb"], "count", "--kallisto", cfg["kallisto"], "--aa", "--parity", "single",
        "-x", cfg["technology"], "-t", str(threads),
        *(["-w", cfg["whitelist"]] if cfg.get("whitelist") else []),
        *(["--num"] if num else []),
        "--tmp", str(out.with_name(out.name + "_tmp")),
        "-i", str(index), "-g", str(t2g), "--h5ad", "-o", str(out), *cfg["samples"][sample],
    ]  # fmt: skip


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--t2g", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True, help="new or empty directory")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--num", action="store_true", help="pass kb count --num")
    args = ap.parse_args()
    cfg = yaml.safe_load(args.config.read_text())
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise SystemExit(f"Output not empty: {args.out_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cmd = command(cfg, args.sample, args.index, args.t2g, args.out_dir, args.threads, args.num)
    start = datetime.now(timezone.utc).isoformat()
    with open(args.out_dir.with_suffix(".log"), "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    end = datetime.now(timezone.utc).isoformat()
    run_info = json.loads((args.out_dir / "run_info.json").read_text())
    receipt = {
        "sample": args.sample,
        "script_sha256": sha256_file(Path(__file__)),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "command": cmd,
        "technology": cfg["technology"],
        "whitelist": cfg.get("whitelist"),
        "index_sha256": sha256_file(args.index),
        "t2g_sha256": sha256_file(args.t2g),
        "kallisto_version": subprocess.check_output(
            [cfg["kallisto"], "version"], text=True
        ).strip(),
        "kb_python_version": subprocess.check_output(
            [
                sys.executable,
                "-c",
                "from importlib.metadata import version; print(version('kb-python'))",
            ],
            text=True,
        ).strip(),
        "fastq": {p: Path(p).stat().st_size for p in cfg["samples"][args.sample]},
        "start_utc": start,
        "end_utc": end,
        "n_processed": run_info.get("n_processed"),
        "n_pseudoaligned": run_info.get("n_pseudoaligned"),
        "p_pseudoaligned": run_info.get("p_pseudoaligned"),
        "matrix_sha256": sha256_file(args.out_dir / "counts_unfiltered" / "adata.h5ad"),
    }
    (args.out_dir / "translated_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
