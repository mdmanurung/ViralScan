"""Per-sample ``results/artifact_inventory.json`` and ``README.txt`` (PLAN ``LAYOUT-V3-01``).

The inventory lists the files of a sample with their sha256, in the shape of
``schemas/v3/artifact_inventory.schema.json``. What it skips is explicit:

* itself, ``posthoc/`` (post-hoc verbs write their own manifests) and ``.snakemake/``;
* ``.done`` sentinels (empty files);
* the bulky machine-only work dirs of ``outputs.SAMPLE_LAYOUT`` (role ``intermediate``),
  except their ``*.h5ad`` matrices, which users load.

Written by ``run_safety.write_run_complete``, so it is refreshed with ``run_complete.json``.
"""

import json
import re
from pathlib import Path
from typing import Any

from viralscan import __version__
from viralscan.outputs import SAMPLE_LAYOUT, V3Outputs, run_readme, sample_readme

INVENTORY_SCHEMA_VERSION = "1.1.0"
COUNTING_VERSION = "3.0.0"  # the ``count_schema_version`` stamped into the h5ad files
INTERMEDIATE_DIRS = frozenset(
    p.rstrip("/") for p, role, _ in SAMPLE_LAYOUT if role == "intermediate"
)
SKIP_DIRS = INTERMEDIATE_DIRS | {"posthoc", ".snakemake"}
KEEP_SUFFIXES = (".h5ad",)  # matrices inside a skipped work dir
INVENTORY_REL = "results/artifact_inventory.json"


def is_ignored(rel: str) -> bool:
    """True when the sample-relative path ``rel`` is deliberately left out of the inventory."""
    top = Path(rel).parts[0]
    if rel == INVENTORY_REL or rel.endswith(".done") or top in ("posthoc", ".snakemake"):
        return True
    return top in INTERMEDIATE_DIRS and not rel.endswith(KEEP_SUFFIXES)


def inventory_paths(sample_dir: Path) -> list[str]:
    """Sorted sample-relative paths of the files the inventory lists."""
    rels = (p.relative_to(sample_dir).as_posix() for p in sample_dir.rglob("*") if p.is_file())
    return sorted(rel for rel in rels if not is_ignored(rel))


def _manifest_options(run_root: Path) -> dict[str, Any]:
    """The run manifest (``{}`` when absent), whose ``options`` gate the README's file list."""
    path = Path(run_root) / "run_manifest.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}


def write_artifact_inventory(sample_dir: Path, options: dict[str, Any] | None = None) -> Path:
    """Write ``README.txt`` and ``results/artifact_inventory.json`` for one sample."""
    from viralscan.run_safety import sha256_file, software_identity

    sample_dir = Path(sample_dir)
    (sample_dir / "README.txt").write_text(
        sample_readme(sample_dir.name, options=options), encoding="utf-8"
    )
    git_sha = software_identity().get("git_commit")
    artifacts: list[dict[str, Any]] = []
    for rel in inventory_paths(sample_dir):
        path = sample_dir / rel
        if rel == "config.yaml":
            kind = "command"
        else:
            kind = "intermediate" if Path(rel).parts[0] in INTERMEDIATE_DIRS else "result"
        artifacts.append(
            {
                "schema_version": INVENTORY_SCHEMA_VERSION,
                "artifact_id": re.sub(r"[^a-z0-9._-]+", "-", rel.lower()).lstrip("-._") or "file",
                "kind": kind,
                "path": rel,
                "scope": "private",
                "software_version": __version__,
                "counting_version": COUNTING_VERSION,
                "git_sha": git_sha,
                "input_sha256s": [],
                "reference_sha256s": [],
                "artifact_sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
                "schema": "unspecified",
                "layer": "unspecified",
                "denominator": "unspecified",
                "generation_command": "viralscan",
                "rebuild_eligibility": "not_applicable" if rel == "config.yaml" else "rebuildable",
                "claim_ids": [],
            }
        )
    target = V3Outputs(sample_dir).artifact_inventory
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = target.with_suffix(".json.tmp")
    payload = {"schema_version": INVENTORY_SCHEMA_VERSION, "artifacts": artifacts}
    staging.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    staging.replace(target)
    return target


def write_layout_docs(run_root: Path) -> None:
    """Run-level ``README.txt`` plus every sample's README and inventory (run finished)."""
    run_root = Path(run_root)
    manifest = _manifest_options(run_root)
    samples = sorted(p.parent for p in run_root.glob("*/config.yaml"))
    for sample_dir in samples:
        write_artifact_inventory(sample_dir, manifest.get("options"))
    (run_root / "README.txt").write_text(
        run_readme([p.name for p in samples], manifest=manifest), encoding="utf-8"
    )


def write_start_docs(run_root: Path, samples: list[str]) -> None:
    """READMEs written before any work, so a failed or running run explains itself.

    Marked NOT FINISHED; ``write_layout_docs`` replaces them when the run completes.
    """
    run_root = Path(run_root)
    manifest = _manifest_options(run_root)
    for sample in samples:
        (run_root / sample).mkdir(parents=True, exist_ok=True)
        (run_root / sample / "README.txt").write_text(
            sample_readme(sample, finished=False, options=manifest.get("options")),
            encoding="utf-8",
        )
    (run_root / "README.txt").write_text(
        run_readme(samples, finished=False, manifest=manifest), encoding="utf-8"
    )
