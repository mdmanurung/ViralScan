#!/usr/bin/env python3
"""Build a gitignored, symlink-only curation view of the repo.

The view under ``curation/`` is a disposable *lens* over the real files — it is
generated from ``scripts/curation_manifest.yaml`` and contains only relative
symlinks. Originals are never moved, renamed, or edited.

Usage:
    python3 scripts/build_curation_view.py          # (re)build curation/
    python3 scripts/build_curation_view.py --clean  # remove curation/

Design guarantees (non-breaking):
  * symlinks only — no copies, no moves;
  * relative links — survive a repo relocation as long as internal layout holds;
  * ``curation/`` is gitignored and fully regenerable;
  * missing manifest targets are reported and skipped (never a dangling link);
  * ``--clean`` and rebuild refuse to touch a ``curation/`` that lacks our marker.
"""

from __future__ import annotations

import argparse
import glob as globmod
import os
import shutil
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VIEW_DIR = os.path.join(REPO_ROOT, "curation")
MANIFEST = os.path.join(REPO_ROOT, "scripts", "curation_manifest.yaml")
MARKER = ".curation-generated"  # guards against clobbering a non-generated dir


def _load_manifest(path: str) -> dict:
    try:
        import yaml
    except ImportError:
        sys.exit("PyYAML is required. Install it, or run inside the viralscan env.")
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def _expand(targets: list[str]) -> list[str]:
    """Expand globs to repo-root-relative paths that exist, preserving order."""
    seen: dict[str, None] = {}
    for t in targets:
        matches = sorted(globmod.glob(os.path.join(REPO_ROOT, t)))
        if matches:
            for m in matches:
                seen[os.path.relpath(m, REPO_ROOT)] = None
        else:
            seen[t] = None  # keep so we can warn about the miss
    return list(seen)


def _symlink(rel_target: str, link_dir: str, warnings: list[str]) -> bool:
    """Create link_dir/<basename> -> rel_target (relative). Return True if made."""
    abs_target = os.path.join(REPO_ROOT, rel_target)
    if not os.path.exists(abs_target):
        warnings.append(f"missing target, skipped: {rel_target}")
        return False
    os.makedirs(link_dir, exist_ok=True)
    link_path = os.path.join(link_dir, os.path.basename(rel_target.rstrip("/")))
    if os.path.lexists(link_path):  # dedup (e.g. globs overlapping explicit paths)
        return False
    rel = os.path.relpath(abs_target, link_dir)
    os.symlink(rel, link_path)
    return True


def _clean() -> None:
    if not os.path.isdir(VIEW_DIR):
        return
    if not os.path.exists(os.path.join(VIEW_DIR, MARKER)):
        sys.exit(f"refusing to remove {VIEW_DIR}: no {MARKER} marker (not our tree)")
    shutil.rmtree(VIEW_DIR)


def build() -> None:
    _clean()
    manifest = _load_manifest(MANIFEST)
    os.makedirs(VIEW_DIR, exist_ok=True)
    open(os.path.join(VIEW_DIR, MARKER), "w").close()

    warnings: list[str] = []
    n_links = 0

    # PRIMARY: by-result
    for result in manifest.get("by_result", []):
        rid = result["id"]
        rdir = os.path.join(VIEW_DIR, "by-result", rid)
        os.makedirs(rdir, exist_ok=True)
        if result.get("section"):
            with open(os.path.join(rdir, "SECTION.txt"), "w") as fh:
                fh.write(result["section"].rstrip() + "\n")
        for tgt in _expand(result.get("links", [])):
            n_links += _symlink(tgt, rdir, warnings)

    # SECONDARY: by-type
    for category, targets in (manifest.get("by_type") or {}).items():
        cdir = os.path.join(VIEW_DIR, "by-type", category)
        os.makedirs(cdir, exist_ok=True)
        for tgt in _expand(targets):
            n_links += _symlink(tgt, cdir, warnings)

    _write_readme(manifest)

    print(f"Built curation/ — {n_links} symlinks "
          f"({len(manifest.get('by_result', []))} results, "
          f"{len(manifest.get('by_type') or {})} type buckets).")
    if warnings:
        print(f"\n{len(warnings)} manifest target(s) missing (edit "
              f"scripts/curation_manifest.yaml):")
        for w in warnings:
            print("  -", w)


def _write_readme(manifest: dict) -> None:
    lines = [
        "# curation/ — generated view (do not edit here)",
        "",
        "This tree is **generated** by `scripts/build_curation_view.py` from",
        "`scripts/curation_manifest.yaml`. It contains only relative symlinks to the",
        "real files; nothing here is canonical. It is gitignored and disposable.",
        "",
        "- Rebuild:   `python3 scripts/build_curation_view.py`",
        "- Tear down: `python3 scripts/build_curation_view.py --clean`",
        "- Curate:    edit `scripts/curation_manifest.yaml`, then rebuild.",
        "",
        "## by-result/ (manuscript narratives)",
        "",
    ]
    for r in manifest.get("by_result", []):
        lines.append(f"- `{r['id']}/` — {r.get('section', '')}")
    lines += ["", "## by-type/ (flat scanning buckets)", ""]
    for cat in (manifest.get("by_type") or {}):
        lines.append(f"- `{cat}/`")
    with open(os.path.join(VIEW_DIR, "README.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="Build the gitignored curation/ symlink view.")
    ap.add_argument("--clean", action="store_true", help="remove curation/ and exit")
    args = ap.parse_args()
    if args.clean:
        _clean()
        print("Removed curation/.")
        return
    build()


if __name__ == "__main__":
    main()
