"""Merge 31-mer-identical repeat copies of a gene into one ``gene_id`` (PLAN ``PROG-13``).

Herpesvirus genomes carry their repeats twice (HSV-1 TRL/IRL and TRS/IRS, VZV TRS/IRS), and
RefSeq annotates each copy as its own gene. A read from a copy cannot be placed on either
one, so kallisto reports it as a two-gene equivalence class and the uniquely-placing layer
drops it: LAT, ICP0 and ICP4 (HSV-1) and ORF62/ORF63 (VZV) can never reach layer 2.

The fix is in the index's ``t2g``: rewrite the copy's *gene* column to the surviving gene, so
both copies' transcripts are one gene to ``kb count`` and to the multimapping correction. The
index itself is untouched; only ``t2g.txt`` / ``panel.t2g`` changes, which is why the same
rewrite can be applied to an index that already exists
(``scripts/apply_repeat_gene_merges.py``).

The table is ``data/repeat_gene_merges.tsv``. It lists each ID in the two spellings the
panels use (the bundled GTFs' prefixed ``HUM_HERP1_HHV1gp00s01`` and the VIRTUS-sourced
indices' bare ``HHV1gp00s01``); an index that uses neither is left alone.
"""

from __future__ import annotations

import csv
import importlib.resources
import os
from pathlib import Path

_TSV_FILENAME = "repeat_gene_merges.tsv"
_PACKAGE_DATA = "viralscan.data"
_COLUMNS = ("gene_id", "merged_into", "virus", "rationale")

#: File names a kb index directory may use for its t2g, in preference order.
T2G_NAMES = ("t2g.txt", "panel.t2g")


def _default_tsv_path() -> Path:
    return Path(str(importlib.resources.files(_PACKAGE_DATA).joinpath(_TSV_FILENAME)))


def load_repeat_merges(path: str | Path | None = None) -> dict[str, str]:
    """``{copy gene_id: surviving gene_id}`` from the merge table.

    Raises :class:`ValueError` for a malformed table. A ``merged_into`` that is itself a
    copy would make a second pass rewrite again, so it is refused: that is what keeps
    :func:`apply_repeat_merges` idempotent.
    """
    tsv_path = Path(path) if path is not None else _default_tsv_path()
    with open(tsv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _COLUMNS:
            raise ValueError(f"{tsv_path}: expected columns {_COLUMNS}, got {reader.fieldnames}")
        merges: dict[str, str] = {}
        for row in reader:
            copy, survivor = row["gene_id"].strip(), row["merged_into"].strip()
            if not copy or not survivor or copy == survivor:
                raise ValueError(f"{tsv_path}: bad row {row}")
            if merges.setdefault(copy, survivor) != survivor:
                raise ValueError(f"{tsv_path}: {copy} is merged into two genes")
    chained = sorted(set(merges) & set(merges.values()))
    if chained:
        raise ValueError(f"{tsv_path}: {chained} are both a copy and a survivor")
    return merges


def apply_repeat_merges(t2g_path: str | Path, merges: dict[str, str] | None = None) -> int:
    """Rewrite the gene column of a kb ``t2g`` in place; return the number of rows changed.

    Only column 2 of a row whose gene is a listed copy changes; every other byte (host
    rows with their empty columns 3-4, line endings, rows of other viruses) is written
    back unchanged. Idempotent: a surviving gene is never a key, so a second call changes
    nothing, and the file is not rewritten when nothing changes.
    """
    if merges is None:
        merges = load_repeat_merges()
    path = Path(t2g_path)
    changed = 0
    out: list[str] = []
    # surrogateescape + newline="" make the round trip byte-exact for unchanged rows.
    with open(path, encoding="utf-8", errors="surrogateescape", newline="") as handle:
        for line in handle:
            body = line.rstrip("\r\n")
            ending = line[len(body) :]
            # read_t2g's split: strict tabs, whitespace only for a tab-less legacy file.
            cols = body.split("\t") if "\t" in body else body.split()
            if len(cols) >= 2 and cols[1] in merges:
                cols[1] = merges[cols[1]]
                line = "\t".join(cols) + ending
                changed += 1
            out.append(line)
    if changed:
        tmp = path.with_name(path.name + ".tmp")
        with open(tmp, "w", encoding="utf-8", errors="surrogateescape", newline="") as handle:
            handle.writelines(out)
        os.replace(tmp, path)
    return changed


def find_t2g(index_dir: str | Path) -> Path:
    """The t2g inside a kb index directory (``t2g.txt``, else ``panel.t2g``)."""
    directory = Path(index_dir)
    for name in T2G_NAMES:
        if (directory / name).is_file():
            return directory / name
    raise FileNotFoundError(f"no {' or '.join(T2G_NAMES)} in {directory}")
