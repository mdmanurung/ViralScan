"""Layer 2 script: infer viral gene programmes for the viruses layer 1 detected.

Runs after ``detection`` and reads two things:

* ``results/viral_summary.tsv`` — the viruses layer 1 actually detected. Layer 2
  is scoped strictly to these; a virus layer 1 missed produces no programme row,
  and that is layer 1's detection-limit problem rather than something layer 2
  can repair.
* the multimap H5AD — for the ``counts_unique_viral`` layer, which is the
  evidence this layer trusts. Only the called cells detection reported over
  (``results/called_cells.tsv``) are scored.

Outputs ``results/gene_program_summary.tsv`` and ``results/gene_program_cells.tsv``.
See :mod:`viralscan.gene_programs` for why the design is what it is; the short
version is that a naive per-gene comparison cannot separate EBV latency from
lytic replication in this data, so breadth is counted over distinct overlap
groups using uniquely-placing molecules.
"""

import csv
import os
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from viralscan.gene_programs import (
    EVIDENCE_LAYER,
    Marker,
    call_cell_programme,
    load_catalogue,
    resolve_markers,
    summarise_programs,
    validate_catalogue,
    write_program_outputs,
)
from viralscan.run_context import RunContext
from viralscan.scripts.cellcalling import load_called_mask
from viralscan.runconfig import RunConfig
from viralscan.utils import setup_script_logging
from viralscan.virus_grouping import load_run_identity

log = setup_script_logging()

#: Fallback used when the layer-1 summary has no recognised molecule column.
#: Older summaries used ``total_umi``; v3 uses ``viral_molecules_total_est``.
_MOLECULE_COLUMNS = ("viral_molecules_total_est", "total_umi")

config: RunConfig = RunConfig()
output: str = ""
panel_form: str = "bundled"


def _detected_viruses(path: str) -> tuple[list[str], str]:
    """Return ``(virus_names, molecule_column)`` from layer 1's summary.

    Raises if the file is missing or has no recognisable molecule column: this
    script's entire scope comes from that file, so guessing would mean emitting
    programme rows for viruses layer 1 never looked at.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"layer-1 summary not found: {path}. Run the detection step before "
            "gene-programme inference."
        )
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"layer-1 summary is empty: {path}")
    column = next((c for c in _MOLECULE_COLUMNS if c in rows[0]), None)
    if column is None:
        raise ValueError(
            f"layer-1 summary {path} has no recognised molecule column "
            f"(looked for {list(_MOLECULE_COLUMNS)}, found {sorted(rows[0])})"
        )
    return [r["virus_name"] for r in rows], column


def _marker_matrix(adata, markers: list[Marker]):
    """Subset the unique and selected layers to marker columns, in marker order.

    ``call_cell_programme`` takes columns positionally aligned with ``markers``,
    so this is the one place the name-to-position mapping happens. A marker
    absent from the matrix becomes an all-zero column, which contributes no
    breadth rather than raising -- an unresolved marker is reported by
    :func:`viralscan.gene_programs.resolve_markers`, not silently padded.

    Column selection goes through scipy's own ``matrix[:, indices]`` rather than
    a hand-rolled ``indptr``/``indices`` traversal. The traversal is easy to get
    subtly wrong, and getting it wrong fails *silently* as an all-zero matrix --
    which is indistinguishable from a genuinely negative sample. Letting scipy
    interpret its own CSR removes that failure mode.
    """
    positions = {name: i for i, name in enumerate(adata.var_names)}
    n_rows, n_cols = adata.n_obs, len(markers)
    present = [
        (target, positions[marker.var_name])
        for target, marker in enumerate(markers)
        if marker.var_name in positions
    ]
    if not present:
        empty = sparse.csr_matrix((n_rows, n_cols), dtype=np.float32)
        return empty, empty.copy()

    targets = [t for t, _ in present]
    sources = [s for _, s in present]
    # `picked.col` indexes positions within the picked matrix, so map position
    # -> target column rather than original-column -> target.
    target_of_position = {i: t for i, t in enumerate(targets)}

    def subset(layer):
        if layer is None:
            return None
        picked = layer[:, sources].tocoo()
        cols = np.fromiter(
            (target_of_position[j] for j in picked.col),
            dtype=np.intp,
            count=picked.nnz,
        )
        return sparse.csr_matrix((picked.data, (picked.row, cols)), shape=(n_rows, n_cols))

    return subset(adata.layers.get(EVIDENCE_LAYER)), subset(adata.X)


def _panel_form(config: RunConfig) -> str:
    """Guess which panel form the index was built from.

    The catalogue carries all three ID forms and resolution tries each in turn,
    so this only affects probe order, never the result. It is derived from the
    configured GTF path where that is informative, and defaults to ``bundled``,
    whose IDs are the most specific.
    """
    gtf = str(getattr(config, "gtf", "") or "")
    if "starsolo" in gtf:
        return "starsolo"
    if "serratus" in gtf or "merged" in gtf:
        return "merged"
    return "bundled"


def run_one(adata, viruses: list[str], catalogue, min_breadth: int, form: str, identity=None):
    """Build the per-cell call table for the detected ``viruses``."""
    facts = {}
    for row in catalogue:
        facts.setdefault(
            row["virus"],
            {
                "latency_observable": row["latency_observable_in_rna"] == "true",
                "completeness": row["panel_completeness"],
            },
        )
    obs_names = list(adata.obs_names)
    records: list[dict[str, Any]] = []

    for virus in viruses:
        info = facts.get(virus)
        if info is None:
            # Layer 1 detected it; the catalogue has no programme model. Emit a
            # marker row so a reader sees "no model" rather than silence.
            records.append(
                {
                    "barcode": "",
                    "virus_name": virus,
                    "state": "not_applicable",
                    "productive_breadth": 0,
                    "latent_breadth": 0,
                    "selected_state": None,
                    "selected_productive_breadth": 0,
                    "selected_latent_breadth": 0,
                    "evidence_layer": EVIDENCE_LAYER,
                    "latency_not_observable": True,
                    "n_cells": 0,
                }
            )
            log.info("%s: no programme model in the catalogue; recorded as not_applicable", virus)
            continue

        resolved, unresolved = resolve_markers(
            virus, adata.var_names, form, catalogue=catalogue, identity=identity
        )
        if not resolved:
            log.warning(
                "%s: %d catalogue marker(s) but none resolved against this index; "
                "the reference and the catalogue disagree. Skipping.",
                virus,
                len(unresolved),
            )
            continue
        if unresolved:
            log.info(
                "%s: %d marker(s) unresolved (%s)", virus, len(unresolved), unresolved[0].reason
            )

        unique, selected = _marker_matrix(adata, resolved)
        calls = call_cell_programme(
            unique,
            resolved,
            min_breadth=min_breadth,
            latency_observable=info["latency_observable"],
            selected_matrix=selected,
        )
        for barcode, call in zip(obs_names, calls):
            if call["state"] == "indeterminate" and not (
                call["productive_breadth"] or call["latent_breadth"]
            ):
                # No marker evidence at all: the virus was called by layer 1 on
                # molecules we cannot attribute to a programme, which is worth
                # recording but not worth a per-cell row for every barcode.
                continue
            records.append({"barcode": barcode, "virus_name": virus, **call})
        states = pd.Series([c["state"] for c in calls]).value_counts().to_dict()
        log.info("%s: %s", virus, states)

    return pd.DataFrame(
        records,
        columns=[
            "barcode",
            "virus_name",
            "state",
            "productive_breadth",
            "latent_breadth",
            "selected_state",
            "selected_productive_breadth",
            "selected_latent_breadth",
            "evidence_layer",
            "latency_not_observable",
            "n_cells",
        ],
    )


def main(adata_path: str, summary_path: str, done_path: str) -> None:
    catalogue = load_catalogue()
    problems = validate_catalogue(catalogue)
    if problems:
        raise ValueError("gene-programme catalogue is invalid:\n  " + "\n  ".join(problems))
    viruses, molecule_column = _detected_viruses(summary_path)
    log.info("layer 1 detected %d virus row(s): %s", len(viruses), ", ".join(viruses))

    adata = sc.read_h5ad(adata_path)
    # Score the called cells only, the denominator layer 1 reports over; the
    # multimap H5AD also holds every empty droplet (PLAN PROG-17).
    called = load_called_mask(adata, config, output)
    adata = adata[called].copy()
    cells = run_one(
        adata,
        viruses,
        catalogue,
        min_breadth=int(getattr(config, "programme_min_breadth", 2)),
        form=_panel_form(config),
        identity=load_run_identity(output),
    )
    summary = summarise_programs(
        cells,
        catalogue,
        min_breadth=int(getattr(config, "programme_min_breadth", 2)),
        viruses=viruses,
    )
    _attach_layer1_totals(summary, summary_path, molecule_column)
    summary["n_called_cells"] = adata.n_obs

    paths = write_program_outputs(cells, summary, output)
    log.info("Wrote %s and %s", os.path.basename(paths[0]), os.path.basename(paths[1]))
    if done_path:
        with open(done_path, "w") as handle:
            handle.write("done\n")


def _attach_layer1_totals(summary, layer1_path: str, molecule_column: str) -> None:
    """Carry layer 1's molecule total onto each programme row.

    A programme call is uninterpretable without the load it was called at, and
    joining it here keeps the two layers from being compared by hand.
    """
    if summary.empty:
        return
    with open(layer1_path, newline="", encoding="utf-8") as handle:
        totals = {
            r["virus_name"]: r.get(molecule_column) for r in csv.DictReader(handle, delimiter="\t")
        }
    summary["layer1_molecules"] = summary["virus_name"].map(totals)
    summary["layer1_molecule_column"] = molecule_column


def run(ctx: RunContext, done_file: str) -> None:
    """Snakemake entry point for one Run."""
    global config, output
    from viralscan.kb_outputs import KbCountOutputs

    config = ctx.config
    output = config.output
    outputs = KbCountOutputs.from_config_output(config.output)
    main(
        str(outputs.adata_multimap),
        os.path.join(config.output, "results", "viral_summary.tsv"),
        done_file,
    )


if "snakemake" in globals():
    run(
        RunContext.from_yaml(snakemake.params.configfile),  # noqa: F821
        snakemake.output[0],  # noqa: F821
    )
