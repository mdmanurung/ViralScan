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
from typing import Any, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from viralscan.gene_programs import (
    EVIDENCE_LAYER,
    Marker,
    Row,
    VetoSetting,
    call_cell_programme,
    catalogue_virus_facts,
    effective_veto,
    load_catalogue,
    load_veto_settings,
    resolve_markers,
    summarise_programs,
    validate_catalogue,
    write_program_outputs,
)
from viralscan.run_context import RunContext
from viralscan.runconfig import RunConfig
from viralscan.scripts.cellcalling import load_called_mask
from viralscan.utils import setup_script_logging
from viralscan.virus_grouping import (
    detected_summary_rows,
    group_genes_by_identity,
    group_genes_by_virus,
    load_run_identity,
)

log = setup_script_logging()

#: Fallback used when the layer-1 summary has no recognised molecule column.
#: Older summaries used ``total_umi``; v3 uses ``viral_molecules_total_est``.
_MOLECULE_COLUMNS = ("viral_molecules_total_est", "total_umi")

config: RunConfig = RunConfig()
output: str = ""
# Injected by Snakemake into the script namespace; declared (not assigned) so the
# `"snakemake" in globals()` guard below still holds when run outside Snakemake.
snakemake: Any
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
    kept = detected_summary_rows(rows)
    if len(kept) < len(rows):
        log.info("skipping %d candidate-only virus row(s) (detected=False)", len(rows) - len(kept))
    return [r["virus_name"] for r in kept], column


def _marker_matrix(
    adata: ad.AnnData, markers: list[Marker]
) -> tuple[Optional[sparse.csr_matrix], Optional[sparse.csr_matrix]]:
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
    if EVIDENCE_LAYER not in adata.layers:
        raise ValueError(
            f"Gene programmes require the {EVIDENCE_LAYER!r} layer. "
            "Rebuild the multimap output; allocated counts cannot replace unique evidence."
        )
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

    def subset(layer: Any) -> Optional[sparse.csr_matrix]:
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


#: Cells per ``call_cell_programme`` block in :func:`run_one` (memory bound only).
_CALL_ROW_BLOCK = 8192


def run_one(
    adata: ad.AnnData,
    viruses: list[str],
    catalogue: list[Row],
    min_breadth: int,
    form: str,
    identity: Any = None,
    latent_min_breadth: int = 1,
    min_umi: float = 0.0,
    productive_min_molecules: int = 1,
    veto_settings: dict[str, VetoSetting] | None = None,
    mixed_veto_override: str = "catalogue",
    mixed_veto_alpha: float = 0.01,
) -> pd.DataFrame:
    """Build the per-cell call table for the detected ``viruses``.

    ``veto_settings`` (default: the packaged ``gene_program_vetoes.tsv``) and
    ``mixed_veto_override`` choose each virus's ``mixed`` veto (PLAN PROG-MIXED-VETO-01);
    under the default every virus uses ``any``, the earlier rule.
    """
    facts = catalogue_virus_facts(catalogue)
    if veto_settings is None:
        veto_settings = load_veto_settings()
    genes_of_virus: dict[str, list[str]] | None = None
    obs_names = list(adata.obs_names)
    records: list[dict[str, Any]] = []
    marker_resolution: dict[str, dict[str, Any]] = {}

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
                }
            )
            marker_resolution[virus] = {"resolved": 0, "unresolved": 0}
            log.info("%s: no programme model in the catalogue; recorded as not_applicable", virus)
            continue

        resolved, unresolved = resolve_markers(
            virus, adata.var_names, form, catalogue=catalogue, identity=identity
        )
        # An unresolved productive marker makes the absence-based states (reactivating,
        # latent) unassessable on this index (audit D-02).
        # A productive marker whose polyA capture is not established (AUD-17) does the same.
        absence_ok = not any(u.programme == "productive" for u in unresolved) and not info.get(
            "capture_uncertain_productive"
        )
        marker_resolution[virus] = {
            "resolved": len(resolved),
            "unresolved": len(unresolved),
            "absence_assessable": absence_ok,
        }
        if not resolved:
            log.warning(
                "%s: %d catalogue marker(s) but none resolved against this index; "
                "the reference and the catalogue disagree. Skipping.",
                virus,
                len(unresolved),
            )
            continue
        if unresolved:
            log.warning(
                "%s: %d of %d catalogue marker(s) unresolved (%s); calls use the reduced "
                "marker set",
                virus,
                len(unresolved),
                len(resolved) + len(unresolved),
                unresolved[0].reason,
            )

        unique, selected = _marker_matrix(adata, resolved)
        veto = effective_veto(virus, veto_settings, mixed_veto_override)
        if (
            mixed_veto_override == "betabinomial"
            and veto.rule != "betabinomial"
            and info["latency_observable"]  # mixed is unreachable otherwise
        ):
            log.warning("%s: no fitted veto, so the beta-binomial veto is not applied (any)", virus)
        totals = sel_totals = None
        if veto.rule == "betabinomial":
            # N per cell: every gene of this virus, not only the markers (PROG-MIXED-VETO-01).
            if genes_of_virus is None:
                genes_of_virus = (
                    group_genes_by_identity(adata.var_names, identity)[0]
                    if identity is not None
                    else group_genes_by_virus(adata.var_names)[0]
                )
            cols = adata.var_names.get_indexer(genes_of_virus.get(virus, []))
            cols = cols[cols >= 0]
            if len(cols) == 0:
                raise ValueError(f"{virus}: no genes of this virus in the matrix for the veto's N")
            totals = np.asarray(adata.layers[EVIDENCE_LAYER][:, cols].sum(axis=1)).ravel()
            sel_totals = np.asarray(adata.X[:, cols].sum(axis=1)).ravel()
        # Call in bounded row blocks so only one block's call dicts are alive at a
        # time; pruned rows are dropped immediately. State counts cover *all* calls
        # (kept or pruned), accumulated in global first-seen order, and are logged once.
        state_counts: dict[str, int] = {}
        # At least one (possibly empty) block, so parameter validation inside
        # call_cell_programme still fires for a zero-cell input.
        for start in range(0, max(adata.n_obs, 1), _CALL_ROW_BLOCK):
            stop = start + _CALL_ROW_BLOCK
            # A None layer passes through unchanged so call_cell_programme validates it as before.
            calls = call_cell_programme(
                unique if unique is None else unique[start:stop],
                resolved,
                min_breadth=min_breadth,
                latency_observable=info["latency_observable"],
                selected_matrix=None if selected is None else selected[start:stop],
                latent_min_breadth=latent_min_breadth,
                min_umi=min_umi,
                productive_min_molecules=productive_min_molecules,
                absence_assessable=absence_ok,
                mixed_veto=veto,
                mixed_veto_alpha=mixed_veto_alpha,
                virus_totals=None if totals is None else totals[start:stop],
                selected_virus_totals=None
                if sel_totals is None or selected is None
                else sel_totals[start:stop],
            )
            for barcode, call in zip(obs_names[start:stop], calls):
                state_counts[call["state"]] = state_counts.get(call["state"], 0) + 1
                if call["state"] == "indeterminate" and not any(
                    call[k] or 0  # selected_* is None when there is no selected matrix
                    for k in (
                        "productive_breadth",
                        "latent_breadth",
                        "selected_productive_breadth",
                        "selected_latent_breadth",
                    )
                ):
                    # No marker evidence on either layer: the virus was called by layer 1 on
                    # molecules we cannot attribute to a programme, which is worth
                    # recording but not worth a per-cell row for every barcode. A cell with
                    # allocated-only evidence is kept: summarise_programs counts it on the
                    # selected layer (the BZLF1 "unique 0 / allocated 9,308" comparison).
                    continue
                records.append({"barcode": barcode, "virus_name": virus, **call})
        # Rebuild the series from the counts in first-seen order so value_counts() sees
        # the same first-seen order and counts, hence the same tie order, as before.
        states = (
            pd.Series([s for s, n in state_counts.items() for _ in range(n)])
            .value_counts()
            .to_dict()
        )
        log.info("%s: %s", virus, states)

    df = pd.DataFrame(
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
        ],
    )
    df.attrs["marker_resolution"] = marker_resolution
    return df


def main(
    adata_path: str, summary_path: str, done_path: str, results_root: str | None = None
) -> None:
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
    veto_settings = load_veto_settings()
    veto_override = str(getattr(config, "programme_mixed_veto", "catalogue"))
    cells = run_one(
        adata,
        viruses,
        catalogue,
        min_breadth=int(getattr(config, "programme_min_breadth", 2)),
        form=_panel_form(config),
        identity=load_run_identity(output),
        latent_min_breadth=int(getattr(config, "programme_latent_min_breadth", 1)),
        min_umi=float(getattr(config, "programme_min_umi", 0.0)),
        productive_min_molecules=int(getattr(config, "programme_productive_min_molecules", 1)),
        veto_settings=veto_settings,
        mixed_veto_override=veto_override,
        mixed_veto_alpha=float(getattr(config, "programme_mixed_veto_alpha", 0.01)),
    )
    summary = summarise_programs(
        cells,
        catalogue,
        min_breadth=int(getattr(config, "programme_min_breadth", 2)),
        latent_min_breadth=int(getattr(config, "programme_latent_min_breadth", 1)),
        min_umi=float(getattr(config, "programme_min_umi", 0.0)),
        productive_min_molecules=int(getattr(config, "programme_productive_min_molecules", 1)),
        viruses=viruses,
        marker_resolution=cells.attrs.get("marker_resolution"),
    )
    _attach_layer1_totals(summary, summary_path, molecule_column)
    summary["n_called_cells"] = adata.n_obs
    # PROG-INDEXKIND-01: which index and host filter produced the unique-layer evidence.
    summary["index_kind"], summary["host_filter"] = index_state(
        adata.uns, getattr(config, "host_filter_aligner", None)
    )
    # PROG-MIXED-VETO-01: which veto each virus was called with (blank: no programme model).
    modelled = set(catalogue_virus_facts(catalogue))
    vetoes = {
        v: effective_veto(v, veto_settings, veto_override)
        for v in summary["virus_name"]
        if v in modelled
    }
    summary["mixed_veto_rule"] = [
        vetoes[v].rule if v in vetoes else "" for v in summary["virus_name"]
    ]
    summary["mixed_veto_mean"] = [
        vetoes[v].mean if v in vetoes else None for v in summary["virus_name"]
    ]
    summary["mixed_veto_rho"] = [
        vetoes[v].rho if v in vetoes else None for v in summary["virus_name"]
    ]
    summary["mixed_veto_alpha"] = [
        float(getattr(config, "programme_mixed_veto_alpha", 0.01)) if r == "betabinomial" else None
        for r in summary["mixed_veto_rule"]
    ]
    cells["lytic_programme"] = lytic_programme_flags(cells["state"])

    paths = write_program_outputs(cells, summary, results_root or output)
    log.info("Wrote %s and %s", os.path.basename(paths[0]), os.path.basename(paths[1]))
    if done_path:
        with open(done_path, "w") as handle:
            handle.write("done\n")


def index_state(uns: Any, host_filter_aligner: str | None) -> tuple[str, str]:
    """``(index_kind, host_filter)`` of a run, for the summary columns (PROG-INDEXKIND-01).

    ``index_kind`` is what multimap stamped on the H5AD (``combined``, ``virus_only``,
    ``host_only``), or ``unknown`` for an H5AD without the stamp. ``host_filter`` is the
    host-filter aligner, or ``none`` when no host filter ran.
    """
    return str(uns.get("index_kind", "unknown")), host_filter_aligner or "none"


def lytic_programme_flags(states: pd.Series) -> pd.Series:
    """Per-cell "lytic programme present": ``productive`` or ``mixed`` (PROG-MIXED-VETO-01).

    The split between the two is load-driven (a heavily lytic cell picks up one latent
    molecule by depth alone), their union is not. Placeholder rows of a virus with no
    programme model (``not_applicable``) are not cells, so they stay blank.
    """
    flags = states.isin(["productive", "mixed"]).astype("boolean")
    return flags.mask(states == "not_applicable", pd.NA)


def _attach_layer1_totals(summary: pd.DataFrame, layer1_path: str, molecule_column: str) -> None:
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
