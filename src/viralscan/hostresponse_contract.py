"""Host-response evaluation rules shared by every entry point.

The integrated pipeline (``RunConfig``), ``viralscan hostresponse`` and
``scripts/hostresponse.run_hostresponse`` must reject the same invalid
combinations, and the pipeline must do so before ``kb count`` rather than inside
the ``hostresponse`` rule hours later. This module is import-light on purpose
(no scanpy/sklearn) so ``runconfig`` can use it.
"""

import argparse

CV_MODES = ("cell", "group")
#: Declared exchangeability units of the structured null (never inferred).
NULL_UNITS = ("cell_within_block", "group")


def validate_null_spec(unit, has_block: bool, n_permutations: int) -> None:
    """Fail fast on an undeclared or inconsistent null specification.

    ``has_block`` says whether a permutation block (column/values) was supplied.
    ``n_permutations == 0`` (the default, null disabled) accepts anything.
    """
    if int(n_permutations) < 0:
        raise ValueError("permutations must be >= 0")
    if int(n_permutations) == 0:
        return
    if unit not in NULL_UNITS:
        raise ValueError(
            f"permutation_unit must be declared explicitly, one of {NULL_UNITS} (got {unit!r}); "
            "exchangeability is never inferred from the group column"
        )
    if unit == "cell_within_block" and not has_block:
        raise ValueError("permutation_unit='cell_within_block' requires a permutation block")
    if unit == "group" and has_block:
        raise ValueError("a block column only applies to permutation_unit='cell_within_block'")


def validate_evaluation(
    *,
    cv_mode,
    groups_column,
    permutations: int,
    permutation_unit,
    permutation_block,
    cell_type_column=None,
    cell_types=None,
) -> None:
    """Reject invalid CV-mode / null / strata combinations (``ValueError`` with the reason)."""
    if cv_mode not in CV_MODES:
        raise ValueError(f"cv_mode must be one of {CV_MODES}, got {cv_mode!r}")
    if cv_mode == "group" and not groups_column:
        raise ValueError("cv_mode='group' requires groups_column (no silent fallback to cell CV)")
    validate_null_spec(permutation_unit, bool(permutation_block), permutations)
    if permutations > 0 and cv_mode != "group":
        raise ValueError("permutations > 0 requires cv_mode='group' (the null is for group CV)")
    if cell_types is not None and (not cell_types or not cell_type_column):
        raise ValueError("cell_types requires a nonempty list and cell_type_column")


def parse_cell_types(tokens):
    """Cell-type names from CLI tokens: space-separated and/or comma-separated.

    ``--cell-types T_cells,B_cells`` and ``--cell-types T_cells B_cells`` are the
    same request on every entry point; names are stripped. A name that itself
    contains a comma cannot be given on the command line (use the config).
    """
    if tokens is None:
        return None
    if isinstance(tokens, str):
        tokens = [tokens]
    names = [name.strip() for token in tokens for name in token.split(",") if name.strip()]
    if not names:
        raise ValueError("at least one cell type name is required")
    return names


class CellTypesAction(argparse.Action):
    """argparse action storing :func:`parse_cell_types` of the given tokens."""

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            setattr(namespace, self.dest, parse_cell_types(values))
        except ValueError as exc:
            parser.error(f"argument {option_string}: {exc}")
