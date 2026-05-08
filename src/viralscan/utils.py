"""Small shared helpers used across the Snakemake worker scripts."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Union, cast

import yaml

# Module-level logger for scripts that import this module.
logger = logging.getLogger("viralscan")


def configure_logging(verbose: bool = False, quiet: bool = False) -> None:
    """Configure the root *viralscan* logger.

    Called once from ``menu.py`` after parsing CLI args.  Snakemake worker
    scripts run in their own processes so they call this themselves from
    the ``setup_script_logging()`` helper below.

    Priority: quiet > verbose > default (INFO).
    """
    level = logging.WARNING if quiet else (logging.DEBUG if verbose else logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(levelname)s  %(message)s"))
    root = logging.getLogger("viralscan")
    root.setLevel(level)
    root.handlers.clear()
    root.addHandler(handler)


def setup_script_logging() -> logging.Logger:
    """Minimal logging setup for Snakemake worker scripts.

    Each script runs in its own interpreter so it must configure its own
    handler.  Returns the ``viralscan`` logger ready to use.
    """
    configure_logging()
    return logger


def load_config(path: Union[str, Path]) -> dict[str, Any]:
    """Read a YAML config file and return it as a plain ``dict``.

    Centralised so every script can use the same loader (and so we have
    one place to evolve the boolean-normalisation work tracked in
    PLAN §1.6).
    """
    with open(path, "r", encoding="utf-8") as f:
        loaded = cast(Any, yaml.safe_load(f))
    if not isinstance(loaded, dict):
        raise ValueError(f"Config file {path} did not contain a YAML mapping.")
    return cast(dict[str, Any], loaded)
