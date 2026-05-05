"""Small shared helpers used across the Snakemake worker scripts."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Union

import yaml


def load_config(path: Union[str, Path]) -> dict[str, Any]:
    """Read a YAML config file and return it as a plain ``dict``.

    Centralised so every script can use the same loader (and so we have
    one place to evolve the boolean-normalisation work tracked in
    PLAN §1.6).
    """
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)
