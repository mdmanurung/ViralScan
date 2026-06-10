"""The passive context handed to each Snakemake worker script's ``run(ctx)``.

A Run Context bundles the two Run-level facts every worker needs — the config
and the kb-python output layout — into one value that is *passed in* rather than
read from Snakemake's magic globals. Moving the seam from the global to a
function parameter is what lets the worker logic be imported and tested directly
(no ``runpy`` / fake-``snakemake`` dance, no mirror re-implementations).

It is deliberately passive: it answers "what / where", it performs no I/O. See
``CONTEXT.md`` ("Run Context").

Note: ``config`` is kept as a plain ``dict`` here (what ``load_config`` returns).
The *typed* :class:`viralscan.runconfig.RunConfig` is the write-side checkpoint
(``createconfig``); the read side still consumes the dict, which is what the
worker helpers expect via ``config.get(...)`` / ``config[...]``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Union

from viralscan.kb_outputs import KbCountOutputs
from viralscan.utils import load_config


@dataclass(frozen=True)
class RunContext:
    config: dict[str, Any]
    outputs: KbCountOutputs

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> "RunContext":
        """Build from an already-loaded config dict (the testable seam)."""
        return cls(config, KbCountOutputs.from_config_output(config["output"]))

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "RunContext":
        """Build from a ``config.yaml`` on disk (the production seam)."""
        return cls.from_config(load_config(path))
