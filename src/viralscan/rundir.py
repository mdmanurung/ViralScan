"""Where a post-hoc verb finds its run: ``--run-dir <run root>`` plus an optional ``--sample``.

One resolver (:func:`resolve_sample`) turns that pair into the sample directory. It reuses
``discovery.resolve_sample_dir``, the rule ``open_run`` uses: the run must be complete, and a
run with several samples needs ``--sample``. Before 3.0.0b1 four verbs (``evidence``,
``rerun-programs``, ``hostresponse``, ``host-loss``) took the sample directory itself as
``--run-dir``. That form still works, with a deprecation warning that names the replacement
(``utils.deprecation_message``, removal 3.1.0).

The freshness guard is ``run_safety.require_fresh_output`` (with ``refuse_primary_tree``). It
is not repeated here: every verb writes under ``<sample>/posthoc/`` by default, and an
``--output`` inside a run's own tree outside ``posthoc/`` is refused, ``--force`` included.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from viralscan.discovery import DiscoveryError, resolve_sample_dir
from viralscan.run_safety import RUN_MANIFEST
from viralscan.utils import deprecation_message

__all__ = ["SAMPLE_VERBS", "is_sample_dir", "resolve_sample", "use_sample_dir", "warn_sample_dir"]

# Verbs that read one sample: ``args.run_dir`` is that sample directory once resolved.
SAMPLE_VERBS = ("evidence", "rerun-programs", "hostresponse", "host-loss")


def is_sample_dir(path: Path) -> bool:
    """A sample directory (holds ``config.yaml``) that is not a run root (no ``run_manifest.json``)."""
    return (path / "config.yaml").is_file() and not (path / RUN_MANIFEST).is_file()


def resolve_sample(
    run_dir: str | Path, sample: str | None = None, *, sample_dir_ok: bool = False
) -> Path:
    """The sample directory for ``--run-dir`` and ``--sample``; raises ``DiscoveryError``.

    ``run_dir`` is the run root. With ``sample_dir_ok`` it may also be one sample directory
    (the form before 3.0.0b1): that prints a deprecation warning and returns it unchanged.
    """
    path = Path(run_dir).resolve()
    if sample_dir_ok and is_sample_dir(path):
        if sample is not None and sample != path.name:
            raise DiscoveryError(
                f"--sample {sample!r} does not match the sample directory {path}; "
                f"use --run-dir {path.parent} --sample {path.name}."
            )
        warn_sample_dir(path)
        return path
    return resolve_sample_dir(path, sample)


def warn_sample_dir(sample_dir: Path) -> None:
    """Print the deprecation warning for a sample directory passed as ``--run-dir``."""
    print(
        "viralscan: warning: "
        + deprecation_message(
            "--run-dir pointing at a sample directory",
            f"--run-dir {sample_dir.parent} --sample {sample_dir.name}",
        ),
        file=sys.stderr,
    )


def use_sample_dir(args: argparse.Namespace) -> Path:
    """Replace ``args.run_dir`` by the resolved sample directory of a ``SAMPLE_VERBS`` verb.

    ``args._posthoc_sample`` keeps the sample name for the post-hoc record. The user's own
    ``--run-dir`` stays in the record's settings, which are taken before this runs.
    """
    sample_dir = resolve_sample(args.run_dir, getattr(args, "sample", None), sample_dir_ok=True)
    args.run_dir = str(sample_dir)
    args._posthoc_sample = sample_dir.name
    return sample_dir
