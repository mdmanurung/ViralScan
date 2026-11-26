"""The base class of every exception ViralScan raises on purpose (PLAN ``ERRORS-01``).

Each package exception inherits ``ViralScanError`` and keeps its older base (``ValueError`` or
``RuntimeError``), so an existing ``except ValueError`` still catches it. ``menu._dispatch`` has
one handler for the base: it prints one ``ERROR`` line and exits 1 (see "Exit codes" in
``docs/cli_reference.md``). A usage error that must exit 2 is caught where it is raised.
"""

__all__ = ["ViralScanError"]


class ViralScanError(Exception):
    """A failure ViralScan reports to the user as one ``ERROR`` line, without a traceback."""
