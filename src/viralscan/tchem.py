"""TusslChemistry (``tchem``) as a second opinion on the chemistry (PLAN CHEM-02).

``tchem`` (https://github.com/Tussl-Co/TusslChemistry, Apache-2.0, a Rust binary)
names the single-cell kit from the first reads of a run. ViralScan's own detector
(:mod:`viralscan.chemistry`) knows 10xv2, 10xv3, 10xv4 and Drop-seq; ``tchem``
knows about 25 layouts, among them Multiome GEX, 5' kits and the non-10x kits.
It is run beside the built-in detector and never replaces it: the two are
compared on the barcode/UMI geometry (the part kb's ``-x`` encodes), and a
disagreement stops the run exactly like a mismatch with ``-x`` does.

``tchem`` is optional. It is looked up by ``--tchem-bin``, ``$VIRALSCAN_TCHEM``
and ``PATH``; without it the check is skipped (``--tchem auto``) or stops the run
(``--tchem require``). Its barcode lists are fetched with the package's
``scripts/fetch_whitelists.sh`` and found through ``--tchem-whitelists`` or
``$VIRALSCAN_TCHEM_WHITELISTS``.

Two calls per sample, because the JSON report (``--format json``, schema 1.x)
carries the verdict and the confidence while the kb ``-x`` string comes from
``--emit kb``, which writes the technology as a ``bc:umi:seq`` geometry such as
``0,0,16:0,16,26:1,0,0``. :func:`viralscan.chemistry.cb_umi_geometry` already
reads that form. The report keys used are the documented ones (``verdict``,
``confidence``, ``layout.name``, ``layout.validation.level``, ``cdna_strand``,
``candidates``, ``hint``); unknown keys are ignored and a major schema version
other than 1 is refused, as the tool's own docs ask of consumers.
"""

from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
from dataclasses import asdict, dataclass
from typing import Any, Optional

from viralscan import chemistry

ENV_BIN = "VIRALSCAN_TCHEM"
ENV_LISTS = "VIRALSCAN_TCHEM_WHITELISTS"
ENV_TCHEM_LISTS = "TUSSL_CHEMISTRY_WHITELISTS"  # the variable tchem itself reads
SCHEMA_MAJOR = 1
#: ``tchem`` exits 0 whenever it read the reads, whatever the verdict.
EXIT_OK = 0
MODES = ("auto", "require", "off")
TIMEOUT_S = 900


class TchemError(RuntimeError):
    """``tchem`` is missing, failed, or wrote something this code cannot read."""


@dataclass
class TchemCall:
    verdict: str  # "match", "ambiguous" or "no_match"
    confidence: Optional[str]  # "high", "medium", "low" for a match
    layout: Optional[str]
    validation: Optional[str]  # tchem's own evidence level for the layout
    technology: Optional[str]  # kb -x from --emit kb, a bc:umi:seq string
    whitelist: Optional[str]  # kb -w from --emit kb
    strand: Optional[str]  # forward / reverse
    reason: str
    version: Optional[str] = None

    def as_block(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class Opinion:
    """How the ``tchem`` call compares with the technology the run will use."""

    status: str  # "agree", "disagree", "no_call", "unavailable"
    reason: str
    call: Optional[TchemCall] = None

    def as_block(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "reason": self.reason,
            "call": self.call.as_block() if self.call else None,
        }


def find_binary(explicit: Optional[str] = None) -> Optional[str]:
    """The ``tchem`` executable: ``--tchem-bin``, then ``$VIRALSCAN_TCHEM``, then ``PATH``."""
    for cand in (explicit, os.environ.get(ENV_BIN)):
        if cand:
            return cand if os.access(cand, os.X_OK) else None
    return shutil.which("tchem")


def find_whitelists(explicit: Optional[str] = None) -> Optional[str]:
    """``--tchem-whitelists``, then ``$VIRALSCAN_TCHEM_WHITELISTS``, then tchem's own variable."""
    return explicit or os.environ.get(ENV_LISTS) or os.environ.get(ENV_TCHEM_LISTS) or None


def _major(version: Any) -> int:
    try:
        return int(str(version).split(".")[0])
    except ValueError as exc:
        raise TchemError(f"tchem schema_version {version!r} is not 'MAJOR.MINOR'") from exc


def parse_emit(stdout: str) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """``(-x, -w, --strand)`` read from the command ``tchem --emit kb`` prints.

    The command may span lines with ``\\`` continuations and carry ``$VAR``
    placeholders for the index; only the three options ViralScan cares about
    are taken from it.
    """
    text = stdout.replace("\\\n", " ")
    try:
        tokens = shlex.split(text, comments=True)
    except ValueError as exc:
        raise TchemError(f"cannot read the kb command tchem printed: {exc}") from exc

    def after(flag: str) -> Optional[str]:
        for i, tok in enumerate(tokens[:-1]):
            if tok == flag:
                return tokens[i + 1]
            if tok.startswith(flag + "="):
                return tok.split("=", 1)[1]
        return None

    return after("-x"), after("-w"), after("--strand")


def parse_report(report: dict[str, Any], emit_stdout: Optional[str] = None) -> TchemCall:
    """Turn the JSON report (and the ``--emit kb`` text) into a :class:`TchemCall`."""
    if report.get("tool") != "tussl-chemistry":
        raise TchemError("the report is not from tussl-chemistry")
    if _major(report.get("schema_version")) != SCHEMA_MAJOR:
        raise TchemError(
            f"tchem report schema {report.get('schema_version')} is not 1.x; update ViralScan "
            "or pin an older tchem"
        )
    verdict = str(report.get("verdict"))
    if verdict not in ("match", "ambiguous", "no_match"):
        raise TchemError(f"unknown tchem verdict {verdict!r}")
    layout = report.get("layout") or {}
    technology = whitelist = emit_strand = None
    if emit_stdout:
        technology, whitelist, emit_strand = parse_emit(emit_stdout)
    if verdict == "match":
        reason = f"{layout.get('name')} ({report.get('confidence')} confidence)"
    elif verdict == "ambiguous":
        reason = "ambiguous between " + ", ".join(report.get("candidates") or ["?"])
    else:
        reason = str(report.get("hint") or "no layout fits these reads")
    return TchemCall(
        verdict=verdict,
        confidence=report.get("confidence"),
        layout=layout.get("name"),
        validation=(layout.get("validation") or {}).get("level"),
        technology=technology,
        whitelist=whitelist,
        strand=report.get("cdna_strand") or emit_strand,
        reason=reason,
        version=report.get("version"),
    )


def _run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=TIMEOUT_S)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise TchemError(f"could not run {cmd[0]}: {exc}") from exc


def call(r1: str, r2: Optional[str], *, binary: str, whitelists: Optional[str]) -> TchemCall:
    """Run ``tchem`` on one sample's reads and parse what it says."""
    base = [binary]
    if whitelists:
        base += ["--whitelists", whitelists]
    base += ["--r1", r1]
    if r2:
        base += ["--r2", r2]
    proc = _run(base + ["--format", "json"])
    if proc.returncode != EXIT_OK:
        raise TchemError(f"tchem exited {proc.returncode}: {proc.stderr.strip()[:300]}")
    try:
        report = json.loads(proc.stdout)
    except ValueError as exc:
        raise TchemError(f"tchem did not write JSON: {exc}") from exc
    emit = None
    if report.get("verdict") == "match":
        proc = _run(base + ["--emit", "kb"])
        # exit 3 = a match with no kb command (e.g. no single barcode + UMI): not an error here
        if proc.returncode == EXIT_OK:
            emit = proc.stdout
        elif proc.returncode != 3:
            raise TchemError(
                f"tchem --emit kb exited {proc.returncode}: {proc.stderr.strip()[:300]}"
            )
    return parse_report(report, emit)


def compare(technology: Optional[str], tc: TchemCall) -> Opinion:
    """Compare ``tchem``'s call with the technology the run will use (``-x``).

    With no technology (``check-chemistry`` on a library neither ``-x`` nor the
    built-in detector can name) the call is reported as ``info`` instead.
    """
    if tc.verdict != "match":
        return Opinion("no_call", f"tchem: {tc.verdict} ({tc.reason})", tc)
    if not tc.technology:
        return Opinion(
            "no_call",
            f"tchem names {tc.layout} but gives no kb technology for it, so the geometry "
            "cannot be compared",
            tc,
        )
    if not technology:
        return Opinion(
            "info",
            f"tchem: {tc.layout} ({tc.confidence}); kb -x {tc.technology}"
            + (f" -w {tc.whitelist}" if tc.whitelist else "")
            + (f" --strand {tc.strand}" if tc.strand else ""),
            tc,
        )
    try:
        theirs = chemistry.cb_umi_geometry(tc.technology)
        ours = chemistry.cb_umi_geometry(technology)
    except ValueError as exc:
        return Opinion("no_call", f"geometry not comparable: {exc}", tc)
    if theirs == ours:
        return Opinion(
            "agree", f"tchem: {tc.layout} ({tc.confidence}); same barcode+UMI geometry", tc
        )
    return Opinion(
        "disagree",
        f"-x {technology} is {ours[0]} bp barcode + {ours[1]} bp UMI but tchem reads "
        f"{tc.layout} ({tc.confidence}) as {theirs[0]} + {theirs[1]} (kb -x {tc.technology})",
        tc,
    )


def second_opinion(
    technology: Optional[str],
    r1: str,
    r2: Optional[str],
    *,
    mode: str = "auto",
    binary: Optional[str] = None,
    whitelists: Optional[str] = None,
) -> Opinion:
    """One sample's :class:`Opinion`; raises :class:`TchemError` only for ``mode == "require"``."""
    if mode not in MODES:
        raise ValueError(f"--tchem must be one of {MODES}")
    if mode == "off":
        return Opinion("unavailable", "tchem check is off")
    exe = find_binary(binary)
    if exe is None:
        if mode == "require":
            raise TchemError(
                f"--tchem require but no tchem binary was found (--tchem-bin, ${ENV_BIN}, PATH)"
            )
        return Opinion("unavailable", "tchem is not installed")
    try:
        return compare(technology, call(r1, r2, binary=exe, whitelists=find_whitelists(whitelists)))
    except TchemError as exc:
        if mode == "require":
            raise
        return Opinion("unavailable", f"tchem could not be used: {exc}")


def enforce(opinions: list[Opinion], force: bool = False) -> None:
    """Stop the run on a disagreement or a ``no_call`` unless the technology is forced.

    ``unavailable`` (not installed, off, or failed in ``auto`` mode) never stops it.
    """
    bad = [o.reason for o in opinions if o.status in ("disagree", "no_call")]
    if bad and not force:
        raise chemistry.ChemistryError(
            "tchem second opinion failed: " + "; ".join(bad) + ". Fix -x/-w, or pass -x with "
            "--force-technology to run anyway, or --tchem off to skip this check."
        )
