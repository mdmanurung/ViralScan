#!/usr/bin/env python3
"""Audit the locked environments against OSV (REL-11). Needs the network, no credentials.

``uv.lock``: every package must come from the PyPI registry with hashes, and each
registry package is queried against OSV at its locked version. ``environment.explicit.txt``:
URL/md5 shape and channel check, then an approximate OSV query by PyPI-matching name
(OSV has no conda ecosystem; the image scan in ``release.yml`` is the authoritative one).
Advisories are reported, not gated: accepting one is a maintainer decision (REL-11).
"""

import argparse
import collections
import json
import re
import sys
import urllib.request
from pathlib import Path

try:
    import tomllib
except ImportError:  # Python < 3.11
    import tomli as tomllib  # type: ignore[no-redef]

OSV_BATCH = "https://api.osv.dev/v1/querybatch"


def osv_hits(pairs):
    """Return {(name, version): [advisory ids]} for the (name, version) pairs OSV flags."""
    hits = {}
    for i in range(0, len(pairs), 500):
        chunk = pairs[i : i + 500]
        body = {"queries": [{"package": {"name": n, "ecosystem": "PyPI"}, "version": v} for n, v in chunk]}
        req = urllib.request.Request(OSV_BATCH, json.dumps(body).encode(), {"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as resp:  # nosec B310 - fixed https URL
            results = json.load(resp)["results"]
        for pair, result in zip(chunk, results):
            if result.get("vulns"):
                hits[pair] = sorted(v["id"] for v in result["vulns"])
    return hits


def audit_uv_lock(path: Path) -> int:
    packages = tomllib.loads(path.read_text())["package"]
    sources = collections.Counter(next(iter(p["source"])) for p in packages)
    missing_hash = [
        p["name"]
        for p in packages
        if "registry" in p["source"]
        and not all("hash" in a for a in [*([p["sdist"]] if "sdist" in p else []), *p.get("wheels", [])])
    ]
    print(f"{path}: {len(packages)} packages, sources {dict(sources)}, without full hashes {missing_hash}")
    registry = [p for p in packages if "registry" in p["source"] and p["source"]["registry"] == "https://pypi.org/simple"]
    hits = osv_hits(sorted({(p["name"], p["version"]) for p in registry}))
    for (name, version), ids in sorted(hits.items()):
        print(f"  {name}=={version}: {len(ids)} advisories ({', '.join(ids[:3])}{'...' if len(ids) > 3 else ''})")
    print(f"  {len(registry)} registry packages queried, {len(hits)} with advisories")
    return int(bool(hits or missing_hash))


def audit_conda_lock(path: Path) -> int:
    pattern = re.compile(r"https://conda\.anaconda\.org/([^/]+)/(?:linux-64|noarch)/(.+)-([^-]+)-([^-]+)\.(?:conda|tar\.bz2)$")
    rows, bad, channels = [], [], collections.Counter()
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#") or line == "@EXPLICIT":
            continue
        url, _, md5 = line.partition("#")
        m = pattern.match(url)
        if m and re.fullmatch(r"[0-9a-f]{32}", md5):
            channels[m.group(1)] += 1
            rows.append((m.group(2), m.group(3)))
        else:
            bad.append(line)
    dup = [n for n, c in collections.Counter(n for n, _ in rows).items() if c > 1]
    print(f"{path}: {len(rows)} packages, channels {dict(channels)}, malformed {bad}, duplicate names {dup}")
    hits = osv_hits(sorted(set(rows)))
    for (name, version), ids in sorted(hits.items()):
        print(f"  {name}=={version}: {len(ids)} advisories (name-matched, approximate)")
    print(f"  {len(hits)} name-matched packages with advisories")
    return int(bool(hits or bad or dup))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uv-lock", type=Path, default=Path("uv.lock"))
    parser.add_argument("--conda-lock", type=Path, default=Path("environment.explicit.txt"))
    args = parser.parse_args()
    sys.exit(audit_uv_lock(args.uv_lock) | audit_conda_lock(args.conda_lock))
