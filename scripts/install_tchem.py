#!/usr/bin/env python3
"""Install the pinned tchem release into a prefix (PLAN TCHEM-SOURCE-01).

tchem (https://github.com/Tussl-Co/TusslChemistry, Apache-2.0) has no conda or PyPI package.
The digests below are the release's published ``.sha256`` files, copied on 2026-11-26. The
archive is checked against them and never against a digest fetched at install time.

    python scripts/install_tchem.py --prefix "$CONDA_PREFIX"   # bin/tchem, opt/tchem/
    pixi run install-tchem                                       # the same, into the pixi env

Barcode lists are not installed: fetch them with ``opt/tchem/scripts/fetch_whitelists.sh DIR`` and
point ``--tchem-whitelists`` or ``$VIRALSCAN_TCHEM_WHITELISTS`` at DIR (docs/installation.md).
"""

import argparse
import hashlib
import platform
import tarfile
import tempfile
import urllib.request
from pathlib import Path

VERSION = "0.1.0"
RELEASE = f"https://github.com/Tussl-Co/TusslChemistry/releases/download/v{VERSION}"
#: Rust target triple -> sha256 of ``tchem-v<VERSION>-<triple>.tar.gz``.
SHA256 = {
    "x86_64-unknown-linux-musl": "9112e5018c26868317fa6547490094b59915927ef1482bb361f0b60ce5fe488e",
    "aarch64-unknown-linux-musl": "df31fd49bfdc3c8ec64cdc54b4cf1c0e62e5078cb9a909234600891d7a19159f",
    "x86_64-apple-darwin": "98478f308a8bdbbc3472fdd625ef3cbc833bf1623fde01dac1975ca3a3802a09",
    "aarch64-apple-darwin": "05e3ed8c74a477b9ea9330cd1b69af1b914ece131b9cac49c01f8c7a432ca639",
}


def host_triple() -> str:
    machine = platform.machine().lower().replace("amd64", "x86_64").replace("arm64", "aarch64")
    system = platform.system()
    if system == "Linux":
        return f"{machine}-unknown-linux-musl"
    if system == "Darwin":
        return f"{machine}-apple-darwin"
    raise SystemExit(f"no tchem release asset for {system} {machine}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--prefix", required=True, type=Path, help="install into PREFIX/opt/tchem")
    parser.add_argument("--target", default=None, choices=sorted(SHA256), help="default: this host")
    args = parser.parse_args()
    target = args.target or host_triple()
    if target not in SHA256:
        raise SystemExit(f"no tchem release asset for {target}")
    asset = f"tchem-v{VERSION}-{target}.tar.gz"
    dest = args.prefix / "opt" / "tchem"
    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp) / asset
        urllib.request.urlretrieve(f"{RELEASE}/{asset}", archive)
        got = hashlib.sha256(archive.read_bytes()).hexdigest()
        if got != SHA256[target]:
            raise SystemExit(f"{asset}: sha256 {got} is not the pinned {SHA256[target]}")
        dest.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive) as tar:
            for member in tar.getmembers():
                member.name = member.name.split("/", 1)[1] if "/" in member.name else ""
                if member.name:
                    tar.extract(member, dest, filter="data")
    link = args.prefix / "bin" / "tchem"
    link.parent.mkdir(parents=True, exist_ok=True)
    link.unlink(missing_ok=True)
    link.symlink_to(dest / "bin" / "tchem")
    print(f"tchem {VERSION} installed: {link}")


if __name__ == "__main__":
    main()
