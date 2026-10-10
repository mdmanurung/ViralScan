#!/usr/bin/env python3
"""Emit an Apptainer definition derived from a published OCI digest (REL-08)."""

import argparse
import re


def definition(oci: str) -> str:
    if not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9./_:-]*@sha256:[a-f0-9]{64}", oci):
        raise ValueError(
            "Pass an OCI reference with a full sha256 digest; mutable tags are refused"
        )
    return (
        f"Bootstrap: docker\nFrom: {oci}\n\n"
        "%runscript\n"
        '    exec /opt/conda/envs/viralscan/bin/viralscan "$@"\n\n'
        "%test\n"
        "    /opt/conda/envs/viralscan/bin/viralscan --help\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("oci", help="ghcr.io/...@sha256:<64 hex digits>")
    try:
        print(definition(parser.parse_args().oci), end="")
    except ValueError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
