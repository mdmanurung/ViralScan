"""CONTRACT-VERSIONS-01: the version constants live in viralscan.contract and the schemas agree."""

import json
import re
from importlib import resources
from pathlib import Path

from viralscan import contract, outputs

SRC = Path(contract.__file__).parent


def _schema(name: str) -> dict:
    return json.loads(resources.files("viralscan.schemas.v3").joinpath(name).read_text("utf-8"))


def test_no_other_module_hard_codes_a_schema_version() -> None:
    pattern = re.compile(
        r"""(schema_version|count_schema_version)["']?\]?\s*(=|!=|:)\s*["']3\.0\.0["']"""
    )
    offenders = [
        str(path.relative_to(SRC))
        for path in SRC.rglob("*.py")
        if path.name != "contract.py" and pattern.search(path.read_text("utf-8"))
    ]
    assert not offenders, f"hard-coded schema version in {offenders}"


def test_outputs_re_exports_the_layout_constants() -> None:
    assert outputs.LAYOUT_VERSION == contract.LAYOUT_VERSION
    assert outputs.SUPPORTED_LAYOUT_VERSIONS == contract.SUPPORTED_LAYOUT_VERSIONS
    assert contract.LAYOUT_VERSION in contract.SUPPORTED_LAYOUT_VERSIONS


def test_schemas_pin_the_contract_version() -> None:
    assert _schema("run_manifest.schema.json")["properties"]["schema_version"]["const"] == (
        contract.SCHEMA_VERSION
    )
    assert _schema("run_complete.schema.json")["properties"]["schema_version"]["const"] == (
        contract.SCHEMA_VERSION
    )
    layouts = _schema("run_manifest.schema.json")["properties"]["layout_version"]["enum"]
    assert tuple(layouts) == contract.SUPPORTED_LAYOUT_VERSIONS


def test_run_manifest_schema_requires_the_fingerprint_and_options() -> None:
    required = set(_schema("run_manifest.schema.json")["required"])
    assert {"run_fingerprint", "options", "schema_version", "viralscan_version"} <= required
