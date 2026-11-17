"""PKG-01: the conda recipe agrees with pyproject.toml where both state a requirement."""

import re
from pathlib import Path

import tomllib

ROOT = Path(__file__).resolve().parents[1]
RECIPE = (ROOT / "conda-recipe" / "meta.yaml").read_text(encoding="utf-8")
PYPROJECT = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))


def test_recipe_entry_point_is_the_pyproject_script() -> None:
    script = PYPROJECT["project"]["scripts"]["viralscan"]
    assert f"- viralscan = {script}" in RECIPE


def test_recipe_python_bound_is_requires_python() -> None:
    requires = PYPROJECT["project"]["requires-python"]
    assert requires == ">=3.10"
    assert f"python {requires}" in RECIPE
    assert "python >=3.9" not in RECIPE


def test_recipe_snakemake_major_matches_the_full_extra() -> None:
    extra = next(
        d for d in PYPROJECT["project"]["optional-dependencies"]["full"] if "snakemake" in d
    )
    assert "snakemake>=9,<10" in extra
    assert re.search(r"snakemake-minimal >=9,<10", RECIPE)


def test_recipe_lists_every_core_dependency() -> None:
    for dependency in PYPROJECT["project"]["dependencies"]:
        name = re.split(r"[<>=\s]", dependency, maxsplit=1)[0].lower()
        assert name in RECIPE.lower() or name == "matplotlib" and "matplotlib-base" in RECIPE, name
