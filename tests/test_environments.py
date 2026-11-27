"""The conda environments of the pipeline's steps that run plexembler.

Each installs plexembler's dependencies from conda-forge, and plexembler alone
with pip, at the release's tag: pip would take a dependency that ~/.local
holds for installed, which the steps, run with PYTHONNOUSERSITE, do not see.
"""

import re
from importlib.metadata import requires, version
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
ENVIRONMENTS = sorted(
    (ROOT / "modules" / "local" / "plexembler").glob("*/environment.yml")
)
#: PyPI's dependencies that the steps do not need: Nextflow's launcher.
NOT_IN_STEPS = {"nextflow"}


def _numbers(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.split("."))


def _floors() -> dict[str, str]:
    """plexembler's dependencies, and their floors."""
    floors = {}
    for requirement in requires("plexembler") or []:
        found = re.fullmatch(r"([A-Za-z0-9._-]+)>=([0-9.]+)", requirement)
        if found and found[1] not in NOT_IN_STEPS:
            floors[found[1].lower()] = found[2]
    return floors


@pytest.mark.skipif(not ENVIRONMENTS, reason="the pipeline is not in the sdist")
@pytest.mark.parametrize("path", ENVIRONMENTS, ids=lambda path: path.parent.name)
def test_plexemblers_dependencies_come_from_conda(path: Path):
    dependencies = yaml.safe_load(path.read_text(encoding="utf-8"))["dependencies"]
    conda = {}
    for dependency in dependencies:
        if isinstance(dependency, str):
            name, _, pinned = dependency.split("::")[-1].partition("=")
            conda[name] = pinned
    floors = _floors()
    assert floors, "plexembler has no dependency"
    for name, floor in floors.items():
        assert name in conda, f"{path}: {name} from conda-forge, not pip"
        assert _numbers(conda[name]) >= _numbers(floor), f"{path}: {name}>={floor}"
    (pip,) = [
        dependency["pip"] for dependency in dependencies if isinstance(dependency, dict)
    ]
    assert pip == [
        f"plexembler @ git+https://github.com/js2264/plexembler@v{version('plexembler')}"
    ]
