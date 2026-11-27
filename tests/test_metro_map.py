"""The README's metro map is its generator's.

docs/assets/metro_map.py writes docs/assets/metro_map_light.svg and
metro_map_dark.svg: a station added to it, or a line moved, without the two
maps written again (``python docs/assets/metro_map.py``) fails here.
"""

import importlib.util
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from types import ModuleType

import pytest

ASSETS = Path(__file__).resolve().parents[1] / "docs" / "assets"
THEMES = ["light", "dark"]


def generator() -> ModuleType:
    """docs/assets/metro_map.py, which is not a module of the package."""
    name = "plexembler_metro_map"
    spec = importlib.util.spec_from_file_location(name, ASSETS / "metro_map.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    # Its dataclasses look their module up while it runs.
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("theme", THEMES)
def test_each_map_is_the_generators(theme: str):
    written = (ASSETS / f"metro_map_{theme}.svg").read_text(encoding="utf-8")
    assert written == generator().svg(theme), "python docs/assets/metro_map.py"


@pytest.mark.parametrize("theme", THEMES)
def test_each_map_is_an_svg_with_its_title(theme: str):
    root = ET.parse(ASSETS / f"metro_map_{theme}.svg").getroot()
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    assert (
        root.findtext("{http://www.w3.org/2000/svg}title") == "plexembler's metro map"
    )


def test_the_generator_writes_both_maps(tmp_path: Path):
    generator().main(tmp_path)
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "metro_map_dark.svg",
        "metro_map_light.svg",
    ]
