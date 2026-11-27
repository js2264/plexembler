"""`plexembler browser`: chromosomix's store of every species' files."""

import importlib.util
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.browser import gather, store, track_label
from plexembler.cli import cli

#: chromosomix's own check of a store, when its repository is beside this one.
CHECKER = (
    Path(__file__).resolve().parents[2] / "chromosomix" / "scripts" / "make_store.py"
)

FILES = [
    ("Scer", "species/Scer/genome/Scer.fa.gz"),
    ("Scer", "species/Scer/hic/Scer.mcool"),
    ("Scer", "species/Scer/tracks/Scer.rna-ypd-1.forward.bigWig"),
    ("Scer", "species/Scer/annotation/Scer.gff3.gz"),
    ("Cgla", "species/Cgla/genome/Cgla.fa.gz"),
]


def test_labels_and_rows():
    assert track_label("a/Scer.gff3.gz", "Scer") == "genes"
    assert track_label("a/Scer.rna-ypd-1.forward.bigWig", "Scer") == "rna-ypd-1 forward"
    species = gather({"Scer": "poolA", "Cgla": "poolA"}, FILES)
    assert [s.name for s in species] == ["Cgla", "Scer"]
    # The genes first, then the coverage.
    assert species[1].tracks[0].endswith("Scer.gff3.gz")
    document = store(species, "test", base="https://example.org/results")
    assert document["base"] == "https://example.org/results/"
    assert [(r["label"], r["kind"]) for r in document["entries"]] == [
        ("Cgla genome", "genome"),
        ("Scer genome", "genome"),
        ("Scer hi-c map", "map"),
    ]
    assert document["entries"][2]["tracks"][1] == {
        "label": "rna-ypd-1 forward",
        "url": "species/Scer/tracks/Scer.rna-ypd-1.forward.bigWig",
    }
    assert document["layout"] == "cards" and document["group"] == "species"


def test_bad_files():
    with pytest.raises(ValueError, match="not one of the species"):
        gather({"Scer": "poolA"}, [("Spar", "species/Spar/hic/Spar.mcool")])
    with pytest.raises(ValueError, match="neither a map"):
        gather({"Scer": "poolA"}, [("Scer", "species/Scer/Scer.report.html")])


def test_the_command(tmp_path: Path, plain):
    arguments = ["browser", "-s", "Scer=poolA", "-s", "Cgla=poolA", "-o", str(tmp_path)]
    arguments += [
        item for species, path in FILES for item in ("-f", f"{species}={path}")
    ]
    result = CliRunner().invoke(cli, arguments)
    assert result.exit_code == 0, result.output
    document = json.loads((tmp_path / "chromosomix.store.json").read_text())
    assert document["kind"] == "chromosomix.store"
    assert len(document["entries"]) == 3
    assert CliRunner().invoke(cli, arguments).exit_code == 1
    bad = CliRunner().invoke(
        cli, ["browser", "-s", "Scer", "-o", str(tmp_path), "--force"]
    )
    assert bad.exit_code == 2
    assert "SPECIES=POOL" in plain(bad.output)
    unknown = CliRunner().invoke(cli, [*arguments, "-f", "Spar=x.mcool", "--force"])
    assert unknown.exit_code == 1
    url = CliRunner().invoke(cli, [*arguments, "--base", "s3://bucket", "--force"])
    assert url.exit_code == 2


@pytest.mark.skipif(not CHECKER.exists(), reason="chromosomix is not beside plexembler")
def test_chromosomix_reads_it(tmp_path: Path):
    spec = importlib.util.spec_from_file_location("make_store", CHECKER)
    assert spec and spec.loader
    make_store = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(make_store)
    path = tmp_path / "chromosomix.store.json"
    path.write_text(
        json.dumps(store(gather({"Scer": "poolA", "Cgla": "poolA"}, FILES), "test"))
    )
    assert make_store.check_store(path) == []
