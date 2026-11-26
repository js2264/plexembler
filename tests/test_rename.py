"""`plexembler rename`: contigs named after their species, and the run's log."""

import gzip
import logging
import re
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler._log import elapsed, log_file_name
from plexembler.cli import cli
from plexembler.rename import check_species, new_name, read_fasta

FASTA = ">ctg1 length=12\nACGTACGTAC\nGT\n>sp1_ctg2\nTTTT\n\n>ctg3\n" + "A" * 85 + "\n"


@pytest.fixture
def fasta(tmp_path: Path) -> Path:
    path = tmp_path / "contigs.fa"
    path.write_text(FASTA, encoding="utf-8")
    return path


def run(*arguments: str):
    return CliRunner().invoke(cli, ["rename", *arguments])


#: A terminal's colour codes, which rich-click writes where it sees a terminal (CI).
ANSI = re.compile(r"\x1b\[[0-9;]*m")


def plain(output: str) -> str:
    """The output without rich-click's panel: its colours, borders and wrapping."""
    text = ANSI.sub("", output)
    return " ".join(text.translate(str.maketrans("", "", "│╭╮╰╯─")).split())


def test_every_contig_starts_with_its_species(fasta: Path, tmp_path: Path):
    output = tmp_path / "out" / "sp1.fa"
    result = run(str(fasta), "-p", "sp1", "-o", str(output))
    assert result.exit_code == 0, result.output
    records = dict(read_fasta(output.open(encoding="utf-8")))
    # A name that already starts with the species stays as it is.
    assert records == {
        "sp1_ctg1": "ACGTACGTACGT",
        "sp1_ctg2": "TTTT",
        "sp1_ctg3": "A" * 85,
    }
    # Sequences are written 80 bases a line.
    assert max(len(line) for line in output.read_text().splitlines()) == 80


def test_the_names_table(fasta: Path, tmp_path: Path):
    output = tmp_path / "sp1.fa.gz"
    assert run(str(fasta), "-p", "sp1", "-o", str(output)).exit_code == 0
    table = (tmp_path / "sp1.names.tsv").read_text().splitlines()
    assert table == [
        "name\toriginal_name\tlength_bp",
        "sp1_ctg1\tctg1\t12",
        "sp1_ctg2\tsp1_ctg2\t4",
        "sp1_ctg3\tctg3\t85",
    ]


def test_gzipped_in_and_out(fasta: Path, tmp_path: Path):
    zipped = tmp_path / "contigs.fa.gz"
    zipped.write_bytes(gzip.compress(fasta.read_bytes()))
    output = tmp_path / "sp1.fasta.gz"
    names = tmp_path / "names.tsv"
    result = run(str(zipped), "-p", "sp1", "-o", str(output), "-n", str(names))
    assert result.exit_code == 0, result.output
    with gzip.open(output, "rt") as handle:
        assert [name for name, _ in read_fasta(handle)] == [
            "sp1_ctg1",
            "sp1_ctg2",
            "sp1_ctg3",
        ]
    assert names.exists()


def test_outputs_are_not_overwritten_without_force(fasta: Path, tmp_path: Path):
    output = tmp_path / "sp1.fa"
    assert run(str(fasta), "-p", "sp1", "-o", str(output)).exit_code == 0
    again = run(str(fasta), "-p", "sp1", "-o", str(output))
    assert again.exit_code == 1
    assert "--force" in plain(again.output)
    assert run(str(fasta), "-p", "sp1", "-o", str(output), "--force").exit_code == 0


@pytest.mark.parametrize("species", ["1sp", "sp_1", "sp 1", ""])
def test_a_species_id_is_checked(fasta: Path, tmp_path: Path, species: str):
    result = run(str(fasta), "-p", species, "-o", str(tmp_path / "x.fa"))
    assert result.exit_code == 2
    with pytest.raises(ValueError, match="not a species ID"):
        check_species(species)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("ACGT\n>ctg1\nACGT\n", "before the first header"),
        (">\nACGT\n", "without a name"),
        (">ctg1\nAC\n>sp1_ctg1\nGT\n", "two contigs would be named"),
    ],
)
def test_malformed_contigs_fail(tmp_path: Path, text: str, message: str):
    path = tmp_path / "bad.fa"
    path.write_text(text, encoding="utf-8")
    result = run(str(path), "-p", "sp1", "-o", str(tmp_path / "out.fa"))
    assert result.exit_code == 1
    assert message in plain(result.output)


def test_the_log_beside_the_outputs(fasta: Path, tmp_path: Path):
    output = tmp_path / "out" / "sp1.fa"
    assert run(str(fasta), "-p", "sp1", "-o", str(output), "-v").exit_code == 0
    (log,) = (tmp_path / "out").glob("sp1.plexembler_*.log")
    assert re.fullmatch(r"sp1\.plexembler_\d{14}\.log", log.name)
    lines = log.read_text().splitlines()
    assert lines[0].startswith("## plexembler: v")
    assert lines[1].startswith("## date: ")
    assert lines[2].startswith("## cmd: ")
    assert "## prefix: sp1" in lines
    assert "---" in lines
    assert any(":: DEBUG :: ctg1 -> sp1_ctg1" in line for line in lines)
    assert re.search(r"Done after \dh \dm \ds$", lines[-1])
    # The handlers leave with the run.
    assert not logging.getLogger("plexembler").handlers


def test_quiet_keeps_the_log(fasta: Path, tmp_path: Path):
    output = tmp_path / "sp1.fa"
    assert run(str(fasta), "-p", "sp1", "-o", str(output), "-q").exit_code == 0
    (log,) = tmp_path.glob("sp1.plexembler_*.log")
    assert "Renamed 3 contigs of sp1" in log.read_text()


def test_log_names_and_durations(tmp_path: Path):
    from datetime import datetime

    started = datetime(2026, 11, 26, 10, 15, 30)
    assert (
        log_file_name(tmp_path, None, started).name == "plexembler_20261126101530.log"
    )
    assert elapsed(1051) == "0h 17m 31s"
    assert new_name("sp1", "sp1_a") == "sp1_a"
