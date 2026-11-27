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
        "name\toriginal_name\tlength_bp\toriginal_length_bp\tcopies\tcircular",
        "sp1_ctg1\tctg1\t12\t12\t1\tNA",
        "sp1_ctg2\tsp1_ctg2\t4\t4\t1\tNA",
        "sp1_ctg3\tctg3\t85\t85\t1\tNA",
    ]


def write_assembly(tmp_path: Path, contigs: dict[str, str], info: dict | None):
    """An assembly, and Flye's table of it: each contig's ``circ.`` and ``mult.``."""
    fasta = tmp_path / "assembly.fasta"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in contigs.items()))
    if info is None:
        return fasta, None
    table = tmp_path / "assembly_info.txt"
    table.write_text(
        "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
        + "".join(
            f"{name}\t{len(contigs[name])}\t40\t{circ}\tN\t{mult}\t*\t{i}\n"
            for i, (name, (circ, mult)) in enumerate(info.items(), start=1)
        )
    )
    return fasta, table


def names_table(path: Path) -> dict[str, list[str]]:
    rows = [line.split("\t") for line in path.read_text().splitlines()]
    return {row[0]: row[1:] for row in rows[1:]}


def test_a_doubled_circle_comes_out_as_one_copy(tmp_path: Path, sequence):
    """Flye's doubled mitochondrion: two copies of a 20 kb circle in one
    contig, cut to one. A single circle, and a linear contig, stay whole."""
    unit = sequence(20_000, seed=1, gc=0.17)
    contigs = {
        "contig_1": unit + unit,
        "contig_2": sequence(20_000, seed=2),
        "contig_3": sequence(30_000, seed=3),
    }
    fasta, info = write_assembly(
        tmp_path, contigs, {"contig_1": ("Y", 2), "contig_2": ("Y", 1), "contig_3": ("N", 1)}
    )  # fmt: skip
    output = tmp_path / "out" / "Cgla.fa"
    result = run(str(fasta), "-p", "Cgla", "-o", str(output), "-a", str(info))
    assert result.exit_code == 0, result.output
    records = dict(read_fasta(output.open(encoding="utf-8")))
    assert records == {
        "Cgla_contig_1": unit,
        "Cgla_contig_2": contigs["contig_2"],
        "Cgla_contig_3": contigs["contig_3"],
    }
    assert names_table(tmp_path / "out" / "Cgla.names.tsv") == {
        "Cgla_contig_1": ["contig_1", "20000", "40000", "2", "yes"],
        "Cgla_contig_2": ["contig_2", "20000", "20000", "1", "yes"],
        "Cgla_contig_3": ["contig_3", "30000", "30000", "1", "no"],
    }
    log = next((tmp_path / "out").glob("Cgla.plexembler_*.log")).read_text()
    assert "## assembly_info: " in log
    assert (
        "Cgla_contig_1: 2 copies of a 20000-bp unit, 100.00 % alike (circular in "
        "the assembly info, multiplicity 2); one kept, 20000 of its 40000 bp" in log
    )
    assert "1 of them cut to one copy of their unit: Cgla_contig_1" in log


def test_a_linear_contig_is_never_cut(tmp_path: Path, sequence):
    """Copies end to end in a contig its assembler says is linear: a tandem
    repeat of a chromosome, such as the rDNA's, which stays whole."""
    unit = sequence(9_100, seed=4)
    fasta, info = write_assembly(
        tmp_path, {"contig_1": unit * 2}, {"contig_1": ("N", 3)}
    )
    output = tmp_path / "sp1.fa"
    assert (
        run(str(fasta), "-p", "sp1", "-o", str(output), "-a", str(info)).exit_code == 0
    )
    assert dict(read_fasta(output.open(encoding="utf-8"))) == {"sp1_contig_1": unit * 2}
    assert names_table(tmp_path / "sp1.names.tsv")["sp1_contig_1"][2:] == [
        "18200",
        "1",
        "no",
    ]


def test_without_assembly_info_the_sequence_decides(tmp_path: Path, sequence):
    """hifiasm gives no table: copies end to end are cut all the same, three
    of them here, and a contig without copies stays whole, its circularity
    unknown."""
    unit = sequence(6_300, seed=5)
    other = sequence(25_000, seed=6)
    fasta, _ = write_assembly(tmp_path, {"ptg1c": unit * 3, "ptg2l": other}, None)
    output = tmp_path / "sp1.fa"
    result = run(str(fasta), "-p", "sp1", "-o", str(output), "-v")
    assert result.exit_code == 0, result.output
    assert dict(read_fasta(output.open(encoding="utf-8"))) == {
        "sp1_ptg1c": unit,
        "sp1_ptg2l": other,
    }
    table = names_table(tmp_path / "sp1.names.tsv")
    assert table["sp1_ptg1c"][2:] == ["18900", "3", "yes"]
    assert table["sp1_ptg2l"][2:] == ["25000", "1", "NA"]
    log = next(tmp_path.glob("sp1.plexembler_*.log")).read_text()
    assert "(by its sequence alone)" in log


def test_a_circle_of_one_copy_with_a_high_multiplicity(tmp_path: Path, sequence):
    """A plasmid deeper than the genome, written once: whole, and the debug log
    says it was looked at."""
    plasmid = sequence(6_300, seed=7)
    fasta, info = write_assembly(
        tmp_path, {"contig_1": plasmid}, {"contig_1": ("Y", 7)}
    )
    output = tmp_path / "sp1.fa"
    result = run(str(fasta), "-p", "sp1", "-o", str(output), "-a", str(info), "-v")
    assert result.exit_code == 0, result.output
    assert dict(read_fasta(output.open(encoding="utf-8"))) == {"sp1_contig_1": plasmid}
    log = next(tmp_path.glob("sp1.plexembler_*.log")).read_text()
    assert (
        "sp1_contig_1: circular, multiplicity 7 in the assembly info; one copy" in log
    )


def test_no_collapse_and_the_thresholds(tmp_path: Path, sequence, plain):
    unit = sequence(2_000, seed=8)
    fasta, info = write_assembly(
        tmp_path, {"contig_1": unit * 2}, {"contig_1": ("Y", 2)}
    )
    output = tmp_path / "sp1.fa"
    arguments = [str(fasta), "-p", "sp1", "-o", str(output), "-a", str(info)]
    assert run(*arguments, "--no-collapse").exit_code == 0
    assert names_table(tmp_path / "sp1.names.tsv")["sp1_contig_1"][2:] == [
        "4000",
        "1",
        "yes",
    ]
    # A unit shorter than --min-unit is a satellite's: left whole.
    assert run(*arguments, "--min-unit", "3000", "--force").exit_code == 0
    assert names_table(tmp_path / "sp1.names.tsv")["sp1_contig_1"][1] == "4000"
    assert run(*arguments, "--min-unit", "1500", "--force").exit_code == 0
    assert names_table(tmp_path / "sp1.names.tsv")["sp1_contig_1"][1] == "2000"
    bad = run(*arguments, "--min-identity", "0.5", "--force")
    assert bad.exit_code == 2
    assert "--min-identity" in plain(bad.output)


def test_a_malformed_assembly_info_fails(fasta: Path, tmp_path: Path, plain):
    info = tmp_path / "assembly_info.txt"
    info.write_text("contig_1\t100\n")
    result = run(
        str(fasta), "-p", "sp1", "-o", str(tmp_path / "sp1.fa"), "-a", str(info)
    )
    assert result.exit_code == 1
    assert "not Flye's assembly_info.txt" in plain(result.output)


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


def test_outputs_are_not_overwritten_without_force(fasta: Path, tmp_path: Path, plain):
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
def test_malformed_contigs_fail(tmp_path: Path, text: str, message: str, plain):
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
