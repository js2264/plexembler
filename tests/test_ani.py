"""`plexembler ani`: the identity between the species of a pool."""

import logging
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.ani import Pair, read_skani, species_of, verdict
from plexembler.cli import cli

HEADER = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name"


def skani(tmp_path: Path, rows: list[str]) -> Path:
    path = tmp_path / "pool.tsv"
    path.write_text("\n".join([HEADER, *rows]) + "\n", encoding="utf-8")
    return path


ROWS = [
    "Scer.fa.gz\tSpar.fa.gz\t88.40\t80.10\t79.00\tScer_ctg1\tSpar_ctg7 length=1000",
    "Scer.fa.gz\tSmik.fa.gz\t96.20\t90.00\t91.50\tScer_ctg1\tSmik_ctg2",
]


def run(*arguments: str):
    return CliRunner().invoke(cli, ["ani", *arguments])


def test_species_from_contigs_or_files():
    assert species_of("Scer_ctg1 length=1000", "x.fa") == "Scer"
    assert species_of("ctg1", "data/Kmar.contigs.fa.gz") == "Kmar"
    assert species_of("", "Kmar.fa") == "Kmar"


def test_verdicts():
    assert verdict(None, 95, 85) == "distinct"
    assert verdict(84.9, 95, 85) == "distinct"
    assert verdict(85.0, 95, 85) == "same genus"
    assert verdict(95.0, 95, 85) == "too close"


def test_read_skani(tmp_path: Path):
    pairs = read_skani(skani(tmp_path, ROWS))
    assert pairs[0] == Pair("Scer", "Spar", 88.4, 80.1, 79.0)
    with pytest.raises(ValueError, match="no column ANI"):
        bad = tmp_path / "bad.tsv"
        bad.write_text("Ref_file\tQuery_file\n", encoding="utf-8")
        read_skani(bad)
    with pytest.raises(ValueError, match="two genomes of 'Scer'"):
        read_skani(skani(tmp_path, ["a.fa\tb.fa\t99\t90\t90\tScer_1\tScer_2"]))


def test_the_matrix_pairs_and_warnings(tmp_path: Path, caplog):
    outdir = tmp_path / "out"
    result = run(
        str(skani(tmp_path, ROWS)), "-p", "poolA", "-o", str(outdir), "-s", "Kmar"
    )
    assert result.exit_code == 0, result.output
    matrix = (outdir / "poolA.ani.tsv").read_text().splitlines()
    assert matrix[0] == "species\tKmar\tScer\tSmik\tSpar"
    assert matrix[2] == "Scer\tNA\t100.00\t96.20\t88.40"
    pairs = (outdir / "poolA.ani_pairs.tsv").read_text().splitlines()
    assert "Scer\tSmik\t96.20\t90.00\t91.50\ttoo close" in pairs
    # A pair is written with its first species' aligned fraction first.
    assert "Scer\tSpar\t88.40\t80.10\t79.00\tsame genus" in pairs
    assert "Kmar\tScer\tNA\tNA\tNA\tdistinct" in pairs
    multiqc = (outdir / "poolA.ani_mqc.tsv").read_text()
    assert "# id: 'plexembler_ani'" in multiqc
    assert "Scer_vs_Smik\tpoolA\t96.20" in multiqc
    (log,) = outdir.glob("poolA.plexembler_*.log")
    text = log.read_text()
    assert ":: WARNING :: Scer and Smik are 96.20 % identical" in text
    assert "Scer and Spar are 88.40 % identical: a pair of the same genus" in text


def test_the_second_species_first(tmp_path: Path):
    rows = ["Spar.fa\tScer.fa\t88.40\t70.00\t60.00\tSpar_1\tScer_1"]
    outdir = tmp_path / "out"
    assert run(str(skani(tmp_path, rows)), "-p", "p", "-o", str(outdir)).exit_code == 0
    pairs = (outdir / "p.ani_pairs.tsv").read_text()
    assert "Scer\tSpar\t88.40\t60.00\t70.00\tsame genus" in pairs


def test_force_and_checks(tmp_path: Path):
    table = skani(tmp_path, ROWS)
    outdir = tmp_path / "out"
    assert run(str(table), "-p", "poolA", "-o", str(outdir), "-q").exit_code == 0
    again = run(str(table), "-p", "poolA", "-o", str(outdir))
    assert again.exit_code == 1
    assert "--force" in again.output
    assert run(str(table), "-p", "poolA", "-o", str(outdir), "--force").exit_code == 0
    assert run(str(table), "-p", "pool_A", "-o", str(outdir)).exit_code == 2
    bad = tmp_path / "bad.tsv"
    bad.write_text("x\n", encoding="utf-8")
    failed = run(str(bad), "-p", "poolA", "-o", str(tmp_path / "o2"))
    assert failed.exit_code == 1
    assert "no column" in failed.output
    assert not logging.getLogger("plexembler").handlers
