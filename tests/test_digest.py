"""`plexembler digest`: Hi-C reads cut at their ligation junctions."""

import gzip
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.digest import (
    Junction,
    cut_points,
    ends,
    junctions,
    parse_enzymes,
    read_fastq,
)

LEFT = "ACGTTGCAAC" * 6  # 60 bases, no site
RIGHT = "TTGGCCAAGG" * 6


def fastq(path: Path, records: list[tuple[str, str]]) -> Path:
    text = "".join(
        f"@{name} comment\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in records
    )
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    return path


def test_enzymes_ends_and_junctions():
    assert ends("DpnII") == ("GATC", "GATC")
    assert ends("HinfI") == ("GANT", "ANTC")
    assert ends("HindIII") == ("AAGCT", "AGCTT")
    assert ends("NlaIII") == ("CATG", "CATG")
    assert parse_enzymes("Arima") == ["DpnII", "HinfI"]
    assert parse_enzymes("DpnII, DpnII") == ["DpnII"]
    found = {j.pattern: j.cut for j in junctions(["DpnII", "HinfI"])}
    assert found == {
        "GATCGATC": 4,
        "GATCA[ACGT]TC": 4,
        "GA[ACGT]TGATC": 4,
        "GA[ACGT]TA[ACGT]TC": 4,
    }
    for bad, message in (("none", "no enzyme"), ("XyzI", "not an enzyme")):
        with pytest.raises(ValueError, match=message):
            parse_enzymes(bad)


def test_cut_points():
    patterns = [Junction("GATCGATC", 4)]
    assert cut_points(LEFT + "GATCGATC" + RIGHT, patterns) == [64]
    assert cut_points(LEFT, patterns) == []
    # Two junctions, non-overlapping; lower case is read too.
    read = LEFT + "gatcgatc" + LEFT + "GATCGATC" + RIGHT
    assert cut_points(read, patterns) == [64, 132]


def test_digest_writes_every_piece_once(tmp_path: Path):
    r1 = fastq(
        tmp_path / "r1.fq.gz",
        [
            ("a", LEFT + "GATCGATC" + RIGHT),
            ("b", LEFT + RIGHT),
            ("c", "GATCGATC" + LEFT),
        ],
    )
    r2 = fastq(tmp_path / "r2.fq.gz", [("a", LEFT + "GAATAGTC" + RIGHT)])
    outdir = tmp_path / "out"
    result = CliRunner().invoke(
        cli,
        ["digest", str(r1), str(r2), "-e", "arima", "-p", "Scer", "-o", str(outdir)],
    )
    assert result.exit_code == 0, result.output
    with gzip.open(outdir / "Scer.digested.fastq.gz", "rt") as handle:
        pieces = {name: (seq, qual) for name, seq, qual in read_fastq(handle)}
    # Each piece keeps its half of the site.
    assert pieces["a_1.0"][0] == LEFT + "GATC"
    assert pieces["a_1.1"][0] == "GATC" + RIGHT
    assert pieces["b_1.0"][0] == LEFT + RIGHT
    # HinfI's junction, GANTANTC, in the second read.
    assert pieces["a_2.0"][0] == LEFT + "GAAT"
    # A piece under 40 bases is dropped: c's first.
    assert "c_1.0" not in pieces
    assert pieces["c_1.1"][0] == "GATC" + LEFT
    assert len(pieces["a_1.0"][1]) == len(pieces["a_1.0"][0])
    numbers = json.loads((outdir / "Scer.digest.json").read_text())
    assert numbers["reads"] == 4
    assert (numbers["no_junction"], numbers["one_junction"]) == (1, 3)
    assert numbers["short_pieces"] == 1
    assert "plexembler_digest" in (outdir / "Scer.digest_mqc.tsv").read_text()
    assert (
        "pieces of 40 bases or more"
        in next(outdir.glob("Scer.plexembler_*.log")).read_text()
    )


def test_digest_checks(tmp_path: Path):
    r1 = fastq(tmp_path / "r1.fq.gz", [("a", LEFT)])
    outdir = tmp_path / "out"
    arguments = ["digest", str(r1), "-e", "DpnII", "-p", "Scer", "-o", str(outdir)]
    assert CliRunner().invoke(cli, arguments).exit_code == 0
    assert CliRunner().invoke(cli, arguments).exit_code == 1
    assert CliRunner().invoke(cli, [*arguments, "--force"]).exit_code == 0
    none = CliRunner().invoke(
        cli,
        ["digest", str(r1), "-e", "none", "-p", "Scer", "-o", str(outdir), "--force"],
    )
    assert none.exit_code == 2
    plain = tmp_path / "bad.fq"
    plain.write_text("not a fastq\nACGT\n-\nIIII\n", encoding="ascii")
    bad = CliRunner().invoke(
        cli,
        ["digest", str(plain), "-e", "DpnII", "-p", "Scer", "-o", str(tmp_path / "b")],
    )
    assert bad.exit_code == 1
    assert "not a FASTQ record" in bad.output
