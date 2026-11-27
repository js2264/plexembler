"""`plexembler dotplot`: a species' genome against another, drawn and tabulated."""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.dotplot import (
    Sequence,
    _unique,
    label_of,
    read_headers,
    read_paf,
    reference_title,
)

#: The other genome: two chromosomes, a mitochondrion and an unlocalized scaffold.
REFERENCE = (
    ">refA Org strain chromosome I, complete sequence\n" + "A" * 1000 + "\n"
    ">refB Org strain chromosome II, complete sequence\n" + "C" * 600 + "\n"
    ">refM Org strain mitochondrion, complete genome\n" + "G" * 100 + "\n"
    ">refC Org strain chromosome I unlocalized genomic scaffold\n" + "T" * 50 + "\n"
)
#: The species' genome: its chr01 is refB turned round, its chr02 refA with an
#: inversion; its unplaced scaffold aligns nowhere.
GENOME = (
    ">Sp_chr01\n" + "C" * 600 + "\n"
    ">Sp_chr02\n" + "A" * 1000 + "\n"
    ">Sp_unplaced_01\n" + "T" * 80 + "\n"
)
PAF = [
    "Sp_chr02\t1000\t0\t600\t+\trefA\t1000\t0\t600\t590\t600\t60\ttp:A:P",
    "Sp_chr02\t1000\t600\t700\t-\trefA\t1000\t600\t700\t95\t100\t60\ttp:A:P",
    "Sp_chr02\t1000\t700\t1000\t+\trefA\t1000\t700\t1000\t297\t300\t60\ttp:A:P",
    "Sp_chr01\t600\t0\t600\t-\trefB\t600\t0\t600\t594\t600\t60\ttp:A:P",
    "Sp_chr01\t600\t100\t200\t+\trefA\t1000\t800\t900\t90\t100\t0\ttp:A:S",
]
#: The sign that marks a sequence turned round on the plot.
MINUS = "\u2212"
CLOSEST = (
    "rank\taccession\tname\tani\tspecies_containment\tmatch_containment\tshared_bp\tlineage\tsequence\n"
    "1\tGCA_000146045.2\tGCA_000146045.2 Saccharomyces cerevisiae S288C\t99.30\t80.00\t24.00\t1\tNA\tdownloaded\n"
)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def inputs(tmp_path: Path) -> dict[str, Path]:
    return {
        "paf": write(tmp_path / "Sp.paf", "\n".join(PAF) + "\n"),
        "genome": write(tmp_path / "Sp.fa", GENOME),
        "reference": write(tmp_path / "Sp.closest.fa", REFERENCE),
        "closest": write(tmp_path / "Sp.closest.tsv", CLOSEST),
    }


def run(*arguments: str):
    return CliRunner().invoke(cli, ["dotplot", *arguments])


def rows(path: Path) -> dict[str, list[str]]:
    lines = path.read_text().splitlines()
    assert lines[0] == (
        "genome\tsequence\tlabel\tlength_bp\taligned_bp\taligned_fraction\tidentity"
        "\tbest_match\tbest_match_share\tstrand"
    )
    return {line.split("\t")[1]: line.split("\t") for line in lines[1:]}


def test_labels_from_headers():
    assert (
        label_of(
            "BK006935.2", "TPA_inf: S. cerevisiae S288C chromosome I, complete sequence"
        )
        == "chrI"
    )
    assert label_of("NC_1", "Candida albicans SC5314 chromosome R") == "chrR"
    assert (
        label_of("KP263414.1", "S. cerevisiae S288C mitochondrion, complete genome")
        == "chrM"
    )
    assert label_of("X1", "Some alga chloroplast, complete genome") == "chrC"
    assert (
        label_of("X2", "S. cerevisiae plasmid 2-micron, complete sequence")
        == "plasmid 2-micron"
    )
    assert (
        label_of("JAA01.1", "S. cerevisiae contig_12, whole genome shotgun")
        == "JAA01.1"
    )
    assert label_of("X3", "chromosome I unlocalized genomic scaffold") == "X3"
    assert label_of("X4", "") == "X4"
    same = _unique(
        [Sequence("a", 1, "chrI"), Sequence("b", 1, "chrI"), Sequence("c", 1, "chrII")]
    )
    assert [s.label for s in same] == ["a", "b", "chrII"]  # fmt: skip


def test_read_headers_and_paf(inputs, tmp_path: Path):
    assert read_headers(inputs["reference"])[2] == (
        "refM",
        100,
        "Org strain mitochondrion, complete genome",
    )
    alignments = read_paf(inputs["paf"])
    assert len(alignments) == 5
    assert alignments[1].reverse
    assert not alignments[4].primary
    with pytest.raises(ValueError, match="a header without a name"):
        read_headers(write(tmp_path / "a.fa", ">\nACGT\n"))
    with pytest.raises(ValueError, match="a sequence before the first header"):
        read_headers(write(tmp_path / "b.fa", "ACGT\n>x\nA\n"))
    with pytest.raises(ValueError, match="line 1: not a PAF line"):
        read_paf(write(tmp_path / "c.paf", "a\tb\tc\n"))
    with pytest.raises(ValueError, match="line 2: invalid literal"):
        read_paf(
            write(
                tmp_path / "d.paf", PAF[0] + "\n" + PAF[1].replace("600", "x", 1) + "\n"
            )
        )


def test_reference_title(tmp_path: Path):
    assert reference_title(write(tmp_path / "a.tsv", CLOSEST)) == (
        "Saccharomyces cerevisiae S288C (GCA_000146045.2)"
    )
    custom = CLOSEST.replace(
        "GCA_000146045.2\tGCA_000146045.2 Saccharomyces", "NA\tmy own"
    )
    assert (
        reference_title(write(tmp_path / "b.tsv", custom)) == "my own cerevisiae S288C"
    )
    assert (
        reference_title(write(tmp_path / "c.tsv", CLOSEST.splitlines()[0] + "\n"))
        is None
    )


def test_the_dotplot_its_figure_and_table(inputs, tmp_path: Path):
    outdir = tmp_path / "out"
    result = run(
        str(inputs["paf"]),
        "-p", "Sp", "-o", str(outdir),
        "-g", str(inputs["genome"]), "-r", str(inputs["reference"]),
        "-c", str(inputs["closest"]), "--min-length", "50",
    )  # fmt: skip
    assert result.exit_code == 0, result.output

    table = rows(outdir / "Sp.dotplot.tsv")
    # The species' sequences in its order, then the reference's in its order.
    assert list(table) == ["Sp_chr01", "Sp_chr02", "Sp_unplaced_01", "refA", "refB", "refM", "refC"]  # fmt: skip
    assert table["Sp_chr01"] == ["species", "Sp_chr01", "chr01", "600", "600", "100.00", "99.00", "refB", "100.00", "-"]  # fmt: skip
    assert table["Sp_chr02"][4:] == ["1000", "100.00", "98.20", "refA", "100.00", "+"]
    assert table["Sp_unplaced_01"][4:] == ["0", "0.00", "NA", "NA", "NA", "NA"]
    assert table["refA"][:3] == ["reference", "refA", "chrI"]
    assert table["refA"][4:8] == ["1000", "100.00", "98.20", "Sp_chr02"]
    assert table["refB"][-1] == "-"
    assert table["refM"][2] == "chrM"
    assert table["refC"][2] == "refC"

    svg = (outdir / "Sp.dotplot.svg").read_text()
    assert svg.startswith('<svg xmlns="http://www.w3.org/2000/svg"')
    assert "Sp against Saccharomyces cerevisiae S288C (GCA_000146045.2)" in svg
    assert f"chr01 ({MINUS})" in svg  # turned round: it aligns to refB's other strand
    assert 'stroke="#0072B2"' in svg and 'stroke="#D55E00"' in svg
    assert ">chrI<" in svg and ">chrII<" in svg

    figure = json.loads((outdir / "Sp.dotplot.json").read_text())
    same, other = figure["data"]
    assert (same["name"], other["name"]) == ("same strand", "opposite strand")
    # Five alignments drawn: three on the same strand (chr01 turned round,
    # chr02's two outer parts), two on the other (the inversion, the secondary).
    assert same["x"].count(None) == 3
    assert other["x"].count(None) == 2
    layout = figure["layout"]
    assert layout["xaxis"]["ticktext"] == ["chrI", "chrII", "chrM", "refC"]
    # chr02 first, beside chrI; chr01 beside chrII; the unplaced scaffold last.
    assert layout["yaxis"]["ticktext"] == ["chr02", f"chr01 ({MINUS})", "unplaced_01"]
    assert layout["yaxis"]["range"] == [0, 1680]
    assert "chr01:600 ↔ chrII:0" in same["text"]

    (log,) = outdir.glob("Sp.plexembler_*.log")
    text = log.read_text()
    assert "Sp: turned round on the plot: chr01" in text
    assert "Drew 5 alignments of 5" in text


def test_the_shortest_alignments_are_not_drawn(inputs, tmp_path: Path):
    outdir = tmp_path / "out"
    result = run(
        str(inputs["paf"]), "-p", "Sp", "-o", str(outdir), "--min-length", "200"
    )
    assert result.exit_code == 0, result.output
    figure = json.loads((outdir / "Sp.dotplot.json").read_text())
    assert [trace["x"].count(None) for trace in figure["data"]] == [3, 0]
    # Without the genomes' files: the sequences the alignments name, the
    # longest first, and a reference named after nothing.
    assert figure["layout"]["xaxis"]["title"]["text"] == "reference"
    assert figure["layout"]["xaxis"]["ticktext"] == ["refA", "refB"]
    table = rows(outdir / "Sp.dotplot.tsv")
    assert list(table) == ["Sp_chr02", "Sp_chr01", "refA", "refB"]


def test_an_empty_dotplot(inputs, tmp_path: Path, plain):
    outdir = tmp_path / "out"
    empty = write(tmp_path / "empty.paf", "")
    result = run(
        str(empty), "-p", "Sp", "-o", str(outdir),
        "-g", str(inputs["genome"]), "-r", str(inputs["reference"]),
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    assert (
        "no alignment of 1000 bases or more to Sp.closest: an empty dotplot"
        in plain(result.output)
    )
    assert (
        "No alignment of 1,000 bases or more" in (outdir / "Sp.dotplot.svg").read_text()
    )
    table = rows(outdir / "Sp.dotplot.tsv")
    assert all(row[4] == "0" for row in table.values())


def test_inputs_that_do_not_agree(inputs, tmp_path: Path, plain):
    other = write(tmp_path / "other.fa", ">Sp_chr09\nACGT\n")
    result = run(
        str(inputs["paf"]), "-p", "Sp", "-o", str(tmp_path / "a"), "-g", str(other)
    )
    assert result.exit_code != 0
    assert "Sp_chr02 is not a sequence of" in plain(result.output)
    result = run(
        str(inputs["paf"]), "-p", "Sp", "-o", str(tmp_path / "b"), "-r", str(other)
    )
    assert result.exit_code != 0
    assert "refA is not a sequence of" in plain(result.output)
    outdir = tmp_path / "c"
    assert run(str(inputs["paf"]), "-p", "Sp", "-o", str(outdir)).exit_code == 0
    result = run(str(inputs["paf"]), "-p", "Sp", "-o", str(outdir))
    assert result.exit_code != 0
    assert "pass --force" in plain(result.output)


def test_a_draft_reference_keeps_its_figure_light(tmp_path: Path):
    # A reference of a large chromosome and 300 small contigs: the figure
    # names and bounds the chromosome alone.
    contigs = "".join(f">ctg{n}\n{'A' * 10}\n" for n in range(300))
    reference = write(
        tmp_path / "draft.fa", ">big Org chromosome I\n" + "C" * 100000 + "\n" + contigs
    )
    paf = write(
        tmp_path / "Sp.paf",
        "Sp_chr01\t100000\t0\t100000\t+\tbig\t100000\t0\t100000\t99000\t100000\t60\ttp:A:P\n",
    )
    outdir = tmp_path / "out"
    result = run(str(paf), "-p", "Sp", "-o", str(outdir), "-r", str(reference))
    assert result.exit_code == 0, result.output
    layout = json.loads((outdir / "Sp.dotplot.json").read_text())["layout"]
    assert layout["xaxis"]["ticktext"] == ["chrI"]
    assert [shape["x0"] for shape in layout["shapes"]] == []
    # The SVG draws every bound all the same.
    assert (outdir / "Sp.dotplot.svg").read_text().count("V804") == 300
