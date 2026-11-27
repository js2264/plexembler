"""`plexembler finalize`: a species' scaffolds named, its genome written."""

import gzip
from pathlib import Path

from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.finalize import name_scaffolds
from plexembler.rename import read_fasta

SCAFFOLDS = {
    "scaffold_1": "A" * 300,
    "scaffold_2": "C" * 500,
    "scaffold_3": "G" * 50,
    "scaffold_4": "T" * 80,
    "scaffold_5": "A" * 20,
}
AGP = """##agp-version\t2.1
scaffold_1\t1\t150\t1\tW\tScer_contig_1\t1\t150\t+
scaffold_1\t151\t250\t2\tN\t100\tscaffold\tyes\tproximity_ligation
scaffold_1\t251\t300\t3\tW\tScer_contig_2\t1\t50\t-
scaffold_2\t1\t500\t1\tW\tScer_contig_3\t1\t500\t+
scaffold_3\t1\t50\t1\tW\tScer_contig_4\t1\t50\t+
scaffold_4\t1\t80\t1\tW\tScer_chrM\t1\t80\t+
scaffold_5\t1\t20\t1\tW\tScer_contig_5\t1\t20\t+
"""


def files(tmp_path: Path) -> tuple[Path, Path]:
    fasta = tmp_path / "scaffolds.fa"
    fasta.write_text(
        "".join(f">{n}\n{s}\n" for n, s in SCAFFOLDS.items()), encoding="ascii"
    )
    agp = tmp_path / "scaffolds.agp"
    agp.write_text(AGP, encoding="ascii")
    return fasta, agp


def test_names_by_role_and_size():
    lengths = {n: len(s) for n, s in SCAFFOLDS.items()}
    named = {
        s.original_name: (s.name, s.role)
        for s in name_scaffolds("Scer", lengths, n_chromosomes=2)
    }
    assert named["scaffold_2"] == ("Scer_chr01", "chromosome")
    assert named["scaffold_1"] == ("Scer_chr02", "chromosome")
    assert named["scaffold_4"] == ("Scer_unplaced_01", "unplaced")
    # Without a count: the scaffolds over the minimum are chromosomes.
    by_size = name_scaffolds(
        "Scer", lengths, min_chromosome=60, organelles={"scaffold_4"}
    )
    assert [(s.name, s.original_name) for s in by_size] == [
        ("Scer_chr01", "scaffold_2"),
        ("Scer_chr02", "scaffold_1"),
        ("Scer_chrM", "scaffold_4"),
        ("Scer_unplaced_01", "scaffold_3"),
        ("Scer_unplaced_02", "scaffold_5"),
    ]
    assert (
        name_scaffolds("Scer", {"Scer_chrM": 10, "Scer_chrM_b": 5})[1].name
        == "Scer_chrM2"
    )


def test_the_genome_and_its_files(tmp_path: Path):
    fasta, agp = files(tmp_path)
    outdir = tmp_path / "out"
    result = CliRunner().invoke(
        cli,
        [
            "finalize",
            str(fasta),
            "-p",
            "Scer",
            "-o",
            str(outdir),
            "-a",
            str(agp),
            "-n",
            "2",
        ],
    )
    assert result.exit_code == 0, result.output
    with gzip.open(outdir / "Scer.fa.gz", "rt") as handle:
        genome = dict(read_fasta(handle))
    # The scaffold of the mitochondrion's contig is the organelle, found through the AGP.
    assert list(genome) == [
        "Scer_chr01",
        "Scer_chr02",
        "Scer_chrM",
        "Scer_unplaced_01",
        "Scer_unplaced_02",
    ]
    assert genome["Scer_chr01"] == "C" * 500
    sizes = (outdir / "Scer.chrom.sizes").read_text().splitlines()
    assert sizes[:2] == ["Scer_chr01\t500", "Scer_chr02\t300"]
    names = (outdir / "Scer.names.tsv").read_text().splitlines()
    assert names[3] == "Scer_chrM\tscaffold_4\t80\torganelle"
    renamed = (outdir / "Scer.agp").read_text().splitlines()
    assert renamed[0] == "##agp-version\t2.1"
    assert renamed[1].startswith("Scer_chr01\t1\t500\t1\tW\tScer_contig_3")
    assert renamed[2].startswith("Scer_chr02\t1\t150")
    assert "chromosomes" in next(outdir.glob("Scer.plexembler_*.log")).read_text()
    again = CliRunner().invoke(
        cli, ["finalize", str(fasta), "-p", "Scer", "-o", str(outdir)]
    )
    assert again.exit_code == 1


def test_an_agp_that_names_another_scaffold(tmp_path: Path, plain):
    fasta, agp = files(tmp_path)
    agp.write_text(
        AGP + "scaffold_9\t1\t10\t1\tW\tScer_contig_9\t1\t10\t+\n", encoding="ascii"
    )
    result = CliRunner().invoke(
        cli,
        [
            "finalize",
            str(fasta),
            "-p",
            "Scer",
            "-o",
            str(tmp_path / "o"),
            "-a",
            str(agp),
        ],
    )
    assert result.exit_code == 1
    assert "not in the scaffolds" in plain(result.output)


#: A scaffold of three contigs, the second reversed, and one of a contig alone.
CUT_AGP = """##agp-version\t2.1
scaffold_1\t1\t150\t1\tW\tScer_contig_1\t1\t150\t+
scaffold_1\t151\t250\t2\tU\t100\tscaffold\tyes\tproximity_ligation
scaffold_1\t251\t300\t3\tW\tScer_contig_2\t1\t50\t-
scaffold_1\t301\t400\t4\tU\t100\tscaffold\tyes\tproximity_ligation
scaffold_1\t401\t600\t5\tW\tScer_contig_3\t1\t200\t+
scaffold_2\t1\t500\t1\tW\tScer_contig_4\t1\t500\t+
"""
JOINS = "scaffold\tposition\tleft_contig\tright_contig\tobserved\texpected\tratio\tverdict\n"


def cut_files(tmp_path: Path, joins: str) -> list[str]:
    fasta = tmp_path / "scaffolds.fa"
    scaffold_1 = "A" * 150 + "N" * 100 + "C" * 50 + "N" * 100 + "G" * 200
    fasta.write_text(f">scaffold_1\n{scaffold_1}\n>scaffold_2\n{'T' * 500}\n")
    agp = tmp_path / "scaffolds.agp"
    agp.write_text(CUT_AGP, encoding="ascii")
    table = tmp_path / "scaffolds.junctions.tsv"
    table.write_text(JOINS + joins, encoding="ascii")
    return ["finalize", str(fasta), "-p", "Scer", "-a", str(agp), "--cut", str(table)]


def test_a_join_cut(tmp_path: Path):
    """The unsupported join is cut, the thin one kept: two scaffolds of the
    first, named as any other, their AGP lines from 1."""
    arguments = cut_files(
        tmp_path,
        "scaffold_1\t150\tScer_contig_1\tScer_contig_2\t0\t40.0\t0.000\tunsupported\n"
        "scaffold_1\t300\tScer_contig_2\tScer_contig_3\t3\t5.0\t0.600\tthin\n",
    )
    outdir = tmp_path / "out"
    result = CliRunner().invoke(cli, [*arguments, "-o", str(outdir), "-n", "2"])
    assert result.exit_code == 0, result.output
    with gzip.open(outdir / "Scer.fa.gz", "rt") as handle:
        genome = dict(read_fasta(handle))
    assert genome == {
        "Scer_chr01": "T" * 500,
        "Scer_chr02": "C" * 50 + "N" * 100 + "G" * 200,
        "Scer_unplaced_01": "A" * 150,
    }
    names = (outdir / "Scer.names.tsv").read_text().splitlines()
    assert names[2:] == [
        "Scer_chr02\tscaffold_1.2\t350\tchromosome",
        "Scer_unplaced_01\tscaffold_1.1\t150\tunplaced",
    ]
    assert (outdir / "Scer.agp").read_text().splitlines()[1:] == [
        "Scer_chr01\t1\t500\t1\tW\tScer_contig_4\t1\t500\t+",
        "Scer_chr02\t1\t50\t1\tW\tScer_contig_2\t1\t50\t-",
        "Scer_chr02\t51\t150\t2\tU\t100\tscaffold\tyes\tproximity_ligation",
        "Scer_chr02\t151\t350\t3\tW\tScer_contig_3\t1\t200\t+",
        "Scer_unplaced_01\t1\t150\t1\tW\tScer_contig_1\t1\t150\t+",
    ]
    log = next(outdir.glob("Scer.plexembler_*.log")).read_text()
    assert (
        "cut scaffold_1 after 150 bases, between Scer_contig_1 and Scer_contig_2" in log
    )
    assert "(0 pairs across, 40.0 expected)" in log


def test_both_joins_cut(tmp_path: Path):
    arguments = cut_files(
        tmp_path,
        "scaffold_1\t300\tScer_contig_2\tScer_contig_3\t0\t40.0\t0.000\tunsupported\n"
        "scaffold_1\t150\tScer_contig_1\tScer_contig_2\t1\t40.0\t0.025\tunsupported\n",
    )
    outdir = tmp_path / "out"
    result = CliRunner().invoke(cli, [*arguments, "-o", str(outdir), "-n", "2"])
    assert result.exit_code == 0, result.output
    names = [
        line.split("\t")[1:3]
        for line in (outdir / "Scer.names.tsv").read_text().splitlines()[1:]
    ]
    assert names == [
        ["scaffold_2", "500"],
        ["scaffold_1.3", "200"],
        ["scaffold_1.1", "150"],
        ["scaffold_1.2", "50"],
    ]


def test_no_join_to_cut(tmp_path: Path):
    """Every join supported or thin: the genome as without --cut."""
    arguments = cut_files(
        tmp_path,
        "scaffold_1\t150\tScer_contig_1\tScer_contig_2\t30\t40.0\t0.750\tsupported\n",
    )
    result = CliRunner().invoke(cli, [*arguments, "-o", str(tmp_path / "out")])
    assert result.exit_code == 0, result.output
    sizes = (tmp_path / "out" / "Scer.chrom.sizes").read_text().splitlines()
    assert sizes == ["Scer_unplaced_01\t600", "Scer_unplaced_02\t500"]


def test_cuts_that_are_not_joins(tmp_path: Path, plain):
    bad = {
        "scaffold_1\t151\tScer_contig_1\tScer_contig_2\t0\t40.0\t0.000\tunsupported\n": "no join of Scer_contig_1 and Scer_contig_2 at 151 in scaffold_1",
        "scaffold_1\t600\tScer_contig_3\tScer_contig_4\t0\t40.0\t0.000\tunsupported\n": "no join of Scer_contig_3 and Scer_contig_4",
        "scaffold_9\t150\tScer_contig_1\tScer_contig_2\t0\t40.0\t0.000\tunsupported\n": "no scaffold scaffold_9 to cut",
        "scaffold_1\tten\tScer_contig_1\tScer_contig_2\t0\t40.0\t0.000\tunsupported\n": "not a position",
    }  # fmt: skip
    for number, (joins, message) in enumerate(bad.items()):
        arguments = cut_files(tmp_path, joins)
        result = CliRunner().invoke(
            cli, [*arguments, "-o", str(tmp_path / f"o{number}")]
        )
        assert result.exit_code == 1, joins
        assert message in plain(result.output), plain(result.output)
    arguments = cut_files(tmp_path, "")
    (tmp_path / "scaffolds.junctions.tsv").write_text("scaffold\tposition\n")
    result = CliRunner().invoke(cli, [*arguments, "-o", str(tmp_path / "o5")])
    assert "not a table of plexembler junctions" in plain(result.output)
    without_agp = [*arguments[:4], *arguments[6:], "-o", str(tmp_path / "o6")]
    (tmp_path / "scaffolds.junctions.tsv").write_text(JOINS)
    result = CliRunner().invoke(cli, without_agp)
    assert result.exit_code == 1
    assert "needs its AGP" in plain(result.output)
    # The AGP longer than the scaffold it lays out.
    (tmp_path / "scaffolds.junctions.tsv").write_text(
        JOINS
        + "scaffold_1\t150\tScer_contig_1\tScer_contig_2\t0\t40.0\t0.000\tunsupported\n"
    )
    (tmp_path / "scaffolds.fa").write_text(">scaffold_1\nACGT\n>scaffold_2\nACGT\n")
    result = CliRunner().invoke(cli, [*arguments, "-o", str(tmp_path / "o7")])
    assert "600 bases long in its AGP, 4 in its FASTA" in plain(result.output)
