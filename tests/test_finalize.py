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
