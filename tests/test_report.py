"""`plexembler report`: a species' report, from the files the pipeline wrote."""

from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.report import build, sort_files

BUSCO = "\tC:95.5%[S:94.0%,D:1.5%],F:2.0%,M:2.5%,n:2137\n"


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def files(tmp_path: Path) -> list[Path]:
    return [
        write(
            tmp_path / "Scer.assembly.selection.tsv",
            "candidate\tkept\tbusco_complete\tbusco_single\tbusco_duplicated\tbusco_fragmented\t"
            "busco_missing\tbusco_n\tqv\tcontigs\tcontig_n50\ttotal_length\treason\n"
            "flye\tyes\t95.5\t94.0\t1.5\t2.0\t2.5\t2137\t45.2\t18\t800000\t12000000\tthe highest QV\n"
            "hifiasm\tno\t95.0\t93.0\t2.0\t2.0\t3.0\t2137\tinf\t25\t600000\t12100000\t\n",
        ),
        write(
            tmp_path / "Scer.polish.selection.tsv",
            "candidate\tkept\tqv\treason\ndraft\tno\t40.0\t\npypolca\tyes\t48.1\tthe highest QV\n",
        ),
        write(
            tmp_path / "poolA.poolA-wgs.split.tsv",
            "species\tpairs\tshare\tholdout\tduplicates\nScer\t9000\t0.6000\t1000\t0\nCgla\t6000\t0.4000\t600\t0\n",
        ),
        write(
            tmp_path / "poolA.screen.tsv",
            "contig\tspecies\tverdict\tnew_name\n"
            "Scer_c1\tScer\tkeep\tScer_c1\nScer_c2\tScer\tmove\tCgla_Scer_c2\nScer_m\tScer\torganelle\tScer_chrM\n"
            "Cgla_c1\tCgla\tkeep\tCgla_c1\n",
        ),
        write(
            tmp_path / "Scer.names.tsv",
            "name\toriginal_name\tlength_bp\trole\n"
            "Scer_chr01\tscaffold_1\t1500000\tchromosome\nScer_chr02\tscaffold_2\t1000000\tchromosome\n"
            "Scer_chrM\tscaffold_3\t85000\torganelle\nScer_unplaced_01\tscaffold_4\t15000\tunplaced\n",
        ),
        write(tmp_path / "short_summary.specific.saccharomycetes_odb12.Scer.fa.txt", BUSCO),
        write(tmp_path / "Scer.genome.wgs.qv", "Scer.fa\t100\t12000000\t48.12\t1.5e-05\n"),
        write(
            tmp_path / "Scer.junctions.tsv",
            "scaffold\tposition\tleft_contig\tright_contig\tobserved\texpected\tratio\tverdict\n"
            "Scer_chr01\t700000\ta\tb\t500\t480.0\t1.042\tsupported\n",
        ),
        write(
            tmp_path / "Scer.genes.tsv",
            "gene\ttranscripts\tchromosome\tstart\tend\tstrand\texons\tprotein_length\toriginal_name\n"
            "Scer_000010\t1\tScer_chr01\t11\t1010\t+\t2\t300\tg1\nScer_000020\t2\tScer_chr01\t2001\t2500\t-\t1\t150\tg2\n",
        ),
        write(tmp_path / "short_summary.specific.saccharomycetes_odb12.Scer.proteins.faa.txt", BUSCO.replace("95.5", "93.0")),
        write(tmp_path / "Scer.counts.tsv", "gene\tlength_bp\tYPD_1\tGAL_1\nScer_000010\t1000\t10\t20\nScer_000020\t500\t5\t0\n"),
        write(tmp_path / "unknown.txt", "nothing\n"),
    ]  # fmt: skip


def test_files_sorted(files):
    found = sort_files(files, "Scer")
    assert found.busco.name.endswith("Scer.fa.txt")
    assert found.proteins_busco.name.endswith("Scer.proteins.faa.txt")
    assert found.qv.name == "Scer.genome.wgs.qv"
    assert [p.name for p in found.splits] == ["poolA.poolA-wgs.split.tsv"]


def test_the_report(files):
    report = build("Scer", "poolA", files)
    assert report.sections["assembly"]["kept"] == "flye"
    assert report.sections["polishing"]["qv"] == 48.1
    assert report.sections["screen"] == {
        "keep": 1,
        "move": 1,
        "drop": 0,
        "organelle": 1,
    }
    assert report.sections["genome"]["chromosomes"] == 2
    assert report.sections["genome"]["length_mb"] == 2.6
    assert report.sections["quality"]["busco_complete"] == 95.5
    assert report.sections["quality"]["qv"] == 48.12
    assert report.sections["scaffolding"]["supported"] == 1
    assert report.sections["annotation"] == {
        "genes": 2,
        "proteins_busco_complete": 93.0,
        "transcripts": 3,
        "genes_with_introns": 1,
    }
    assert report.sections["expression"] == {"libraries": 2, "pairs_on_genes": 35}
    assert report.tables["libraries"].rows == [("poolA-wgs", 9000, 0.6)]
    # An infinite QV is written as missing.
    assert report.tables["candidates"].rows[1][3] is None


def test_a_species_with_its_genome_alone(tmp_path: Path, files):
    names = next(f for f in files if f.name == "Scer.names.tsv")
    report = build("Scer", None, [names])
    assert report.sections["quality"] == {"busco_complete": None, "qv": None}
    assert "expression" not in report.sections


def test_the_command(files, tmp_path: Path):
    outdir = tmp_path / "out"
    arguments = [
        "report",
        *map(str, files),
        "-p",
        "Scer",
        "--pool",
        "poolA",
        "-o",
        str(outdir),
    ]
    result = CliRunner().invoke(cli, arguments)
    assert result.exit_code == 0, result.output
    text = (outdir / "Scer.report.txt").read_text()
    assert "[genome]" in text and "chromosomes" in text
    assert "plexembler_report" in (outdir / "Scer.report_mqc.yaml").read_text()
    assert "<html" in (outdir / "Scer.report.html").read_text().lower()
    log = next(outdir.glob("Scer.plexembler_*.log")).read_text()
    assert "2600000 bases in 2 chromosomes; BUSCO 95.5 %, QV 48.12" in log
    assert CliRunner().invoke(cli, arguments).exit_code == 1
    bad = tmp_path / "Scer.genes.tsv"
    bad.write_text("gene\n1\n")
    failed = CliRunner().invoke(
        cli, ["report", str(bad), "-p", "Scer", "-o", str(outdir), "--force"]
    )
    assert failed.exit_code == 1
