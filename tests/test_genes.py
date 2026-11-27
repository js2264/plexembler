"""`plexembler genes`: a species' gene models, named, translated and written."""

import gzip
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.genes import annotate, code, read_models, revcomp, translate
from plexembler.rename import read_fasta

#: Two genes: g1 on the plus strand, two exons (ATG AAA CTG | TGG GGG TAA), its
#: CUG a leucine in table 1, a serine in table 12; g2 on the minus strand, one
#: exon (ATG CCC AAA TTT TAG).
GENOME = (
    "C" * 10
    + "ATGAAACTG"
    + "GTAAGTCCCCCCCCCCCCAG"
    + "TGGGGGTAA"
    + "C" * 51
    + revcomp("ATGCCCAAATTTTAG")
    + "C" * 86
)


def rows(*lines):
    return "".join("\t".join(str(c) for c in line) + "\n" for line in lines)


def ids(transcript, gene):
    return f'transcript_id "{transcript}"; gene_id "{gene}";'


#: BRAKER's GTF: g2 first, its stop codon apart from its CDS.
GTF = rows(
    ("Scer_chr01", "AUGUSTUS", "gene", 100, 114, ".", "-", ".", "g2"),
    ("Scer_chr01", "AUGUSTUS", "transcript", 100, 114, ".", "-", ".", "g2.t1"),
    ("Scer_chr01", "AUGUSTUS", "CDS", 103, 114, ".", "-", 0, ids("g2.t1", "g2")),
    ("Scer_chr01", "AUGUSTUS", "stop_codon", 100, 102, ".", "-", 0, ids("g2.t1", "g2")),
    ("Scer_chr01", "AUGUSTUS", "gene", 11, 48, ".", "+", ".", "g1"),
    ("Scer_chr01", "AUGUSTUS", "transcript", 11, 48, ".", "+", ".", "g1.t1"),
    ("Scer_chr01", "AUGUSTUS", "CDS", 11, 19, ".", "+", 0, ids("g1.t1", "g1")),
    ("Scer_chr01", "AUGUSTUS", "exon", 11, 19, ".", "+", ".", ids("g1.t1", "g1")),
    ("Scer_chr01", "AUGUSTUS", "intron", 20, 39, ".", "+", ".", ids("g1.t1", "g1")),
    ("Scer_chr01", "AUGUSTUS", "CDS", 40, 45, ".", "+", 0, ids("g1.t1", "g1")),
    ("Scer_chr01", "AUGUSTUS", "exon", 40, 48, ".", "+", ".", ids("g1.t1", "g1")),
    ("Scer_chr01", "AUGUSTUS", "stop_codon", 46, 48, ".", "+", 0, ids("g1.t1", "g1")),
)  # fmt: skip

#: AUGUSTUS's GFF3 of g1, its header naming its authors.
GFF3 = "##gff-version 3\n# O. Keller, S. König\n" + rows(
    ("Scer_chr01", "AUGUSTUS", "gene", 11, 48, ".", "+", ".", "ID=g1;"),
    ("Scer_chr01", "AUGUSTUS", "transcript", 11, 48, ".", "+", ".", "ID=g1.t1;Parent=g1"),
    ("Scer_chr01", "AUGUSTUS", "CDS", 11, 19, ".", "+", 0, "ID=g1.t1.cds;Parent=g1.t1"),
    ("Scer_chr01", "AUGUSTUS", "CDS", 40, 48, ".", "+", 0, "ID=g1.t1.cds;Parent=g1.t1"),
)  # fmt: skip


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    path = tmp_path / "Scer.fa.gz"
    with gzip.open(path, "wt") as sink:
        sink.write(f">Scer_chr01\n{GENOME}\n")
    return path


def fasta(path: Path) -> dict[str, str]:
    with gzip.open(path, "rt") as handle:
        return dict(read_fasta(handle))


def test_codes():
    assert translate("ATGCTGTAA", code(1)) == "ML*"
    assert translate("ATGCTGTAA", code(12)) == "MS*"
    assert translate("ATGCTGTAA", code(26)) == "MA*"
    assert translate("ATGTGANNN", code(4)) == "MWX"
    with pytest.raises(ValueError, match="genetic code 7"):
        code(7)


def test_braker_gtf(genome, tmp_path: Path):
    models = tmp_path / "braker.gtf"
    models.write_text(GTF)
    genes = annotate(genome, models, "Scer", tmp_path / "out", genetic_code=12)
    assert [(g.new_name, g.name, g.strand) for g in genes] == [
        ("Scer_000010", "g1", "+"),
        ("Scer_000020", "g2", "-"),
    ]
    out = tmp_path / "out"
    assert fasta(out / "Scer.proteins.faa.gz") == {
        "Scer_000010-T1": "MKSWG",
        "Scer_000020-T1": "MPKF",
    }
    assert fasta(out / "Scer.cds.fna.gz")["Scer_000020-T1"] == "ATGCCCAAATTTTAG"
    with gzip.open(out / "Scer.gff3.gz", "rt") as handle:
        gff3 = [line.split("\t") for line in handle if not line.startswith("#")]
    cds = [(r[3], r[4], r[7]) for r in gff3 if r[2] == "CDS" and "000010" in r[8]]
    assert cds == [("11", "19", "0"), ("40", "48", "0")]
    assert [r[2] for r in gff3[:5]] == ["gene", "mRNA", "exon", "exon", "CDS"]
    with gzip.open(out / "Scer.gtf.gz", "rt") as handle:
        assert 'gene_id "Scer_000020"; transcript_id "Scer_000020-T1";' in handle.read()
    table = (out / "Scer.genes.tsv").read_text().splitlines()
    assert table[1].split("\t") == [
        "Scer_000010", "1", "Scer_chr01", "11", "48", "+", "2", "5", "g1",
    ]  # fmt: skip
    assert "Scer\t2\t2\t12\t0" in (out / "Scer.genes_mqc.tsv").read_text()


def test_augustus_gff3(genome, tmp_path: Path):
    models = tmp_path / "augustus.gff3"
    models.write_text(GFF3)
    (gene,) = read_models(models)
    assert gene.name == "g1"
    genes = annotate(genome, models, "Scer", tmp_path / "out")
    assert fasta(tmp_path / "out" / "Scer.proteins.faa.gz") == {
        "Scer_000010-T1": "MKLWG"
    }
    assert genes[0].transcripts[0].exons == [(11, 19), (40, 48)]


def test_a_wrong_code_warns(genome, tmp_path: Path, caplog):
    """A table that reads a codon of the genes as a stop: proteins with internal stops."""
    models = tmp_path / "braker.gtf"
    # g1 read as if its CTG were TGA: table 1 against a gene of table 4.
    models.write_text(GTF)
    stop = tmp_path / "stop.fa.gz"
    with gzip.open(stop, "wt") as sink:
        sink.write(f">Scer_chr01\n{GENOME.replace('ATGAAACTG', 'ATGAAATGA')}\n")
    annotate(stop, models, "Scer", tmp_path / "out")
    assert "internal stop codon" in caplog.text


def test_models_off_the_genome(genome, tmp_path: Path):
    models = tmp_path / "braker.gtf"
    models.write_text(GTF.replace("Scer_chr01", "Scer_chr09"))
    with pytest.raises(ValueError, match="not in the genome"):
        annotate(genome, models, "Scer", tmp_path / "out")
    models.write_text("Scer_chr01\tAUGUSTUS\tCDS\n")
    with pytest.raises(ValueError, match="3 columns"):
        annotate(genome, models, "Scer", tmp_path / "out")


def test_the_command(genome, tmp_path: Path, plain):
    models = tmp_path / "braker.gtf"
    models.write_text(GTF)
    outdir = tmp_path / "out"
    arguments = ["genes", str(genome), str(models), "-p", "Scer", "-o", str(outdir), "-g", "12", "-s", "BRAKER3"]  # fmt: skip
    result = CliRunner().invoke(cli, arguments)
    assert result.exit_code == 0, result.output
    with gzip.open(outdir / "Scer.gff3.gz", "rt") as handle:
        assert "\tBRAKER3\tgene\t" in handle.read()
    log = next(outdir.glob("Scer.plexembler_*.log")).read_text()
    assert "2 genes, 2 transcripts, genetic code 12" in log
    assert CliRunner().invoke(cli, arguments).exit_code == 1
    models.write_text(GTF.replace("Scer_chr01", "Scer_chr09"))
    bad = CliRunner().invoke(cli, [*arguments, "--force"])
    assert bad.exit_code == 1
    assert "not in the genome" in plain(bad.output)
