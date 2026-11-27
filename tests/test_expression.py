"""`plexembler expression`: a species' counts and TPM, genes by libraries."""

from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.expression import library_of, read_featurecounts

#: featureCounts' table of two genes in two libraries.
TABLE = (
    "# Program:featureCounts v2.0.6; Command:...\n"
    "Geneid\tChr\tStart\tEnd\tStrand\tLength\t"
    "split/Scer.rna-ypd-1.bam\tScer.rna-gal-1.bam\n"
    "Scer_000010\tScer_chr01\t11\t1010\t+\t1000\t100\t0\n"
    "Scer_000020\tScer_chr01\t2001\t2500\t-\t500\t100\t40\n"
)


@pytest.fixture
def counts(tmp_path: Path) -> Path:
    path = tmp_path / "Scer.featureCounts.tsv"
    path.write_text(TABLE)
    return path


def test_libraries_and_tpm(counts):
    assert library_of("a/b/Scer.rna-ypd-1.bam", "Scer") == "rna-ypd-1"
    assert library_of("other", "Scer") == "other"
    table = read_featurecounts(counts, "Scer", {"rna-ypd-1": "YPD_1"})
    assert table.libraries == ["YPD_1", "rna-gal-1"]
    tpm = table.tpm()
    # Library 1: rates 0.1 and 0.2, a third and two thirds of a million.
    assert tpm[0][0] == pytest.approx(1e6 / 3)
    assert tpm[1][0] == pytest.approx(2e6 / 3)
    # Library 2: all of it on the second gene.
    assert [row[1] for row in tpm] == [0.0, pytest.approx(1e6)]


def test_bad_tables(tmp_path: Path, counts):
    other = tmp_path / "other.tsv"
    other.write_text("gene\tcount\n")
    with pytest.raises(ValueError, match="not a featureCounts table"):
        read_featurecounts(other, "Scer", {})
    with pytest.raises(ValueError, match="named the same"):
        read_featurecounts(counts, "Scer", {"rna-ypd-1": "x", "rna-gal-1": "x"})


def test_the_command(counts, tmp_path: Path):
    outdir = tmp_path / "out"
    arguments = ["expression", str(counts), "-p", "Scer", "-o", str(outdir), "-l", "rna-ypd-1=YPD_1", "-l", "rna-gal-1=GAL_1"]  # fmt: skip
    result = CliRunner().invoke(cli, arguments)
    assert result.exit_code == 0, result.output
    assert (outdir / "Scer.counts.tsv").read_text().splitlines() == [
        "gene\tlength_bp\tYPD_1\tGAL_1",
        "Scer_000010\t1000\t100\t0",
        "Scer_000020\t500\t100\t40",
    ]
    tpm = (outdir / "Scer.tpm.tsv").read_text().splitlines()
    assert tpm[1] == "Scer_000010\t333333.333\t0.000"
    multiqc = (outdir / "Scer.expression_mqc.tsv").read_text()
    assert "Scer_YPD_1\t200\t2\t2" in multiqc
    assert "Scer_GAL_1\t40\t1\t2" in multiqc
    assert CliRunner().invoke(cli, arguments).exit_code == 1
    bad = CliRunner().invoke(cli, [*arguments[:-1], "nolabel", "--force"])
    assert bad.exit_code == 2
    assert "LIBRARY=LABEL" in bad.output
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    failed = CliRunner().invoke(
        cli, ["expression", str(empty), "-p", "Scer", "-o", str(outdir), "--force"]
    )
    assert failed.exit_code == 1
