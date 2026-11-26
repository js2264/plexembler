"""`plexembler select`: a step's candidates, their checks, and the one kept."""

import gzip
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.select import (
    Busco,
    Candidate,
    choose_assembly,
    copy_sequence,
    read_busco,
    read_candidates,
    read_gfastats,
    read_qv,
)

BUSCO_TXT = """# BUSCO version is: 6.1.0
# The lineage dataset is: saccharomycetes_odb12 (Creation date: 2025-07-01, number of genomes: 120, number of BUSCOs: 2137)

    ***** Results: *****

    C:{c}%[S:{s}%,D:{d}%],F:0.3%,M:1.2%,n:2137
    2105    Complete BUSCOs (C)
"""

COMPLEASM = """## lineage: saccharomycetes_odb10
S:97.90%, 2092
D:0.61%, 13
F:0.20%, 4
I:0.00%, 0
M:1.29%, 28
N:2137
"""

GFASTATS_TABULAR = """# scaffolds\t17
Total scaffold length\t12071326
# contigs\t{contigs}
Contig N50\t{n50}
Total contig length\t12071326
Largest contig\t1531933
GC content %\t38.15
"""


def summary(tmp_path: Path, name: str, c=98.5, s=97.9, d=0.6) -> Path:
    path = tmp_path / f"short_summary.{name}.txt"
    path.write_text(BUSCO_TXT.format(c=c, s=s, d=d), encoding="utf-8")
    return path


def merqury(tmp_path: Path, name: str, qv: str) -> Path:
    path = tmp_path / f"{name}.qv"
    path.write_text(f"{name}\t120\t12000000\t{qv}\t3.4e-05\n", encoding="utf-8")
    return path


def gfastats(tmp_path: Path, name: str, contigs=17, n50=924431) -> Path:
    path = tmp_path / f"{name}.assembly_summary"
    path.write_text(GFASTATS_TABULAR.format(contigs=contigs, n50=n50), "utf-8")
    return path


def fasta(tmp_path: Path, name: str) -> Path:
    path = tmp_path / f"{name}.fa.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(f">Scer_{name}_1\nACGT\n")
    return path


def table(tmp_path: Path, rows: list[dict[str, str]]) -> Path:
    path = tmp_path / "candidates.tsv"
    header = ["candidate", "fasta", "busco", "qv", "stats"]
    lines = ["\t".join(header)]
    lines += ["\t".join(row.get(key, "") for key in header) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def candidate_row(tmp_path, name, s=97.9, d=0.6, qv="45.2", n50=924431):
    return {
        "candidate": name,
        "fasta": fasta(tmp_path, name).name,
        "busco": summary(tmp_path, name, c=s + d, s=s, d=d).name,
        "qv": merqury(tmp_path, name, qv).name,
        "stats": gfastats(tmp_path, name, n50=n50).name,
    }


def run(*arguments: str):
    return CliRunner().invoke(cli, ["select", *arguments])


def test_busco_txt_json_and_compleasm(tmp_path: Path):
    busco = read_busco(summary(tmp_path, "flye"))
    assert busco == Busco(98.5, 97.9, 0.6, 0.3, 1.2, 2137)
    as_json = tmp_path / "short_summary.json"
    as_json.write_text(
        json.dumps(
            {
                "results": {
                    "one_line_summary": "C:90.0%[S:89.0%,D:1.0%],F:2.0%,M:8.0%,n:100"
                }
            }
        ),
        encoding="utf-8",
    )
    assert read_busco(as_json).single == 89.0
    compleasm = tmp_path / "summary.txt"
    compleasm.write_text(COMPLEASM, encoding="utf-8")
    assert read_busco(compleasm) == Busco(
        pytest.approx(98.51), 97.9, 0.61, 0.2, 1.29, 2137
    )
    nothing = tmp_path / "nothing.txt"
    nothing.write_text("no summary here\n", encoding="utf-8")
    with pytest.raises(ValueError, match="neither a BUSCO nor a compleasm"):
        read_busco(nothing)


def test_qv_and_gfastats(tmp_path: Path):
    assert read_qv(merqury(tmp_path, "a", "45.67")) == 45.67
    assert read_qv(merqury(tmp_path, "b", "inf")) == float("inf")
    empty = tmp_path / "empty.qv"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="no QV"):
        read_qv(empty)
    stats = read_gfastats(gfastats(tmp_path, "a", contigs=17, n50=924431))
    assert stats["contigs"] == 17
    assert stats["contig_n50"] == 924431
    assert stats["gc_percent"] == 38.15
    # gfastats' default, human-readable summary reads as well.
    plain = tmp_path / "plain.assembly_summary"
    plain.write_text("+++Assembly summary+++: \n# contigs: 3\nContig N50: n/a\n")
    assert read_gfastats(plain) == {"contigs": 3}


def test_the_rules_in_order(tmp_path: Path):
    def make(name, single=None, qv=None, n50=None):
        busco = Busco(single + 1, single, 1, 0, 0, 100) if single is not None else None
        stats = {"contig_n50": n50} if n50 is not None else {}
        return Candidate(name, tmp_path / name, busco, qv, stats)

    # BUSCO decides first, then QV, then N50; ties within the tolerances.
    assert (
        choose_assembly([make("a", 97, 40, 1), make("b", 98, 30, 1)]).kept.name == "b"
    )
    tie = choose_assembly([make("a", 97.8, 40, 1), make("b", 98, 30, 1)])
    assert (tie.kept.name, tie.reason) == ("a", "the highest QV")
    n50 = choose_assembly([make("a", 98, 40.5, 5), make("b", 98, 40, 9)])
    assert (n50.kept.name, n50.reason) == ("b", "the longest contig N50")
    # A number one candidate lacks leaves its rule out.
    lacking = choose_assembly([make("a", None, 40, 1), make("b", 98, 30, 1)])
    assert (lacking.kept.name, lacking.reason) == ("a", "the highest QV")
    # An infinite QV (no error found) beats every finite one.
    assert (
        choose_assembly([make("a", 98, float("inf")), make("b", 98, 60)]).kept.name
        == "a"
    )
    # A tie on every rule keeps the first.
    same = choose_assembly([make("a", 98, 40, 1), make("b", 98, 40, 1)])
    assert same.reason.startswith("a tie on every rule")
    assert choose_assembly([make("a")]).reason == "the only candidate"
    named = choose_assembly([make("a", 98), make("b", 10)], keep="b")
    assert (named.kept.name, named.reason) == ("b", "named by the samplesheet")
    with pytest.raises(ValueError, match="no candidate 'c'"):
        choose_assembly([make("a"), make("b")], keep="c")


def test_select_writes_the_kept_sequence_and_the_table(tmp_path: Path):
    rows = [
        candidate_row(tmp_path, "flye", s=97.9, qv="45.2", n50=900_000),
        candidate_row(tmp_path, "hifiasm", s=96.0, qv="50.1", n50=1_200_000),
    ]
    outdir = tmp_path / "out"
    result = run(str(table(tmp_path, rows)), "-p", "Scer", "-o", str(outdir))
    assert result.exit_code == 0, result.output
    with gzip.open(outdir / "Scer.assembly.fa.gz", "rt") as handle:
        assert handle.read().startswith(">Scer_flye_1")
    lines = (outdir / "Scer.assembly.selection.tsv").read_text().splitlines()
    assert lines[0].split("\t")[:3] == ["candidate", "kept", "busco_complete"]
    flye = dict(zip(lines[0].split("\t"), lines[1].split("\t"), strict=True))
    assert flye["kept"] == "yes"
    assert flye["busco_single"] == "97.9"
    assert flye["qv"] == "45.20"
    assert flye["contig_n50"] == "900000"
    assert flye["reason"] == "the most single-copy complete BUSCOs"
    assert lines[2].split("\t")[1] == "no"
    multiqc = (outdir / "Scer.assembly_selection_mqc.tsv").read_text()
    assert "# id: 'plexembler_assembly_selection'" in multiqc
    assert "Scer_hifiasm\tno\t96.0" in multiqc
    (log,) = outdir.glob("Scer.plexembler_*.log")
    assert "Kept flye for Scer" in log.read_text()


def test_keep_and_force(tmp_path: Path):
    rows = [candidate_row(tmp_path, "flye"), candidate_row(tmp_path, "hifiasm", s=50)]
    candidates = table(tmp_path, rows)
    outdir = tmp_path / "out"
    result = run(str(candidates), "-p", "Scer", "-o", str(outdir), "-k", "hifiasm")
    assert result.exit_code == 0, result.output
    with gzip.open(outdir / "Scer.assembly.fa.gz", "rt") as handle:
        assert ">Scer_hifiasm_1" in handle.read()
    again = run(str(candidates), "-p", "Scer", "-o", str(outdir))
    assert again.exit_code == 1
    assert "--force" in again.output
    missing = run(
        str(candidates), "-p", "Scer", "-o", str(outdir), "-k", "x", "--force"
    )
    assert missing.exit_code == 1
    assert "no candidate 'x'" in missing.output


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("name\tfasta\nflye\ta.fa\n", "no column candidate"),
        ("candidate\tfasta\n", "no candidate"),
        ("candidate\tfasta\nflye\ta.fa\nflye\tb.fa\n", "there twice"),
        ("candidate\tfasta\nflye\tNA\n", "has no sequence"),
    ],
)
def test_malformed_tables_fail(tmp_path: Path, text: str, message: str):
    path = tmp_path / "candidates.tsv"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        read_candidates(path)
    result = run(str(path), "-p", "Scer", "-o", str(tmp_path / "out"))
    assert result.exit_code == 1
    assert message in result.output


def test_missing_checks_are_empty_cells(tmp_path: Path):
    rows = [{"candidate": "flye", "fasta": fasta(tmp_path, "flye").name, "qv": "NA"}]
    (candidate,) = read_candidates(table(tmp_path, rows))
    assert (candidate.busco, candidate.qv, candidate.stats) == (None, None, {})


def test_copy_sequence_compresses_as_named(tmp_path: Path):
    plain = tmp_path / "a.fa"
    plain.write_text(">a\nAC\n", encoding="utf-8")
    copy_sequence(plain, tmp_path / "b.fa.gz")
    with gzip.open(tmp_path / "b.fa.gz", "rt") as handle:
        assert handle.read() == ">a\nAC\n"
    copy_sequence(tmp_path / "b.fa.gz", tmp_path / "c.fa")
    assert (tmp_path / "c.fa").read_text() == ">a\nAC\n"
    copy_sequence(plain, tmp_path / "d.fa")
    assert (tmp_path / "d.fa").read_text() == ">a\nAC\n"


def test_a_species_id_is_checked(tmp_path: Path):
    rows = [candidate_row(tmp_path, "flye")]
    result = run(str(table(tmp_path, rows)), "-p", "S_cer", "-o", str(tmp_path / "o"))
    assert result.exit_code == 2
