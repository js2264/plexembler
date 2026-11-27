"""`plexembler split`: a pooled library's pairs, given to the species of their mates."""

import gzip
import io
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.split import (
    Tally,
    epsilon,
    held_out,
    parse_record,
    read_sam,
    species_of,
    split_pairs,
    summary,
)

HEADER = "\n".join(
    [
        "@HD\tVN:1.6\tSO:unsorted",
        "@SQ\tSN:Scer_ctg1\tLN:10000",
        "@SQ\tSN:Scer_ctg2\tLN:10000",
        "@SQ\tSN:Spar_ctg1\tLN:10000",
        "@SQ\tSN:Cgla_ctg1\tLN:10000",
        "@PG\tID:bwa-mem2",
    ]
)
SEQ = "ACGTACGTAC" * 10
QUAL = "I" * 100


def record(name, flag, contig, pos, mapq=60, cigar="100M", nm=0, xa=None, seq=SEQ):
    tags = [f"NM:i:{nm}"]
    if xa:
        tags.append(f"XA:Z:{xa}")
    fields = [
        name,
        str(flag),
        contig,
        str(pos),
        str(mapq),
        cigar,
        "=",
        "1",
        "0",
        seq,
        QUAL,
    ]
    return "\t".join(fields + tags)


def pair(name, first, second):
    """Two mates: each (contig, pos, options) as record() takes them."""
    one = record(name, 65 | first.pop("flag", 0), *first.pop("at"), **first)
    two = record(name, 129 | second.pop("flag", 0), *second.pop("at"), **second)
    return [one, two]


def sam(lines):
    return io.StringIO(HEADER + "\n" + "\n".join(lines) + "\n")


def names(path):
    with gzip.open(path, "rt") as handle:
        return [line[1:].strip() for i, line in enumerate(handle) if i % 4 == 0]


def free(prefix: str) -> str:
    """A pair's name that the hash does not hold out."""
    return next(
        f"{prefix}{i}" for i in range(1000) if not held_out(f"{prefix}{i}", 0.1)
    )


#: A name that the hash keeps apart, and one it does not.
HELD = next(f"held{i}" for i in range(1000) if held_out(f"held{i}", 0.1))
KEPT = free("kept")
CGLA = free("cgla")


def wgs_lines():
    return [
        *pair(KEPT, {"at": ("Scer_ctg1", 100)}, {"at": ("Scer_ctg1", 400), "flag": 16}),
        *pair(HELD, {"at": ("Scer_ctg1", 1000)}, {"at": ("Scer_ctg1", 1300)}),
        # A mate below the MAPQ, whose alternative hit is in Spar.
        *pair(
            "amb",
            {"at": ("Scer_ctg1", 2000), "mapq": 3, "xa": "Spar_ctg1,+2000,100M,1;"},
            {"at": ("Scer_ctg1", 2300)},
        ),
        # Ambiguous between Scer and Spar, its mate in Cgla: three species.
        *pair(
            "amb3",
            {"at": ("Scer_ctg1", 3000), "mapq": 3, "xa": "Spar_ctg1,+3000,100M,1;"},
            {"at": ("Cgla_ctg1", 3300)},
        ),
        # A repeat of Scer: an ambiguous pair within one species.
        *pair("rep", {"at": ("Scer_ctg2", 10), "mapq": 0}, {"at": ("Scer_ctg2", 300)}),
        # 10 mismatches in 100 columns: a read of another species.
        *pair("div", {"at": ("Spar_ctg1", 10), "nm": 10}, {"at": ("Spar_ctg1", 300)}),
        *pair("btw", {"at": ("Scer_ctg1", 10)}, {"at": ("Cgla_ctg1", 300)}),
        *pair(
            "unm",
            {"at": ("*", 0), "flag": 4, "mapq": 0, "cigar": "*"},
            {"at": ("Cgla_ctg1", 1)},
        ),
        *pair(CGLA, {"at": ("Cgla_ctg1", 500)}, {"at": ("Cgla_ctg1", 800), "flag": 16}),
    ]


def test_species_and_records():
    assert species_of("Scer_ctg_1") == "Scer"
    mate = parse_record(record("r", 81, "Scer_ctg1", 100, cigar="5S90M2D5M", nm=4))
    assert mate.reverse and mate.mapped
    assert mate.identity() == pytest.approx(1 - 4 / 97)
    # The 5' end of a reversed read is its alignment's last base.
    assert mate.five_prime() == 100 + 97 - 1
    assert (
        mate.fastq("r").splitlines()[1]
        == SEQ.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    )
    assert parse_record(record("r", 4, "*", 0, cigar="*")).identity() == 0.0


def test_a_wgs_library(tmp_path: Path):
    species, tally = split_pairs(sam(wgs_lines()), tmp_path, "wgs", "wgs")
    assert species == ["Scer", "Spar", "Cgla"]
    assert tally.pairs == 9
    assert tally.kept == {"Scer": 1, "Cgla": 1}
    assert tally.holdout == {"Scer": 1}
    assert tally.ambiguous_pairs == 3
    assert tally.ambiguous == {
        ("Scer", "Spar"): 2,
        ("Scer", "Scer"): 1,
        ("Cgla", "Scer"): 1,
        ("Cgla", "Spar"): 1,
    }
    assert (tally.divergent, tally.unmapped) == (1, 1)
    assert tally.between == {("Cgla", "Scer"): 1}
    assert names(tmp_path / "Scer.wgs_R1.fastq.gz") == [KEPT]
    assert names(tmp_path / "Scer.wgs.holdout_R2.fastq.gz") == [HELD]
    # The second mate aligned reversed: written as it was sequenced.
    with gzip.open(tmp_path / "Cgla.wgs_R2.fastq.gz", "rt") as handle:
        assert (
            handle.read().splitlines()[1]
            == SEQ.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        )
    numbers = summary(species, tally, "wgs")
    assert numbers["species"]["Scer"]["share"] == 0.5
    assert numbers["species"]["Spar"]["pairs"] == 0
    assert "epsilon" not in numbers


def test_a_hic_library_its_duplicates_and_epsilon(tmp_path: Path):
    lines = [
        *pair("a", {"at": ("Scer_ctg1", 100)}, {"at": ("Scer_ctg2", 900)}),
        # The same pair again, its mates swapped: a duplicate.
        *pair("b", {"at": ("Scer_ctg2", 900)}, {"at": ("Scer_ctg1", 100)}),
        *pair("c", {"at": ("Spar_ctg1", 100)}, {"at": ("Spar_ctg1", 5000)}),
        *pair("d", {"at": ("Cgla_ctg1", 100)}, {"at": ("Cgla_ctg1", 3000)}),
        *pair("e", {"at": ("Scer_ctg1", 4000)}, {"at": ("Cgla_ctg1", 3000)}),
    ]
    species, tally = split_pairs(sam(lines), tmp_path, "hic", "hic")
    assert tally.kept == {"Scer": 1, "Spar": 1, "Cgla": 1}
    assert tally.duplicates == {"Scer": 1}
    # Mates: Scer 5, Spar 2, Cgla 3; 1 pair between species of 5.
    expected = (1 / 5) / (1 - (0.5**2 + 0.2**2 + 0.3**2))
    assert epsilon(tally) == pytest.approx(expected)
    # The unique pairs between two contigs, the duplicate left out.
    assert tally.contacts == {
        ("Scer_ctg1", "Scer_ctg2"): 1,
        ("Cgla_ctg1", "Scer_ctg1"): 1,
    }
    numbers = summary(species, tally, "hic")
    assert numbers["between_pairs"] == [
        {"species_a": "Cgla", "species_b": "Scer", "pairs": 1}
    ]
    kept, _ = split_pairs(sam(lines), tmp_path / "all", "hic", "hic", dedup=False)
    assert (
        split_pairs(sam(lines), tmp_path / "all2", "hic", "hic", dedup=False)[1].kept[
            "Scer"
        ]
        == 2
    )
    assert kept == species


def test_epsilon_needs_pairs():
    assert epsilon(Tally()) is None


@pytest.mark.parametrize(
    ("lines", "message"),
    [
        (
            pair("a", {"at": ("Scer_ctg1", 1)}, {"at": ("Scer_ctg1", 9)})[:1],
            "without its mate",
        ),
        (
            [
                record("a", 65, "Scer_ctg1", 1),
                record("b", 129, "Scer_ctg1", 9),
            ],
            "not together",
        ),
    ],
)
def test_malformed_alignments(tmp_path: Path, lines, message):
    _, pairs = read_sam(sam(lines))
    with pytest.raises(ValueError, match=message):
        list(pairs)


def test_supplementary_records_are_skipped(tmp_path: Path):
    lines = pair("a", {"at": ("Scer_ctg1", 100)}, {"at": ("Scer_ctg1", 400)})
    lines.insert(1, record("a", 65 | 2048, "Spar_ctg1", 1, cigar="50H50M"))
    _, tally = split_pairs(sam(lines), tmp_path, "hic", "hic")
    assert tally.kept == {"Scer": 1}


def test_the_command(tmp_path: Path, plain):
    path = tmp_path / "wgs.sam"
    path.write_text(HEADER + "\n" + "\n".join(wgs_lines()) + "\n", encoding="utf-8")
    outdir = tmp_path / "out"
    arguments = [
        "split",
        str(path),
        "-p",
        "poolA",
        "-o",
        str(outdir),
        "-l",
        "wgs",
        "-t",
        "wgs",
    ]
    result = CliRunner().invoke(cli, arguments)
    assert result.exit_code == 0, result.output
    numbers = json.loads((outdir / "poolA.wgs.split.json").read_text())
    assert numbers["kept"] == 2
    assert (outdir / "poolA.wgs.split.tsv").read_text().splitlines()[
        1
    ] == "Scer\t1\t0.5000\t1\t0"
    assert (
        "# id: 'plexembler_split'" in (outdir / "poolA.wgs.split_mqc.tsv").read_text()
    )
    (log,) = outdir.glob("poolA.plexembler_*.log")
    assert "Scer: 1 pairs (50.0 %)" in log.read_text()
    again = CliRunner().invoke(cli, arguments)
    assert again.exit_code == 1
    assert "--force" in plain(again.output)
    # From stdin, and gzipped.
    zipped = tmp_path / "wgs.sam.gz"
    zipped.write_bytes(gzip.compress(path.read_bytes()))
    from_gzip = ["split", str(zipped), *arguments[2:], "--force"]
    assert CliRunner().invoke(cli, from_gzip).exit_code == 0
    piped = CliRunner().invoke(
        cli,
        [
            "split",
            "-",
            "-p",
            "poolA",
            "-o",
            str(tmp_path / "piped"),
            "-l",
            "hic",
            "-t",
            "hic",
        ],
        input=path.read_text(),
    )
    assert piped.exit_code == 0, piped.output
    assert (
        "Random ligation rate" in next((tmp_path / "piped").glob("*.log")).read_text()
    )
    broken = tmp_path / "broken.sam"
    broken.write_text(
        HEADER + "\n" + record("a", 65, "Scer_ctg1", 1) + "\n", encoding="utf-8"
    )
    failed = CliRunner().invoke(
        cli,
        [
            "split",
            str(broken),
            "-p",
            "poolA",
            "-o",
            str(tmp_path / "b"),
            "-l",
            "x",
            "-t",
            "wgs",
        ],
    )
    assert failed.exit_code == 1
    assert "without its mate" in plain(failed.output)


def test_an_rnaseq_library(tmp_path: Path):
    """STAR's alignments: each species' pairs written as SAM of its own
    sequences; a multimapper's secondary records say between which species."""
    lines = [
        *pair("r1", {"at": ("Scer_ctg1", 100), "cigar": "40M500N60M"}, {"at": ("Scer_ctg1", 900)}),
        # Unique in Cgla.
        *pair("r2", {"at": ("Cgla_ctg1", 100)}, {"at": ("Cgla_ctg1", 400)}),
        # A multimapper: its primary in Scer, its secondary in Spar.
        *pair("r3", {"at": ("Scer_ctg2", 100), "mapq": 3}, {"at": ("Scer_ctg2", 400), "mapq": 3}),
        *pair(
            "r3",
            {"at": ("Spar_ctg1", 100), "mapq": 3, "flag": 256},
            {"at": ("Spar_ctg1", 400), "mapq": 3, "flag": 256},
        ),
        # The same pair twice: RNA-seq keeps it.
        *pair("r4", {"at": ("Scer_ctg1", 100), "cigar": "40M500N60M"}, {"at": ("Scer_ctg1", 900)}),
    ]  # fmt: skip
    _, tally = split_pairs(sam(lines), tmp_path, "rna", "rnaseq", holdout=0.5)
    assert tally.kept == {"Scer": 2, "Cgla": 1}
    assert tally.holdout == {}
    assert tally.ambiguous == {("Scer", "Spar"): 1}
    with gzip.open(tmp_path / "Scer.rna.sam.gz", "rt") as handle:
        written = handle.read().splitlines()
    header = [line for line in written if line.startswith("@")]
    assert [line.split("\t")[1] for line in header if line.startswith("@SQ")] == [
        "SN:Scer_ctg1",
        "SN:Scer_ctg2",
    ]
    records = [line.split("\t") for line in written if not line.startswith("@")]
    assert [(r[0], r[5]) for r in records][:2] == [("r1", "40M500N60M"), ("r1", "100M")]
    assert len(records) == 4
    assert not list(tmp_path.glob("*.fastq.gz"))


def test_a_read_with_three_primary_records(tmp_path: Path):
    lines = pair("a", {"at": ("Scer_ctg1", 1)}, {"at": ("Scer_ctg1", 9)})
    lines.append(lines[0])
    _, pairs = read_sam(sam(lines))
    with pytest.raises(ValueError, match="more than two"):
        list(pairs)
