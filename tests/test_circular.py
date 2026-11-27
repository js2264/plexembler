"""Circular contigs that an assembler wrote as copies of themselves."""

import gzip
import random
from pathlib import Path

import pytest

from plexembler.circular import (
    AssemblyInfo,
    Copies,
    find_copies,
    informative,
    read_assembly_info,
)

#: Flye 2.9's assembly_info.txt; 2.8's has no alt_group.
FLYE_HEADER = "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"


def mutated(sequence: str, seed: int, substitutions: int, indels: int) -> str:
    """The sequence with differences, none in its first kilobase."""
    rng = random.Random(seed)
    bases = list(sequence)
    for _ in range(substitutions):
        at = rng.randrange(1000, len(bases))
        bases[at] = "A" if bases[at] != "A" else "C"
    for _ in range(indels):
        at = rng.randrange(1000, len(bases))
        if rng.random() < 0.5:
            del bases[at]
        else:
            bases.insert(at, "G")
    return "".join(bases)


def test_two_copies(sequence):
    unit = sequence(20_000, seed=1)
    assert find_copies(unit + unit) == Copies(2, 20_000, 1.0)


def test_copies_that_differ_a_little(sequence):
    """Each copy of a circle polished apart: a few differences, which leave
    them copies; the first is kept, whole."""
    unit = sequence(20_063, seed=2, gc=0.17)
    found = find_copies(unit + mutated(unit, 3, substitutions=10, indels=4))
    assert found is not None
    assert (found.copies, found.unit_bp) == (2, 20_063)
    assert 0.99 <= found.identity < 1


def test_the_most_copies_that_fit(sequence):
    """Four copies of a unit are not two of a unit twice as long."""
    unit = sequence(5_000, seed=4)
    assert find_copies(unit * 4) == Copies(4, 5_000, 1.0)
    assert find_copies(unit * 3) == Copies(3, 5_000, 1.0)


def test_a_repeat_at_its_start_is_no_period(sequence):
    """Its first anchor found again 300 bases on, but nothing else: no
    satellite, and its copies cut all the same."""
    unit = sequence(20_000, seed=15)
    unit = unit[:300] + unit[:21] + unit[321:]
    assert find_copies(unit + unit) == Copies(2, 20_000, 1.0)


def test_soft_masked_bases(sequence):
    unit = sequence(3_000, seed=5)
    unit = unit[:1000] + unit[1000:2000].lower() + unit[2000:]
    assert find_copies(unit + unit) == Copies(2, 3_000, 1.0)


@pytest.mark.parametrize(
    "build",
    [
        pytest.param(lambda s: s(20_000, seed=6), id="one copy"),
        pytest.param(lambda s: s(1_500, seed=7), id="shorter than two units"),
        pytest.param(
            lambda s: (u := s(10_000, seed=8)) + u + u[:5_000], id="2.5 copies"
        ),
        pytest.param(
            lambda s: (u := s(10_000, seed=9)) + mutated(u, 10, 300, 200),
            id="95 % alike",
        ),
        pytest.param(lambda s: s(800, seed=10) * 4, id="a unit under 1 kb"),
        pytest.param(lambda _: "AT" * 5_000, id="a short repeat"),
        pytest.param(lambda _: "N" * 10_000, id="gaps"),
        # Informative at its start, then (AT)n: too few anchors to judge.
        pytest.param(
            lambda s: (u := s(40, seed=11) + "AT" * 980) + u, id="too few anchors"
        ),
    ],
)
def test_no_copies(build, sequence):
    assert find_copies(build(sequence)) is None


def test_the_thresholds(sequence):
    unit = sequence(10_000, seed=12)
    diverged = unit + mutated(unit, 13, 300, 200)
    assert find_copies(diverged, min_identity=0.9) is not None
    short = sequence(800, seed=14) * 4
    assert find_copies(short, min_unit=500) == Copies(4, 800, 1.0)


def test_informative():
    assert informative("ACGTTGCAAGCTTCGATCCGA")
    assert not informative("ACGT")  # too short
    assert not informative("ACGTTGCAAGCNTCGATCCGA")  # a gap
    assert not informative("ATATATATATATATATATATA")  # a short repeat
    # A-and-T alone, but no short repeat: it tells where it is.
    assert informative("AATTTATAAATTATTTAAATA")


def test_flye_assembly_info(tmp_path: Path):
    path = tmp_path / "Cgla.flye.assembly_info.txt"
    path.write_text(
        FLYE_HEADER
        + "contig_1\t40128\t90\tY\tY\t2\t*\t*,1,*\n"
        + "contig_2\t230000\t40\tN\tN\t1\t*\t2,-3\n\n",
        encoding="ascii",
    )
    assert read_assembly_info(path) == {
        "contig_1": AssemblyInfo(circular=True, multiplicity=2),
        "contig_2": AssemblyInfo(circular=False, multiplicity=1),
    }


def test_flye_28_and_gzip(tmp_path: Path):
    """Flye 2.8 wrote no alt_group; a gzipped table is read as well."""
    path = tmp_path / "assembly_info.txt.gz"
    path.write_bytes(
        gzip.compress(
            b"#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\tgraph_path\n"
            b"contig_1\t6300\t352\tY\tY\t*\t1\n"
        )
    )
    assert read_assembly_info(path) == {"contig_1": AssemblyInfo(True, None)}


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("", "empty"),
        ("contig_1\n", "not Flye's assembly_info.txt"),
        ("#seq_name\tlength\n", "not Flye's assembly_info.txt"),
        (FLYE_HEADER + "contig_1\t100\t9\tX\tN\t1\t*\t1\n", "circ. is 'X'"),
        (FLYE_HEADER + "contig_1\t100\n", "circ. is None"),
    ],
)
def test_a_malformed_assembly_info(tmp_path: Path, text: str, message: str):
    path = tmp_path / "assembly_info.txt"
    path.write_text(text, encoding="ascii")
    with pytest.raises(ValueError, match=message):
        read_assembly_info(path)
