"""`plexembler junctions`: each join of a genome's scaffolds, scored by its Hi-C pairs."""

import math
import random
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.junctions import (
    Decay,
    Expectation,
    Lift,
    bin_edges,
    bin_of,
    position_pairs_within,
    read_agp,
    score,
)

LENGTH = 40_000
GAP = 100

#: One scaffold of three contigs: A and B truly neighbours, B and C not (a
#: misjoin); and an unplaced contig alone.
AGP = "\n".join(
    "\t".join(str(c) for c in row)
    for row in [
        ("Scer_chr01", 1, LENGTH, 1, "W", "Scer_a", 1, LENGTH, "+"),
        ("Scer_chr01", LENGTH + 1, LENGTH + GAP, 2, "N", GAP, "scaffold", "yes", "proximity_ligation"),
        ("Scer_chr01", LENGTH + GAP + 1, 2 * LENGTH + GAP, 3, "W", "Scer_b", 1, LENGTH, "+"),
        ("Scer_chr01", 2 * LENGTH + GAP + 1, 2 * LENGTH + 2 * GAP, 4, "N", GAP, "scaffold", "yes", "proximity_ligation"),
        ("Scer_chr01", 2 * LENGTH + 2 * GAP + 1, 3 * LENGTH + 2 * GAP, 5, "W", "Scer_c", 1, LENGTH, "-"),
        ("Scer_unplaced_01", 1, 5_000, 1, "W", "Scer_d", 1, 5_000, "+"),
    ]
)  # fmt: skip


def molecule_pairs(rng, length, n):
    """Pairs on one molecule ``length`` long, their distance log-uniform."""
    for _ in range(n):
        distance = int(10 ** rng.uniform(1, math.log10(length)))
        start = rng.randrange(0, length - distance)
        yield start, start + distance


def pairs_file(
    path: Path,
    seed: int = 1,
    n: int = 30_000,
    hidden: tuple[int, int] | None = None,
    extra: list[str] | None = None,
) -> Path:
    """The pairs of AGP's layout: A and B one molecule, C another. Without the
    pairs with a mate in ``hidden`` (1-based, end excluded), which no mate is
    placed in, as in a repeat; and with the ``extra`` lines."""
    rng = random.Random(seed)
    lines = [
        "## pairs format v1.0",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
    ]

    def to_genome(position):
        """A and B, one molecule: B after A's gap."""
        return position + 1 if position < LENGTH else position + GAP + 1

    for i, (x, y) in enumerate(molecule_pairs(rng, 2 * LENGTH, 2 * n)):
        one, two = to_genome(x), to_genome(y)
        if hidden and (hidden[0] <= one < hidden[1] or hidden[0] <= two < hidden[1]):
            continue
        lines.append(f"r{i}\tScer_chr01\t{one}\tScer_chr01\t{two}\t+\t-\tUU")
    offset = 2 * LENGTH + 2 * GAP + 1
    for i, (x, y) in enumerate(molecule_pairs(rng, LENGTH, n)):
        lines.append(
            f"c{i}\tScer_chr01\t{x + offset}\tScer_chr01\t{y + offset}\t+\t-\tUU"
        )
    # A pair between chromosomes, and one in a gap: neither counts.
    lines.append("t\tScer_chr01\t10\tScer_unplaced_01\t10\t+\t-\tUU")
    lines.append(f"g\tScer_chr01\t{LENGTH + 5}\tScer_chr01\t{LENGTH + 500}\t+\t-\tUU")
    lines += extra or []
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def on_contigs(pairs: Path, agp: Path, path: Path) -> Path:
    """The pairs of ``pairs_file``, each mate on its contig: where an aligner
    would have put it on the contigs the AGP lays out."""
    layout = read_agp(agp)

    def where(scaffold, position):
        for c in layout.components.get(scaffold, []):
            if c.start < position <= c.end:
                offset = position - 1 - c.start
                if c.orientation == "-":
                    return c.name, c.contig_end - offset
                return c.name, c.contig_start + offset + 1
        return None

    lines = []
    for line in pairs.read_text().splitlines():
        if line.startswith("#"):
            lines.append(line)
            continue
        fields = line.split("\t")
        one, two = where(fields[1], int(fields[2])), where(fields[3], int(fields[4]))
        if one and two:
            fields[1], fields[2] = one[0], str(one[1])
            fields[3], fields[4] = two[0], str(two[1])
            lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


@pytest.fixture
def agp(tmp_path: Path) -> Path:
    path = tmp_path / "Scer.agp"
    path.write_text("##agp-version\t2.1\n" + AGP + "\n", encoding="ascii")
    return path


def test_bins_and_position_pairs():
    assert bin_of(5) == 0
    assert bin_of(10) == 1
    assert bin_of(100) == 11
    low, high = bin_edges(11)
    assert bin_of(low) == 11 and bin_of(high) == 12
    # Every distance 1 to 4 apart in 5 bases: 4 + 3 + 2 + 1.
    assert position_pairs_within(5, 1, 5) == 10
    assert position_pairs_within(5, 5, 9) == 0


def test_expected_from_a_flat_decay():
    """A decay of one pair per pair of positions expects every pair across:
    each pair of bins (1 kb, the last one shorter) its pairs of positions, none
    for a bin whose pairs do not count."""
    decay = Decay()
    for index in range(60):
        decay.pairs[index] = 1
        decay.positions[index] = 1
    expectation = Expectation(decay, 6_000)
    assert expectation.across(30, 40, [1.0], [1.0]) == pytest.approx(30 * 40)
    assert expectation.across(2_500, 3_000, [1.0] * 3, [1.0] * 3) == pytest.approx(
        2_500 * 3_000
    )
    assert expectation.across(2_500, 3_000, [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]) == (
        pytest.approx(1_500 * 2_000)
    )


def test_a_join_and_a_misjoin(agp, tmp_path: Path):
    layout = read_agp(agp)
    assert [c.name for c in layout.components["Scer_chr01"]] == [
        "Scer_a",
        "Scer_b",
        "Scer_c",
    ]
    found = score(agp, pairs_file(tmp_path / "Scer.pairs"), window=10_000)
    assert [(j.left.name, j.right.name) for j in found] == [
        ("Scer_a", "Scer_b"),
        ("Scer_b", "Scer_c"),
    ]
    join, misjoin = found
    assert join.verdict == "supported"
    assert 0.7 < join.ratio < 1.3
    assert misjoin.observed == 0
    assert misjoin.verdict == "unsupported"


def test_the_lift(tmp_path: Path):
    """A contig's position on its scaffold: either way round, and on each piece
    of a contig the scaffolder broke in two."""
    path = tmp_path / "broken.agp"
    rows = [
        ("scaffold_1", 1, 100, 1, "W", "ctg1", 1, 100, "+"),
        ("scaffold_1", 101, 200, 2, "U", 100, "scaffold", "yes", "proximity_ligation"),
        ("scaffold_1", 201, 250, 3, "W", "ctg2", 51, 100, "-"),
        ("scaffold_2", 1, 50, 1, "W", "ctg2", 1, 50, "+"),
    ]  # fmt: skip
    path.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")
    lift = Lift(read_agp(path))
    assert lift("ctg1", 0) == ("scaffold_1", 0)
    assert lift("ctg1", 99) == ("scaffold_1", 99)
    assert lift("ctg1", 100) is None  # past the contig's part
    assert lift("ctg2", 0) == ("scaffold_2", 0)
    assert lift("ctg2", 49) == ("scaffold_2", 49)
    assert lift("ctg2", 50) == (
        "scaffold_1",
        249,
    )  # reversed: its part's first base last
    assert lift("ctg2", 99) == ("scaffold_1", 200)
    assert lift("ctg3", 10) is None


def test_the_joins_scored_from_pairs_on_the_contigs(agp, tmp_path: Path):
    """The pairs on the contigs, placed through the AGP, score the joins as
    the same pairs on the scaffolds do."""
    pairs = pairs_file(tmp_path / "Scer.pairs")
    contigs = on_contigs(pairs, agp, tmp_path / "Scer.contigs.pairs")
    assert "Scer_c\t" in contigs.read_text()
    on_scaffolds = score(agp, pairs, window=10_000)
    lifted = score(agp, contigs, window=10_000, on_contigs=True)
    assert [(j.observed, j.verdict) for j in lifted] == [
        (j.observed, j.verdict) for j in on_scaffolds
    ]
    # The same, but for the one mate of the pair in a gap, which only the
    # scaffolds' pairs hold (it counts in the genome's mean coverage).
    assert [j.expected for j in lifted] == pytest.approx(
        [j.expected for j in on_scaffolds], rel=1e-4
    )
    # Pairs on the scaffolds are not on the contigs: none is placed.
    assert {j.verdict for j in score(agp, pairs, on_contigs=True)} == {"thin"}


def test_a_join_between_repeats(agp, tmp_path: Path):
    """A true join whose contigs end in 6 kb where no mate is placed (a repeat
    whose copies all assembled): judged on the rest of the window, it is
    supported, as without the repeats."""
    hidden = (LENGTH - 6_000 + 1, LENGTH + GAP + 6_000 + 1)
    pairs = pairs_file(tmp_path / "Scer.pairs", hidden=hidden)
    join, misjoin = score(agp, pairs, window=10_000)
    assert join.verdict == "supported"
    assert 0.7 < join.ratio < 1.3
    assert join.left_seen == pytest.approx(0.4) and join.right_seen == pytest.approx(
        0.4
    )
    assert misjoin.verdict == "unsupported"


def test_a_misjoin_through_a_collapsed_repeat(agp, tmp_path: Path):
    """A misjoin whose right contig starts in a collapsed repeat: its 2 kb hold
    four times the mates, some of them across the misjoin (another copy's
    contacts). They do not count: the misjoin stays unsupported."""
    start = 2 * LENGTH + 2 * GAP + 1  # C's first base
    rng = random.Random(7)
    extra = [
        f"x{i}\tScer_chr01\t{start + rng.randrange(2_000)}\tScer_unplaced_01\t{1 + rng.randrange(5_000)}\t+\t-\tUU"
        for i in range(9_000)
    ]
    extra += [
        f"y{i}\tScer_chr01\t{start - GAP - 1 - rng.randrange(2_000)}\tScer_chr01\t{start + rng.randrange(2_000)}\t+\t-\tUU"
        for i in range(600)
    ]  # fmt: skip
    pairs = pairs_file(tmp_path / "Scer.pairs", extra=extra)
    _, misjoin = score(agp, pairs, window=10_000)
    assert misjoin.verdict == "unsupported"
    assert misjoin.observed < 20
    assert misjoin.right_seen == pytest.approx(0.8)


def test_too_few_pairs_to_judge(agp, tmp_path: Path):
    found = score(agp, pairs_file(tmp_path / "Scer.pairs", n=20), window=1_000)
    assert {j.verdict for j in found} == {"thin"}


def test_a_bad_agp(tmp_path: Path):
    bad = tmp_path / "bad.agp"
    bad.write_text("Scer_chr01\t1\t10\t1\tW\n", encoding="ascii")
    with pytest.raises(ValueError, match="5 columns"):
        read_agp(bad)


def test_the_command(agp, tmp_path: Path, plain):
    pairs = pairs_file(tmp_path / "Scer.pairs")
    outdir = tmp_path / "out"
    arguments = ["junctions", str(pairs), "-a", str(agp), "-p", "Scer", "-o", str(outdir), "-w", "10000"]  # fmt: skip
    result = CliRunner().invoke(cli, arguments)
    assert result.exit_code == 0, result.output
    table = (outdir / "Scer.junctions.tsv").read_text().splitlines()
    assert table[0].split("\t") == [
        "scaffold",
        "position",
        "left_contig",
        "right_contig",
        "observed",
        "expected",
        "ratio",
        "verdict",
        "left_seen",
        "right_seen",
    ]
    assert table[2].split("\t")[7:] == ["unsupported", "1.00", "1.00"]
    multiqc = (outdir / "Scer.junctions_mqc.tsv").read_text()
    assert "Scer\t2\t1\t1\t0" in multiqc
    log = next(outdir.glob("Scer.plexembler_*.log")).read_text()
    assert "the join of Scer_b and Scer_c" in log
    assert CliRunner().invoke(cli, arguments).exit_code == 1
    broken = tmp_path / "broken.agp"
    broken.write_text("Scer_chr01\t1\n", encoding="ascii")
    bad = CliRunner().invoke(
        cli, [*arguments[:3], str(broken), *arguments[4:], "--force"]
    )
    assert bad.exit_code == 1
    assert "columns" in plain(bad.output)
    contigs = on_contigs(pairs, agp, tmp_path / "Scer.contigs.pairs")
    lifted = CliRunner().invoke(
        cli, [arguments[0], str(contigs), *arguments[2:], "--contigs", "--force"]
    )
    assert lifted.exit_code == 0, lifted.output
    assert (outdir / "Scer.junctions.tsv").read_text().splitlines() == table
