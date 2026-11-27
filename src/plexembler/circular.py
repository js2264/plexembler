"""Circular contigs that an assembler wrote as copies of themselves.

An assembler can write a circular genome, a mitochondrion or a plasmid, as two
or more copies of it joined end to end: on the benchmark pool, Flye wrote
*N. glabratus*' 20,063-bp mitochondrion as one 40,128-bp contig, which QUAST
then read as a misassembly. It happens to high-copy circular genomes, and from
one run to the next. ``plexembler rename`` cuts such a contig to one copy.

A contig is ``k`` copies of a unit when each of its bases, but the last copy's,
is found again one unit on, and when it is a whole number of units long: the
copies of a circle run end to end, the last one closing on the first. A tandem
repeat of a linear chromosome seldom ends on a whole number of units; and a
contig that its assembler says is linear (Flye's ``assembly_info.txt``,
``circ.`` N) is not looked at. Flye's table is the only hint: other assemblers
give none, and the sequence decides alone.

The sequence is read through anchors, 21-mers of it. Those at its start,
looked for where a second copy would start, give the unit's length; then up to
a thousand anchors along it, each looked for one unit on (within as many bases
as the copies may differ by), give the share of it found again, and so the
copies' identity: an anchor holds 21 bases, so the share is about the
identity to the 21st power.
"""

from __future__ import annotations

import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO, cast

#: An anchor's length: a 21-mer of a genome is found once in it, but by chance.
ANCHOR = 21
#: The shortest unit looked for, in bases: shorter repeats are satellites.
MIN_UNIT = 1_000
#: The lowest identity between a contig's copies.
MIN_IDENTITY = 0.99
#: The anchors along a contig, at most, and the fewest that can judge it.
ANCHORS = 1_000
MIN_ANCHORS = 50
#: The anchors of a contig's start that look for its second copy.
START_ANCHORS = 8
#: How far an anchor's copy may lie from one unit on, in bases, at least.
MIN_SLACK = 20
#: The columns of Flye's table that say whether a contig is circular.
INFO_COLUMNS = ("seq_name", "circ.")
BASES = frozenset("ACGT")


@dataclass(frozen=True)
class Copies:
    """A contig found to be copies of a unit, end to end.

    ``unit_bp`` is where its second copy starts: the first copy's length.
    ``identity`` is the copies' identity, from the anchors found one unit on.
    """

    copies: int
    unit_bp: int
    identity: float


@dataclass(frozen=True)
class AssemblyInfo:
    """What the assembler says of a contig: whether it is circular, and its
    multiplicity (Flye's ``mult.``: its depth over the genome's)."""

    circular: bool
    multiplicity: int | None = None


def informative(anchor: str) -> bool:
    """Whether an anchor tells where it is: A, C, G and T alone, and no short
    repeat (eight different 4-mers at least, where (AT)n has two)."""
    return (
        len(anchor) == ANCHOR
        and set(anchor) <= BASES
        and len({anchor[i : i + 4] for i in range(ANCHOR - 3)}) >= 8
    )


def slack(unit: float, min_identity: float) -> int:
    """How many bases two copies of a unit may differ in length by, at the
    identity asked: as many as their differences, ``MIN_SLACK`` at least."""
    return max(MIN_SLACK, math.ceil((1 - min_identity) * unit))


def found_near(sequence: str, start: int, offset: int, within: int) -> int:
    """Where the anchor at ``start`` is found again, ``offset`` bases on, give
    or take ``within``; -1 where it is not."""
    near = start + offset
    anchor = sequence[start : start + ANCHOR]
    return sequence.find(anchor, max(start + 1, near - within), near + within + ANCHOR)


def second_copy(
    sequence: str, starts: list[int], expected: float, within: int
) -> int | None:
    """Where a second copy starts, ``expected`` bases on give or take
    ``within``: from the first of the anchors at ``starts`` found again there.
    None when none is."""
    for start in starts:
        found = found_near(sequence, start, round(expected), within)
        if found != -1:
            return found - start
    return None


def copy_identity(sequence: str, unit: int, within: int) -> float | None:
    """The identity between each copy of a sequence and the next, from the
    anchors found again one unit on, give or take ``within``; None when the
    sequence has too few anchors to tell."""
    span = len(sequence) - unit - ANCHOR
    tried = found = 0
    for start in range(0, span + 1, max(1, span // ANCHORS)):
        if not informative(sequence[start : start + ANCHOR]):
            continue
        tried += 1
        found += found_near(sequence, start, unit, within) != -1
    if tried < MIN_ANCHORS:
        return None
    return (found / tried) ** (1 / ANCHOR)


def short_period(
    sequence: str, start: int, below: int, min_identity: float
) -> int | None:
    """A period of under ``below`` bases that a sequence has throughout, from
    where the anchor at ``start`` is found again: a tandem repeat of a short
    unit, a satellite, not a circle. None when it has none."""
    anchor = sequence[start : start + ANCHOR]
    end = start + below + ANCHOR
    found = sequence.find(anchor, start + 1, end)
    while found != -1:
        period = found - start
        identity = copy_identity(sequence, period, slack(period, min_identity))
        if identity is not None and identity >= min_identity:
            return period
        found = sequence.find(anchor, found + 1, end)
    return None


def find_copies(
    sequence: str, min_unit: int = MIN_UNIT, min_identity: float = MIN_IDENTITY
) -> Copies | None:
    """The copies of a unit that a sequence is, end to end, if it is some.

    Parameters
    ----------
    sequence
        The contig's sequence.
    min_unit
        The shortest unit looked for, in bases.
    min_identity
        The lowest identity between its copies.

    Returns
    -------
    Copies or None
        Its copies, the most of them that fit (four copies of a unit are not
        two of a unit twice as long); None when it is not two or more copies
        of a unit of ``min_unit`` bases or more, whole and near-identical, or
        when it is a tandem repeat of a shorter unit throughout.
    """
    upper = sequence.upper()
    length = len(upper)
    if length < 2 * min_unit:
        return None
    starts = [
        start
        for start in range(0, max(ANCHOR, min_unit // 2), ANCHOR)
        if informative(upper[start : start + ANCHOR])
    ][:START_ANCHORS]
    if not starts:
        return None
    for copies in range(length // min_unit, 1, -1):
        expected = length / copies
        within = slack(expected, min_identity)
        unit = second_copy(upper, starts, expected, within)
        # The copies of a circle are whole: each as long as the first, give
        # or take their differences.
        if unit is None or abs(length - copies * unit) > (copies - 1) * within:
            continue
        identity = copy_identity(upper, unit, within)
        if identity is not None and identity >= min_identity:
            # Copies of a unit under min_unit tile it too, two at a time.
            if short_period(upper, starts[0], min_unit, min_identity):
                return None
            return Copies(copies, unit, identity)
    return None


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, "rt", encoding="utf-8"))
    return cast(TextIO, path.open(encoding="utf-8"))


def read_assembly_info(path: Path) -> dict[str, AssemblyInfo]:
    """Flye's ``assembly_info.txt``: whether each contig is circular, and its
    multiplicity, by its name.

    Raises
    ------
    ValueError
        When the table is empty, has no header naming its ``seq_name`` and
        ``circ.`` columns, or a row's ``circ.`` is neither Y nor N.
    """
    found: dict[str, AssemblyInfo] = {}
    header: list[str] = []
    with open_text(path) as handle:
        for number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if not header:
                header = [field.lstrip("#") for field in fields]
                if not line.startswith("#") or not set(INFO_COLUMNS) <= set(header):
                    raise ValueError(
                        f"line {number}: not Flye's assembly_info.txt, whose header "
                        f"names {' and '.join(INFO_COLUMNS)}"
                    )
                continue
            row = dict(zip(header, fields, strict=False))
            circular = row.get("circ.")
            if circular not in ("Y", "N"):
                raise ValueError(f"line {number}: circ. is {circular!r}, not Y or N")
            multiplicity = row.get("mult.", "")
            found[row["seq_name"]] = AssemblyInfo(
                circular == "Y", int(multiplicity) if multiplicity.isdigit() else None
            )
    if not header:
        raise ValueError("empty: not Flye's assembly_info.txt")
    return found
