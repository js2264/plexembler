"""Each join of a species' scaffolds, scored by the Hi-C contacts across it.

A scaffolder joins two contigs when Hi-C pairs link their ends; a join it made
wrongly (a misjoin) has few pairs across it, fewer than the species' contact
decay says two neighbouring loci share. ``plexembler junctions`` scores each
join of an AGP (docs/plan.md, M3):

- the species' contact decay P(s), from the pairs whose mates lie in one
  contig, where no join can mislead it: the pairs at each distance (in
  logarithmic bins), over the pairs of positions that distance apart;
- for each join, the pairs across it, a mate within ``--window`` bases on each
  side, against those P(s) expects there;
- the join ``supported`` when their ratio is ``--min-ratio`` or more,
  ``unsupported`` below, ``thin`` when P(s) expects too few pairs to tell
  (``--min-expected``).

The pairs are the species' deduplicated pairs on its genome (4DN pairs, from
pairtools), and the AGP the genome's (``plexembler finalize``).
"""

from __future__ import annotations

import bisect
import gzip
import itertools
import logging
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO, cast

logger = logging.getLogger(__name__)

#: The decay's bins: ten a decade, from 10 bases.
BINS_PER_DECADE = 10
SMALLEST = 10


@dataclass
class Component:
    """A contig of a scaffold: where it lies in the scaffold (0-based, end
    excluded), and its name."""

    start: int
    end: int
    name: str


@dataclass
class Junction:
    """A join between two consecutive contigs of a scaffold."""

    scaffold: str
    left: Component
    right: Component
    observed: int = 0
    expected: float = 0.0
    verdict: str = ""

    @property
    def ratio(self) -> float:
        return self.observed / self.expected if self.expected else math.nan


@dataclass
class Layout:
    """A genome's scaffolds, as their contigs."""

    components: dict[str, list[Component]] = field(default_factory=dict)

    def starts(self, scaffold: str) -> list[int]:
        return [c.start for c in self.components[scaffold]]


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, "rt", encoding="ascii"))
    return cast(TextIO, path.open(encoding="ascii"))


def read_agp(path: Path) -> Layout:
    """An AGP's scaffolds and their contigs, gaps left out."""
    layout = Layout()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"{path}: an AGP line of {len(fields)} columns")
            if fields[4] in ("N", "U"):
                continue
            layout.components.setdefault(fields[0], []).append(
                Component(int(fields[1]) - 1, int(fields[2]), fields[5])
            )
    for components in layout.components.values():
        components.sort(key=lambda c: c.start)
    return layout


def junctions(layout: Layout) -> list[Junction]:
    """Every join of the layout, scaffold by scaffold, in order."""
    found = []
    for scaffold, components in layout.components.items():
        for left, right in itertools.pairwise(components):
            found.append(Junction(scaffold, left, right))
    return found


def bin_of(distance: int) -> int:
    """A distance's bin: 0 below ``SMALLEST``, then ten a decade."""
    if distance < SMALLEST:
        return 0
    return 1 + int(BINS_PER_DECADE * math.log10(distance / SMALLEST))


def bin_edges(index: int) -> tuple[int, int]:
    """The distances of a bin, its first and one past its last."""
    if index == 0:
        return 0, SMALLEST
    low = math.ceil(SMALLEST * 10 ** ((index - 1) / BINS_PER_DECADE))
    high = math.ceil(SMALLEST * 10 ** (index / BINS_PER_DECADE))
    return low, max(high, low + 1)


def position_pairs_within(length: int, low: int, high: int) -> int:
    """The pairs of positions of a sequence ``length`` long, ``low`` to
    ``high - 1`` apart: the sum of ``length - s`` over those distances ``s``."""
    high = min(high, length)
    if high <= low:
        return 0
    count = high - low
    return count * length - (low + high - 1) * count // 2


def position_pairs_across(left: int, right: int, distance: int) -> int:
    """The pairs of positions ``distance`` apart, one in the last ``left`` bases
    before a join and one in the first ``right`` after it."""
    if distance < 1 or distance > left + right - 1:
        return 0
    return min(distance, left, right, left + right - distance)


@dataclass
class Decay:
    """The species' contact decay: pairs, and pairs of positions, by bin."""

    pairs: dict[int, int] = field(default_factory=lambda: defaultdict(int))
    positions: dict[int, int] = field(default_factory=lambda: defaultdict(int))

    def add_contig(self, length: int, longest: int) -> None:
        for index in range(bin_of(longest) + 1):
            low, high = bin_edges(index)
            self.positions[index] += position_pairs_within(length, low, high)

    def probability(self, distance: int) -> float:
        index = bin_of(distance)
        positions = self.positions.get(index, 0)
        return self.pairs.get(index, 0) / positions if positions else 0.0


def expected_across(decay: Decay, left: int, right: int) -> float:
    """The pairs the decay expects across a join, ``left`` and ``right`` bases
    of window on each side."""
    expected = 0.0
    longest = left + right - 1
    index = 0
    while True:
        low, high = bin_edges(index)
        if low > longest:
            break
        probability = decay.probability(max(low, 1))
        if probability:
            expected += probability * sum(
                position_pairs_across(left, right, s)
                for s in range(max(low, 1), min(high, longest + 1))
            )
        index += 1
    return expected


def read_pairs(
    path: Path, layout: Layout, found: list[Junction], window: int
) -> tuple[Decay, int]:
    """The decay from the pairs within contigs; the pairs across each join.

    Returns the decay and the number of pairs read.
    """
    decay = Decay()
    across: dict[tuple[str, int], Junction] = {}
    for junction in found:
        index = layout.components[junction.scaffold].index(junction.left)
        across[(junction.scaffold, index)] = junction
    starts = {name: layout.starts(name) for name in layout.components}
    read = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 6)
            if len(fields) < 5:
                continue
            chrom1, chrom2 = fields[1], fields[3]
            if chrom1 != chrom2 or chrom1 not in starts:
                continue
            read += 1
            pos1, pos2 = int(fields[2]) - 1, int(fields[4]) - 1
            if pos2 < pos1:
                pos1, pos2 = pos2, pos1
            components = layout.components[chrom1]
            one = bisect.bisect_right(starts[chrom1], pos1) - 1
            two = bisect.bisect_right(starts[chrom1], pos2) - 1
            if one < 0 or pos1 >= components[one].end:
                continue  # in a gap
            if two < 0 or pos2 >= components[two].end:
                continue
            if one == two:
                decay.pairs[bin_of(pos2 - pos1)] += 1
            elif two == one + 1:
                junction = across[(chrom1, one)]
                if (
                    pos1 >= junction.left.end - window
                    and pos2 < junction.right.start + window
                ):
                    junction.observed += 1
    return decay, read


def score(
    agp: Path,
    pairs: Path,
    window: int = 20_000,
    min_ratio: float = 0.2,
    min_expected: float = 10.0,
) -> list[Junction]:
    """Each join of the AGP, its pairs across, those expected, and its verdict."""
    layout = read_agp(agp)
    found = junctions(layout)
    decay, read = read_pairs(pairs, layout, found, window)
    for components in layout.components.values():
        for component in components:
            length = component.end - component.start
            decay.add_contig(length, length - 1)
    for junction in found:
        left = min(window, junction.left.end - junction.left.start)
        right = min(window, junction.right.end - junction.right.start)
        junction.expected = expected_across(decay, left, right)
        if junction.expected < min_expected:
            junction.verdict = "thin"
        elif junction.ratio >= min_ratio:
            junction.verdict = "supported"
        else:
            junction.verdict = "unsupported"
    logger.info(
        "%d pairs within scaffolds; %d joins: %d supported, %d unsupported, %d thin",
        read,
        len(found),
        sum(j.verdict == "supported" for j in found),
        sum(j.verdict == "unsupported" for j in found),
        sum(j.verdict == "thin" for j in found),
    )
    return found


def write_outputs(found: list[Junction], outdir: Path, prefix: str) -> None:
    """``<prefix>.junctions.tsv``, a row per join, and its MultiQC table."""
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / f"{prefix}.junctions.tsv").open("w", encoding="ascii") as sink:
        sink.write(
            "scaffold\tposition\tleft_contig\tright_contig\tobserved\texpected\tratio\tverdict\n"
        )
        for j in found:
            ratio = "NA" if math.isnan(j.ratio) else f"{j.ratio:.3f}"
            sink.write(
                f"{j.scaffold}\t{j.left.end}\t{j.left.name}\t{j.right.name}\t"
                f"{j.observed}\t{j.expected:.1f}\t{ratio}\t{j.verdict}\n"
            )
    counts = {
        verdict: sum(j.verdict == verdict for j in found)
        for verdict in ("supported", "unsupported", "thin")
    }
    with (outdir / f"{prefix}.junctions_mqc.tsv").open("w", encoding="ascii") as sink:
        sink.write(
            "# id: plexembler_junctions\n"
            "# section_name: Scaffold joins\n"
            "# description: each join of the genome's scaffolds, scored by the Hi-C pairs "
            "across it against the species' contact decay (plexembler junctions)\n"
            "# plot_type: table\n"
            "Sample\tJoins\tSupported\tUnsupported\tThin\n"
            f"{prefix}\t{len(found)}\t{counts['supported']}\t"
            f"{counts['unsupported']}\t{counts['thin']}\n"
        )
