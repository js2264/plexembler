"""A species' genome against another: a dotplot, and how the two align.

``plexembler dotplot`` reads minimap2's alignments (PAF) of a species' genome,
the query, to another genome, the reference: its closest known genome
(``plexembler closest``). It draws them, the reference's sequences along x in
their order, the species' along y, each placed beside the reference sequence it
aligns to most and turned round when most of it aligns to that one's other
strand: the chromosomes the two genomes share make a diagonal, and an
inversion, a translocation or a fusion stands out of it. Alignments on the same
strand are blue, on the other vermilion (the plot's strands, once the turned
sequences are turned).

It writes the dotplot as SVG; the same as a plotly figure, which the species'
report draws (``plexembler report``); and a table of every sequence of both
genomes: how much of it aligns, how well, and to which sequence of the other
genome most. The table counts every primary alignment, the plot draws every
alignment of ``--min-length`` bases or more, secondary ones included.
"""

from __future__ import annotations

import csv
import json
import logging
import re
from collections import defaultdict
from dataclasses import dataclass
from html import escape
from pathlib import Path

from .rename import open_text
from .screen import merged_length

logger = logging.getLogger(__name__)

#: The colours of the alignments on the same strand and on the other: Okabe
#: and Ito's blue and vermilion, which colour-blind eyes tell apart.
SAME = "#0072B2"
OTHER = "#D55E00"
#: The shortest alignment drawn, in bases, by default.
MIN_LENGTH = 1000
#: The plot's side, and its margins, in the SVG's pixels.
_PLOT = 720
_LEFT, _RIGHT, _TOP, _BOTTOM = 130, 30, 84, 150
#: The narrowest sequence, in pixels, whose name is written on its axis.
_MIN_LABEL = 9.0

_CHROMOSOME = re.compile(r"\bchromosome[\s:]+([^\s,;]+)", re.IGNORECASE)
_PLASMID = re.compile(r"\bplasmid[\s:]+([^\s,;]+)", re.IGNORECASE)
_MITOCHONDRION = re.compile(r"\bmitochondri", re.IGNORECASE)
_PLASTID = re.compile(r"\b(chloroplast|plastid)\b", re.IGNORECASE)
_PIECE = re.compile(r"\b(unlocalized|unplaced|scaffold|contig)\b", re.IGNORECASE)


@dataclass(frozen=True)
class Alignment:
    """A line of minimap2's PAF: a stretch of the query on one of the target."""

    query: str
    query_length: int
    query_start: int
    query_end: int
    reverse: bool
    target: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    block: int
    primary: bool


@dataclass
class Sequence:
    """A sequence on an axis: its name, length, label and place."""

    name: str
    length: int
    label: str
    offset: int = 0
    flipped: bool = False


@dataclass(frozen=True)
class Segment:
    """An alignment as drawn, in the axes' bases, and where it starts and ends
    in its two sequences."""

    x0: int
    y0: int
    x1: int
    y1: int
    reverse: bool
    start: tuple[int, int]
    end: tuple[int, int]
    query: Sequence
    target: Sequence


@dataclass(frozen=True)
class Row:
    """A sequence of either genome, and how it aligns to the other."""

    genome: str
    name: str
    label: str
    length: int
    aligned: int
    identity: float | None
    best_match: str | None
    share: float | None
    strand: str | None


def read_paf(path: Path) -> list[Alignment]:
    """minimap2's alignments, as they come.

    Raises
    ------
    ValueError
        When a line is not PAF.
    """
    alignments: list[Alignment] = []
    with open_text(path) as handle:
        for number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                if line.strip():
                    raise ValueError(f"{path}, line {number}: not a PAF line")
                continue
            try:
                alignments.append(
                    Alignment(
                        query=fields[0],
                        query_length=int(fields[1]),
                        query_start=int(fields[2]),
                        query_end=int(fields[3]),
                        reverse=fields[4] == "-",
                        target=fields[5],
                        target_length=int(fields[6]),
                        target_start=int(fields[7]),
                        target_end=int(fields[8]),
                        matches=int(fields[9]),
                        block=int(fields[10]),
                        primary="tp:A:S" not in fields[12:],
                    )
                )
            except ValueError as error:
                raise ValueError(f"{path}, line {number}: {error}") from error
    return alignments


def read_headers(path: Path) -> list[tuple[str, int, str]]:
    """Each sequence of a FASTA file: its name, length and the rest of its header.

    Raises
    ------
    ValueError
        When a header has no name, or a sequence comes before the first header.
    """
    found: list[tuple[str, int, str]] = []
    name: str | None = None
    description, length = "", 0
    with open_text(path) as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if name is not None:
                    found.append((name, length, description))
                words = line[1:].strip().split(maxsplit=1)
                if not words:
                    raise ValueError(f"{path}, line {number}: a header without a name")
                name, length = words[0], 0
                description = words[1] if len(words) > 1 else ""
            elif line.strip():
                if name is None:
                    raise ValueError(
                        f"{path}, line {number}: a sequence before the first header"
                    )
                length += len(line.strip())
    if name is not None:
        found.append((name, length, description))
    return found


def label_of(name: str, description: str) -> str:
    """A sequence's short name on an axis, from its FASTA header.

    NCBI's ``Saccharomyces cerevisiae S288C chromosome IV, complete sequence``
    is ``chrIV``, a mitochondrion ``chrM``, a plastid ``chrC``, a plasmid
    ``plasmid <name>``; a scaffold or a contig keeps its name.
    """
    if _PIECE.search(description):
        return name
    if found := _CHROMOSOME.search(description):
        return f"chr{found.group(1)}"
    if _MITOCHONDRION.search(description):
        return "chrM"
    if _PLASTID.search(description):
        return "chrC"
    if found := _PLASMID.search(description):
        return f"plasmid {found.group(1)}"
    return name


def _unique(sequences: list[Sequence]) -> list[Sequence]:
    """The sequences, a label two of them share replaced by their names."""
    counts: dict[str, int] = defaultdict(int)
    for sequence in sequences:
        counts[sequence.label] += 1
    for sequence in sequences:
        if counts[sequence.label] > 1:
            sequence.label = sequence.name
    return sequences


def species_sequences(
    headers: list[tuple[str, int, str]], prefix: str
) -> list[Sequence]:
    """The species' sequences, each labelled by its name without ``<species>_``."""
    start = f"{prefix}_"
    return _unique(
        [
            Sequence(name, length, name.removeprefix(start) or name)
            for name, length, _description in headers
        ]
    )


def reference_sequences(headers: list[tuple[str, int, str]]) -> list[Sequence]:
    """The reference's sequences, each labelled from its header."""
    return _unique(
        [
            Sequence(name, length, label_of(name, description))
            for name, length, description in headers
        ]
    )


def from_alignments(
    alignments: list[Alignment], side: str
) -> list[tuple[str, int, str]]:
    """The sequences that alignments name, longest first: for a genome given
    without its FASTA file."""
    lengths: dict[str, int] = {}
    for alignment in alignments:
        if side == "query":
            lengths[alignment.query] = alignment.query_length
        else:
            lengths[alignment.target] = alignment.target_length
    return [
        (name, length, "")
        for name, length in sorted(
            lengths.items(), key=lambda item: (-item[1], item[0])
        )
    ]


def lay_out(
    queries: list[Sequence], targets: list[Sequence], alignments: list[Alignment]
) -> list[Sequence]:
    """The queries in their order along y, each turned round when most of it
    aligns to the other strand of its target; every sequence's offset set.

    A query's target is the one its primary alignments cover most of; queries
    follow their targets' order, then the middle of their alignments on them.
    The queries without an alignment come last, in their order.
    """
    index = {target.name: number for number, target in enumerate(targets)}
    bases: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    position: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    strand: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for alignment in alignments:
        if not alignment.primary:
            continue
        aligned = alignment.query_end - alignment.query_start
        middle = (alignment.target_start + alignment.target_end) / 2
        bases[alignment.query][alignment.target] += aligned
        position[alignment.query][alignment.target] += aligned * middle
        strand[alignment.query][alignment.target] += (
            -aligned if alignment.reverse else aligned
        )

    def home(name: str) -> str:
        found = bases[name]
        return max(found, key=lambda target: (found[target], -index[target]))

    placed = []
    for number, query in enumerate(queries):
        if bases.get(query.name):
            target = home(query.name)
            query.flipped = strand[query.name][target] < 0
            where = position[query.name][target] / bases[query.name][target]
            placed.append(((index[target], where, number), query))
        else:
            placed.append(((len(targets), 0.0, number), query))
    ordered = [query for _key, query in sorted(placed, key=lambda item: item[0])]
    for sequences in (ordered, targets):
        offset = 0
        for sequence in sequences:
            sequence.offset = offset
            offset += sequence.length
    return ordered


def segments(
    alignments: list[Alignment],
    queries: dict[str, Sequence],
    targets: dict[str, Sequence],
    min_length: int,
) -> list[Segment]:
    """The alignments drawn: those of min_length bases or more, in the axes'
    bases, the turned queries turned."""
    drawn = []
    for alignment in alignments:
        if alignment.block < min_length:
            continue
        query, target = queries[alignment.query], targets[alignment.target]
        start, end = alignment.query_start, alignment.query_end
        reverse = alignment.reverse
        if query.flipped:
            start, end = query.length - end, query.length - start
            reverse = not reverse
        x0 = target.offset + alignment.target_start
        x1 = target.offset + alignment.target_end
        y0, y1 = query.offset + start, query.offset + end
        if reverse:
            y0, y1 = y1, y0
        # Where the alignment starts and ends in its sequences, as they are.
        first, last = (
            (alignment.query_end, alignment.query_start)
            if alignment.reverse
            else (alignment.query_start, alignment.query_end)
        )
        drawn.append(
            Segment(
                x0,
                y0,
                x1,
                y1,
                reverse,
                (alignment.target_start, first),
                (alignment.target_end, last),
                query,
                target,
            )
        )
    return drawn


def alignment_table(
    queries: list[Sequence], targets: list[Sequence], alignments: list[Alignment]
) -> list[Row]:
    """Every sequence of both genomes, and how it aligns to the other: its
    primary alignments' bases (each base once), identity, and the sequence of
    the other genome it aligns to most, with that one's share and strand."""
    intervals: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    matches: dict[tuple[str, str], int] = defaultdict(int)
    blocks: dict[tuple[str, str], int] = defaultdict(int)
    shared: dict[tuple[str, str], dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    strands: dict[tuple[str, str], dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    for alignment in alignments:
        if not alignment.primary:
            continue
        sides = (
            ("species", alignment.query, alignment.query_start, alignment.query_end, alignment.target),
            ("reference", alignment.target, alignment.target_start, alignment.target_end, alignment.query),
        )  # fmt: skip
        for genome, name, start, end, other in sides:
            key = (genome, name)
            intervals[key].append((start, end))
            matches[key] += alignment.matches
            blocks[key] += alignment.block
            shared[key][other] += end - start
            strands[key][other] += -(end - start) if alignment.reverse else end - start

    rows = []
    for genome, sequences in (("species", queries), ("reference", targets)):
        for sequence in sequences:
            key = (genome, sequence.name)
            if key not in intervals:
                rows.append(
                    Row(genome, sequence.name, sequence.label, sequence.length, 0, None, None, None, None)
                )  # fmt: skip
                continue
            others = shared[key]
            best = max(others, key=lambda other: (others[other], other))
            rows.append(
                Row(
                    genome=genome,
                    name=sequence.name,
                    label=sequence.label,
                    length=sequence.length,
                    aligned=min(merged_length(intervals[key]), sequence.length),
                    identity=100 * matches[key] / blocks[key] if blocks[key] else None,
                    best_match=best,
                    share=100 * others[best] / sum(others.values()),
                    strand="-" if strands[key][best] < 0 else "+",
                )
            )
    return rows


def _number(value: float | None) -> str:
    return f"{value:.2f}" if value is not None else "NA"


def write_table(rows: list[Row], path: Path) -> None:
    """The table of both genomes' sequences (docs/formats.md)."""
    columns = (
        "genome",
        "sequence",
        "label",
        "length_bp",
        "aligned_bp",
        "aligned_fraction",
        "identity",
        "best_match",
        "best_match_share",
        "strand",
    )
    with path.open("w", encoding="utf-8", newline="") as sink:
        writer = csv.writer(sink, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            fraction = 100 * row.aligned / row.length if row.length else None
            writer.writerow(
                [
                    row.genome,
                    row.name,
                    row.label,
                    row.length,
                    row.aligned,
                    _number(fraction),
                    _number(row.identity),
                    row.best_match or "NA",
                    _number(row.share),
                    row.strand or "NA",
                ]
            )


def _axis_label(sequence: Sequence) -> str:
    """A sequence's label on its axis; a turned one is marked with a minus sign."""
    return f"{sequence.label} (\u2212)" if sequence.flipped else sequence.label


def _coordinate(value: float) -> str:
    return f"{value:.1f}".removesuffix(".0")


def svg(
    queries: list[Sequence],
    targets: list[Sequence],
    drawn: list[Segment],
    species: str,
    reference: str,
    min_length: int,
) -> str:
    """The dotplot, as a standalone SVG."""
    width = _LEFT + _PLOT + _RIGHT
    height = _TOP + _PLOT + _BOTTOM
    x_total = max(sum(target.length for target in targets), 1)
    y_total = max(sum(query.length for query in queries), 1)
    x_scale, y_scale = _PLOT / x_total, _PLOT / y_total

    def x(value: float) -> str:
        return _coordinate(_LEFT + value * x_scale)

    def y(value: float) -> str:
        return _coordinate(_TOP + _PLOT - value * y_scale)

    bottom, right = _TOP + _PLOT, _LEFT + _PLOT
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Helvetica, Arial, sans-serif" '
        'font-size="11" fill="#222222">',
        f"<title>{escape(species)} against {escape(reference)}</title>",
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
        f'<text x="{width / 2:g}" y="28" text-anchor="middle" font-size="16" '
        f'font-weight="bold">{escape(species)} against {escape(reference)}</text>',
    ]
    # The legend, under the title.
    legend_y = 54
    for number, (colour, text) in enumerate(
        ((SAME, "same strand"), (OTHER, "opposite strand"))
    ):
        left = _LEFT + number * 150
        lines.append(
            f'<line x1="{left}" y1="{legend_y}" x2="{left + 24}" y2="{legend_y}" '
            f'stroke="{colour}" stroke-width="3"/>'
            f'<text x="{left + 30}" y="{legend_y + 4}">{text}</text>'
        )
    lines.append(
        f'<text x="{right}" y="{legend_y + 4}" text-anchor="end" fill="#666666">'
        f"alignments of {min_length:,} bases or more</text>"
    )
    # Each sequence's bounds, and the plot's frame.
    bounds = [f"M{x(t.offset)} {_TOP}V{bottom}" for t in targets[1:]]
    bounds += [f"M{_LEFT} {y(q.offset)}H{right}" for q in queries[1:]]
    if bounds:
        lines.append(
            f'<path d="{"".join(bounds)}" stroke="#cccccc" stroke-width="0.6" fill="none"/>'
        )
    lines.append(
        f'<rect x="{_LEFT}" y="{_TOP}" width="{_PLOT}" height="{_PLOT}" '
        'fill="none" stroke="#888888" stroke-width="1"/>'
    )
    # The alignments: a path per colour.
    for reverse, colour in ((False, SAME), (True, OTHER)):
        path = "".join(
            f"M{x(segment.x0)} {y(segment.y0)}L{x(segment.x1)} {y(segment.y1)}"
            for segment in drawn
            if segment.reverse == reverse
        )
        if path:
            lines.append(
                f'<path d="{path}" stroke="{colour}" stroke-width="1.6" '
                'stroke-linecap="round" fill="none"/>'
            )
    if not drawn:
        lines.append(
            f'<text x="{_LEFT + _PLOT / 2:g}" y="{_TOP + _PLOT / 2:g}" '
            f'text-anchor="middle" fill="#666666">No alignment of {min_length:,} '
            "bases or more</text>"
        )
    # The sequences' names, wide enough to be written.
    for target in targets:
        if target.length * x_scale < _MIN_LABEL:
            continue
        middle = x(target.offset + target.length / 2)
        lines.append(
            f'<text x="{middle}" y="{bottom + 10}" text-anchor="end" '
            f'transform="rotate(-60 {middle} {bottom + 10})">'
            f"{escape(_axis_label(target))}</text>"
        )
    for query in queries:
        if query.length * y_scale < _MIN_LABEL:
            continue
        middle = y(query.offset + query.length / 2)
        lines.append(
            f'<text x="{_LEFT - 6}" y="{middle}" text-anchor="end" '
            f'dominant-baseline="middle">{escape(_axis_label(query))}</text>'
        )
    # The axes' titles.
    lines.append(
        f'<text x="{_LEFT + _PLOT / 2:g}" y="{height - 14}" text-anchor="middle" '
        f'font-size="13">{escape(reference)}</text>'
    )
    lines.append(
        f'<text x="18" y="{_TOP + _PLOT / 2:g}" text-anchor="middle" font-size="13" '
        f'transform="rotate(-90 18 {_TOP + _PLOT / 2:g})">{escape(species)}</text>'
    )
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def figure(
    queries: list[Sequence],
    targets: list[Sequence],
    drawn: list[Segment],
    species: str,
    reference: str,
) -> dict:
    """The dotplot as a plotly figure, for the species' report."""
    traces = []
    for reverse, name, colour in ((False, "same strand", SAME), (True, "opposite strand", OTHER)):  # fmt: skip
        xs: list[int | None] = []
        ys: list[int | None] = []
        texts: list[str | None] = []
        for segment in drawn:
            if segment.reverse != reverse:
                continue
            xs += [segment.x0, segment.x1, None]
            ys += [segment.y0, segment.y1, None]
            for target_at, query_at in (segment.start, segment.end):
                texts.append(
                    f"{segment.query.label}:{query_at:,} ↔ "
                    f"{segment.target.label}:{target_at:,}"
                )
            texts.append(None)
        traces.append(
            {
                "type": "scatter",
                "mode": "lines",
                "name": name,
                "x": xs,
                "y": ys,
                "text": texts,
                "hoverinfo": "text",
                "line": {"color": colour, "width": 2},
            }
        )

    def large(sequences: list[Sequence]) -> list[Sequence]:
        """The sequences wide enough to be named, and bounded, on their axis:
        a draft's thousands of contigs would make the figure slow to draw."""
        total = sum(sequence.length for sequence in sequences)
        return [s for s in sequences if total and s.length / total >= 1 / 150]

    def axis(sequences: list[Sequence], title: str) -> dict:
        total = sum(sequence.length for sequence in sequences)
        shown = large(sequences)
        return {
            "title": {"text": title},
            "range": [0, total],
            "tickmode": "array",
            "tickvals": [s.offset + s.length / 2 for s in shown],
            "ticktext": [_axis_label(s) for s in shown],
            "showgrid": False,
            "zeroline": False,
        }

    def bound(sequence: Sequence, along: str) -> dict:
        line = {"color": "rgba(128, 128, 128, 0.45)", "width": 1}
        if along == "x":
            return {"type": "line", "xref": "x", "yref": "paper", "x0": sequence.offset, "x1": sequence.offset, "y0": 0, "y1": 1, "line": line}  # fmt: skip
        return {"type": "line", "xref": "paper", "yref": "y", "x0": 0, "x1": 1, "y0": sequence.offset, "y1": sequence.offset, "line": line}  # fmt: skip

    x_axis = axis(targets, reference)
    x_axis["tickangle"] = -60
    return {
        "data": traces,
        "layout": {
            "xaxis": x_axis,
            "yaxis": axis(queries, species),
            "shapes": [bound(t, "x") for t in large(targets) if t.offset]
            + [bound(q, "y") for q in large(queries) if q.offset],
            "height": 640,
            "hovermode": "closest",
            "legend": {"orientation": "h", "x": 0, "y": 1.08},
            "margin": {"l": 90, "r": 20, "t": 40, "b": 110},
        },
    }


def reference_title(closest: Path) -> str | None:
    """The closest genome's name and accession, from ``plexembler closest``'s
    table: its first match, or None when it has none."""
    with closest.open(encoding="utf-8", newline="") as handle:
        first = next(csv.DictReader(handle, delimiter="\t"), None)
    if not first or not first.get("name"):
        return None
    name, accession = first["name"], first.get("accession") or "NA"
    if accession == "NA":
        return name
    return f"{name.removeprefix(accession).strip() or name} ({accession})"


def dotplot(
    paf: Path,
    prefix: str,
    outdir: Path,
    genome: Path | None = None,
    reference: Path | None = None,
    reference_name: str | None = None,
    min_length: int = MIN_LENGTH,
) -> dict[str, Path]:
    """Draw a species' genome against another, and tabulate how they align.

    Parameters
    ----------
    paf
        minimap2's alignments of the species' genome (query) to the other (target).
    prefix
        The species' ID: it names the outputs and the y axis, and is taken off
        its sequences' names.
    outdir
        Where the outputs go.
    genome, reference
        The two genomes, FASTA: every sequence, in order, those without an
        alignment included. Without them, the sequences the alignments name,
        longest first.
    reference_name
        The other genome's name, for the x axis; by default its file's.
    min_length
        The shortest alignment drawn, in bases.

    Returns
    -------
    dict
        The paths written: ``svg``, ``figure`` and ``table``.

    Raises
    ------
    ValueError
        When a file cannot be read, or an alignment names a sequence that its
        genome's FASTA file has not.
    """
    alignments = read_paf(paf)
    queries = species_sequences(
        read_headers(genome) if genome else from_alignments(alignments, "query"), prefix
    )
    targets = reference_sequences(
        read_headers(reference) if reference else from_alignments(alignments, "target")
    )
    by_query = {query.name: query for query in queries}
    by_target = {target.name: target for target in targets}
    for alignment in alignments:
        if alignment.query not in by_query:
            raise ValueError(f"{paf}: {alignment.query} is not a sequence of {genome}")
        if alignment.target not in by_target:
            raise ValueError(
                f"{paf}: {alignment.target} is not a sequence of {reference}"
            )
    if reference_name is None:
        reference_name = "reference"
        if reference is not None:
            reference_name = reference.name.removesuffix(".gz")
            for extension in (".fasta", ".fa", ".fna"):
                reference_name = reference_name.removesuffix(extension)

    # The table keeps the species' order; the plot, the reference's.
    rows = alignment_table(queries, targets, alignments)
    queries = lay_out(queries, targets, alignments)
    drawn = segments(alignments, by_query, by_target, min_length)

    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "svg": outdir / f"{prefix}.dotplot.svg",
        "figure": outdir / f"{prefix}.dotplot.json",
        "table": outdir / f"{prefix}.dotplot.tsv",
    }
    paths["svg"].write_text(
        svg(queries, targets, drawn, prefix, reference_name, min_length),
        encoding="utf-8",
    )
    paths["figure"].write_text(
        json.dumps(figure(queries, targets, drawn, prefix, reference_name)) + "\n",
        encoding="utf-8",
    )
    write_table(rows, paths["table"])

    if not drawn:
        logger.warning(
            "%s: no alignment of %d bases or more to %s: an empty dotplot",
            prefix,
            min_length,
            reference_name,
        )
    for genome_name, label in (("species", prefix), ("reference", reference_name)):
        mine = [row for row in rows if row.genome == genome_name]
        length = sum(row.length for row in mine)
        aligned = sum(row.aligned for row in mine)
        logger.info(
            "%s: %.2f %% of its %d bases aligned to the other genome",
            label,
            100 * aligned / length if length else 0.0,
            length,
        )
    turned = [query.label for query in queries if query.flipped]
    if turned:
        logger.info("%s: turned round on the plot: %s", prefix, ", ".join(turned))
    logger.info("Drew %d alignments of %d", len(drawn), len(alignments))
    return paths
