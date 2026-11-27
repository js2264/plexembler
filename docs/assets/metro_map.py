"""plexembler's metro map, in nf-core's style: light and dark SVGs, written
from one description of its lines, stations, sections and legend.

    python docs/assets/metro_map.py    # docs/assets/metro_map_{light,dark}.svg

Each line is a polyline whose corners are rounded; parallel lines keep one
spacing (S), and turn together on concentric arcs. Coordinates are pixels of
the SVG's view box. A step added to the pipeline is a station added here, and
the two maps written again (tests/test_metro_map.py checks that they are).
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

W, H = 1500, 940
S = 22  # the spacing of parallel tracks
LW = 10  # a line's width
HAIRLINE = 1.4  # its dark edge, on each side
R = 72  # the trunk's U-turn radius
FONT = "Inter, 'Helvetica Neue', Helvetica, Arial, sans-serif"

# The trunk, each species' reads then its genome, runs right on row 1, left on
# row 2 and right on row 3; the pooled libraries join it where they are used.
Y_RNA_TOP = 76
Y_HIC_1, Y_WGS_1 = 141, 206
Y1, Y2, Y3 = 296, 508, 730
Y_BRANCH = Y3 + 80  # the RNA-seq's branch, for the expression
Y_RNA_BOTTOM = Y3 + 178
XR, XL = 1392, 72  # the trunk's vertical, in the right and the left U-turn
X_LOOP = 1468  # the RNA-seq loop's right side

PALETTES = {
    "light": {
        "ont": "#7C3AED",
        "wgs": "#2563EB",
        "hic": "#EA580C",
        "rnaseq": "#059669",
        "pooled": "#1E293B",
        "ont_edge": "#4C1D95",
        "wgs_edge": "#1E3A8A",
        "hic_edge": "#7C2D12",
        "rnaseq_edge": "#064E3B",
        "pooled_edge": "#0F172A",
        "station_top": "#FFFFFF",
        "station_bottom": "#E6EBF1",
        "edge": "#1E293B",
        "text": "#0F172A",
        "muted": "#64748B",
        "card_top": "#FFFFFF",
        "card_bottom": "#F3F5F8",
        "card_edge": "#DDE3EA",
        "chip_top": "#475569",
        "chip_bottom": "#1E293B",
        "chip_text": "#FFFFFF",
        "paper": "#FFFFFF",
        "shadow": "#0F172A",
        "lift": "0.22",
        "card_lift": "0.07",
    },
    "dark": {
        "ont": "#A78BFA",
        "wgs": "#60A5FA",
        "hic": "#FB923C",
        "rnaseq": "#34D399",
        "pooled": "#C9D1D9",
        "ont_edge": "#3B0F7A",
        "wgs_edge": "#172E6E",
        "hic_edge": "#6B260E",
        "rnaseq_edge": "#053D2F",
        "pooled_edge": "#0D1117",
        "station_top": "#262C36",
        "station_bottom": "#0D1117",
        "edge": "#E6EDF3",
        "text": "#E6EDF3",
        "muted": "#8B949E",
        "card_top": "#1A2029",
        "card_bottom": "#12171E",
        "card_edge": "#30363D",
        "chip_top": "#F0F6FC",
        "chip_bottom": "#B1BAC4",
        "chip_text": "#0D1117",
        "paper": "#0D1117",
        "shadow": "#000000",
        "lift": "0.55",
        "card_lift": "0.45",
    },
}


def rounded(points: list[tuple[float, float]], radii: float | list[float]) -> str:
    """An SVG path through points, each corner an arc of its radius."""
    corners = len(points) - 2
    if not isinstance(radii, list):
        radii = [radii] * corners
    path = [f"M {points[0][0]:.1f} {points[0][1]:.1f}"]
    for i in range(1, len(points) - 1):
        (x0, y0), (x1, y1), (x2, y2) = points[i - 1], points[i], points[i + 1]
        ux, uy = x0 - x1, y0 - y1
        vx, vy = x2 - x1, y2 - y1
        lu, lv = math.hypot(ux, uy), math.hypot(vx, vy)
        ux, uy, vx, vy = ux / lu, uy / lu, vx / lv, vy / lv
        theta = math.acos(max(-1.0, min(1.0, ux * vx + uy * vy)))
        if math.isclose(theta, math.pi):  # straight on: no corner
            path.append(f"L {x1:.1f} {y1:.1f}")
            continue
        # A corner shares each segment with its neighbour: half of it each,
        # all of it at the line's ends.
        d = min(
            radii[i - 1] / math.tan(theta / 2),
            lu if i == 1 else lu / 2,
            lv if i == len(points) - 2 else lv / 2,
        )
        r = d * math.tan(theta / 2)
        ax, ay = x1 + ux * d, y1 + uy * d
        bx, by = x1 + vx * d, y1 + vy * d
        sweep = 1 if (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1) > 0 else 0
        path.append(f"L {ax:.1f} {ay:.1f}")
        path.append(f"A {r:.1f} {r:.1f} 0 0 {sweep} {bx:.1f} {by:.1f}")
    path.append(f"L {points[-1][0]:.1f} {points[-1][1]:.1f}")
    return " ".join(path)


Label = list[list[tuple[str, str]]]  # lines of (text, style) runs


@dataclass
class Station:
    """A stop: on one track, a circle; on several, a capsule."""

    x: float
    ys: tuple[float, ...]
    label: Label
    where: str = "below"  # the label's side: above or below
    optional: bool = False


def tool(*lines: str) -> Label:
    """A label of third-party tools, a line each."""
    return [[(line, "")] for line in lines]


def own(*lines: str) -> Label:
    """A label of plexembler's own commands, a line each."""
    return [[(line, "own")] for line in lines]


@dataclass
class Section:
    number: int
    title: str
    x0: float
    y0: float
    x1: float
    y1: float


@dataclass
class Icon:
    """A file, at a line's start or end."""

    x: float
    y: float
    kind: str
    line: str
    caption: list[str] = field(default_factory=list)


# The lines ------------------------------------------------------------------

X0 = 92  # where the lines start, right of their input files
LINES: dict[str, tuple[list[tuple[float, float]], float | list[float]]] = {
    # Each species' nanopore reads, then its contigs, then its genome.
    "ont": ([(X0, Y1), (XR, Y1), (XR, Y2), (XL, Y2), (XL, Y3), (1400, Y3)], R),
    # The pooled WGS: its k-mers; the pool's splits; the short-read polishing.
    "wgs": (
        [
            (X0, Y_HIC_1),
            (108, Y_HIC_1),
            (173, Y_WGS_1),
            (1236, Y_WGS_1),
            (1302, Y1 - S),
            (XR + S, Y1 - S),
            (XR + S, Y2 + S),
            (768, Y2 + S),
        ],
        [26, 26, 26, 26, R + S, R + S],
    ),
    # The pooled Hi-C: the splits, the polishing without WGS, the scaffolds,
    # the maps.
    "hic": (
        [
            (X0, Y_HIC_1),
            (1202, Y_HIC_1),
            (1311, Y1 - 2 * S),
            (XR + 2 * S, Y1 - 2 * S),
            (XR + 2 * S, Y2 + 2 * S),
            (750, Y2 + 2 * S),
            (728, Y2 + S),
            (296, Y2 + S),
        ],
        [26, 26, R + 2 * S, R + 2 * S, 20, 20],
    ),
    # The pooled RNA-seq: round the map, to the annotation and the expression.
    "rnaseq": (
        [
            (X0, Y_HIC_1),
            (108, Y_HIC_1),
            (173, Y_RNA_TOP),
            (X_LOOP, Y_RNA_TOP),
            (X_LOOP, Y_RNA_BOTTOM),
            (162, Y_RNA_BOTTOM),
            (162, Y3 + S),
            (510, Y3 + S),
            (568, Y_BRANCH),
            (900, Y_BRANCH),
        ],
        [26, 26, 56, 56, 40, 40, 26, 26],
    ),
}

# The stations ---------------------------------------------------------------

POOLED_1 = (Y_RNA_TOP, Y_HIC_1, Y_WGS_1)
POOLED_2 = (Y2, Y2 + S, Y2 + 2 * S)

STATIONS = [
    # 1. Reads: each species' nanopore reads; the pooled libraries.
    Station(222, POOLED_1, tool("FastQC")),
    Station(322, (Y_WGS_1,), tool("meryl"), "above"),
    Station(142, (Y1,), tool("samtools"), optional=True),
    Station(232, (Y1,), tool("NanoPlot")),
    Station(322, (Y1,), tool("chopper")),
    Station(412, (Y1,), tool("LRGE"), optional=True),
    Station(502, (Y1,), tool("Filtlong")),
    Station(602, (Y1,), tool("meryl", "GenomeScope2")),
    # 2. Each species' contigs; the identity between the species of a pool.
    Station(712, (Y1,), [[("Flye", "")], [("hifiasm", "muted")]]),
    Station(802, (Y1,), own("rename")),
    Station(892, (Y1,), [[("Medaka", "")], [("Dorado", "muted")]]),
    Station(986, (Y1,), tool("gfastats", "BUSCO", "Merqury")),
    Station(1078, (Y1,), own("select")),
    Station(1166, (Y1,), [[("skani", "")], [("ani", "own")]]),
    # 3. The pool, split; the contigs screened, split again, polished.
    Station(1300, POOLED_2, [[("bwa-mem2", "")], [("split", "own")]], "above"),
    Station(1206, (Y2,), [[("minimap2", "")], [("screen", "own")]], "above"),
    Station(1112, POOLED_2, [[("split", "own")], [("screened", "muted")]], "above"),
    Station(1028, (Y2 + 2 * S,), own("digest"), optional=True),
    Station(944, POOLED_2, tool("BWA", "Polypolish"), "above"),
    Station(856, POOLED_2, tool("pypolca"), "above"),
    Station(
        768,
        POOLED_2,
        [[("meryl", "")], [("Merqury", "")], [("select", "own")]],
        "above",
    ),
    # 4. The chromosomes; each species' Hi-C map.
    Station(
        664, (Y2, Y2 + S), tool("bwa-mem2", "YaHS \N{MULTIPLICATION SIGN}2"), "above"
    ),
    Station(574, (Y2,), own("finalize"), "above"),
    Station(484, (Y2,), tool("gfastats", "BUSCO", "Merqury"), "above"),
    Station(526, (Y2 + S,), tool("bwa-mem2", "pairtools")),
    Station(436, (Y2 + S,), tool("cooler")),
    Station(346, (Y2 + S,), own("junctions")),
    # 5. The annotation; the expression and the tracks.
    Station(250, (Y3, Y3 + S), [[("STAR", "")], [("split", "own")]], "above"),
    Station(365, (Y3,), tool("RepeatModeler", "RepeatMasker"), "above"),
    Station(480, (Y3, Y3 + S), [[("BRAKER3", "")], [("AUGUSTUS", "muted")]], "above"),
    Station(590, (Y3,), own("genes"), "above"),
    Station(695, (Y3,), tool("BUSCO"), "above"),
    Station(800, (Y3,), tool("tRNAscan-SE", "barrnap"), "above"),
    Station(625, (Y_BRANCH,), tool("featureCounts")),
    Station(740, (Y_BRANCH,), own("expression")),
    Station(855, (Y_BRANCH,), tool("deepTools")),
    # 6. The reports.
    Station(
        1075,
        (Y3,),
        [[("sourmash", "")], [("closest", "own")], [("dotplot", "own")]],
        "above",
        True,
    ),
    Station(1170, (Y3,), own("report"), "above"),
    Station(1268, (Y3,), [[("browser", "own")], [("chromosomix", "")]], "above"),
    Station(1365, (Y3,), tool("MultiQC"), "above"),
]

SECTIONS = [
    Section(1, "Reads", 112, 30, 648, Y1 + 80),
    Section(2, "Contigs", 664, 30, 1214, Y1 + 80),
    Section(3, "The pool: split, screened, polished", 726, Y2 - 118, 1350, Y2 + 90),
    Section(4, "Chromosomes and maps", 232, Y2 - 118, 712, Y2 + 90),
    Section(5, "Annotation and expression", 196, Y3 - 118, 955, Y_BRANCH + 46),
    Section(6, "Reports", 967, Y3 - 118, 1452, Y3 + 34),
]

ICONS = [
    Icon(48, Y_HIC_1, "FASTQ", "pooled", ["pooled", "libraries"]),
    Icon(48, Y1, "FASTQ", "ont", ["nanopore,", "per species"]),
    Icon(266, Y2 + S + 14, "MCOOL", "hic"),
    Icon(928, Y_BRANCH, "BIGWIG", "rnaseq"),
    Icon(1426, Y3, "HTML", "ont"),
]

LEGEND = (967, Y3 + 46, 485, 112)  # x, y, width, height


# Drawing --------------------------------------------------------------------


def draw_label(station: Station) -> str:
    leading = 19.5
    lines = station.label
    if station.where == "above":
        top = min(station.ys) - 18 - leading * (len(lines) - 1)
    else:
        top = max(station.ys) + 28
    out = []
    for i, line in enumerate(lines):
        spans = "".join(
            f'<tspan class="{style}">{text}</tspan>'
            if style
            else f"<tspan>{text}</tspan>"
            for text, style in line
        )
        cls = "label muted" if station.optional else "label"
        out.append(
            f'<text class="{cls}" x="{station.x:.1f}" y="{top + i * leading:.1f}">{spans}</text>'
        )
    return "\n".join(out)


def draw_station(station: Station) -> str:
    dash = ' stroke-dasharray="4 3.2"' if station.optional else ""
    if len(station.ys) == 1:
        return f'<circle class="station" cx="{station.x:.1f}" cy="{station.ys[0]:.1f}" r="9"{dash}/>'
    top, bottom = min(station.ys) - 11, max(station.ys) + 11
    return (
        f'<rect class="station" x="{station.x - 11:.1f}" y="{top:.1f}" width="22" '
        f'height="{bottom - top:.1f}" rx="11"{dash}/>'
    )


def draw_icon(icon: Icon, colors: dict[str, str]) -> str:
    color = colors[icon.line]
    out = [
        f'<g transform="translate({icon.x - 16:.1f} {icon.y - 21:.1f})">',
        '<path d="M 2 0 H 22 L 32 10 V 40 Q 32 42 30 42 H 2 Q 0 42 0 40 V 2 Q 0 0 2 0 Z" '
        f'fill="{colors["paper"]}" stroke="{color}" stroke-width="2.5" stroke-linejoin="round"/>',
        f'<path d="M 22 0 V 8 Q 22 10 24 10 H 32" fill="none" stroke="{color}" '
        'stroke-width="2.5" stroke-linejoin="round"/>',
        f'<rect x="-6" y="22" width="44" height="14" rx="4" fill="{color}" '
        f'stroke="{colors[icon.line + "_edge"]}" stroke-width="1"/>',
        f'<text class="icon" x="16" y="32.6">{icon.kind}</text>',
        "</g>",
    ]
    for i, line in enumerate(icon.caption):
        out.append(
            f'<text class="caption" x="{icon.x:.1f}" y="{icon.y + 41 + i * 16:.1f}">{line}</text>'
        )
    return "\n".join(out)


def draw_card(x0: float, y0: float, width: float, height: float) -> str:
    return (
        f'<rect class="card" x="{x0:.1f}" y="{y0:.1f}" width="{width:.1f}" '
        f'height="{height:.1f}" rx="16"/>'
    )


def draw_section(section: Section) -> str:
    x0, y0 = section.x0, section.y0
    return "\n".join(
        [
            draw_card(x0, y0, section.x1 - x0, section.y1 - y0),
            f'<circle class="chip" cx="{x0 + 22:.1f}" cy="{y0 + 22:.1f}" r="11.5"/>',
            f'<text class="chip-number" x="{x0 + 22:.1f}" y="{y0 + 26.8:.1f}">{section.number}</text>',
            f'<text class="section-title" x="{x0 + 41:.1f}" y="{y0 + 27.5:.1f}">{section.title}</text>',
        ]
    )


def lane(line: str, d: str, colors: dict[str, str]) -> tuple[str, str]:
    """A line's two strokes: its dark edge, wider, then its colour."""
    return (
        f'<path class="lane-edge" stroke="{colors[line + "_edge"]}" d="{d}"/>',
        f'<path class="lane" id="{line}" stroke="{colors[line]}" d="{d}"/>',
    )


def draw_legend(colors: dict[str, str]) -> str:
    x0, y0, width, height = LEGEND
    edges, lanes, rest = [], [], []
    names = [
        ("ont", "nanopore, per species"),
        ("wgs", "WGS, pooled"),
        ("hic", "Hi-C, pooled"),
        ("rnaseq", "RNA-seq, pooled"),
    ]
    for i, (line, name) in enumerate(names):
        y = y0 + 22 + i * 23
        edge, stroke = lane(
            line, f"M {x0 + 22:.1f} {y:.1f} L {x0 + 54:.1f} {y:.1f}", colors
        )
        edges.append(edge)
        lanes.append(stroke.replace(f' id="{line}"', ""))
        rest.append(
            f'<text class="legend-text" x="{x0 + 68:.1f}" y="{y + 5:.1f}">{name}</text>'
        )
    x = x0 + 262
    shapes = [
        f'<circle class="station" cx="{x:.1f}" cy="{y0 + 22:.1f}" r="9"/>',
        f'<rect class="station" x="{x - 11:.1f}" y="{y0 + 33:.1f}" width="22" height="27" rx="11"/>',
        f'<circle class="station" cx="{x:.1f}" cy="{y0 + 68:.1f}" r="9" stroke-dasharray="4 3.2"/>',
    ]
    for i, name in enumerate(["a step", "a step on several lines", "optional"]):
        rest.append(
            f'<text class="legend-text" x="{x + 22:.1f}" y="{y0 + 27 + i * 23:.1f}">{name}</text>'
        )
    rest.append(
        f'<text class="legend-text" x="{x - 11:.1f}" y="{y0 + 96:.1f}">'
        '<tspan class="own">bold</tspan>: plexembler\'s own</text>'
    )
    return "\n".join(
        [
            f'<g filter="url(#card-lift)">{draw_card(x0, y0, width, height)}</g>',
            '<g filter="url(#lift)">',
            *edges,
            "</g>",
            *lanes,
            '<g filter="url(#lift-small)">',
            *shapes,
            "</g>",
            *rest,
        ]
    )


def shadow(name: str, blur: float, dy: float, opacity: str, color: str) -> str:
    """A soft shadow beneath what the filter is set on."""
    return "\n".join(
        [
            f'<filter id="{name}" filterUnits="userSpaceOnUse" x="-20" y="-20" '
            f'width="{W + 40}" height="{H + 40}" color-interpolation-filters="sRGB">',
            f'<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/>',
            f'<feOffset dy="{dy}" result="offset"/>',
            f'<feFlood flood-color="{color}" flood-opacity="{opacity}"/>',
            '<feComposite in2="offset" operator="in"/>',
            '<feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
            "</filter>",
        ]
    )


def gradient(name: str, top: str, bottom: str) -> str:
    """A vertical gradient, from top to bottom of each shape it fills."""
    return (
        f'<linearGradient id="{name}" x1="0" y1="0" x2="0" y2="1">'
        f'<stop offset="0" stop-color="{top}"/><stop offset="1" stop-color="{bottom}"/>'
        "</linearGradient>"
    )


def svg(theme: str) -> str:
    colors = PALETTES[theme]
    rules = [
        f"text {{ font-family: {FONT}; fill: {colors['text']}; }}",
        ".label { font-size: 17px; font-weight: 500; text-anchor: middle; }",
        ".own { font-weight: 750; }",
        f".muted {{ fill: {colors['muted']}; }}",
        ".caption { font-size: 13.5px; font-weight: 500; text-anchor: middle; "
        f"fill: {colors['muted']}; }}",
        ".icon { font-size: 10.5px; font-weight: 750; text-anchor: middle; "
        f"fill: {colors['paper']}; letter-spacing: 0.4px; }}",
        f".station {{ fill: url(#station-fill); stroke: {colors['edge']}; stroke-width: 3.2; }}",
        f".card {{ fill: url(#card-fill); stroke: {colors['card_edge']}; stroke-width: 1.2; }}",
        ".chip { fill: url(#chip-fill); }",
        ".chip-number { font-size: 13.5px; font-weight: 750; text-anchor: middle; "
        f"fill: {colors['chip_text']}; }}",
        ".section-title { font-size: 15.5px; font-weight: 700; }",
        ".legend-text { font-size: 14px; font-weight: 500; }",
        ".lane, .lane-edge { fill: none; stroke-linecap: round; stroke-linejoin: round; }",
        f".lane {{ stroke-width: {LW}; }}",
        f".lane-edge {{ stroke-width: {LW + 2 * HAIRLINE}; }}",
    ]
    edges, lanes = [], []
    for line, (points, radii) in LINES.items():
        edge, stroke = lane(line, rounded(points, radii), colors)
        edges.append(edge)
        lanes.append(stroke)
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
        'role="img" aria-labelledby="title desc">',
        '<title id="title">plexembler\'s metro map</title>',
        '<desc id="desc">Each species\' nanopore reads are assembled into its contigs; the '
        "pooled WGS, Hi-C and RNA-seq libraries are split between the species, and polish, "
        "scaffold, map and annotate each genome; plexembler reports each species.</desc>",
        "<defs>",
        shadow("lift", 1.4, 1.5, colors["lift"], colors["shadow"]),
        shadow("lift-small", 1.0, 1.0, colors["lift"], colors["shadow"]),
        shadow("card-lift", 5.0, 2.0, colors["card_lift"], colors["shadow"]),
        gradient("station-fill", colors["station_top"], colors["station_bottom"]),
        gradient("card-fill", colors["card_top"], colors["card_bottom"]),
        gradient("chip-fill", colors["chip_top"], colors["chip_bottom"]),
        "</defs>",
        "<style>",
        *rules,
        "</style>",
        '<g id="sections">',
        '<g filter="url(#card-lift)">',
        *(draw_card(s.x0, s.y0, s.x1 - s.x0, s.y1 - s.y0) for s in SECTIONS),
        "</g>",
        *(line for s in SECTIONS for line in draw_section(s).split("\n")[1:]),
        "</g>",
        '<g id="lines">',
        # Every line's edge before any colour: where lines meet, they merge.
        '<g filter="url(#lift)">',
        *edges,
        "</g>",
        *lanes,
        "</g>",
        '<g id="stations" filter="url(#lift-small)">',
        *(draw_station(station) for station in STATIONS),
        "</g>",
        '<g id="labels">',
        *(draw_label(station) for station in STATIONS),
        "</g>",
        '<g id="files" filter="url(#lift-small)">',
        *(draw_icon(icon, colors) for icon in ICONS),
        "</g>",
        '<g id="legend">',
        draw_legend(colors),
        "</g>",
        "</svg>",
    ]
    return "\n".join(parts) + "\n"


def main(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for theme in PALETTES:
        (outdir / f"metro_map_{theme}.svg").write_text(svg(theme), encoding="utf-8")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
