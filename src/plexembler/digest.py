"""Hi-C reads cut at their ligation junctions: genome, to polish with.

A Hi-C read is chimeric at a ligation junction: the filled-in end of one
restriction fragment, ligated to another's. Between junctions, a read is
genome, sequenced with the same chemistry as a WGS read. ``plexembler digest``
cuts each read where a junction's left end stops, as hicstuff's ``cutsite``
does, so that each piece keeps its half of the site; drops the pieces shorter
than ``--min-piece``; and writes every piece once, as a single read: a Hi-C
pair is two loci, not a fragment, and a polisher would read its pairing as an
insert size (docs/plan.md, M2).

A junction, for two enzymes of the library (an enzyme with itself included):
the first's site up to the end of its overhang, filled in, then the second's
from the start of its overhang. DpnII's (^GATC) with itself is GATCGATC; HinfI's
(G^ANTC) with itself GANTANTC; Arima's kit has both enzymes, and four junctions.
"""

from __future__ import annotations

import gzip
import itertools
import json
import logging
import re
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO, cast

logger = logging.getLogger(__name__)

#: Each enzyme's site and where it cuts its top strand (an index into the site).
ENZYMES: dict[str, tuple[str, int]] = {
    "DpnII": ("GATC", 0),
    "MboI": ("GATC", 0),
    "Sau3AI": ("GATC", 0),
    "HinfI": ("GANTC", 1),
    "DdeI": ("CTNAG", 1),
    "MseI": ("TTAA", 1),
    "Csp6I": ("GTAC", 1),
    "HindIII": ("AAGCTT", 1),
    "BglII": ("AGATCT", 1),
    "EcoRI": ("GAATTC", 1),
    "NlaIII": ("CATG", 4),
}
#: Kits, by the name a samplesheet may give them.
KITS = {"arima": ("DpnII", "HinfI")}


@dataclass(frozen=True)
class Junction:
    """A ligation junction: its pattern, and where it is cut, from its start."""

    pattern: str
    cut: int


def ends(enzyme: str) -> tuple[str, str]:
    """An enzyme's ligated ends, filled in: the left one up to the far end of its
    overhang, the right one from the near end of it."""
    if enzyme not in ENZYMES:
        known = ", ".join(sorted(ENZYMES))
        raise ValueError(f"{enzyme!r} is not an enzyme plexembler knows: {known}")
    site, cut = ENZYMES[enzyme]
    return site[: max(cut, len(site) - cut)], site[min(cut, len(site) - cut) :]


def parse_enzymes(text: str) -> list[str]:
    """A library's enzymes, comma-separated, as the samplesheet names them."""
    names: list[str] = []
    for word in (part.strip() for part in text.split(",") if part.strip()):
        for name in KITS.get(word.lower(), (word,)):
            if name not in names:
                names.append(name)
    if not names or names == ["none"]:
        raise ValueError(
            f"{text!r} names no enzyme: a library made without one (Micro-C) has no junction"
        )
    for name in names:
        ends(name)
    return names


def junctions(enzymes: list[str]) -> list[Junction]:
    """Every junction of a library's enzymes, the longest first."""
    found = {}
    for first, second in itertools.product(enzymes, repeat=2):
        left, _ = ends(first)
        _, right = ends(second)
        pattern = (left + right).replace("N", "[ACGT]")
        found[pattern] = Junction(pattern, len(left))
    return sorted(found.values(), key=lambda j: -len(j.pattern))


def cut_points(sequence: str, patterns: list[Junction]) -> list[int]:
    """Where a read is cut: each junction's cut, the junctions not overlapping."""
    points = []
    regex = re.compile(
        "|".join(f"(?P<j{i}>{j.pattern})" for i, j in enumerate(patterns))
    )
    for match in regex.finditer(sequence.upper()):
        index = int(match.lastgroup[1:]) if match.lastgroup else 0
        points.append(match.start() + patterns[index].cut)
    return points


def read_fastq(handle: TextIO) -> Iterator[tuple[str, str, str]]:
    """Each record of a FASTQ file: its name (the header's first word), sequence, quality.

    Raises
    ------
    ValueError
        When a record is not four lines, its header an @.
    """
    for number, header in enumerate(handle, start=1):
        if not header.strip():
            continue
        sequence, plus, quality = (next(handle, "") for _ in range(3))
        if not header.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"record {number}: not a FASTQ record")
        yield header[1:].split()[0], sequence.strip(), quality.strip()


def open_text(path: Path, mode: str = "r") -> TextIO:
    if path.suffix == ".gz":
        return cast(
            TextIO, gzip.open(path, mode + "t", encoding="ascii", compresslevel=1)
        )
    return cast(TextIO, path.open(mode, encoding="ascii"))


@dataclass
class Counts:
    """What the digestion did: reads by their junctions, and pieces."""

    reads: Counter
    pieces: int = 0
    short: int = 0
    bases: int = 0


def digest(
    fastqs: list[Path], output: Path, enzymes: list[str], min_piece: int = 40
) -> Counts:
    """Cut the reads of ``fastqs`` at their junctions; write the pieces to ``output``.

    Parameters
    ----------
    fastqs
        A Hi-C library's reads: its first and second reads, or one of them.
    output
        The pieces, as single reads (gzipped when the name ends in .gz).
    enzymes
        The library's enzymes.
    min_piece
        The shortest piece kept.

    Returns
    -------
    Counts
        The reads with 0, 1 and 2 or more junctions, and the pieces.
    """
    patterns = junctions(enzymes)
    counts = Counts(Counter())
    with open_text(output, "w") as sink:
        for mate, fastq in enumerate(fastqs, start=1):
            with open_text(fastq) as source:
                for name, sequence, quality in read_fastq(source):
                    points = cut_points(sequence, patterns)
                    counts.reads[min(len(points), 2)] += 1
                    bounds = [0, *points, len(sequence)]
                    for k, (start, stop) in enumerate(itertools.pairwise(bounds)):
                        if stop - start < min_piece:
                            counts.short += 1
                            continue
                        counts.pieces += 1
                        counts.bases += stop - start
                        sink.write(
                            f"@{name}_{mate}.{k}\n{sequence[start:stop]}\n+\n{quality[start:stop]}\n"
                        )
    return counts


def write_outputs(counts: Counts, outdir: Path, prefix: str) -> dict[str, Path]:
    """The digestion's numbers, and MultiQC's table."""
    reads = sum(counts.reads.values())
    numbers = {
        "reads": reads,
        "no_junction": counts.reads[0],
        "one_junction": counts.reads[1],
        "junctions_2_or_more": counts.reads[2],
        "ligation_share": (counts.reads[1] + counts.reads[2]) / reads if reads else 0.0,
        "pieces": counts.pieces,
        "short_pieces": counts.short,
        "bases": counts.bases,
    }
    paths = {
        "json": outdir / f"{prefix}.digest.json",
        "multiqc": outdir / f"{prefix}.digest_mqc.tsv",
    }
    paths["json"].write_text(json.dumps(numbers, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# id: 'plexembler_digest'",
        "# section_name: 'Hi-C reads, cut at their junctions'",
        "# description: 'Each species'' Hi-C reads, by their ligation junctions, and "
        "the pieces cut from them to polish with (plexembler digest).'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'plexembler_digest_table'",
        "#     namespace: 'plexembler'",
        "Sample\tReads\tWith a junction (%)\tPieces\tBases",
        f"{prefix}\t{reads}\t{100 * numbers['ligation_share']:.2f}\t{counts.pieces}\t{counts.bases}",
    ]
    paths["multiqc"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths
