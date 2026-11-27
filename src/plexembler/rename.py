"""Contigs named after their species.

Every contig of the pipeline carries its species in its name,
``<species>_<contig>``, from the moment it is assembled: the pooled libraries
are mapped on the contigs of every species at once, and a read is given to the
species its contig's name starts with (docs/concepts.md).

A contig that is two or more copies of a unit, end to end, a circular genome
its assembler went round more than once, is cut to its first copy on the way
(``circular``): everything after the assembly sees one copy of it.
"""

from __future__ import annotations

import gzip
import logging
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO, cast

from .circular import MIN_IDENTITY, MIN_UNIT, AssemblyInfo, find_copies

logger = logging.getLogger(__name__)

#: A species' ID: a letter, then letters, digits and hyphens. No underscore: the
#: first one of a contig's name ends the species.
SPECIES = re.compile(r"^[A-Za-z][A-Za-z0-9-]*$")
#: The width of a sequence's lines in the FASTA files written.
LINE_WIDTH = 80


@dataclass(frozen=True)
class Renamed:
    """A contig, its new name and its original one; its length as written and
    as assembled; the copies of its unit it was, cut to one; and whether it is
    circular, or None when nothing says."""

    name: str
    original_name: str
    length_bp: int
    original_length_bp: int
    copies: int = 1
    circular: bool | None = None


def check_species(species: str) -> str:
    """Return ``species`` if it can start a contig's name, else raise ValueError."""
    if not SPECIES.match(species):
        raise ValueError(
            f"{species!r} is not a species ID: a letter, then letters, digits and hyphens"
        )
    return species


def open_text(path: Path, mode: str = "r") -> TextIO:
    """Open a text file, compressed with gzip when its name ends in .gz."""
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, mode + "t", encoding="utf-8"))
    return cast(TextIO, path.open(mode, encoding="utf-8"))


def read_fasta(handle: TextIO) -> Iterator[tuple[str, str]]:
    """Each record of a FASTA file: its name (the header's first word) and sequence.

    Raises
    ------
    ValueError
        When the file does not start with a header, or a header has no name.
    """
    name: str | None = None
    chunks: list[str] = []
    for number, line in enumerate(handle, start=1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            words = line[1:].split()
            if not words:
                raise ValueError(f"line {number}: a header without a name")
            name, chunks = words[0], []
        elif name is None:
            raise ValueError(f"line {number}: a sequence before the first header")
        else:
            chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def new_name(species: str, name: str) -> str:
    """The contig's name, starting with its species; a name that already does stays."""
    prefix = f"{species}_"
    return name if name.startswith(prefix) else prefix + name


def one_copy(
    name: str,
    sequence: str,
    hint: AssemblyInfo | None,
    min_unit: int,
    min_identity: float,
) -> tuple[str, int]:
    """A contig's sequence, cut to its first copy when it is copies of a unit;
    and the copies it was. A contig its assembler says is linear stays whole."""
    if hint is not None and not hint.circular:
        return sequence, 1
    found = find_copies(sequence, min_unit, min_identity)
    if found is None:
        if hint is not None and (hint.multiplicity or 0) >= 2:
            logger.debug(
                "%s: circular, multiplicity %d in the assembly info; one copy "
                "of itself",
                name,
                hint.multiplicity,
            )
        return sequence, 1
    said = (
        f"circular in the assembly info, multiplicity {hint.multiplicity or 'NA'}"
        if hint is not None
        else "by its sequence alone"
    )
    logger.info(
        "%s: %d copies of a %d-bp unit, %.2f %% alike (%s); one kept, %d of its %d bp",
        name,
        found.copies,
        found.unit_bp,
        100 * found.identity,
        said,
        found.unit_bp,
        len(sequence),
    )
    return sequence[: found.unit_bp], found.copies


def rename_contigs(
    fasta: Path,
    output: Path,
    species: str,
    info: dict[str, AssemblyInfo] | None = None,
    collapse: bool = True,
    min_unit: int = MIN_UNIT,
    min_identity: float = MIN_IDENTITY,
) -> list[Renamed]:
    """Write ``fasta``'s records to ``output``, each named after its species.

    Parameters
    ----------
    fasta
        The contigs (FASTA, gzipped or not).
    output
        The renamed contigs (gzipped when the name ends in .gz).
    species
        The species' ID, which starts every new name.
    info
        What the assembler says of each contig, by its name in ``fasta``
        (Flye's ``assembly_info.txt``); a contig it says is linear is not cut.
    collapse
        Whether a contig that is copies of a unit, end to end, is cut to one.
    min_unit, min_identity
        The shortest unit looked for, in bases, and the lowest identity
        between copies.

    Returns
    -------
    list of Renamed
        Each contig, in the order of the file.

    Raises
    ------
    ValueError
        When the species' ID is not one, the FASTA is malformed, or two contigs
        get the same name.
    """
    check_species(species)
    info = info or {}
    renamed: list[Renamed] = []
    seen: set[str] = set()
    with open_text(fasta) as source, open_text(output, "w") as sink:
        for name, sequence in read_fasta(source):
            name_now = new_name(species, name)
            if name_now in seen:
                raise ValueError(f"two contigs would be named {name_now!r}")
            seen.add(name_now)
            hint = info.get(name)
            written, copies = sequence, 1
            if collapse:
                written, copies = one_copy(
                    name_now, sequence, hint, min_unit, min_identity
                )
            sink.write(f">{name_now}\n")
            for start in range(0, len(written), LINE_WIDTH):
                sink.write(written[start : start + LINE_WIDTH] + "\n")
            # Copies end to end are a circle's: the last closes on the first.
            circular = hint.circular if hint is not None else (copies > 1 or None)
            renamed.append(
                Renamed(name_now, name, len(written), len(sequence), copies, circular)
            )
            logger.debug("%s -> %s (%d bp)", name, name_now, len(written))
    logger.info("Renamed %d contigs of %s into %s", len(renamed), species, output)
    cut = [contig for contig in renamed if contig.copies > 1]
    if cut:
        logger.info(
            "%d of them cut to one copy of their unit: %s",
            len(cut),
            ", ".join(contig.name for contig in cut),
        )
    return renamed


def write_names(renamed: list[Renamed], table: Path) -> None:
    """Write the table of names: ``name``, ``original_name``, ``length_bp``,
    ``original_length_bp``, ``copies`` and ``circular`` (``yes``, ``no``, or
    ``NA`` when nothing says)."""
    said = {True: "yes", False: "no", None: "NA"}
    with table.open("w", encoding="utf-8") as sink:
        sink.write(
            "name\toriginal_name\tlength_bp\toriginal_length_bp\tcopies\tcircular\n"
        )
        for contig in renamed:
            cells = [
                contig.name,
                contig.original_name,
                str(contig.length_bp),
                str(contig.original_length_bp),
                str(contig.copies),
                said[contig.circular],
            ]
            sink.write("\t".join(cells) + "\n")
