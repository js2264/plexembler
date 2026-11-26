"""Contigs named after their species.

Every contig of the pipeline carries its species in its name,
``<species>_<contig>``, from the moment it is assembled: the pooled libraries
are mapped on the contigs of every species at once, and a read is given to the
species its contig's name starts with (docs/concepts.md).
"""

from __future__ import annotations

import gzip
import logging
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO, cast

logger = logging.getLogger(__name__)

#: A species' ID: a letter, then letters, digits and hyphens. No underscore: the
#: first one of a contig's name ends the species.
SPECIES = re.compile(r"^[A-Za-z][A-Za-z0-9-]*$")
#: The width of a sequence's lines in the FASTA files written.
LINE_WIDTH = 80


@dataclass(frozen=True)
class Renamed:
    """A contig, its new name and its original one."""

    name: str
    original_name: str
    length_bp: int


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


def rename_contigs(fasta: Path, output: Path, species: str) -> list[Renamed]:
    """Write ``fasta``'s records to ``output``, each named after its species.

    Parameters
    ----------
    fasta
        The contigs (FASTA, gzipped or not).
    output
        The renamed contigs (gzipped when the name ends in .gz).
    species
        The species' ID, which starts every new name.

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
    renamed: list[Renamed] = []
    seen: set[str] = set()
    with open_text(fasta) as source, open_text(output, "w") as sink:
        for name, sequence in read_fasta(source):
            name_now = new_name(species, name)
            if name_now in seen:
                raise ValueError(f"two contigs would be named {name_now!r}")
            seen.add(name_now)
            sink.write(f">{name_now}\n")
            for start in range(0, len(sequence), LINE_WIDTH):
                sink.write(sequence[start : start + LINE_WIDTH] + "\n")
            renamed.append(Renamed(name_now, name, len(sequence)))
            logger.debug("%s -> %s (%d bp)", name, name_now, len(sequence))
    logger.info("Renamed %d contigs of %s into %s", len(renamed), species, output)
    return renamed


def write_names(renamed: list[Renamed], table: Path) -> None:
    """Write the table of names: ``name``, ``original_name``, ``length_bp``."""
    with table.open("w", encoding="utf-8") as sink:
        sink.write("name\toriginal_name\tlength_bp\n")
        for contig in renamed:
            sink.write(f"{contig.name}\t{contig.original_name}\t{contig.length_bp}\n")
