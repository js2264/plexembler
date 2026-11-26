"""A species' genome, finished: its scaffolds named, ordered and described.

The scaffolder's sequences come out named after the scaffolder (``scaffold_1``,
…). ``plexembler finalize`` names them as the species' genome is read
(docs/plan.md, M3):

- its chromosomes, longest first: ``<species>_chr01``, ``<species>_chr02``, …:
  the samplesheet's ``n_chromosomes`` longest scaffolds when it gives a number,
  else every scaffold of ``--min-chromosome`` bases or more;
- its mitochondrion, named by the contig screen, as it was: ``<species>_chrM``;
- the others, longest first: ``<species>_unplaced_01``, …

It writes the genome (FASTA, gzipped), its chromosome sizes, the table of
names, and the scaffolder's AGP with every scaffold renamed.
"""

from __future__ import annotations

import gzip
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO, cast

from .rename import LINE_WIDTH, read_fasta

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Named:
    """A scaffold, its new name, its old one and its role."""

    name: str
    original_name: str
    length_bp: int
    role: str


def open_text(path: Path, mode: str = "r") -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, mode + "t", encoding="ascii"))
    return cast(TextIO, path.open(mode, encoding="ascii"))


def agp_components(agp: Path) -> dict[str, list[str]]:
    """Each AGP object's components (its sequence lines' ``W`` parts)."""
    components: dict[str, list[str]] = {}
    with open_text(agp) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 5 and fields[4] not in ("N", "U"):
                components.setdefault(fields[0], []).append(fields[5])
    return components


def name_scaffolds(
    species: str,
    lengths: dict[str, int],
    n_chromosomes: int | None = None,
    min_chromosome: int = 100_000,
    organelles: set[str] | None = None,
) -> list[Named]:
    """Each scaffold's new name and role, chromosomes first, longest first.

    ``organelles`` names the scaffolds made of the species' organelle contigs
    (a scaffolder renames them); by default, those already named ``<species>_chrM``.
    """
    found = organelles or set()
    organelle_names = sorted(
        name for name in lengths if name.startswith(f"{species}_chrM") or name in found
    )
    nuclear = sorted(
        (name for name in lengths if name not in organelle_names),
        key=lambda name: (-lengths[name], name),
    )
    if n_chromosomes is not None:
        chromosomes = nuclear[:n_chromosomes]
    else:
        chromosomes = [name for name in nuclear if lengths[name] >= min_chromosome]
    unplaced = [name for name in nuclear if name not in chromosomes]
    width = max(2, len(str(len(chromosomes))))
    named = [
        Named(f"{species}_chr{i:0{width}d}", name, lengths[name], "chromosome")
        for i, name in enumerate(chromosomes, start=1)
    ]
    named += [
        Named(
            f"{species}_chrM" + ("" if i == 0 else str(i + 1)),
            name,
            lengths[name],
            "organelle",
        )
        for i, name in enumerate(
            sorted(organelle_names, key=lambda n: (-lengths[n], n))
        )
    ]
    width = max(2, len(str(len(unplaced))))
    named += [
        Named(f"{species}_unplaced_{i:0{width}d}", name, lengths[name], "unplaced")
        for i, name in enumerate(unplaced, start=1)
    ]
    return named


def rename_agp(agp: Path, names: dict[str, str], output: Path) -> int:
    """The scaffolder's AGP, every object renamed; its lines, in the genome's order."""
    by_object: dict[str, list[list[str]]] = {}
    with open_text(agp) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] not in names:
                raise ValueError(
                    f"{agp}: {fields[0]} is in the AGP but not in the scaffolds"
                )
            fields[0] = names[fields[0]]
            by_object.setdefault(fields[0], []).append(fields)
    order = {new: index for index, new in enumerate(names.values())}
    written = 0
    with output.open("w", encoding="ascii") as sink:
        sink.write("##agp-version\t2.1\n")
        for name in sorted(by_object, key=lambda n: order[n]):
            for fields in by_object[name]:
                sink.write("\t".join(fields) + "\n")
                written += 1
    return written


def finalize(
    fasta: Path,
    species: str,
    outdir: Path,
    agp: Path | None = None,
    n_chromosomes: int | None = None,
    min_chromosome: int = 100_000,
) -> list[Named]:
    """Name a species' scaffolds and write its genome, sizes, names and AGP.

    Raises
    ------
    ValueError
        When the FASTA is malformed, or the AGP names a scaffold it lacks.
    """
    with open_text(fasta) as handle:
        sequences = dict(read_fasta(handle))
    organelles = set()
    if agp is not None:
        prefix = f"{species}_chrM"
        organelles = {
            name
            for name, parts in agp_components(agp).items()
            if parts and all(part.startswith(prefix) for part in parts)
        }
    named = name_scaffolds(
        species,
        {n: len(s) for n, s in sequences.items()},
        n_chromosomes,
        min_chromosome,
        organelles,
    )
    outdir.mkdir(parents=True, exist_ok=True)
    with gzip.open(outdir / f"{species}.fa.gz", "wt", encoding="ascii") as sink:
        for scaffold in named:
            sequence = sequences[scaffold.original_name]
            sink.write(f">{scaffold.name}\n")
            for start in range(0, len(sequence), LINE_WIDTH):
                sink.write(sequence[start : start + LINE_WIDTH] + "\n")
    with (outdir / f"{species}.chrom.sizes").open("w", encoding="ascii") as sink:
        for scaffold in named:
            sink.write(f"{scaffold.name}\t{scaffold.length_bp}\n")
    with (outdir / f"{species}.names.tsv").open("w", encoding="ascii") as sink:
        sink.write("name\toriginal_name\tlength_bp\trole\n")
        for scaffold in named:
            sink.write(
                f"{scaffold.name}\t{scaffold.original_name}\t{scaffold.length_bp}\t{scaffold.role}\n"
            )
    if agp is not None:
        rename_agp(
            agp, {s.original_name: s.name for s in named}, outdir / f"{species}.agp"
        )
    roles = {
        role: sum(s.role == role for s in named)
        for role in ("chromosome", "organelle", "unplaced")
    }
    logger.info(
        "%s: %d chromosomes, %d organelle, %d unplaced; %d bases",
        species,
        roles["chromosome"],
        roles["organelle"],
        roles["unplaced"],
        sum(s.length_bp for s in named),
    )
    return named
