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

With ``--cut``, the table of ``plexembler junctions --contigs`` on the
scaffolder's joins, it first cuts each join that table finds unsupported: the
scaffold's two sides become two scaffolds, ``<scaffold>.1`` and
``<scaffold>.2`` (more for more cuts), the gap between them dropped, before
they are named. Joins too thin to judge are kept.
"""

from __future__ import annotations

import gzip
import logging
import tempfile
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


@dataclass(frozen=True)
class Cut:
    """A join to cut: its scaffold, where it is (its left contig's last base,
    in the scaffold), its two contigs, and its pairs across, seen and expected."""

    scaffold: str
    position: int
    left: str
    right: str
    observed: str = ""
    expected: str = ""


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


def read_cuts(path: Path) -> list[Cut]:
    """The joins ``plexembler junctions`` found unsupported, from its table."""
    needed = ["scaffold", "position", "left_contig", "right_contig", "verdict"]
    cuts = []
    with open_text(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        missing = [column for column in needed if column not in header]
        if missing:
            raise ValueError(
                f"{path}: no column {', '.join(missing)}, not a table of plexembler junctions"
            )
        for line in handle:
            if not line.strip():
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t"), strict=False))
            if row.get("verdict") != "unsupported":
                continue
            try:
                position = int(row["position"])
            except (KeyError, ValueError):
                raise ValueError(
                    f"{path}: a join at {row.get('position')!r}, not a position"
                ) from None
            cuts.append(
                Cut(
                    row["scaffold"],
                    position,
                    row.get("left_contig", ""),
                    row.get("right_contig", ""),
                    row.get("observed", ""),
                    row.get("expected", ""),
                )
            )
    return cuts


def cut_agp(
    agp: Path, cuts: list[Cut], output: Path
) -> dict[str, list[tuple[str, int, int]]]:
    """The AGP cut at each of ``cuts``, written to ``output``: a scaffold's two
    sides of a join two scaffolds, ``<scaffold>.1``, ``<scaffold>.2``, …, the gap
    between them dropped, their lines numbered again from 1.

    Returns each cut scaffold's pieces: their names, and where each lies in the
    scaffold (0-based, end excluded).

    Raises
    ------
    ValueError
        When a cut is not a join of the AGP.
    """
    rows: dict[str, list[list[str]]] = {}
    with open_text(agp) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"{agp}: an AGP line of {len(fields)} columns")
            rows.setdefault(fields[0], []).append(fields)
    at: dict[str, list[Cut]] = {}
    for cut in cuts:
        if cut.scaffold not in rows:
            raise ValueError(f"{agp}: no scaffold {cut.scaffold} to cut")
        at.setdefault(cut.scaffold, []).append(cut)
    pieces: dict[str, list[tuple[str, int, int]]] = {}
    with output.open("w", encoding="ascii") as sink:
        for scaffold, lines in rows.items():
            groups = [lines]
            if scaffold in at:
                parts = [k for k, f in enumerate(lines) if f[4] not in ("N", "U")]
                bounds = []
                for cut in sorted(at[scaffold], key=lambda c: c.position):
                    index = next(
                        (
                            n
                            for n, k in enumerate(parts[:-1])
                            if int(lines[k][2]) == cut.position
                            and lines[k][5] == cut.left
                            and lines[parts[n + 1]][5] == cut.right
                        ),
                        None,
                    )
                    if index is None:
                        raise ValueError(
                            f"{agp}: no join of {cut.left} and {cut.right} "
                            f"at {cut.position:,} in {scaffold}"
                        )
                    bounds.append((parts[index], parts[index + 1]))
                groups, start = [], 0
                for last, first in bounds:
                    groups.append(lines[start : last + 1])
                    start = first
                groups.append(lines[start:])
                pieces[scaffold] = []
            for number, group in enumerate(groups, start=1):
                name = f"{scaffold}.{number}" if scaffold in at else scaffold
                offset = int(group[0][1]) - 1
                for part, fields in enumerate(group, start=1):
                    moved = [
                        name,
                        str(int(fields[1]) - offset),
                        str(int(fields[2]) - offset),
                        str(part),
                        *fields[4:],
                    ]
                    sink.write("\t".join(moved) + "\n")
                if scaffold in at:
                    pieces[scaffold].append((name, offset, int(group[-1][2])))
    return pieces


def cut_sequences(
    sequences: dict[str, str], pieces: dict[str, list[tuple[str, int, int]]]
) -> dict[str, str]:
    """The scaffolds, each cut one replaced by its pieces' sequences."""
    cut = {}
    for name, sequence in sequences.items():
        if name not in pieces:
            cut[name] = sequence
            continue
        if pieces[name][-1][2] != len(sequence):
            raise ValueError(
                f"{name} is {pieces[name][-1][2]:,} bases long in its AGP, "
                f"{len(sequence):,} in its FASTA"
            )
        for piece, start, end in pieces[name]:
            cut[piece] = sequence[start:end]
    missing = sorted(set(pieces) - set(sequences))
    if missing:
        raise ValueError(f"{missing[0]} is cut but not in the scaffolds")
    return cut


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
    cuts: Path | None = None,
) -> list[Named]:
    """Name a species' scaffolds and write its genome, sizes, names and AGP;
    first, with ``cuts``, cut the joins that table finds unsupported.

    Raises
    ------
    ValueError
        When the FASTA is malformed, the AGP names a scaffold it lacks, or a
        join to cut is not one of the AGP's.
    """
    with open_text(fasta) as handle:
        sequences = dict(read_fasta(handle))
    if cuts is None:
        return write_genome(
            sequences, species, outdir, agp, n_chromosomes, min_chromosome
        )
    if agp is None:
        raise ValueError("cutting the scaffolder's joins needs its AGP (--agp)")
    found = read_cuts(cuts)
    with tempfile.TemporaryDirectory() as scratch:
        cut = Path(scratch) / "cut.agp"
        pieces = cut_agp(agp, found, cut)
        sequences = cut_sequences(sequences, pieces)
        for join in found:
            logger.info(
                "%s: cut %s after %s bases, between %s and %s "
                "(%s pairs across, %s expected)",
                species,
                join.scaffold,
                f"{join.position:,}",
                join.left,
                join.right,
                join.observed or "?",
                join.expected or "?",
            )
        logger.info(
            "%s: %d joins cut, of %d scaffolds",
            species,
            len(found),
            len(pieces),
        )
        return write_genome(
            sequences, species, outdir, cut, n_chromosomes, min_chromosome
        )


def write_genome(
    sequences: dict[str, str],
    species: str,
    outdir: Path,
    agp: Path | None,
    n_chromosomes: int | None,
    min_chromosome: int,
) -> list[Named]:
    """The scaffolds named, and the genome, its sizes, names and AGP written."""
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
