"""A pooled library split between its species (docs/concepts.md, deconvolution).

Every read of a pooled library is aligned to the contigs of every species of
its pool at once (bwa-mem2; ``-5SP`` for Hi-C); a contig's name starts with its
species, ``<species>_<contig>``. ``plexembler split`` reads those alignments as
bwa-mem2 writes them, a pair's records together, and gives each pair to the
species of its two mates' primary alignments when both are:

- **unique**: a mapping quality of ``--min-mapq`` or more. A mate below it is
  ambiguous; bwa-mem2's alternative hits (its ``XA`` tag) say between which
  species, and the pool's table of ambiguous pairs counts it there;
- **near-identical**: an identity of ``--min-identity`` or more, one minus the
  edit distance (``NM``) over the alignment's columns. A read of a region
  missing from its own species' contigs aligns, unique, to a sister species'
  copy, 85 to 95 % identical: the identity, not the uniqueness, catches it;
- **of one species**: a pair whose mates land in two species is counted, not
  kept. For Hi-C, those pairs are the pool's random ligations: their share,
  over the share two random partners would be of two species (1 - sum of the
  species' squared shares), estimates the random-ligation rate, epsilon.

Each library keeps a share of its pairs apart (``--holdout``, by a hash of the
read's name, the same at every run): the k-mers the polishing's QV guard counts
against, which the polishing never sees. A Hi-C library loses its duplicates:
pairs whose mates start at the same places, on the same strands, as a pair
before them; and its unique pairs between two contigs are counted, the contig
screen's Hi-C signal.

Each species' pairs are written as FASTQ, their reads as sequenced; the pool's
numbers as a table, a JSON file and MultiQC's table.

An RNA-seq library (``--type rnaseq``) is aligned by STAR to the pool's
genomes, whose sequences are each species' own: its pairs are written as the
alignments they are, each species' in a SAM file of its genome's sequences,
ready to sort and count. STAR has no ``XA`` tag: a multimapper's secondary
records say between which species it is ambiguous. RNA-seq keeps its
duplicates and holds no pair out.
"""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import logging
import re
import sys
import zlib
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, TextIO

logger = logging.getLogger(__name__)

#: A CIGAR's operations, and those that make an alignment's columns.
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_COLUMNS = set("MID=X")
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
#: The flags of a SAM record this module reads.
UNMAPPED, REVERSE, FIRST, SECOND, SECONDARY, SUPPLEMENTARY = 4, 16, 64, 128, 256, 2048


@dataclass
class Mate:
    """A mate's primary alignment: where it lands, and the read as sequenced."""

    name: str
    flag: int
    contig: str
    position: int
    mapq: int
    cigar: str
    sequence: str
    quality: str
    tags: dict[str, str]
    #: The record, as it came.
    line: str = ""
    #: The species of the read's secondary alignments (STAR's multimappers).
    alternatives: frozenset[str] = frozenset()

    @property
    def mapped(self) -> bool:
        return not self.flag & UNMAPPED and self.contig != "*"

    @property
    def reverse(self) -> bool:
        return bool(self.flag & REVERSE)

    @property
    def species(self) -> str:
        return species_of(self.contig)

    def identity(self) -> float:
        """One minus the edit distance over the alignment's columns."""
        columns = sum(int(n) for n, op in _CIGAR.findall(self.cigar) if op in _COLUMNS)
        if not columns:
            return 0.0
        return 1 - int(self.tags.get("NM", "0")) / columns

    def five_prime(self) -> int:
        """The reference position of the read's first base."""
        if not self.reverse:
            return self.position
        span = sum(int(n) for n, op in _CIGAR.findall(self.cigar) if op in "MDN=X")
        return self.position + span - 1

    def alternative_species(self) -> set[str]:
        """The species of its other hits: bwa-mem2's alternative hits
        (``XA:Z:contig,±pos,CIGAR,NM;``), and its secondary records'."""
        hits = self.tags.get("XA", "")
        found = {species_of(hit.split(",")[0]) for hit in hits.split(";") if hit}
        return found | self.alternatives

    def fastq(self, name: str) -> str:
        """The read as it was sequenced: reverse-complemented back when it aligned reversed."""
        sequence, quality = self.sequence, self.quality
        if self.reverse:
            sequence = sequence.translate(_COMPLEMENT)[::-1]
            quality = quality[::-1]
        if quality == "*":
            quality = "I" * len(sequence)
        return f"@{name}\n{sequence}\n+\n{quality}\n"


def species_of(contig: str) -> str:
    """A contig's species: its name up to the first underscore."""
    return contig.split("_", 1)[0]


def parse_record(line: str) -> Mate:
    fields = line.rstrip("\n").split("\t")
    tags = {}
    for tag in fields[11:]:
        key, _, value = tag.split(":", 2)
        tags[key] = value
    return Mate(
        name=fields[0],
        flag=int(fields[1]),
        contig=fields[2],
        position=int(fields[3]),
        mapq=int(fields[4]),
        cigar=fields[5],
        sequence=fields[9],
        quality=fields[10],
        tags=tags,
        line=line if line.endswith("\n") else line + "\n",
    )


def read_sam(
    handle: TextIO, header: list[str] | None = None
) -> tuple[list[str], Iterator[tuple[Mate, Mate]]]:
    """The species of the header's contigs, and each pair's two primary records.

    A read's records come together, as bwa-mem2 and STAR write them; its
    secondary records give the species of its other hits. ``header``, when
    given, gets the header's lines.

    Raises
    ------
    ValueError
        When a pair lacks a mate's primary record, or its records are not
        together.
    """
    species: list[str] = []
    first = None
    for line in handle:
        if line.startswith("@"):
            if header is not None:
                header.append(line)
            if line.startswith("@SQ"):
                name = next(f[3:] for f in line.split("\t") if f.startswith("SN:"))
                if (found := species_of(name.strip())) not in species:
                    species.append(found)
            continue
        first = line
        break

    def pair_of(primaries: list[Mate], others: set[str]) -> tuple[Mate, Mate]:
        if len(primaries) > 2:
            raise ValueError(f"{primaries[0].name}: more than two primary records")
        one, two = primaries
        alternatives = frozenset(others)
        one.alternatives = two.alternatives = alternatives
        return (one, two) if one.flag & FIRST else (two, one)

    def pairs() -> Iterator[tuple[Mate, Mate]]:
        name = None
        primaries: list[Mate] = []
        others: set[str] = set()
        lines = itertools.chain([first], handle) if first is not None else iter(())
        for line in lines:
            if not line.strip():
                continue
            record = parse_record(line)
            if record.name != name:
                if name is not None:
                    if len(primaries) < 2:
                        raise ValueError(
                            f"{name} and {record.name}: a pair's records are not together"
                        )
                    yield pair_of(primaries, others)
                name, primaries, others = record.name, [], set()
            if record.flag & SECONDARY:
                if record.mapped:
                    others.add(record.species)
            elif not record.flag & SUPPLEMENTARY:
                primaries.append(record)
        if name is not None:
            if len(primaries) < 2:
                raise ValueError(f"{name}: a read without its mate")
            yield pair_of(primaries, others)

    return species, pairs()


@dataclass
class Tally:
    """What the split did with a library's pairs."""

    pairs: int = 0
    kept: Counter = field(default_factory=Counter)
    holdout: Counter = field(default_factory=Counter)
    duplicates: Counter = field(default_factory=Counter)
    between: Counter = field(default_factory=Counter)
    #: Each pair of species an ambiguous pair involves; a pair between three
    #: species counts for each two of them.
    ambiguous: Counter = field(default_factory=Counter)
    ambiguous_pairs: int = 0
    unmapped: int = 0
    divergent: int = 0
    mates: Counter = field(default_factory=Counter)
    #: Hi-C: the unique pairs between two contigs, of a species or two.
    contacts: Counter = field(default_factory=Counter)


def verdict(mate: Mate, min_mapq: int, min_identity: float) -> str:
    """A mate's verdict: its species, or why it has none."""
    if not mate.mapped:
        return "unmapped"
    if mate.mapq < min_mapq:
        return "ambiguous"
    if mate.identity() < min_identity:
        return "divergent"
    return "assigned"


def held_out(name: str, share: float) -> bool:
    """Whether a pair is held out: a hash of its name, the same at every run."""
    return zlib.crc32(name.encode()) % 10_000 < share * 10_000


def duplicate_key(one: Mate, two: Mate) -> bytes:
    """A pair's place: both mates' contig, 5' position and strand, in order."""
    ends = sorted(
        (
            (one.contig, one.five_prime(), one.reverse),
            (two.contig, two.five_prime(), two.reverse),
        )
    )
    return hashlib.blake2b(repr(ends).encode(), digest_size=8).digest()


class Writers:
    """Each species' FASTQ files, opened when its first pair comes."""

    def __init__(self, outdir: Path, library: str, level: int) -> None:
        self.outdir, self.library, self.level = outdir, library, level
        self.handles: dict[tuple[str, str], IO[str]] = {}

    def path(self, species: str, part: str, read: int) -> Path:
        middle = f".{part}" if part else ""
        return self.outdir / f"{species}.{self.library}{middle}_R{read}.fastq.gz"

    def write(self, species: str, part: str, one: Mate, two: Mate) -> None:
        for read, mate in ((1, one), (2, two)):
            key = (f"{species}:{part}", str(read))
            if key not in self.handles:
                self.handles[key] = gzip.open(  # noqa: SIM115
                    self.path(species, part, read), "wt", compresslevel=self.level
                )
            self.handles[key].write(mate.fastq(one.name))

    def close(self) -> None:
        for handle in self.handles.values():
            handle.close()


class SamWriters:
    """Each species' alignments, as SAM of its own sequences, gzipped."""

    def __init__(self, outdir: Path, library: str, header: list[str], level: int):
        self.outdir, self.library, self.header, self.level = (
            outdir,
            library,
            header,
            level,
        )
        self.handles: dict[str, IO[str]] = {}

    def path(self, species: str) -> Path:
        return self.outdir / f"{species}.{self.library}.sam.gz"

    def write(self, species: str, one: Mate, two: Mate) -> None:
        if species not in self.handles:
            handle = gzip.open(  # noqa: SIM115
                self.path(species), "wt", compresslevel=self.level
            )
            for line in self.header:
                if line.startswith("@SQ"):
                    name = next(
                        f[3:] for f in line.split("\t") if f.startswith("SN:")
                    ).strip()
                    if species_of(name) != species:
                        continue
                handle.write(line)
            self.handles[species] = handle
        self.handles[species].write(one.line)
        self.handles[species].write(two.line)

    def close(self) -> None:
        for handle in self.handles.values():
            handle.close()


def split_pairs(
    sam: TextIO,
    outdir: Path,
    library: str,
    kind: str,
    min_mapq: int = 20,
    min_identity: float = 0.98,
    holdout: float = 0.1,
    dedup: bool = True,
    level: int = 1,
) -> tuple[list[str], Tally]:
    """Split a library's pairs between its pool's species; write each one's FASTQ.

    Parameters
    ----------
    sam
        bwa-mem2's alignments of the library to the pool's contigs, as SAM.
    outdir
        Where each species' FASTQ files go.
    library
        The library's ID, which names the files.
    kind
        ``wgs``, ``hic`` or ``rnaseq`` (STAR's alignments, written as SAM:
        no duplicates dropped, no pair held out).
    min_mapq, min_identity
        What makes a mate's alignment unique and near-identical.
    holdout
        The share of a WGS library's pairs kept apart, for the QV guard.
    dedup
        Drop a Hi-C library's duplicate pairs.
    level
        The gzip level of the FASTQ files.

    Returns
    -------
    tuple
        The pool's species, from the header, and the tally.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    header: list[str] = []
    species, pairs = read_sam(sam, header)
    tally = Tally()
    seen: set[bytes] = set()
    writers = Writers(outdir, library, level)
    alignments = SamWriters(outdir, library, header, level)
    if kind == "rnaseq":
        holdout, dedup = 0.0, False
    try:
        for one, two in pairs:
            tally.pairs += 1
            verdicts = (
                verdict(one, min_mapq, min_identity),
                verdict(two, min_mapq, min_identity),
            )
            if "unmapped" in verdicts:
                tally.unmapped += 1
                continue
            if "ambiguous" in verdicts:
                involved = set()
                for mate, said in zip((one, two), verdicts, strict=True):
                    involved.add(mate.species)
                    if said == "ambiguous":
                        involved |= mate.alternative_species()
                tally.ambiguous_pairs += 1
                names = sorted(involved)
                for key in (
                    itertools.combinations(names, 2)
                    if len(names) > 1
                    else [(names[0], names[0])]
                ):
                    tally.ambiguous[key] += 1
                continue
            if "divergent" in verdicts:
                tally.divergent += 1
                continue
            tally.mates[one.species] += 1
            tally.mates[two.species] += 1
            duplicate = False
            if kind == "hic" and dedup:
                key = duplicate_key(one, two)
                duplicate = key in seen
                seen.add(key)
            between = one.species != two.species
            if kind == "hic" and not duplicate and one.contig != two.contig:
                tally.contacts[tuple(sorted((one.contig, two.contig)))] += 1
            if between:
                tally.between[tuple(sorted((one.species, two.species)))] += 1
                continue
            name = one.species
            if duplicate:
                tally.duplicates[name] += 1
                continue
            if kind == "rnaseq":
                tally.kept[name] += 1
                alignments.write(name, one, two)
                continue
            if held_out(one.name, holdout):
                tally.holdout[name] += 1
                writers.write(name, "holdout", one, two)
                continue
            tally.kept[name] += 1
            writers.write(name, "", one, two)
    finally:
        writers.close()
        alignments.close()
    return species, tally


def epsilon(tally: Tally) -> float | None:
    """The random-ligation rate: the pairs between species over the share of
    random pairs that land in two species, 1 - sum of the squared shares."""
    within = (
        sum(tally.kept.values())
        + sum(tally.holdout.values())
        + sum(tally.duplicates.values())
    )
    between = sum(tally.between.values())
    mates = sum(tally.mates.values())
    if not mates or not within + between:
        return None
    shares = [count / mates for count in tally.mates.values()]
    random_between = 1 - sum(share**2 for share in shares)
    if random_between <= 0:
        return None
    return between / (within + between) / random_between


def summary(species: list[str], tally: Tally, kind: str) -> dict:
    """Every number of the split, for the JSON file and the tables."""
    kept = sum(tally.kept.values())
    numbers = {
        "pairs": tally.pairs,
        "kept": kept,
        "holdout": sum(tally.holdout.values()),
        "duplicates": sum(tally.duplicates.values()),
        "between_species": sum(tally.between.values()),
        "ambiguous": tally.ambiguous_pairs,
        "divergent": tally.divergent,
        "unmapped": tally.unmapped,
        "species": {
            name: {
                "pairs": tally.kept[name],
                "share": tally.kept[name] / kept if kept else 0.0,
                "holdout": tally.holdout[name],
                "duplicates": tally.duplicates[name],
            }
            for name in species
        },
        "ambiguous_pairs": [
            {"species_a": a, "species_b": b, "pairs": count}
            for (a, b), count in sorted(tally.ambiguous.items())
        ],
        "between_pairs": [
            {"species_a": a, "species_b": b, "pairs": count}
            for (a, b), count in sorted(tally.between.items())
        ],
    }
    if kind == "hic":
        numbers["epsilon"] = epsilon(tally)
    return numbers


def write_outputs(
    numbers: dict, outdir: Path, prefix: str, library: str
) -> dict[str, Path]:
    """The split's table, its JSON file and MultiQC's table."""
    paths = {
        "table": outdir / f"{prefix}.{library}.split.tsv",
        "json": outdir / f"{prefix}.{library}.split.json",
        "multiqc": outdir / f"{prefix}.{library}.split_mqc.tsv",
    }
    paths["json"].write_text(json.dumps(numbers, indent=2) + "\n", encoding="utf-8")
    with paths["table"].open("w", encoding="utf-8") as sink:
        sink.write("species\tpairs\tshare\tholdout\tduplicates\n")
        for name, row in numbers["species"].items():
            sink.write(
                f"{name}\t{row['pairs']}\t{row['share']:.4f}\t{row['holdout']}\t{row['duplicates']}\n"
            )
    lines = [
        "# id: 'plexembler_split'",
        "# section_name: 'Pooled libraries, split'",
        "# description: 'Each species'' pairs of each pooled library, and the pairs "
        "no species gets: ambiguous, divergent, between species, unmapped "
        "(plexembler split).'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'plexembler_split_table'",
        "#     namespace: 'plexembler'",
        "Sample\tPairs\tShare (%)\tHeld out\tDuplicates",
    ]
    for name, row in numbers["species"].items():
        lines.append(
            f"{library}_{name}\t{row['pairs']}\t{100 * row['share']:.2f}\t{row['holdout']}\t{row['duplicates']}"
        )
    paths["multiqc"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths


def open_sam(path: Path) -> TextIO:
    """A SAM file, gzipped or not; ``-`` for stdin."""
    if str(path) == "-":
        return sys.stdin
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def write_contacts(tally: Tally, outdir: Path, prefix: str, library: str) -> Path:
    """A Hi-C library's unique pairs between two contigs: ``contig_a``,
    ``contig_b``, ``pairs``, each pair of contigs once (the contig screen's
    Hi-C signal)."""
    path = outdir / f"{prefix}.{library}.contacts.tsv"
    with path.open("w", encoding="utf-8") as sink:
        sink.write("contig_a\tcontig_b\tpairs\n")
        for (a, b), count in sorted(tally.contacts.items()):
            sink.write(f"{a}\t{b}\t{count}\n")
    return path
