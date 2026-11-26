"""The contig screen: each contig of a pool kept, moved to another species, or dropped.

A species' nanopore run can hold reads of another species of the pool (a
barcode misassigned: bleed-through); assembled, they make contigs of the wrong
species. ``plexembler screen`` weighs three signals for each contig of a pool
(docs/plan.md, M2):

- **the nanopore depth of each species' reads** on it (minimap2's map-ont
  alignments of every species' reads to every species' contigs), against the
  depth of each species' reads on its own contigs;
- **its identity to the other species' contigs** (minimap2's asm alignments of
  the pool's contigs to each other): the share of it found, 99 % identical or
  more, in another species' contigs;
- **its Hi-C contacts** with each species' other contigs (``plexembler split``'s
  unique pairs between contigs).

And gives a verdict:

- ``organelle``: much deeper than its species' contigs (``--organelle-depth``
  times) and AT-rich (GC under ``--organelle-gc``): its mitochondrion, kept as
  ``<species>_chrM``;
- ``drop``: most of it (``--duplicate-share``) found in another species' contigs,
  and its own species' reads thin on it (under ``--own-depth`` of their depth):
  that species' DNA, which it holds already;
- ``move``: its Hi-C contacts with another species' contigs more than
  ``--hic-ratio`` times those with its own (and ``--min-contacts`` of them); or
  its own species' reads thin on it while another species' reads cover it
  (``--other-depth`` of their depth): that species' DNA, which it lacks. It is
  renamed ``<that species>_<contig>``;
- ``keep``: the others.
"""

from __future__ import annotations

import csv
import gzip
import logging
import statistics
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO, cast

from .rename import LINE_WIDTH, read_fasta
from .split import species_of

logger = logging.getLogger(__name__)


@dataclass
class Contig:
    """A contig of the pool and its signals."""

    name: str
    sequence: str
    depth: dict[str, float] = field(default_factory=dict)
    found_in: dict[str, float] = field(default_factory=dict)
    contacts: dict[str, int] = field(default_factory=dict)
    verdict: str = "keep"
    new_name: str = ""
    reason: str = ""

    @property
    def species(self) -> str:
        return species_of(self.name)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc(self) -> float:
        bases = sum(self.sequence.count(b) for b in "ACGTacgt")
        if not bases:
            return 0.0
        return sum(self.sequence.count(b) for b in "GCgc") / bases


@dataclass(frozen=True)
class Thresholds:
    """What each verdict needs."""

    organelle_depth: float = 3.0
    organelle_gc: float = 0.3
    duplicate_share: float = 0.9
    duplicate_identity: float = 0.99
    own_depth: float = 0.2
    other_depth: float = 0.5
    hic_ratio: float = 2.0
    min_contacts: int = 20


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, "rt", encoding="utf-8"))
    return cast(TextIO, path.open(encoding="utf-8"))


def read_contigs(path: Path) -> dict[str, Contig]:
    with open_text(path) as handle:
        return {name: Contig(name, sequence) for name, sequence in read_fasta(handle)}


def merged_length(intervals: list[tuple[int, int]]) -> int:
    """The bases covered by intervals, each base once."""
    total, end = 0, -1
    for start, stop in sorted(intervals):
        if stop <= end:
            continue
        total += stop - max(start, end)
        end = stop
    return total


def read_depths(paf: Path, species: str, contigs: dict[str, Contig]) -> None:
    """A species' reads' depth on each contig: their primary alignments' bases
    over the contig's length."""
    bases: dict[str, int] = defaultdict(int)
    with open_text(paf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or "tp:A:S" in fields[12:]:
                continue
            target = fields[5]
            if target in contigs:
                bases[target] += int(fields[8]) - int(fields[7])
    for name, contig in contigs.items():
        contig.depth[species] = bases[name] / contig.length if contig.length else 0.0


def read_identity(paf: Path, contigs: dict[str, Contig], identity: float) -> None:
    """The share of each contig found in each other species' contigs, at
    ``identity`` or more: minimap2's alignments of the pool's contigs to each
    other (``-X``: no self-alignments)."""
    intervals: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with open_text(paf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query, target = fields[0], fields[5]
            if query not in contigs or species_of(query) == species_of(target):
                continue
            block = int(fields[10])
            if not block or int(fields[9]) / block < identity:
                continue
            intervals[(query, species_of(target))].append(
                (int(fields[2]), int(fields[3]))
            )
    for (query, other), found in intervals.items():
        contig = contigs[query]
        contig.found_in[other] = (
            merged_length(found) / contig.length if contig.length else 0.0
        )


def read_contacts(table: Path, contigs: dict[str, Contig]) -> None:
    """Each contig's Hi-C contacts with each species' contigs (its own others)."""
    with open_text(table) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            a, b, pairs = row["contig_a"], row["contig_b"], int(row["pairs"])
            for this, other in ((a, b), (b, a)):
                if this in contigs:
                    species = species_of(other)
                    contacts = contigs[this].contacts
                    contacts[species] = contacts.get(species, 0) + pairs


def typical_depths(contigs: dict[str, Contig]) -> dict[str, float]:
    """Each species' reads' depth on its own contigs: the median, by length."""
    typical = {}
    by_species: dict[str, list[Contig]] = defaultdict(list)
    for contig in contigs.values():
        by_species[contig.species].append(contig)
    for species, own in by_species.items():
        depths = []
        for contig in own:
            depths += [contig.depth.get(species, 0.0)] * max(1, contig.length // 1000)
        typical[species] = statistics.median(depths) if depths else 0.0
    return typical


def judge(contigs: dict[str, Contig], thresholds: Thresholds) -> None:
    """Each contig's verdict, its new name and the reason."""
    typical = typical_depths(contigs)
    for contig in contigs.values():
        own = contig.species
        own_typical = typical.get(own, 0.0)
        own_ratio = contig.depth.get(own, 0.0) / own_typical if own_typical else 1.0
        contig.new_name = contig.name
        if (
            own_ratio >= thresholds.organelle_depth
            and contig.gc < thresholds.organelle_gc
        ):
            contig.verdict = "organelle"
            contig.reason = (
                f"{own_ratio:.1f} times its species' depth, {100 * contig.gc:.0f} % GC"
            )
            continue
        duplicates = {
            other: share
            for other, share in contig.found_in.items()
            if share >= thresholds.duplicate_share
        }
        if duplicates and own_ratio < thresholds.own_depth:
            other = max(duplicates, key=lambda name: duplicates[name])
            contig.verdict = "drop"
            contig.reason = (
                f"{100 * duplicates[other]:.0f} % of it in {other}'s contigs, "
                f"{own}'s reads at {own_ratio:.2f} of their depth"
            )
            continue
        own_contacts = contig.contacts.get(own, 0)
        others = {k: v for k, v in contig.contacts.items() if k != own}
        if others:
            other = max(others, key=lambda name: others[name])
            count = others[other]
            if (
                count >= thresholds.min_contacts
                and count > thresholds.hic_ratio * own_contacts
            ):
                move(
                    contig,
                    other,
                    f"{count} Hi-C contacts with {other}, {own_contacts} with {own}",
                )
                continue
        if own_ratio < thresholds.own_depth:
            covered = {
                other: contig.depth.get(other, 0.0) / typical[other]
                for other in typical
                if other != own and typical[other]
            }
            if covered:
                other = max(covered, key=lambda name: covered[name])
                if covered[other] >= thresholds.other_depth:
                    move(
                        contig,
                        other,
                        f"{own}'s reads at {own_ratio:.2f} of their depth, {other}'s at {covered[other]:.2f}",
                    )
                    continue
        contig.verdict = "keep"


def move(contig: Contig, other: str, reason: str) -> None:
    contig.verdict = "move"
    contig.new_name = f"{other}_{contig.name.split('_', 1)[1]}"
    contig.reason = reason


def name_organelles(contigs: dict[str, Contig]) -> None:
    """A species' organelle contigs, named ``<species>_chrM`` (then ``chrM2``, …),
    the longest first."""
    by_species: dict[str, list[Contig]] = defaultdict(list)
    for contig in contigs.values():
        if contig.verdict == "organelle":
            by_species[contig.species].append(contig)
    for species, organelles in by_species.items():
        for index, contig in enumerate(sorted(organelles, key=lambda c: -c.length)):
            suffix = "" if index == 0 else str(index + 1)
            contig.new_name = f"{species}_chrM{suffix}"


def write_outputs(
    contigs: dict[str, Contig], species: list[str], outdir: Path, prefix: str
) -> dict[str, Path]:
    """Each species' screened contigs, the verdicts' table and MultiQC's."""
    outdir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    every = sorted(set(species) | {c.species for c in contigs.values()})
    for name in every:
        path = outdir / f"{name}.screened.fa.gz"
        paths[name] = path
        with gzip.open(path, "wt", encoding="ascii") as sink:
            for contig in sorted(contigs.values(), key=lambda c: c.new_name):
                if contig.verdict == "drop" or species_of(contig.new_name) != name:
                    continue
                sink.write(f">{contig.new_name}\n")
                for start in range(0, contig.length, LINE_WIDTH):
                    sink.write(contig.sequence[start : start + LINE_WIDTH] + "\n")
    paths["table"] = outdir / f"{prefix}.screen.tsv"
    with paths["table"].open("w", encoding="utf-8") as sink:
        header = ["contig", "species", "length_bp", "gc_percent"]
        header += [f"depth_{s}" for s in every] + [f"found_in_{s}" for s in every]
        header += [f"contacts_{s}" for s in every] + ["verdict", "new_name", "reason"]
        sink.write("\t".join(header) + "\n")
        for contig in sorted(contigs.values(), key=lambda c: c.name):
            cells = [
                contig.name,
                contig.species,
                str(contig.length),
                f"{100 * contig.gc:.1f}",
            ]
            cells += [f"{contig.depth.get(s, 0.0):.1f}" for s in every]
            cells += [f"{contig.found_in.get(s, 0.0):.2f}" for s in every]
            cells += [str(contig.contacts.get(s, 0)) for s in every]
            cells += [contig.verdict, contig.new_name, contig.reason]
            sink.write("\t".join(cells) + "\n")
    paths["multiqc"] = outdir / f"{prefix}.screen_mqc.tsv"
    lines = [
        "# id: 'plexembler_screen'",
        "# section_name: 'Contig screen'",
        "# description: 'Each species'' contigs: kept, moved to another species, "
        "dropped, or its organelle (plexembler screen).'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'plexembler_screen_table'",
        "#     namespace: 'plexembler'",
        "Sample\tKept\tMoved in\tMoved out\tDropped\tOrganelle",
    ]
    for name in every:
        own = [c for c in contigs.values() if c.species == name]
        moved_in = [
            c
            for c in contigs.values()
            if c.verdict == "move" and species_of(c.new_name) == name
        ]
        lines.append(
            "\t".join(
                [
                    name,
                    str(sum(c.verdict == "keep" for c in own)),
                    str(len(moved_in)),
                    str(sum(c.verdict == "move" for c in own)),
                    str(sum(c.verdict == "drop" for c in own)),
                    str(sum(c.verdict == "organelle" for c in own)),
                ]
            )
        )
    paths["multiqc"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths


def screen(
    fasta: Path,
    reads: dict[str, Path],
    identity_paf: Path | None,
    contacts: list[Path],
    outdir: Path,
    prefix: str,
    thresholds: Thresholds | None = None,
) -> dict[str, Contig]:
    """Screen a pool's contigs; write each species' and the verdicts.

    Parameters
    ----------
    fasta
        Every species' contigs of the pool, named ``<species>_<contig>``.
    reads
        Each species' nanopore reads aligned to them (PAF), by species.
    identity_paf
        The contigs aligned to each other (PAF), or None.
    contacts
        ``plexembler split``'s contacts between contigs, one per Hi-C library.
    outdir, prefix
        Where the outputs go, and the pool's ID, which names the table.
    thresholds
        What each verdict needs.
    """
    thresholds = thresholds or Thresholds()
    contigs = read_contigs(fasta)
    for species, paf in reads.items():
        read_depths(paf, species, contigs)
    if identity_paf is not None:
        read_identity(identity_paf, contigs, thresholds.duplicate_identity)
    for table in contacts:
        read_contacts(table, contigs)
    judge(contigs, thresholds)
    name_organelles(contigs)
    write_outputs(contigs, list(reads), outdir, prefix)
    for contig in contigs.values():
        if contig.verdict != "keep":
            logger.info(
                "%s: %s (%s) %s",
                contig.name,
                contig.verdict,
                contig.new_name,
                contig.reason,
            )
    return contigs
