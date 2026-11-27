"""How alike the species of a pool are: skani's identities, as a matrix.

The pooled short reads are split between the species by mapping
(docs/concepts.md): two species whose genomes are too alike share reads that
no alignment can give to one of them. ``plexembler ani`` reads skani's pairs
(``skani triangle --sparse``), names each genome by its contigs' species, and
writes the pool's matrix and its verdict on each pair:

- above ``--warn-above`` (95 % by default), the pair is too close for the pooled
  reads to be split well: a warning;
- from ``--genus-above`` (85 %) to that, a pair of the same genus: its reads are
  split by identity as well as uniqueness (docs/plan.md, M2);
- below, or below skani's floor (about 80 %, where it reports no pair), the two
  species are told apart by their reads' alignments alone.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

#: The columns of skani's sparse output that this module reads.
_SKANI_COLUMNS = (
    "Ref_file",
    "Query_file",
    "ANI",
    "Align_fraction_ref",
    "Align_fraction_query",
    "Ref_name",
    "Query_name",
)


@dataclass(frozen=True)
class Pair:
    """Two species of a pool, and skani's identity between their genomes."""

    first: str
    second: str
    ani: float
    fraction_first: float
    fraction_second: float


def species_of(name: str, path: str) -> str:
    """A genome's species: the start of its contigs' names (``<species>_…``).

    A name without an underscore falls back on the file's name, its extensions
    removed (``Scer.contigs.fa.gz`` is ``Scer``).
    """
    first = name.split()[0] if name.split() else ""
    if "_" in first:
        return first.split("_", 1)[0]
    return Path(path).name.split(".", 1)[0]


def read_skani(table: Path) -> list[Pair]:
    """The pairs of ``skani triangle --sparse``.

    Raises
    ------
    ValueError
        When a column is missing, or a pair is of one species with itself.
    """
    pairs: list[Pair] = []
    with table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(_SKANI_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{table}: no column {', '.join(sorted(missing))}")
        for row in reader:
            first = species_of(row["Ref_name"], row["Ref_file"])
            second = species_of(row["Query_name"], row["Query_file"])
            if first == second:
                raise ValueError(f"{table}: two genomes of {first!r}")
            pairs.append(
                Pair(
                    first,
                    second,
                    float(row["ANI"]),
                    float(row["Align_fraction_ref"]),
                    float(row["Align_fraction_query"]),
                )
            )
    return pairs


def verdict(ani: float | None, warn_above: float, genus_above: float) -> str:
    """What a pair's identity means for the split of the pool's reads."""
    if ani is None:
        return "distinct"
    if ani >= warn_above:
        return "too close"
    if ani >= genus_above:
        return "same genus"
    return "distinct"


def ani_matrix(
    pool: str,
    pairs: list[Pair],
    species: list[str],
    outdir: Path,
    warn_above: float = 95.0,
    genus_above: float = 85.0,
) -> dict[str, Path]:
    """Write the pool's matrix, its pairs and MultiQC's heatmap; log the verdicts.

    Parameters
    ----------
    pool
        The pool's ID, which names the outputs.
    pairs
        skani's pairs.
    species
        Every species of the pool, whether skani paired it or not.
    outdir
        Where the outputs go.
    warn_above, genus_above
        The identities, in percent, of a pair too close and of a pair of the
        same genus.

    Returns
    -------
    dict
        The paths written: ``matrix``, ``pairs`` and ``multiqc``.
    """
    every = sorted(set(species) | {p.first for p in pairs} | {p.second for p in pairs})
    identity: dict[tuple[str, str], Pair] = {}
    for pair in pairs:
        identity[(pair.first, pair.second)] = pair
        identity[(pair.second, pair.first)] = pair

    def cell(a: str, b: str) -> str:
        if a == b:
            return "100.00"
        pair = identity.get((a, b))
        return f"{pair.ani:.2f}" if pair else "NA"

    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "matrix": outdir / f"{pool}.ani.tsv",
        "pairs": outdir / f"{pool}.ani_pairs.tsv",
        "multiqc": outdir / f"{pool}.ani_mqc.tsv",
    }
    with paths["matrix"].open("w", encoding="utf-8") as sink:
        sink.write("\t".join(["species", *every]) + "\n")
        for a in every:
            sink.write("\t".join([a, *(cell(a, b) for b in every)]) + "\n")

    rows = []
    for index, a in enumerate(every):
        for b in every[index + 1 :]:
            pair = identity.get((a, b))
            if pair is None:
                fractions: tuple[float | None, float | None] = (None, None)
            elif pair.first == a:
                fractions = (pair.fraction_first, pair.fraction_second)
            else:
                fractions = (pair.fraction_second, pair.fraction_first)
            ani = pair.ani if pair else None
            rows.append((a, b, ani, *fractions, verdict(ani, warn_above, genus_above)))

    def number(value: float | None) -> str:
        return f"{value:.2f}" if value is not None else "NA"

    with paths["pairs"].open("w", encoding="utf-8") as sink:
        sink.write(
            "species_a\tspecies_b\tani\taligned_fraction_a\taligned_fraction_b\tverdict\n"
        )
        for a, b, ani, fraction_a, fraction_b, judged in rows:
            cells = [a, b, number(ani), number(fraction_a), number(fraction_b), judged]
            sink.write("\t".join(cells) + "\n")
            if judged == "too close":
                logger.warning(
                    "%s and %s are %.2f %% identical: above %.0f %%, the pool's "
                    "short reads cannot be split well between them",
                    a,
                    b,
                    ani,
                    warn_above,
                )
            elif judged == "same genus":
                logger.info(
                    "%s and %s are %.2f %% identical: a pair of the same genus",
                    a,
                    b,
                    ani,
                )

    # One table of every pool's pairs: each pool's file has the same section ID.
    lines = [
        "# id: 'plexembler_ani'",
        "# section_name: 'Identity between species'",
        "# description: 'The average nucleotide identity (skani) of each pair of "
        "species of a pool, in percent; NA below skani''s floor, about 80 %. Above "
        f"{warn_above:.0f} %, the pool''s short reads cannot be split well "
        "between the two species (plexembler ani).'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'plexembler_ani_table'",
        "#     namespace: 'plexembler'",
        "Sample\tPool\tANI (%)\tAligned fraction, first (%)\tAligned fraction, second (%)\tVerdict",
    ]
    for a, b, ani, fraction_a, fraction_b, judged in rows:
        cells = [f"{a}_vs_{b}", pool, number(ani), number(fraction_a)]
        cells += [number(fraction_b), judged]
        lines.append("\t".join(cells))
    paths["multiqc"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info("Wrote the identities of %d species of %s", len(every), pool)
    return paths
