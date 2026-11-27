"""A species' expression tables: counts and TPM, genes by RNA-seq libraries.

featureCounts counts each species' read pairs of each pooled RNA-seq library
on its genes (``plexembler split --type rnaseq`` gave it only its own pairs).
``plexembler expression`` reads its table and writes the species' tables
(docs/plan.md, M5):

- ``<species>.counts.tsv``: each gene's pairs in each library, a column per
  library named ``<condition>_<replicate>`` (``--label``), after the gene's
  length (its exons', without overlaps, as featureCounts measures it);
- ``<species>.tpm.tsv``: transcripts per million, within the species: each
  gene's pairs over its length, as a share of the library's sum, times a
  million;
- a MultiQC table: each library's pairs on genes, and its genes detected.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class Table:
    """A species' counts: genes, their lengths, and a column per library."""

    genes: list[str]
    lengths: list[int]
    libraries: list[str]
    counts: list[list[int]]

    def tpm(self) -> list[list[float]]:
        """Each gene's TPM in each library."""
        columns = []
        for index in range(len(self.libraries)):
            rates = [
                row[index] / length if length else 0.0
                for row, length in zip(self.counts, self.lengths, strict=True)
            ]
            total = sum(rates)
            columns.append([1e6 * rate / total if total else 0.0 for rate in rates])
        return [list(row) for row in zip(*columns, strict=True)] if columns else []


def library_of(column: str, species: str) -> str:
    """A featureCounts column's library: ``<species>.<library>.bam``, its path
    left out."""
    name = column.rsplit("/", 1)[-1]
    if name.endswith(".bam"):
        name = name[: -len(".bam")]
    if name.startswith(f"{species}."):
        name = name[len(species) + 1 :]
    return name


def read_featurecounts(path: Path, species: str, labels: dict[str, str]) -> Table:
    """featureCounts' table: its genes, their lengths and each library's counts.

    Raises
    ------
    ValueError
        When the table is not featureCounts', or two libraries share a label.
    """
    with path.open(encoding="utf-8") as handle:
        lines = [
            line.rstrip("\n").split("\t") for line in handle if not line.startswith("#")
        ]
    if not lines or lines[0][:6] != [
        "Geneid",
        "Chr",
        "Start",
        "End",
        "Strand",
        "Length",
    ]:
        raise ValueError(f"{path}: not a featureCounts table")
    libraries = [library_of(column, species) for column in lines[0][6:]]
    named = [labels.get(library, library) for library in libraries]
    if len(set(named)) != len(named):
        raise ValueError(f"{path}: two libraries named the same: {', '.join(named)}")
    genes, lengths, counts = [], [], []
    for row in lines[1:]:
        genes.append(row[0])
        lengths.append(int(row[5]))
        counts.append([int(value) for value in row[6:]])
    return Table(genes, lengths, named, counts)


def write_outputs(
    table: Table, outdir: Path, species: str, min_count: int = 10
) -> None:
    """The counts, the TPM and the MultiQC table."""
    outdir.mkdir(parents=True, exist_ok=True)
    header = "\t".join(table.libraries)
    with (outdir / f"{species}.counts.tsv").open("w", encoding="utf-8") as sink:
        sink.write(f"gene\tlength_bp\t{header}\n")
        for gene, length, row in zip(
            table.genes, table.lengths, table.counts, strict=True
        ):
            sink.write(f"{gene}\t{length}\t" + "\t".join(map(str, row)) + "\n")
    with (outdir / f"{species}.tpm.tsv").open("w", encoding="utf-8") as sink:
        sink.write(f"gene\t{header}\n")
        for gene, row in zip(table.genes, table.tpm(), strict=True):
            sink.write(f"{gene}\t" + "\t".join(f"{v:.3f}" for v in row) + "\n")
    lines = [
        "# id: plexembler_expression",
        "# section_name: Expression",
        "# description: each species' read pairs on its genes in each pooled RNA-seq "
        f"library, and its genes with {min_count} pairs or more (plexembler expression)",
        "# plot_type: table",
        "Sample\tPairs on genes\tGenes detected\tGenes",
    ]
    for index, library in enumerate(table.libraries):
        column = [row[index] for row in table.counts]
        detected = sum(value >= min_count for value in column)
        lines.append(f"{species}_{library}\t{sum(column)}\t{detected}\t{len(column)}")
    (outdir / f"{species}.expression_mqc.tsv").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    for index, library in enumerate(table.libraries):
        logger.info(
            "%s: %d pairs on %d genes",
            library,
            sum(row[index] for row in table.counts),
            len(table.genes),
        )
