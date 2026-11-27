"""A species' report: what the pipeline made of it, on chromokit's report kit.

``plexembler report`` reads the files the pipeline wrote for a species and
writes its report (docs/plan.md, M6): ``<species>.report.txt``, every number,
for scripts; ``<species>.report_mqc.yaml``, its row of MultiQC's General
Statistics table; and ``<species>.report.html``, its page. Each file is
known by its name, as the pipeline writes it:

========================================  =====================================
File                                      What the report takes from it
========================================  =====================================
``<species>.assembly.selection.tsv``      the assembly's candidates and choice
``<species>.polish.selection.tsv``        the short-read polishing's choice
``<pool>.<library>.split.tsv``            the species' share of each library
``<pool>.screen.tsv``                     the screen's verdicts on its contigs
``<species>.names.tsv``                   its chromosomes, organelle, unplaced
``short_summary.*.<species>.fa.txt``      BUSCO on its genome
``<species>.genome.<k-mers>.qv``          Merqury's QV of its genome
``<species>.junctions.tsv``               its scaffolds' joins
``<species>.genes.tsv``                   its genes
``short_summary.*.<species>.proteins*``   BUSCO on its proteins
``<species>.counts.tsv``                  its pairs on genes, by library
``<species>.closest.tsv``                 its closest known genomes
``<species>.dotplot.tsv``                 its alignment to the closest one
``<species>.dotplot.json``                the dotplot, a plotly figure
========================================  =====================================

A file it is not given leaves its numbers out (``NA``): a species of a pool
without Hi-C has no joins, without RNA-seq no expression.
"""

from __future__ import annotations

import csv
import json
import logging
import math
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

from . import __version__
from ._report import Report, Tool, report_section, write_report
from .select import read_busco, read_qv

logger = logging.getLogger(__name__)

TOOL = Tool(
    name="plexembler",
    version=__version__,
    format_version=1,
    columns=(
        ("length_mb", "genome", {"title": "Genome (Mb)", "description": "The genome's length", "scale": "Blues", "format": "{:,.2f}"}),
        ("chromosomes", "genome", {"title": "Chromosomes", "description": "Its chromosomes", "format": "{:,.0f}"}),
        ("busco_complete", "quality", {"title": "BUSCO (%)", "description": "BUSCO's complete genes in the genome", "max": 100, "suffix": "%", "format": "{:,.1f}"}),
        ("qv", "quality", {"title": "QV", "description": "Merqury's QV of the genome", "format": "{:,.1f}"}),
        ("genes", "annotation", {"title": "Genes", "description": "Its genes", "format": "{:,.0f}"}),
        ("proteins_busco_complete", "annotation", {"title": "BUSCO, proteins (%)", "description": "BUSCO's complete genes in its proteins", "max": 100, "suffix": "%", "format": "{:,.1f}"}),
    ),
    id="plexembler_report",
    description="each species' genome, its quality, genes and expression (plexembler report)",
    docs="https://js2264.github.io/plexembler/formats/#the-report-of-a-species",
)  # fmt: skip


@dataclass
class Rows:
    """A table, as the report kit reads one: its columns and its rows."""

    columns: list[str]
    rows: list[tuple] = field(default_factory=list)

    def itertuples(self, index: bool = False) -> Iterator[tuple]:
        return iter(self.rows)

    def to_dict(self, orient: str = "records") -> list[dict]:
        return [dict(zip(self.columns, row, strict=True)) for row in self.rows]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str | None) -> float | None:
    """A cell as a number: ``None`` for ``NA`` or an empty one."""
    if value in (None, "", "NA"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def finite(value: float | None) -> float | None:
    """A number the report can write: ``None`` for a missing or infinite one."""
    return None if value is None or math.isinf(value) else value


@dataclass
class Inputs:
    """A species' files, by what they hold."""

    assembly: Path | None = None
    polish: Path | None = None
    splits: list[Path] = field(default_factory=list)
    screen: Path | None = None
    names: Path | None = None
    busco: Path | None = None
    qv: Path | None = None
    junctions: Path | None = None
    genes: Path | None = None
    proteins_busco: Path | None = None
    counts: Path | None = None
    closest: Path | None = None
    dotplot: Path | None = None
    dotplot_figure: Path | None = None


def sort_files(files: list[Path], species: str) -> Inputs:
    """A species' files, each known by its name; the others, left out."""
    found = Inputs()
    escaped = re.escape(species)
    for path in files:
        name = path.name
        if name == f"{species}.assembly.selection.tsv":
            found.assembly = path
        elif name == f"{species}.polish.selection.tsv":
            found.polish = path
        elif name.endswith(".split.tsv"):
            found.splits.append(path)
        elif name.endswith(".screen.tsv"):
            found.screen = path
        elif name == f"{species}.names.tsv":
            found.names = path
        elif re.fullmatch(rf"short_summary\..*\.{escaped}\.fa\.txt", name):
            found.busco = path
        elif re.fullmatch(rf"short_summary\..*\.{escaped}\.proteins.*\.txt", name):
            found.proteins_busco = path
        elif re.fullmatch(rf"{escaped}\.genome\.\w+\.qv", name):
            found.qv = path
        elif name == f"{species}.junctions.tsv":
            found.junctions = path
        elif name == f"{species}.genes.tsv":
            found.genes = path
        elif name == f"{species}.counts.tsv":
            found.counts = path
        elif name == f"{species}.closest.tsv":
            found.closest = path
        elif name == f"{species}.dotplot.tsv":
            found.dotplot = path
        elif name == f"{species}.dotplot.json":
            found.dotplot_figure = path
        else:
            logger.debug("%s: not a file the report reads", name)
    found.splits.sort()
    return found


def build(species: str, pool: str | None, files: list[Path]) -> Report:
    """A species' report, from its files."""
    inputs = sort_files(files, species)
    report = Report()
    report.sections["report"] = report_section(
        TOOL,
        species,
        {"pool": pool or "NA", "files": len(files)},
    )

    if inputs.assembly is not None:
        rows = read_tsv(inputs.assembly)
        kept = next(
            (r for r in rows if r.get("kept") == "yes"), rows[0] if rows else {}
        )
        report.sections["assembly"] = {
            "candidates": len(rows),
            "kept": kept.get("candidate"),
            "reason": kept.get("reason"),
            "contigs": number(kept.get("contigs")),
            "contig_n50": number(kept.get("contig_n50")),
        }
        report.tables["candidates"] = Rows(
            [
                "candidate",
                "kept",
                "busco_complete",
                "qv",
                "contigs",
                "contig_n50",
                "total_length",
            ],
            [
                (
                    r["candidate"],
                    r["kept"],
                    number(r.get("busco_complete")),
                    finite(number(r.get("qv"))),
                    number(r.get("contigs")),
                    number(r.get("contig_n50")),
                    number(r.get("total_length")),
                )
                for r in rows
            ],
        )
    if inputs.polish is not None:
        rows = read_tsv(inputs.polish)
        kept = next((r for r in rows if r.get("kept") == "yes"), {})
        report.sections["polishing"] = {
            "kept": kept.get("candidate"),
            "qv": finite(number(kept.get("qv"))),
            "reason": kept.get("reason"),
        }

    shares = []
    for path in inputs.splits:
        library = path.name.removesuffix(".split.tsv")
        if pool and library.startswith(f"{pool}."):
            library = library[len(pool) + 1 :]
        for row in read_tsv(path):
            if row.get("species") == species:
                shares.append((library, int(row["pairs"]), number(row["share"])))
    if shares:
        report.tables["libraries"] = Rows(["library", "pairs", "share"], shares)

    if inputs.screen is not None:
        verdicts: dict[str, int] = {}
        for row in read_tsv(inputs.screen):
            if row.get("species") == species:
                verdicts[row["verdict"]] = verdicts.get(row["verdict"], 0) + 1
        report.sections["screen"] = {
            verdict: verdicts.get(verdict, 0)
            for verdict in ("keep", "move", "drop", "organelle")
        }

    genome: dict = {
        "length_bp": None,
        "length_mb": None,
        "chromosomes": None,
        "organelles": None,
        "unplaced": None,
    }
    if inputs.names is not None:
        rows = read_tsv(inputs.names)
        length = sum(int(r["length_bp"]) for r in rows)
        genome = {
            "length_bp": length,
            "length_mb": round(length / 1e6, 3),
            "chromosomes": sum(r["role"] == "chromosome" for r in rows),
            "chromosomes_bp": sum(
                int(r["length_bp"]) for r in rows if r["role"] == "chromosome"
            ),
            "organelles": sum(r["role"] == "organelle" for r in rows),
            "unplaced": sum(r["role"] == "unplaced" for r in rows),
        }
        report.tables["chromosomes"] = Rows(
            ["name", "length_bp", "role", "original_name"],
            [
                (r["name"], int(r["length_bp"]), r["role"], r["original_name"])
                for r in rows
            ],
        )
    report.sections["genome"] = genome

    quality: dict = {"busco_complete": None, "qv": None}
    if inputs.busco is not None:
        busco = read_busco(inputs.busco)
        quality.update(
            busco_complete=busco.complete,
            busco_single=busco.single,
            busco_duplicated=busco.duplicated,
            busco_fragmented=busco.fragmented,
            busco_missing=busco.missing,
            busco_n=busco.n,
        )
    if inputs.qv is not None:
        quality["qv"] = finite(read_qv(inputs.qv))
    report.sections["quality"] = quality

    if inputs.junctions is not None:
        joins = [r["verdict"] for r in read_tsv(inputs.junctions)]
        report.sections["scaffolding"] = {
            "joins": len(joins),
            "supported": joins.count("supported"),
            "unsupported": joins.count("unsupported"),
            "thin": joins.count("thin"),
        }

    annotation: dict = {"genes": None, "proteins_busco_complete": None}
    if inputs.genes is not None:
        rows = read_tsv(inputs.genes)
        annotation["genes"] = len(rows)
        annotation["transcripts"] = sum(int(r["transcripts"]) for r in rows)
        annotation["genes_with_introns"] = sum(int(r["exons"]) > 1 for r in rows)
    if inputs.proteins_busco is not None:
        busco = read_busco(inputs.proteins_busco)
        annotation["proteins_busco_complete"] = busco.complete
    report.sections["annotation"] = annotation

    if inputs.counts is not None:
        with inputs.counts.open(encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            totals = [0] * (len(header) - 2)
            for line in handle:
                cells = line.rstrip("\n").split("\t")
                for index, value in enumerate(cells[2:]):
                    totals[index] += int(value)
        report.sections["expression"] = {
            "libraries": len(totals),
            "pairs_on_genes": sum(totals),
        }
        report.tables["expression"] = Rows(
            ["library", "pairs_on_genes"], list(zip(header[2:], totals, strict=True))
        )

    if inputs.closest is not None:
        closest_genome(report, inputs.closest, inputs.dotplot)
    if inputs.dotplot_figure is not None:
        report.plots["plot-dotplot"] = json.loads(
            inputs.dotplot_figure.read_text(encoding="utf-8")
        )
    return report


def closest_genome(report: Report, closest: Path, dotplot: Path | None) -> None:
    """The species' closest known genome (plexembler closest) and, when its
    sequence could be downloaded, how the two genomes align (plexembler dotplot)."""
    rows = read_tsv(closest)
    best = rows[0] if rows else {}
    accession = best.get("accession")
    name = best.get("name") or None
    if name and accession and accession != "NA":
        name = name.removeprefix(accession).strip() or name
    section: dict = {
        "matches": len(rows),
        "accession": None if accession in (None, "NA") else accession,
        "name": name,
        "ani": number(best.get("ani")),
        "species_containment": number(best.get("species_containment")),
        "match_containment": number(best.get("match_containment")),
        "lineage": None if best.get("lineage") in (None, "", "NA") else best["lineage"],
        "sequence": best.get("sequence"),
        "species_aligned": None,
        "reference_aligned": None,
        "identity": None,
    }
    report.tables["closest_genomes"] = Rows(
        ["rank", "accession", "name", "ani", "species_containment", "match_containment", "lineage"],
        [
            (
                int(row["rank"]),
                row["accession"],
                row["name"],
                number(row.get("ani")),
                number(row.get("species_containment")),
                number(row.get("match_containment")),
                row.get("lineage"),
            )
            for row in rows
        ],
    )  # fmt: skip
    if dotplot is not None:
        aligned = read_tsv(dotplot)
        for genome, key in (("species", "species_aligned"), ("reference", "reference_aligned")):  # fmt: skip
            mine = [row for row in aligned if row["genome"] == genome]
            length = sum(int(row["length_bp"]) for row in mine)
            if length:
                section[key] = round(
                    100 * sum(int(row["aligned_bp"]) for row in mine) / length, 2
                )
        weighed = [
            (int(row["aligned_bp"]), number(row["identity"]))
            for row in aligned
            if row["genome"] == "species" and number(row["identity"]) is not None
        ]
        bases = sum(count for count, _identity in weighed)
        if bases:
            section["identity"] = round(
                sum(count * (identity or 0.0) for count, identity in weighed) / bases, 2
            )
        report.tables["alignment_to_closest"] = Rows(
            # report.txt writes a _fraction as a fraction: this one is in percent.
            ["genome", "sequence", "label", "length_bp", "aligned_pct", "identity", "best_match", "best_match_share", "strand"],
            [
                (
                    row["genome"],
                    row["sequence"],
                    row["label"],
                    int(row["length_bp"]),
                    number(row["aligned_fraction"]),
                    number(row["identity"]),
                    row["best_match"],
                    number(row["best_match_share"]),
                    row["strand"],
                )
                for row in aligned
            ],
        )  # fmt: skip
    report.sections["closest_genome"] = section


def write(species: str, pool: str | None, files: list[Path], outdir: Path) -> list[str]:
    """Write a species' report in outdir: its text, MultiQC row and page."""
    report = build(species, pool, files)
    paths = write_report(report, str(outdir), TOOL, figures=report.plots or None)
    genome = report.sections["genome"]
    quality = report.sections["quality"]
    logger.info(
        "%s: %s bases in %s chromosomes; BUSCO %s %%, QV %s",
        species,
        genome.get("length_bp"),
        genome.get("chromosomes"),
        quality.get("busco_complete"),
        quality.get("qv"),
    )
    return paths
