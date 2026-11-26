"""One candidate kept for a step: the candidates' numbers, the rules, the choice.

A step of the pipeline can make several candidates for a species, such as the
contigs of each assembler (docs/plan.md, M1). ``plexembler select`` reads each
candidate's checks, keeps one by the step's rules, copies its sequence, and
writes a table of every candidate with its numbers and the reason for the
choice, which the species' report shows.

The rules of the assembly step, in order, each breaking the ties of the one
before (a tie is a difference within the rule's tolerance):

1. the most complete, single-copy BUSCOs (complete minus duplicated);
2. the highest consensus quality (Merqury's QV);
3. the longest contig N50.

A number that a candidate lacks leaves its rule out, for every candidate.
"""

from __future__ import annotations

import csv
import gzip
import logging
import math
import re
import shutil
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

#: The steps that have rules, and what each one's kept sequence is called.
STEPS = ("assembly",)

#: BUSCO's one-line summary: C:98.5%[S:97.9%,D:0.6%],F:0.3%,M:1.2%,n:2137
_BUSCO_LINE = re.compile(
    r"C:(?P<complete>[\d.]+)%\[S:(?P<single>[\d.]+)%,D:(?P<duplicated>[\d.]+)%\],"
    r"F:(?P<fragmented>[\d.]+)%,M:(?P<missing>[\d.]+)%,n:(?P<n>\d+)"
)
#: compleasm's summary: one line per class, "S:97.90%, 2092".
_COMPLEASM_LINE = re.compile(r"^(?P<key>[SDFIMN]):(?P<value>[\d.]+)%?", re.MULTILINE)
#: A cell that says nothing: no file for this number.
_EMPTY = {"", "NA", "na", "-"}


@dataclass(frozen=True)
class Busco:
    """A BUSCO run's summary, in percent of its lineage's ``n`` genes."""

    complete: float
    single: float
    duplicated: float
    fragmented: float
    missing: float
    n: int


@dataclass
class Candidate:
    """A candidate of a step: its name, its sequence, and its numbers."""

    name: str
    fasta: Path
    busco: Busco | None = None
    qv: float | None = None
    stats: dict[str, float] = field(default_factory=dict)

    @property
    def contigs(self) -> float | None:
        return self.stats.get("contigs")

    @property
    def contig_n50(self) -> float | None:
        return self.stats.get("contig_n50")

    @property
    def total_length(self) -> float | None:
        return self.stats.get("total_length")


@dataclass(frozen=True)
class Choice:
    """The candidate kept, and why."""

    kept: Candidate
    reason: str


def read_busco(path: Path) -> Busco:
    """A BUSCO short summary (its .txt or .json), or compleasm's summary.txt.

    Raises
    ------
    ValueError
        When the file holds neither summary.
    """
    text = path.read_text(encoding="utf-8")
    found = _BUSCO_LINE.search(text)
    if found:
        values = found.groupdict()
        return Busco(
            complete=float(values["complete"]),
            single=float(values["single"]),
            duplicated=float(values["duplicated"]),
            fragmented=float(values["fragmented"]),
            missing=float(values["missing"]),
            n=int(values["n"]),
        )
    classes = {m["key"]: float(m["value"]) for m in _COMPLEASM_LINE.finditer(text)}
    if {"S", "D", "F", "M", "N"} <= classes.keys():
        return Busco(
            complete=classes["S"] + classes["D"],
            single=classes["S"],
            duplicated=classes["D"],
            fragmented=classes["F"] + classes.get("I", 0.0),
            missing=classes["M"],
            n=int(classes["N"]),
        )
    raise ValueError(f"{path}: neither a BUSCO nor a compleasm summary")


def read_qv(path: Path) -> float:
    """The QV of the first assembly of Merqury's ``.qv`` file.

    Its columns: the assembly, its k-mers found in it alone, its k-mers, the QV
    and the error rate. A QV of ``inf`` (no error found) is kept as infinite.

    Raises
    ------
    ValueError
        When the file has no line with a QV.
    """
    for line in path.read_text(encoding="utf-8").splitlines():
        cells = line.split("\t")
        if len(cells) >= 4:
            try:
                return float(cells[3])
            except ValueError:
                continue
    raise ValueError(f"{path}: no QV in this Merqury file")


#: gfastats' names for the numbers the rules and the table use.
_GFASTATS_KEYS = {
    "# contigs": "contigs",
    "contig n50": "contig_n50",
    "total contig length": "total_length",
    "largest contig": "largest_contig",
    "gc content %": "gc_percent",
    "# gaps": "gaps",
    "# scaffolds": "scaffolds",
    "scaffold n50": "scaffold_n50",
}


def read_gfastats(path: Path) -> dict[str, float]:
    """gfastats' summary, tabular or not: the numbers the table uses, by name."""
    stats: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition("\t")
        if not separator:
            key, separator, value = line.partition(":")
        name = _GFASTATS_KEYS.get(key.strip().lower())
        if name is None:
            continue
        try:
            stats[name] = float(value.strip())
        except ValueError:
            logger.debug("%s: %r is not a number", path, value)
    return stats


def _cell_path(cell: str, folder: Path) -> Path | None:
    cell = cell.strip()
    if cell in _EMPTY:
        return None
    path = Path(cell)
    return path if path.is_absolute() else folder / path


def read_candidates(table: Path) -> list[Candidate]:
    """The candidates of a table: ``candidate``, ``fasta``, and their checks.

    The columns ``busco``, ``qv`` and ``stats`` name a candidate's BUSCO
    summary, Merqury QV and gfastats summary, and may be empty or ``NA``. A
    relative path is read from the table's folder.

    Raises
    ------
    ValueError
        When the table lacks a column it needs, names a candidate twice, or has
        none.
    """
    folder = table.parent
    candidates: list[Candidate] = []
    with table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"candidate", "fasta"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{table}: no column {', '.join(sorted(missing))}")
        for row in reader:
            name = row["candidate"].strip()
            if any(c.name == name for c in candidates):
                raise ValueError(f"{table}: the candidate {name!r} is there twice")
            fasta = _cell_path(row["fasta"], folder)
            if fasta is None:
                raise ValueError(f"{table}: the candidate {name!r} has no sequence")
            candidate = Candidate(name, fasta)
            if busco := _cell_path(row.get("busco") or "", folder):
                candidate.busco = read_busco(busco)
            if qv := _cell_path(row.get("qv") or "", folder):
                candidate.qv = read_qv(qv)
            if stats := _cell_path(row.get("stats") or "", folder):
                candidate.stats = read_gfastats(stats)
            candidates.append(candidate)
    if not candidates:
        raise ValueError(f"{table}: no candidate")
    return candidates


def _narrow(
    candidates: list[Candidate],
    value: Callable[[Candidate], float | None],
    tolerance: float,
) -> list[Candidate]:
    """The candidates within ``tolerance`` of the best by ``value``; all of them
    when one lacks the value, which then decides nothing."""
    values = [value(c) for c in candidates]
    if any(v is None for v in values):
        return candidates
    best = max(v for v in values if v is not None)
    if math.isinf(best):
        return [c for c, v in zip(candidates, values, strict=True) if v == best]
    return [
        c
        for c, v in zip(candidates, values, strict=True)
        if v is not None and v >= best - tolerance
    ]


def choose_assembly(
    candidates: list[Candidate],
    keep: str | None = None,
    busco_tolerance: float = 0.5,
    qv_tolerance: float = 1.0,
) -> Choice:
    """The assembly kept, by the rules of the module's docstring.

    Parameters
    ----------
    candidates
        The species' candidates.
    keep
        A candidate named by the samplesheet, kept whatever its numbers.
    busco_tolerance
        Single-copy BUSCOs within this many points of the best are a tie.
    qv_tolerance
        QVs within this many units of the best are a tie.

    Raises
    ------
    ValueError
        When ``keep`` names no candidate.
    """
    if keep is not None:
        for candidate in candidates:
            if candidate.name == keep:
                return Choice(candidate, "named by the samplesheet")
        names = ", ".join(c.name for c in candidates)
        raise ValueError(f"no candidate {keep!r} to keep: the candidates are {names}")
    if len(candidates) == 1:
        return Choice(candidates[0], "the only candidate")

    rules: list[tuple[str, Callable[[Candidate], float | None], float]] = [
        (
            "the most single-copy complete BUSCOs",
            lambda c: c.busco.single if c.busco else None,
            busco_tolerance,
        ),
        ("the highest QV", lambda c: c.qv, qv_tolerance),
        ("the longest contig N50", lambda c: c.contig_n50, 0.0),
    ]
    left = candidates
    for reason, value, tolerance in rules:
        narrowed = _narrow(left, value, tolerance)
        if len(narrowed) == 1:
            return Choice(narrowed[0], reason)
        left = narrowed
    return Choice(left[0], "a tie on every rule: the first candidate")


def copy_sequence(source: Path, output: Path) -> None:
    """Copy a FASTA to ``output``, gzipped when its name ends in .gz."""
    output.parent.mkdir(parents=True, exist_ok=True)
    gzipped_in = source.suffix == ".gz"
    gzipped_out = output.suffix == ".gz"
    if gzipped_in == gzipped_out:
        shutil.copyfile(source, output)
        return
    opener_in = gzip.open if gzipped_in else open
    opener_out = gzip.open if gzipped_out else open
    with opener_in(source, "rb") as reading, opener_out(output, "wb") as writing:
        shutil.copyfileobj(reading, writing)


def _number(value: float | int | None, digits: int = 2) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float) and math.isinf(value):
        return "inf"
    if isinstance(value, float) and value.is_integer() and digits == 0:
        return str(int(value))
    return f"{value:.{digits}f}" if digits else str(int(value))


#: The table's columns: its header, and how a candidate fills it.
COLUMNS: list[tuple[str, Callable[[Candidate], str]]] = [
    ("busco_complete", lambda c: _number(c.busco.complete if c.busco else None, 1)),
    ("busco_single", lambda c: _number(c.busco.single if c.busco else None, 1)),
    (
        "busco_duplicated",
        lambda c: _number(c.busco.duplicated if c.busco else None, 1),
    ),
    (
        "busco_fragmented",
        lambda c: _number(c.busco.fragmented if c.busco else None, 1),
    ),
    ("busco_missing", lambda c: _number(c.busco.missing if c.busco else None, 1)),
    ("busco_n", lambda c: _number(c.busco.n if c.busco else None, 0)),
    ("qv", lambda c: _number(c.qv, 2)),
    ("contigs", lambda c: _number(c.contigs, 0)),
    ("contig_n50", lambda c: _number(c.contig_n50, 0)),
    ("total_length", lambda c: _number(c.total_length, 0)),
]


def write_selection(
    candidates: Iterable[Candidate], choice: Choice, table: Path
) -> None:
    """The table of candidates: their numbers, the one kept and why."""
    with table.open("w", encoding="utf-8") as sink:
        header = ["candidate", "kept", *(name for name, _ in COLUMNS), "reason"]
        sink.write("\t".join(header) + "\n")
        for candidate in candidates:
            kept = candidate is choice.kept
            cells = [candidate.name, "yes" if kept else "no"]
            cells += [fill(candidate) for _, fill in COLUMNS]
            cells.append(choice.reason if kept else "")
            sink.write("\t".join(cells) + "\n")


def write_multiqc(
    candidates: Iterable[Candidate],
    choice: Choice,
    species: str,
    step: str,
    table: Path,
) -> None:
    """MultiQC's custom content: a row per candidate, ``<species>_<candidate>``.

    Every species' file has the same section ID, so that MultiQC makes one
    table of every species' candidates.
    """
    lines = [
        f"# id: 'plexembler_{step}_selection'",
        f"# section_name: 'Candidates: {step}'",
        f"# description: 'The {step} candidates of each species, their checks, "
        "and the one plexembler kept (plexembler select).'",
        "# plot_type: 'table'",
        "# pconfig:",
        f"#     id: 'plexembler_{step}_selection_table'",
        "#     namespace: 'plexembler'",
        "Sample\tKept\tBUSCO single (%)\tBUSCO duplicated (%)\tQV\tContigs\tContig N50",
    ]
    for c in candidates:
        lines.append(
            "\t".join(
                [
                    f"{species}_{c.name}",
                    "yes" if c is choice.kept else "no",
                    _number(c.busco.single if c.busco else None, 1),
                    _number(c.busco.duplicated if c.busco else None, 1),
                    _number(c.qv, 2),
                    _number(c.contigs, 0),
                    _number(c.contig_n50, 0),
                ]
            )
        )
    table.write_text("\n".join(lines) + "\n", encoding="utf-8")


def output_paths(outdir: Path, prefix: str, step: str) -> dict[str, Path]:
    """The paths ``plexembler select`` writes, by what they hold."""
    return {
        "fasta": outdir / f"{prefix}.{step}.fa.gz",
        "table": outdir / f"{prefix}.{step}.selection.tsv",
        "multiqc": outdir / f"{prefix}.{step}_selection_mqc.tsv",
    }


def select(
    table: Path,
    outdir: Path,
    prefix: str,
    step: str = "assembly",
    keep: str | None = None,
    busco_tolerance: float = 0.5,
    qv_tolerance: float = 1.0,
) -> Choice:
    """Read a step's candidates, keep one, and write the outputs.

    Returns
    -------
    Choice
        The candidate kept, and why.

    Raises
    ------
    ValueError
        When the step has no rules, or the table or ``keep`` is wrong.
    """
    if step not in STEPS:
        raise ValueError(f"no rules for the step {step!r}: {', '.join(STEPS)}")
    candidates = read_candidates(table)
    choice = choose_assembly(candidates, keep, busco_tolerance, qv_tolerance)
    paths = output_paths(outdir, prefix, step)
    copy_sequence(choice.kept.fasta, paths["fasta"])
    write_selection(candidates, choice, paths["table"])
    write_multiqc(candidates, choice, prefix, step, paths["multiqc"])
    for candidate in candidates:
        logger.info(
            "%s: BUSCO single %s %%, QV %s, contig N50 %s",
            candidate.name,
            _number(candidate.busco.single if candidate.busco else None, 1),
            _number(candidate.qv, 2),
            _number(candidate.contig_n50, 0),
        )
    logger.info("Kept %s for %s: %s", choice.kept.name, prefix, choice.reason)
    return choice
