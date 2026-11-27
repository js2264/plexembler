"""A species' gene models, named and written as its annotation.

A gene predictor names its genes after itself (BRAKER's ``g1``, AUGUSTUS's
``g1.t1``). ``plexembler genes`` reads its models, GTF (BRAKER) or GFF3
(AUGUSTUS), and writes the species' annotation (docs/plan.md, M4):

- the genes in the genome's order, each named ``<species>_000010``,
  ``<species>_000020``, … (a step of 10, room for genes found later), its
  transcripts ``<species>_000010-T1``, …;
- GFF3 and GTF, each gene, transcript, exon and CDS, the stop codon in the CDS;
- the proteins, translated with the species' genetic code (``--genetic-code``:
  NCBI's table; CUG is a sense codon in tables 1, 12 and 26 alike, so the
  gene structures do not change, the proteins do), and the CDS;
- a table of the genes: their transcripts, places, lengths and old names.
"""

from __future__ import annotations

import gzip
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO, cast

from .rename import LINE_WIDTH, read_fasta

logger = logging.getLogger(__name__)

_BASES = "TCAG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
#: NCBI's standard code, its codons in TCAG order.
_STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
#: How each table differs from the standard one.
_DIFFERENCES: dict[int, dict[str, str]] = {
    1: {},
    3: {"ATA": "M", "CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T", "TGA": "W"},
    4: {"TGA": "W"},
    11: {},
    12: {"CTG": "S"},
    26: {"CTG": "A"},
}
GENETIC_CODES = sorted(_DIFFERENCES)
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_GTF_ATTRIBUTE = re.compile(r'(\S+)\s+"([^"]*)"')


def code(table: int) -> dict[str, str]:
    """A genetic code, codon to amino acid."""
    if table not in _DIFFERENCES:
        raise ValueError(
            f"genetic code {table}: not one of {', '.join(map(str, GENETIC_CODES))}"
        )
    standard = dict(zip(_CODONS, _STANDARD, strict=True))
    return standard | _DIFFERENCES[table]


def translate(sequence: str, table: dict[str, str]) -> str:
    """A CDS's protein; ``X`` for a codon with an ambiguous base."""
    sequence = sequence.upper()
    return "".join(
        table.get(sequence[i : i + 3], "X") for i in range(0, len(sequence) - 2, 3)
    )


def revcomp(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass
class Transcript:
    """A transcript: its exons and its CDS, 1-based and inclusive."""

    name: str
    gene: str
    chromosome: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)
    #: The phase of its 5'-most CDS piece, as the predictor gave it.
    phase: int = 0
    new_name: str = ""

    @property
    def start(self) -> int:
        return min(start for start, _ in self.exons or self.cds)

    @property
    def end(self) -> int:
        return max(end for _, end in self.exons or self.cds)


@dataclass
class Gene:
    """A gene and its transcripts."""

    name: str
    chromosome: str
    strand: str
    transcripts: list[Transcript] = field(default_factory=list)
    new_name: str = ""

    @property
    def start(self) -> int:
        return min(t.start for t in self.transcripts)

    @property
    def end(self) -> int:
        return max(t.end for t in self.transcripts)


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Intervals, sorted; those that overlap or touch, joined."""
    joined: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if joined and start <= joined[-1][1] + 1:
            joined[-1] = (joined[-1][0], max(joined[-1][1], end))
        else:
            joined.append((start, end))
    return joined


def open_text(path: Path, mode: str = "r") -> TextIO:
    """A text file, gzipped or not: UTF-8, as predictors write their comments
    (AUGUSTUS's name its authors)."""
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, mode + "t", encoding="utf-8"))
    return cast(TextIO, path.open(mode, encoding="utf-8"))


def attributes(column: str) -> tuple[dict[str, str], bool]:
    """A ninth column's attributes, and whether they are GTF's."""
    if '"' in column:
        return dict(_GTF_ATTRIBUTE.findall(column)), True
    found = {}
    for part in column.strip().strip(";").split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            found[key.strip()] = value.strip()
    return found, False


def read_models(path: Path) -> list[Gene]:
    """A predictor's genes, GTF or GFF3: each transcript's exons and CDS (the
    stop codon in it)."""
    transcripts: dict[str, Transcript] = {}
    gene_of: dict[str, str] = {}
    stops: dict[str, list[tuple[int, int]]] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"{path}: a line of {len(fields)} columns")
            chromosome, feature = fields[0], fields[2]
            start, end, strand = int(fields[3]), int(fields[4]), fields[6]
            found, gtf = attributes(fields[8])
            if feature in ("mRNA", "transcript"):
                # GFF3's say their gene; GTF's (BRAKER's: a bare name) say nothing more.
                if not gtf and "ID" in found:
                    gene_of[found["ID"]] = found.get("Parent", found["ID"])
                continue
            if feature not in ("exon", "CDS", "stop_codon"):
                continue
            if gtf:
                parents = [found["transcript_id"]]
                gene_of.setdefault(parents[0], found.get("gene_id", parents[0]))
            else:
                parents = found.get("Parent", "").split(",")
            for parent in parents:
                transcript = transcripts.setdefault(
                    parent, Transcript(parent, "", chromosome, strand)
                )
                if feature == "exon":
                    transcript.exons.append((start, end))
                elif feature == "CDS":
                    transcript.cds.append((start, end))
                    five_prime = (
                        start == min(s for s, _ in transcript.cds)
                        if strand == "+"
                        else end == max(e for _, e in transcript.cds)
                    )
                    if five_prime and fields[7] in "012":
                        transcript.phase = int(fields[7])
                else:
                    stops.setdefault(parent, []).append((start, end))
    genes: dict[str, Gene] = {}
    for name, transcript in transcripts.items():
        if not transcript.cds:
            continue
        transcript.cds = merge(transcript.cds + stops.get(name, []))
        transcript.exons = merge(transcript.exons + transcript.cds)
        transcript.gene = gene_of.get(name, name)
        gene = genes.setdefault(
            transcript.gene,
            Gene(transcript.gene, transcript.chromosome, transcript.strand),
        )
        gene.transcripts.append(transcript)
    return list(genes.values())


def name_genes(genes: list[Gene], species: str, order: list[str]) -> None:
    """Each gene's and transcript's new name, the genes in the genome's order."""
    rank = {name: index for index, name in enumerate(order)}
    genes.sort(key=lambda g: (rank.get(g.chromosome, len(rank)), g.start, g.end))
    for index, gene in enumerate(genes, start=1):
        gene.new_name = f"{species}_{10 * index:06d}"
        gene.transcripts.sort(key=lambda t: (t.start, t.end, t.name))
        for number, transcript in enumerate(gene.transcripts, start=1):
            transcript.new_name = f"{gene.new_name}-T{number}"


def cds_sequence(transcript: Transcript, genome: dict[str, str]) -> str:
    """A transcript's CDS, 5' to 3', its 5' phase skipped."""
    chromosome = genome[transcript.chromosome]
    pieces = [chromosome[start - 1 : end] for start, end in transcript.cds]
    sequence = "".join(pieces)
    if transcript.strand == "-":
        sequence = revcomp(sequence)
    return sequence[transcript.phase :]


def phases(transcript: Transcript) -> list[int]:
    """Each CDS piece's phase, in the genome's order."""
    ordered = (
        transcript.cds if transcript.strand == "+" else list(reversed(transcript.cds))
    )
    found, done = [], -transcript.phase
    for start, end in ordered:
        found.append((-done) % 3)
        done += end - start + 1
    return found if transcript.strand == "+" else list(reversed(found))


def write_gff3(genes: list[Gene], path: Path, source: str) -> None:
    with open_text(path, "w") as sink:
        sink.write("##gff-version 3\n")
        for gene in genes:
            sink.write(
                f"{gene.chromosome}\t{source}\tgene\t{gene.start}\t{gene.end}\t.\t"
                f"{gene.strand}\t.\tID={gene.new_name};Name={gene.new_name}\n"
            )
            for t in gene.transcripts:
                name = t.new_name
                sink.write(
                    f"{t.chromosome}\t{source}\tmRNA\t{t.start}\t{t.end}\t.\t{t.strand}\t.\t"
                    f"ID={name};Parent={gene.new_name}\n"
                )
                for number, (start, end) in enumerate(t.exons, start=1):
                    sink.write(
                        f"{t.chromosome}\t{source}\texon\t{start}\t{end}\t.\t{t.strand}\t.\t"
                        f"ID={name}.exon{number};Parent={name}\n"
                    )
                for (start, end), phase in zip(t.cds, phases(t), strict=True):
                    sink.write(
                        f"{t.chromosome}\t{source}\tCDS\t{start}\t{end}\t.\t{t.strand}\t{phase}\t"
                        f"ID={name}.cds;Parent={name}\n"
                    )


def write_gtf(genes: list[Gene], path: Path, source: str) -> None:
    with open_text(path, "w") as sink:
        for gene in genes:
            ids = f'gene_id "{gene.new_name}";'
            sink.write(
                f"{gene.chromosome}\t{source}\tgene\t{gene.start}\t{gene.end}\t.\t"
                f"{gene.strand}\t.\t{ids}\n"
            )
            for t in gene.transcripts:
                both = f'{ids} transcript_id "{t.new_name}";'
                sink.write(
                    f"{t.chromosome}\t{source}\ttranscript\t{t.start}\t{t.end}\t.\t"
                    f"{t.strand}\t.\t{both}\n"
                )
                for start, end in t.exons:
                    sink.write(
                        f"{t.chromosome}\t{source}\texon\t{start}\t{end}\t.\t{t.strand}\t.\t{both}\n"
                    )
                for (start, end), phase in zip(t.cds, phases(t), strict=True):
                    sink.write(
                        f"{t.chromosome}\t{source}\tCDS\t{start}\t{end}\t.\t{t.strand}\t{phase}\t{both}\n"
                    )


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with open_text(path, "w") as sink:
        for name, sequence in records:
            sink.write(f">{name}\n")
            for start in range(0, len(sequence), LINE_WIDTH):
                sink.write(sequence[start : start + LINE_WIDTH] + "\n")


def annotate(
    genome_path: Path,
    models: Path,
    species: str,
    outdir: Path,
    genetic_code: int = 1,
    source: str = "plexembler",
) -> list[Gene]:
    """Name a species' gene models and write its annotation.

    Raises
    ------
    ValueError
        When the models are malformed, lie on a sequence the genome lacks, or the
        genetic code is not one the package knows.
    """
    table = code(genetic_code)
    with open_text(genome_path) as handle:
        genome = dict(read_fasta(handle))
    genes = read_models(models)
    for gene in genes:
        if gene.chromosome not in genome:
            raise ValueError(
                f"{models}: {gene.name} lies on {gene.chromosome}, not in the genome"
            )
    name_genes(genes, species, list(genome))
    outdir.mkdir(parents=True, exist_ok=True)
    proteins, coding, internal = [], [], 0
    for gene in genes:
        for t in gene.transcripts:
            sequence = cds_sequence(t, genome)
            protein = translate(sequence, table)
            if protein.endswith("*"):
                protein = protein[:-1]
            internal += "*" in protein
            proteins.append((t.new_name, protein))
            coding.append((t.new_name, sequence))
    write_gff3(genes, outdir / f"{species}.gff3.gz", source)
    write_gtf(genes, outdir / f"{species}.gtf.gz", source)
    write_fasta(proteins, outdir / f"{species}.proteins.faa.gz")
    write_fasta(coding, outdir / f"{species}.cds.fna.gz")
    lengths = dict(proteins)
    with (outdir / f"{species}.genes.tsv").open("w", encoding="ascii") as sink:
        sink.write(
            "gene\ttranscripts\tchromosome\tstart\tend\tstrand\texons\tprotein_length\toriginal_name\n"
        )
        for gene in genes:
            first = gene.transcripts[0]
            sink.write(
                f"{gene.new_name}\t{len(gene.transcripts)}\t{gene.chromosome}\t{gene.start}\t"
                f"{gene.end}\t{gene.strand}\t{len(first.exons)}\t"
                f"{len(lengths[first.new_name])}\t{gene.name}\n"
            )
    transcripts = sum(len(g.transcripts) for g in genes)
    with (outdir / f"{species}.genes_mqc.tsv").open("w", encoding="ascii") as sink:
        sink.write(
            "# id: plexembler_genes\n"
            "# section_name: Genes\n"
            "# description: each species' gene models, named and translated with its "
            "genetic code (plexembler genes)\n"
            "# plot_type: table\n"
            "Sample\tGenes\tTranscripts\tGenetic code\tInternal stops\n"
            f"{species}\t{len(genes)}\t{transcripts}\t{genetic_code}\t{internal}\n"
        )
    if internal:
        logger.warning(
            "%s: %d proteins with an internal stop codon: is the genetic code %d right?",
            species,
            internal,
            genetic_code,
        )
    logger.info(
        "%s: %d genes, %d transcripts, genetic code %d",
        species,
        len(genes),
        transcripts,
        genetic_code,
    )
    return genes
