"""The benchmark pool's small set: three yeasts, and the reads of their pool.

Two or three chromosomes and the mitochondrion of three yeasts with complete
references (docs/plan.md, T2): *Saccharomyces cerevisiae* S288C and its sister
species *S. paradoxus* CBS432, about 88 % identical, the pair of a genus whose
reads are hard to split; and *Nakaseomyces glabratus* (*Candida glabrata*)
CBS138, distant from both. Their genes are the references' annotations (NCBI).

Each species' genome is its reference with a repeat planted in it: a 9 kb
sequence of its own, a copy for every 200 kb of each chromosome, between its
genes. Its nanopore reads (Badread, R10's error model), none longer than 7 kb,
cannot cross a copy, so its assembly breaks there and the Hi-C scaffolding has
joins to make. The reads share errors (a homopolymer a base too long or too
short, a substitution), as a nanopore consensus keeps them: one every 20 kb,
which the polishing with short reads has to fix. And a planted share of
another species' reads (bleed-through).

For the pool: a WGS library (InSilicoSeq, NovaSeq's errors), deep enough that
the polishing's held-out pairs see past a genome's QV of 40; an Arima Hi-C
library (DpnII and HinfI) simulated here, each pair a ligation of two
restriction fragments' ends, its read through the junction when the molecule
is short, its contacts decaying with distance, some between chromosomes and a
few between species (random ligation), PCR duplicates and unligated molecules;
and four RNA-seq libraries, two conditions of two replicates, reverse-stranded
(dUTP), each gene's fragments drawn at its known level.

Every read's name says where it comes from, which the benchmark's checks
compare against; the levels, the planted repeats and errors are in ``truth/``
and ``truth.json``. The set is deterministic: the same GENERATOR_VERSION writes
the same reads. It sits in ``<outdir>/small/``, with its samplesheets and the
references, and is published as the release ``benchmark-data-v<VERSION>``
(.github/workflows/benchmarks.yml).

    python benchmarks/make_benchmark_data.py --outdir benchmarks/data --threads 4

Needs badread (0.4.1) and iss (InSilicoSeq, 2.0.1) on PATH, Python 3.12 or
later, and the network to download the references from NCBI once.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import itertools
import json
import math
import random
import re
import shutil
import subprocess
import tempfile
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from pathlib import Path

#: Raised when what this script writes changes: it names the data's release.
GENERATOR_VERSION = 2
NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class Species:
    """A species of the pool, its reference, the sequences kept and their names."""

    name: str
    organism: str
    accession: str
    ftp_path: str
    sequences: dict[str, str]
    ont_depth: float
    wgs_depth: float
    hic_share: float
    rnaseq_share: float
    n_chromosomes: int


SPECIES = (
    Species(
        "Scer",
        "Saccharomyces cerevisiae S288C",
        "GCF_000146045.2",
        "GCF/000/146/045/GCF_000146045.2_R64",
        {
            "NC_001133.9": "chrI",
            "NC_001135.5": "chrIII",
            "NC_001138.5": "chrVI",
            "NC_001224.1": "chrM",
        },
        ont_depth=40,
        wgs_depth=100,
        hic_share=0.4,
        rnaseq_share=0.4,
        n_chromosomes=3,
    ),
    Species(
        "Spar",
        "Saccharomyces paradoxus CBS432",
        "GCA_002079055.1",
        "GCA/002/079/055/GCA_002079055.1_ASM207905v1",
        {
            "CP020242.1": "chrI",
            "CP020244.1": "chrIII",
            "CP020247.1": "chrVI",
            "CP020258.1": "chrM",
        },
        ont_depth=40,
        wgs_depth=100,
        hic_share=0.35,
        rnaseq_share=0.35,
        n_chromosomes=3,
    ),
    Species(
        "Cgla",
        "Nakaseomyces glabratus CBS138",
        "GCF_000002545.3",
        "GCF/000/002/545/GCF_000002545.3_ASM254v2",
        {"NC_005967.2": "chrA", "NC_005968.1": "chrB", "NC_004691.1": "chrM"},
        ont_depth=40,
        wgs_depth=100,
        hic_share=0.25,
        rnaseq_share=0.25,
        n_chromosomes=2,
    ),
)

#: The mitochondrion's copies per nuclear genome, in every library.
MITO_COPIES = 5
#: Bleed-through: the share of a species' nanopore reads that come from another
#: species of the pool (barcodes misassigned between runs). Spar's run holds
#: enough of Cgla's reads to assemble some: contigs the screen must move or
#: drop. The others, at 1 %, are too thin to assemble.
BLEED_THROUGH = {"Scer": ("Spar", 0.01), "Spar": ("Cgla", 0.15), "Cgla": ("Scer", 0.01)}

#: The planted repeat: a species' own sequence, a copy for every
#: REPEAT_SPACING bases of each nuclear chromosome (one at least), between
#: genes, REPEAT_MARGIN from the chromosome's ends.
REPEAT_LENGTH = 9_000
REPEAT_SPACING = 200_000
REPEAT_MARGIN = 40_000
REPEAT_GENE_GAP = 300
#: Nanopore reads: their lengths (mean, standard deviation), none longer than
#: ONT_MAX_LENGTH, so that no read crosses a planted repeat; simulated deeper
#: by ONT_OVERSAMPLE, which the cap takes back.
ONT_LENGTH = (4_000, 1_800)
ONT_MAX_LENGTH = 7_000
ONT_OVERSAMPLE = 1.15
#: Errors every nanopore read of a place carries, as a consensus keeps them:
#: one every ONT_ERROR_SPACING bases, ONT_HOMOPOLYMER_SHARE of them a
#: homopolymer (4 bases or more) a base too long or too short, the others a
#: substitution. Never in a planted repeat, whose copies stay identical.
ONT_ERROR_SPACING = 20_000
ONT_HOMOPOLYMER_SHARE = 0.7
#: Hi-C: pairs, read length, and how the contacts fall.
HIC_PAIRS = 400_000
READ_LENGTH = 150
HIC_TRANS = 0.12  # between chromosomes of a species
HIC_RANDOM_LIGATION = 0.03  # a partner from anywhere in the pool, by abundance: epsilon
HIC_UNLIGATED = 0.05  # a molecule without a junction
HIC_DUPLICATES = 0.05  # PCR duplicates
#: Arima's enzymes: the site, and where the filled-in end stops (GATC|, GANT|).
SITES = {"DpnII": "GATC", "HinfI": "GA[ACGT]TC"}
#: RNA-seq: the libraries (condition, replicate), each one's pairs, reverse-
#: stranded (dUTP: the first read on the transcript's other strand). A gene's
#: level is log-normal; in heat, RNA_CHANGED of the genes change by a log2
#: fold drawn from N(0, RNA_LOG2FC_SD); each replicate adds its own noise.
RNA_LIBRARIES = (("YPD", 1), ("YPD", 2), ("heat", 1), ("heat", 2))
RNA_PAIRS = 500_000
RNA_READ_LENGTH = 100
RNA_FRAGMENT = (250, 40)
RNA_MIN_LENGTH = 300
RNA_LEVEL_SD = 1.5
RNA_CHANGED = 0.15
RNA_LOG2FC_SD = 1.5
RNA_REPLICATE_SD = 0.15
#: Short reads' substitutions, per base (Hi-C and RNA-seq; InSilicoSeq models
#: the WGS' own).
SHORT_READ_ERRORS = 0.002


def revcomp(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def gzip_writer(path: Path):
    """A gzip text file without a timestamp: the same bytes each time."""
    raw = path.open("wb")
    return raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0)


def write_gzip(path: Path, text_chunks) -> None:
    raw, handle = gzip_writer(path)
    with raw, handle:
        for chunk in text_chunks:
            handle.write(chunk.encode("ascii"))


def read_fasta(text: str) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name is not None:
            records[name].append(line.strip())
    return {name: "".join(chunks).upper() for name, chunks in records.items()}


def download(species: Species, cache: Path, suffix: str) -> Path:
    """A file of the species' assembly at NCBI, downloaded once."""
    cache.mkdir(parents=True, exist_ok=True)
    name = species.ftp_path.rsplit("/", 1)[1]
    path = cache / f"{name}_{suffix}"
    if not path.exists():
        url = f"{NCBI}/{species.ftp_path}/{name}_{suffix}"
        with urllib.request.urlopen(url) as response, path.open("wb") as sink:
            shutil.copyfileobj(response, sink)
    return path


def reference(species: Species, cache: Path) -> dict[str, str]:
    """The species' kept sequences, named ``<species>_<chromosome>``."""
    with gzip.open(download(species, cache, "genomic.fna.gz"), "rt") as handle:
        records = read_fasta(handle.read())
    return {
        f"{species.name}_{new}": records[old] for old, new in species.sequences.items()
    }


@dataclass
class Gene:
    """A gene of a reference's annotation, its first transcript: 1-based,
    inclusive coordinates."""

    gene_id: str
    chromosome: str
    strand: str
    biotype: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)

    @property
    def start(self) -> int:
        return min(start for start, _ in self.exons)

    @property
    def end(self) -> int:
        return max(end for _, end in self.exons)

    def shifted(self, after: int, by: int) -> Gene:
        """The gene, its coordinates past ``after`` moved by ``by``."""

        def move(intervals):
            return [
                (s + by if s > after else s, e + by if e > after else e)
                for s, e in intervals
            ]

        return replace(self, exons=move(self.exons), cds=move(self.cds))


def annotation(species: Species, cache: Path) -> list[Gene]:
    """The genes of the species' kept sequences, from its assembly's GTF."""
    names = {old: f"{species.name}_{new}" for old, new in species.sequences.items()}
    biotypes: dict[str, str] = {}
    first: dict[str, str] = {}
    genes: dict[str, Gene] = {}
    with gzip.open(download(species, cache, "genomic.gtf.gz"), "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] not in names:
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            gene_id = attributes["gene_id"]
            if fields[2] == "gene":
                biotypes[gene_id] = attributes.get("gene_biotype", "unknown")
                continue
            transcript = attributes.get("transcript_id", "")
            if fields[2] not in ("exon", "CDS", "stop_codon") or not transcript:
                continue
            # A gene's first transcript only.
            if first.setdefault(gene_id, transcript) != transcript:
                continue
            gene = genes.setdefault(
                gene_id,
                Gene(gene_id, names[fields[0]], fields[6], biotypes.get(gene_id, "")),
            )
            interval = (int(fields[3]), int(fields[4]))
            # The CDS with its stop codon, as the pipeline's genes write theirs.
            (gene.exons if fields[2] == "exon" else gene.cds).append(interval)
    for gene in genes.values():
        gene.biotype = biotypes.get(gene.gene_id, gene.biotype)
        gene.exons.sort()
        gene.cds = merge_intervals(gene.cds)
    return sorted(
        (g for g in genes.values() if g.exons), key=lambda g: (g.chromosome, g.start)
    )


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Intervals sorted, those touching or overlapping as one."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def spliced(sequence: str, intervals: list[tuple[int, int]], strand: str) -> str:
    """The sequence of 1-based inclusive intervals, joined, on the gene's strand."""
    joined = "".join(sequence[start - 1 : end] for start, end in sorted(intervals))
    return revcomp(joined) if strand == "-" else joined


def random_sequence(length: int, gc: float, rng: random.Random) -> str:
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=weights, k=length))


def plant_repeats(
    genome: dict[str, str], genes: list[Gene], rng: random.Random
) -> tuple[dict[str, str], list[Gene], list[dict]]:
    """The genome with its planted repeat, the genes moved past each copy, and
    where each copy sits (1-based, inclusive)."""
    nuclear = {n: s for n, s in genome.items() if not n.endswith("_chrM")}
    gc = sum(s.count("G") + s.count("C") for s in nuclear.values()) / sum(
        len(s) for s in nuclear.values()
    )
    repeat = random_sequence(REPEAT_LENGTH, gc, rng)
    planted = dict(genome)
    copies: list[dict] = []
    for name, sequence in nuclear.items():
        spans = merge_intervals(
            [(g.start, g.end) for g in genes if g.chromosome == name]
        )
        gaps = [
            (left_end + REPEAT_GENE_GAP, right_start - REPEAT_GENE_GAP)
            for (_, left_end), (right_start, _) in itertools.pairwise(spans)
            if right_start - left_end > 2 * REPEAT_GENE_GAP
        ]
        n_copies = max(1, len(sequence) // REPEAT_SPACING)
        fits = [
            (a + b) // 2
            for a, b in gaps
            if REPEAT_MARGIN < (a + b) // 2 < len(sequence) - REPEAT_MARGIN
        ]
        points = []
        for k in range(1, n_copies + 1):
            target = len(sequence) * k // (n_copies + 1)
            # The gap between genes nearest the copy's place, each copy its own.
            point = min(fits, key=lambda p: abs(p - target))
            fits.remove(point)
            points.append(point)
        # Planted from the chromosome's end, so that each point stands.
        for point in sorted(points, reverse=True):
            sequence = sequence[:point] + repeat + sequence[point:]
            genes = [
                g.shifted(point, REPEAT_LENGTH) if g.chromosome == name else g
                for g in genes
            ]
        planted[name] = sequence
        for index, point in enumerate(sorted(points)):
            start = point + index * REPEAT_LENGTH + 1
            copies.append(
                {"chromosome": name, "start": start, "end": start + REPEAT_LENGTH - 1}
            )
    return planted, genes, copies


def systematic_errors(
    genome: dict[str, str], repeats: list[dict], rng: random.Random
) -> tuple[dict[str, str], list[tuple[str, int, str]]]:
    """The genome as its nanopore reads all see it: one error every
    ONT_ERROR_SPACING bases, outside the planted repeats; and each error (the
    true genome's sequence, 1-based position, kind)."""
    mutated = {}
    planted: list[tuple[str, int, str]] = []
    for name, sequence in genome.items():
        avoid = [(r["start"], r["end"]) for r in repeats if r["chromosome"] == name]

        def clear(position: int, avoid: list[tuple[int, int]] = avoid) -> bool:
            return all(not (start - 50 <= position <= end + 50) for start, end in avoid)

        runs = [
            m
            for m in re.finditer(r"A{4,}|C{4,}|G{4,}|T{4,}", sequence)
            if clear(m.start() + 1)
        ]
        n_errors = round(len(sequence) / ONT_ERROR_SPACING)
        n_homopolymers = round(n_errors * ONT_HOMOPOLYMER_SHARE)
        edits: list[tuple[int, str, str]] = []  # (0-based position, kind, new bases)
        for run in rng.sample(runs, min(n_homopolymers, len(runs))):
            if rng.random() < 0.5:
                edits.append((run.start(), "homopolymer deletion", ""))
            else:
                edits.append((run.start(), "homopolymer insertion", run.group()[0] * 2))
        taken = {position for position, _, _ in edits}
        while len(edits) < n_errors:
            position = rng.randrange(len(sequence))
            if position in taken or not clear(position + 1):
                continue
            base = sequence[position]
            if base not in "ACGT":
                continue
            taken.add(position)
            edits.append(
                (position, "substitution", rng.choice("ACGT".replace(base, "")))
            )
        chars = list(sequence)
        for position, kind, new in sorted(edits, reverse=True):
            # A deletion takes the run's first base; an insertion doubles it.
            chars[position : position + 1] = list(new)
            planted.append((name, position + 1, kind))
        mutated[name] = "".join(chars)
    return mutated, sorted(planted)


def write_reference(genome: dict[str, str], path: Path, depth: bool = False) -> None:
    """A FASTA file; with ``depth``, Badread's relative depth and circularity in
    the headers (the mitochondrion: circular, MITO_COPIES times as deep)."""
    with path.open("w", encoding="ascii") as sink:
        for name, sequence in genome.items():
            header = name
            if depth and name.endswith("_chrM"):
                header += f" depth={MITO_COPIES} circular=true"
            sink.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                sink.write(sequence[start : start + 80] + "\n")


def write_gtf(genes: list[Gene], path: Path) -> None:
    """The genes, as GTF: a gene, a transcript, its exons and CDS."""
    with path.open("w", encoding="ascii") as sink:
        for g in genes:
            ids = f'gene_id "{g.gene_id}"; transcript_id "{g.gene_id}.t1";'
            attributes = f'{ids} gene_biotype "{g.biotype}";'
            sink.write(
                f"{g.chromosome}\treference\tgene\t{g.start}\t{g.end}\t.\t{g.strand}\t.\t{attributes}\n"
            )
            sink.write(
                f"{g.chromosome}\treference\ttranscript\t{g.start}\t{g.end}\t.\t{g.strand}\t.\t{ids}\n"
            )
            for start, end in g.exons:
                sink.write(
                    f"{g.chromosome}\treference\texon\t{start}\t{end}\t.\t{g.strand}\t.\t{ids}\n"
                )
            for start, end in g.cds:
                sink.write(
                    f"{g.chromosome}\treference\tCDS\t{start}\t{end}\t.\t{g.strand}\t.\t{ids}\n"
                )


def nanopore_reads(species: Species, ref: Path, seed: int, workdir: Path) -> list[str]:
    """Badread's reads of a species, R10.4.1's error and quality models, as
    FASTQ records; none longer than ONT_MAX_LENGTH."""
    out = workdir / f"{species.name}.badread.fastq"
    with out.open("w") as sink:
        subprocess.run(
            [
                "badread", "simulate",
                "--reference", str(ref),
                "--quantity", f"{species.ont_depth * ONT_OVERSAMPLE:g}x",
                "--error_model", "nanopore2023",
                "--qscore_model", "nanopore2023",
                "--identity", "98,99.5,1",
                "--length", f"{ONT_LENGTH[0]},{ONT_LENGTH[1]}",
                "--seed", str(seed),
            ],
            check=True, stdout=sink, stderr=subprocess.DEVNULL,
        )  # fmt: skip
    lines = out.read_text().splitlines()
    records = ["\n".join(lines[i : i + 4]) + "\n" for i in range(0, len(lines), 4)]
    return [r for r in records if len(r.split("\n", 2)[1]) <= ONT_MAX_LENGTH]


def wgs_reads(
    genomes: dict[str, dict[str, str]], threads: int, seed: int, workdir: Path
):
    """InSilicoSeq's NovaSeq pairs of the pool, each sequence at its species'
    depth (the mitochondrion MITO_COPIES times deeper)."""
    every = workdir / "pool.fa"
    write_reference({n: s for g in genomes.values() for n, s in g.items()}, every)
    counts = workdir / "readcount.txt"
    with counts.open("w") as sink:
        for species in SPECIES:
            for name, sequence in genomes[species.name].items():
                depth = species.wgs_depth * (
                    MITO_COPIES if name.endswith("_chrM") else 1
                )
                # InSilicoSeq counts reads, the two of a pair each.
                reads = 2 * round(depth * len(sequence) / (2 * READ_LENGTH))
                sink.write(f"{name}\t{reads}\n")
    prefix = workdir / "wgs"
    subprocess.run(
        [
            "iss", "generate",
            "--genomes", str(every),
            "--readcount_file", str(counts),
            "--model", "novaseq",
            "--seed", str(seed),
            "--cpus", str(threads),
            "--output", str(prefix),
            "--quiet",
        ],
        check=True,
    )  # fmt: skip
    return Path(f"{prefix}_R1.fastq"), Path(f"{prefix}_R2.fastq")


@dataclass(frozen=True)
class Sites:
    """A sequence's restriction sites: where each site starts, and its enzyme."""

    starts: list[int]
    enzymes: list[str]


def find_sites(sequence: str) -> Sites:
    found = sorted(
        (m.start(), enzyme)
        for enzyme, pattern in SITES.items()
        for m in re.finditer(f"(?={pattern})", sequence)
    )
    return Sites([start for start, _ in found], [enzyme for _, enzyme in found])


def end_at(sequence: str, start: int, enzyme: str, length: int) -> str:
    """A fragment's end, filled in: the ``length`` bases up to the site's
    filled-in end (GATC for DpnII, GANT for HinfI's G^ANTC)."""
    stop = start + 4
    return sequence[max(0, stop - length) : stop]


def start_at(sequence: str, start: int, enzyme: str, length: int) -> str:
    """The other fragment's end: from its site's filled-in start (GATC, ANTC)."""
    begin = start if enzyme == "DpnII" else start + 1
    return sequence[begin : begin + length]


def mutate(read: str, rng: random.Random, rate: float = SHORT_READ_ERRORS) -> str:
    """A read with its substitutions: their number binomial, their places
    uniform."""
    n_errors = rng.binomialvariate(len(read), rate)
    if not n_errors:
        return read
    chars = list(read)
    for position in rng.sample(range(len(read)), n_errors):
        if chars[position] in "ACGT":
            chars[position] = rng.choice("ACGT".replace(chars[position], ""))
    return "".join(chars)


def hic_pairs(genomes: dict[str, dict[str, str]], seed: int):
    """The pool's Hi-C pairs, as FASTQ records of its two reads."""
    rng = random.Random(seed)
    sites = {
        name: find_sites(sequence)
        for genome in genomes.values()
        for name, sequence in genome.items()
    }
    sequences = {n: s for g in genomes.values() for n, s in g.items()}
    by_species = {
        s.name: [n for n in genomes[s.name] if not n.endswith("_chrM")] for s in SPECIES
    }
    weights = [s.hic_share for s in SPECIES]
    quality = "I" * READ_LENGTH

    def random_site(species: str, chromosome: str | None = None):
        names = [chromosome] if chromosome else by_species[species]
        name = rng.choices(names, weights=[len(sequences[n]) for n in names])[0]
        index = rng.randrange(len(sites[name].starts))
        return name, index

    def nearby_site(name: str, index: int):
        """A site at a distance drawn from P(s) ~ 1/s, from 1 kb to the chromosome.

        A distance that falls off the chromosome is drawn again: held at the
        end, it would pile contacts on the chromosome's last site.
        """
        here = sites[name].starts[index]
        length = len(sequences[name])
        while True:
            distance = math.exp(rng.uniform(math.log(1_000), math.log(length)))
            target = here + rng.choice((-1, 1)) * distance
            if 0 <= target < length:
                break
        return min(
            bisect.bisect_left(sites[name].starts, target), len(sites[name].starts) - 1
        )

    records_1: list[str] = []
    records_2: list[str] = []
    made = 0
    last = None
    while made < HIC_PAIRS:
        if last is not None and rng.random() < HIC_DUPLICATES:
            read_1, read_2, label = last
        else:
            species = rng.choices([s.name for s in SPECIES], weights=weights)[0]
            name_a, index_a = random_site(species)
            insert = max(READ_LENGTH + 20, int(rng.gauss(350, 80)))
            roll = rng.random()
            if roll < HIC_UNLIGATED:
                start = sites[name_a].starts[index_a]
                molecule = sequences[name_a][
                    max(0, start - insert // 2) : start + insert // 2
                ]
                label = f"{name_a}:{start}:unligated"
            else:
                if roll < HIC_UNLIGATED + HIC_RANDOM_LIGATION:
                    # Random ligation: a partner from the whole pool, its own
                    # species included, by abundance. The pairs between species
                    # are then epsilon (1 - sum of p squared) of the ligated.
                    other = rng.choices([s.name for s in SPECIES], weights=weights)[0]
                    name_b, index_b = random_site(other)
                elif (
                    roll < HIC_UNLIGATED + HIC_RANDOM_LIGATION + HIC_TRANS
                    and len(by_species[species]) > 1
                ):
                    name_b, index_b = random_site(
                        species,
                        rng.choice([n for n in by_species[species] if n != name_a]),
                    )
                else:
                    name_b, index_b = name_a, nearby_site(name_a, index_a)
                split = rng.randint(30, insert - 30)
                site_a, enzyme_a = (
                    sites[name_a].starts[index_a],
                    sites[name_a].enzymes[index_a],
                )
                site_b, enzyme_b = (
                    sites[name_b].starts[index_b],
                    sites[name_b].enzymes[index_b],
                )
                part_a = end_at(sequences[name_a], site_a, enzyme_a, split)
                part_b = start_at(sequences[name_b], site_b, enzyme_b, insert - split)
                # Either strand of the other fragment: both enzymes' sites are palindromes.
                if rng.random() < 0.5:
                    part_b = revcomp(
                        end_at(sequences[name_b], site_b, enzyme_b, insert - split)
                    )
                molecule = part_a + part_b
                label = f"{name_a}:{site_a}:{enzyme_a} {name_b}:{site_b}:{enzyme_b}"
            if rng.random() < 0.5:
                molecule = revcomp(molecule)
            if len(molecule) < READ_LENGTH:
                continue
            read_1 = mutate(molecule[:READ_LENGTH], rng)
            read_2 = mutate(revcomp(molecule)[:READ_LENGTH], rng)
        last = (read_1, read_2, label)
        records_1.append(
            f"@hic_{made:07d} {label}\n{read_1}\n+\n{quality[: len(read_1)]}\n"
        )
        records_2.append(
            f"@hic_{made:07d} {label}\n{read_2}\n+\n{quality[: len(read_2)]}\n"
        )
        made += 1
    return records_1, records_2


def library_label(condition: str, replicate: int) -> str:
    """A library's column in the expression tables, as the pipeline names it."""
    return f"{condition}_{replicate}"


def expression_levels(
    transcripts: dict[str, str], rng: random.Random
) -> dict[str, dict[str, float]]:
    """Each RNA-seq library's level of each gene (molecules, relative)."""
    genes = sorted(transcripts)
    base = {g: math.exp(rng.gauss(0, RNA_LEVEL_SD)) for g in genes}
    changed = set(rng.sample(genes, round(RNA_CHANGED * len(genes))))
    folds = {g: 2 ** rng.gauss(0, RNA_LOG2FC_SD) for g in sorted(changed)}
    levels: dict[str, dict[str, float]] = {}
    for condition, replicate in RNA_LIBRARIES:
        levels[library_label(condition, replicate)] = {
            g: base[g]
            * (folds.get(g, 1.0) if condition == "heat" else 1.0)
            * math.exp(rng.gauss(0, RNA_REPLICATE_SD))
            for g in genes
        }
    return levels


def rnaseq_pairs(
    transcripts: dict[str, dict[str, str]],
    levels: dict[str, dict[str, dict[str, float]]],
    label: str,
    seed: int,
):
    """A pooled RNA-seq library's pairs, as FASTQ records of its two reads:
    each species by its share, each of its genes' fragments by its level and
    length; the first read on the transcript's other strand (dUTP)."""
    rng = random.Random(seed)
    quality = "I" * RNA_READ_LENGTH
    species_names = [s.name for s in SPECIES]
    shares = [s.rnaseq_share for s in SPECIES]
    tables = {}
    for name in species_names:
        genes = sorted(transcripts[name])
        weights = [
            levels[name][label][g]
            * max(1, len(transcripts[name][g]) - RNA_FRAGMENT[0] + 1)
            for g in genes
        ]
        tables[name] = (genes, list(itertools.accumulate(weights)))
    records_1: list[str] = []
    records_2: list[str] = []
    for index in range(RNA_PAIRS):
        species = rng.choices(species_names, weights=shares)[0]
        genes, cumulative = tables[species]
        gene = genes[bisect.bisect(cumulative, rng.random() * cumulative[-1])]
        sequence = transcripts[species][gene]
        length = min(
            len(sequence),
            max(RNA_READ_LENGTH + 20, int(rng.gauss(*RNA_FRAGMENT))),
        )
        start = rng.randint(0, len(sequence) - length)
        fragment = sequence[start : start + length]
        read_1 = mutate(revcomp(fragment)[:RNA_READ_LENGTH], rng)
        read_2 = mutate(fragment[:RNA_READ_LENGTH], rng)
        name = f"rna_{label}_{index:07d} {species}:{gene}"
        records_1.append(f"@{name}\n{read_1}\n+\n{quality[: len(read_1)]}\n")
        records_2.append(f"@{name}\n{read_2}\n+\n{quality[: len(read_2)]}\n")
    return records_1, records_2


def main(outdir: Path, threads: int) -> None:
    small = outdir / "small"
    if small.exists():
        shutil.rmtree(small)
    for sub in ("ont", "pool", "references", "truth"):
        (small / sub).mkdir(parents=True, exist_ok=True)
    rng = random.Random(f"plexembler-benchmark-v{GENERATOR_VERSION}")

    genomes: dict[str, dict[str, str]] = {}
    genes: dict[str, list[Gene]] = {}
    repeats: dict[str, list[dict]] = {}
    for species in SPECIES:
        genome = reference(species, outdir / ".cache")
        genome, genes[species.name], repeats[species.name] = plant_repeats(
            genome, annotation(species, outdir / ".cache"), rng
        )
        genomes[species.name] = genome
        write_reference(genome, small / "references" / f"{species.name}.fa")
        write_gtf(genes[species.name], small / "references" / f"{species.name}.gtf")

    # The nuclear protein-coding genes: their CDS (the M4 check finds them among
    # the pipeline's) and their transcripts (the RNA-seq is drawn from them).
    transcripts: dict[str, dict[str, str]] = {}
    for species in SPECIES:
        coding = [
            g
            for g in genes[species.name]
            if g.biotype == "protein_coding"
            and g.cds
            and not g.chromosome.endswith("_chrM")
        ]
        genome = genomes[species.name]
        with (small / "references" / f"{species.name}.cds.fa").open("w") as sink:
            for g in coding:
                sink.write(
                    f">{g.gene_id}\n{spliced(genome[g.chromosome], g.cds, g.strand)}\n"
                )
        transcripts[species.name] = {
            g.gene_id: spliced(genome[g.chromosome], g.exons, g.strand)
            for g in coding
            if sum(e - s + 1 for s, e in g.exons) >= RNA_MIN_LENGTH
        }

    errors: dict[str, list[tuple[str, int, str]]] = {}
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for species in SPECIES:
            seen, errors[species.name] = systematic_errors(
                genomes[species.name], repeats[species.name], rng
            )
            write_reference(seen, work / f"{species.name}.fa", depth=True)
            with (small / "truth" / f"{species.name}.ont_errors.tsv").open("w") as sink:
                sink.write("sequence\tposition\tkind\n")
                for name, position, kind in errors[species.name]:
                    sink.write(f"{name}\t{position}\t{kind}\n")
        with ThreadPoolExecutor(max_workers=threads) as pool:
            jobs = {
                s.name: pool.submit(
                    nanopore_reads, s, work / f"{s.name}.fa", 1000 + i, work
                )
                for i, s in enumerate(SPECIES)
            }
            reads = {name: job.result() for name, job in jobs.items()}
        bleed = {}
        ont_bases = {}
        for species in SPECIES:
            source, share = BLEED_THROUGH[species.name]
            extra = rng.sample(reads[source], round(share * len(reads[species.name])))
            bleed[species.name] = {"from": source, "share": share, "reads": len(extra)}
            ont_bases[species.name] = sum(
                len(r.split("\n", 2)[1]) for r in reads[species.name]
            )
            write_gzip(
                small / "ont" / f"{species.name}.fastq.gz", reads[species.name] + extra
            )

        r1, r2 = wgs_reads(genomes, threads, 2000, work)
        for read, name in ((r1, "R1"), (r2, "R2")):
            with read.open() as handle:
                write_gzip(small / "pool" / f"wgs_{name}.fastq.gz", handle)

    hic_1, hic_2 = hic_pairs(genomes, 3000)
    write_gzip(small / "pool" / "hic_R1.fastq.gz", hic_1)
    write_gzip(small / "pool" / "hic_R2.fastq.gz", hic_2)
    del hic_1, hic_2

    levels = {s.name: expression_levels(transcripts[s.name], rng) for s in SPECIES}
    for index, (condition, replicate) in enumerate(RNA_LIBRARIES):
        label = library_label(condition, replicate)
        rna_1, rna_2 = rnaseq_pairs(transcripts, levels, label, 4000 + index)
        stem = f"rna-{condition.lower()}-{replicate}"
        write_gzip(small / "pool" / f"{stem}_R1.fastq.gz", rna_1)
        write_gzip(small / "pool" / f"{stem}_R2.fastq.gz", rna_2)
    labels = [library_label(c, r) for c, r in RNA_LIBRARIES]
    for species in SPECIES:
        # The true TPM: a gene's share of the species' molecules, per million.
        totals = {label: sum(levels[species.name][label].values()) for label in labels}
        with (small / "truth" / f"{species.name}.tpm.tsv").open("w") as sink:
            sink.write("gene\tlength_bp\t" + "\t".join(labels) + "\n")
            for gene in sorted(transcripts[species.name]):
                cells = [
                    f"{1e6 * levels[species.name][label][gene] / totals[label]:.3f}"
                    for label in labels
                ]
                sink.write(
                    f"{gene}\t{len(transcripts[species.name][gene])}\t"
                    + "\t".join(cells)
                    + "\n"
                )

    (small / "species.csv").write_text(
        "species,pool,ont,chemistry,genome_size,ploidy,n_chromosomes,busco_lineage,genetic_code,telomere_repeat\n"
        + "".join(
            f"{s.name},bench,ont/{s.name}.fastq.gz,r10,"
            f"{sum(len(x) for x in genomes[s.name].values())},1,{s.n_chromosomes},"
            "saccharomycetes_odb12,1,NA\n"
            for s in SPECIES
        )
    )
    (small / "libraries.csv").write_text(
        "library,pool,type,fastq_1,fastq_2,condition,replicate,enzymes\n"
        'bench-hic,bench,hic,pool/hic_R1.fastq.gz,pool/hic_R2.fastq.gz,,,"DpnII,HinfI"\n'
        "bench-wgs,bench,wgs,pool/wgs_R1.fastq.gz,pool/wgs_R2.fastq.gz,,,\n"
        + "".join(
            f"bench-rna-{c.lower()}-{r},bench,rnaseq,pool/rna-{c.lower()}-{r}_R1.fastq.gz,"
            f"pool/rna-{c.lower()}-{r}_R2.fastq.gz,{c},{r},\n"
            for c, r in RNA_LIBRARIES
        )
    )
    truth = {
        "generator_version": GENERATOR_VERSION,
        "species": {
            s.name: {
                "organism": s.organism,
                "accession": s.accession,
                "sequences": {
                    f"{s.name}_{new}": {
                        "accession": old,
                        "length": len(genomes[s.name][f"{s.name}_{new}"]),
                    }
                    for old, new in s.sequences.items()
                },
                "n_chromosomes": s.n_chromosomes,
                "genes": len(genes[s.name]),
                "protein_coding_nuclear": sum(
                    1
                    for g in genes[s.name]
                    if g.biotype == "protein_coding"
                    and g.cds
                    and not g.chromosome.endswith("_chrM")
                ),
                "expressed": len(transcripts[s.name]),
                "repeats": repeats[s.name],
                "ont_depth": s.ont_depth,
                "ont_depth_simulated": round(
                    ont_bases[s.name] / sum(len(x) for x in genomes[s.name].values()), 1
                ),
                "ont_errors": len(errors[s.name]),
                "wgs_depth": s.wgs_depth,
                "hic_share": s.hic_share,
                "rnaseq_share": s.rnaseq_share,
                "bleed_through": bleed[s.name],
            }
            for s in SPECIES
        },
        "mito_copies": MITO_COPIES,
        "repeat_length": REPEAT_LENGTH,
        "ont": {
            "length_mean": ONT_LENGTH[0],
            "length_sd": ONT_LENGTH[1],
            "max_length": ONT_MAX_LENGTH,
            "error_spacing": ONT_ERROR_SPACING,
            "homopolymer_share": ONT_HOMOPOLYMER_SHARE,
        },
        "hic": {
            "pairs": HIC_PAIRS,
            "enzymes": list(SITES),
            "trans": HIC_TRANS,
            "random_ligation": HIC_RANDOM_LIGATION,
            "between_species_expected": HIC_RANDOM_LIGATION
            * (1 - sum(s.hic_share**2 for s in SPECIES)),
            "unligated": HIC_UNLIGATED,
            "duplicates": HIC_DUPLICATES,
        },
        "rnaseq": {
            "libraries": labels,
            "pairs": RNA_PAIRS,
            "read_length": RNA_READ_LENGTH,
            "strandedness": "reverse",
            "changed_in_heat": RNA_CHANGED,
        },
    }
    (small / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "-o", "--outdir", type=Path, default=Path(__file__).parent / "data"
    )
    parser.add_argument("-t", "--threads", type=int, default=4)
    arguments = parser.parse_args()
    main(arguments.outdir, arguments.threads)
