"""The benchmark pool's three sets: three yeasts, and the reads of their pool.

Two or three chromosomes and the mitochondrion of three yeasts with complete
references (docs/plan.md, T2): *Saccharomyces cerevisiae* S288C and its sister
species *S. paradoxus* CBS432, about 88 % identical, the pair of a genus whose
reads are hard to split; and *Nakaseomyces glabratus* (*Candida glabrata*)
CBS138, distant from both. Their genes are the references' annotations (NCBI),
each protein-coding gene characterized or not (``truth/<species>.genes.tsv``).

Three sets, from the easiest to the hardest (LEVELS): each plants more of the
repeats yeast genomes hold, and gives nanopore reads (R10's models throughout)
and Hi-C less able to resolve them.

- easy: a Ty1 element on each chromosome; reads of 99.5 % identity, an N50 of
  about 25 kb, 50x; a million Hi-C pairs, 1 % of them random ligations.
- medium: a Ty1 element upstream of a tRNA gene every 100 kb, and Y' elements,
  one or two in tandem, at half the chromosomes' ends; reads of 98 %, an N50 of
  about 11 kb, 40x; 400,000 Hi-C pairs, 3 % random ligations.
- hard: a Ty1 every 70 kb and three in tandem once per chromosome, two to four
  Y' in tandem at three quarters of the ends, and an rDNA array of twelve
  9.1 kb units on each species' longest chromosome; reads of 96 %, an N50 of
  about 6 kb, 30x; 150,000 Hi-C pairs, 6 % random ligations.

The repeats are S288C's own (FAMILY_SOURCES): its Ty1 element YARCTy1-1, the
Y' element at chrXIV's left end, a unit of its rDNA; each species' family is
S288C's diverged (FAMILY_DIVERGENCE), each copy its family's, diverged a
little more (COPY_DIVERGENCE). The nanopore reads share errors, as a nanopore
consensus keeps them (a homopolymer a base too long or too short, a
substitution), which the polishing with short reads has to fix; and a share of
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
and ``truth.json``. The sets are deterministic: the same GENERATOR_VERSION
writes the same reads. Each sits in ``<outdir>/<level>/``, with its
samplesheets and the references, and is published as an asset of the release
``benchmark-data-v<VERSION>`` (.github/workflows/benchmarks.yml).

    python benchmarks/make_benchmark_data.py --outdir benchmarks/data --threads 4
    python benchmarks/make_benchmark_data.py --level easy

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
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

#: Raised when what this script writes changes: it names the data's release.
GENERATOR_VERSION = 3
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
        wgs_depth=100,
        hic_share=0.25,
        rnaseq_share=0.25,
        n_chromosomes=2,
    ),
)


@dataclass(frozen=True)
class Level:
    """A set's difficulty: the repeats planted in each genome, its nanopore
    reads, its Hi-C."""

    name: str
    #: A Ty1 element upstream of a tRNA gene, one for every ``ty1_spacing``
    #: bases of a chromosome (None: one per chromosome).
    ty1_spacing: int | None
    #: Ty1 elements in tandem, at one more place on each chromosome (0: none).
    ty1_cluster: int
    #: The share of the chromosomes' ends given Y' elements, and how many in
    #: tandem (fewest, most).
    yprime_ends: float
    yprime_copies: tuple[int, int]
    #: The units of an rDNA array on each species' longest chromosome (0: none).
    rdna_units: int
    #: Nanopore reads: their depth, their length (mean, standard deviation:
    #: Badread's gamma), their identity (mean, max, standard deviation: its
    #: beta), and the errors they all share, one every ``ont_error_spacing``
    #: bases.
    ont_depth: float
    ont_length: tuple[int, int]
    ont_identity: tuple[float, float, float]
    ont_error_spacing: int
    #: Hi-C: its pairs, and how they fall: between chromosomes of a species,
    #: random ligations (a partner from anywhere in the pool, by abundance:
    #: epsilon), molecules without a junction, PCR duplicates.
    hic_pairs: int
    hic_trans: float
    hic_random_ligation: float
    hic_unligated: float
    hic_duplicates: float


LEVELS = (
    Level(
        "easy",
        ty1_spacing=None,
        ty1_cluster=0,
        yprime_ends=0.0,
        yprime_copies=(0, 0),
        rdna_units=0,
        ont_depth=50,
        ont_length=(18_000, 14_000),
        ont_identity=(99.5, 99.9, 0.3),
        ont_error_spacing=50_000,
        hic_pairs=1_000_000,
        hic_trans=0.10,
        hic_random_ligation=0.01,
        hic_unligated=0.03,
        hic_duplicates=0.03,
    ),
    Level(
        "medium",
        ty1_spacing=100_000,
        ty1_cluster=0,
        yprime_ends=0.5,
        yprime_copies=(1, 2),
        rdna_units=0,
        ont_depth=40,
        ont_length=(8_000, 6_000),
        ont_identity=(98.0, 99.5, 1.0),
        ont_error_spacing=20_000,
        hic_pairs=400_000,
        hic_trans=0.12,
        hic_random_ligation=0.03,
        hic_unligated=0.05,
        hic_duplicates=0.05,
    ),
    Level(
        "hard",
        ty1_spacing=70_000,
        ty1_cluster=3,
        yprime_ends=0.75,
        yprime_copies=(2, 4),
        rdna_units=12,
        ont_depth=30,
        ont_length=(4_500, 3_000),
        ont_identity=(96.0, 98.5, 2.0),
        ont_error_spacing=10_000,
        hic_pairs=150_000,
        hic_trans=0.15,
        hic_random_ligation=0.06,
        hic_unligated=0.08,
        hic_duplicates=0.10,
    ),
)

#: The repeat families, S288C's own (1-based, inclusive): its Ty1 element
#: YARCTy1-1 from LTR to LTR (5,925 bp, its 337 bp LTRs identical); the Y'
#: element at chrXIV's left end, from the chromosome's start to the telomeric
#: repeat between it and the X element (6,617 bp); and a unit of its rDNA, from
#: one RDN37 to the next (9,137 bp).
FAMILY_SOURCES = {
    "Ty1": ("NC_001133.9", 160_238, 166_162),
    "Y'": ("NC_001146.8", 1, 6_617),
    "rDNA": ("NC_001144.5", 451_575, 460_711),
}
#: Each species' families: S288C's, diverged by this share of substitutions
#: (and a tenth as many 1-base indels): S. paradoxus's own, N. glabratus's
#: distant ones.
FAMILY_DIVERGENCE = {"Scer": 0.0, "Spar": 0.08, "Cgla": 0.25}
#: Each planted copy: its family's sequence, diverged by this share.
COPY_DIVERGENCE = {"Ty1": 0.002, "Y'": 0.005, "rDNA": 0.0005}
#: A Ty1 element lands upstream of a tRNA gene, this far from its 5' end, and
#: duplicates the 5 bases it lands on.
TY1_UPSTREAM = (100, 700)
TY1_TARGET_DUPLICATION = 5
#: Ty1 elements stay this far from a chromosome's ends, where the Y' sit.
TY1_MARGIN = 25_000
#: The telomeric repeat between two tandem Y' elements: its length.
ITS_LENGTH = (50, 130)
#: The rDNA array's place: the gap between genes nearest this share of its
#: chromosome.
RDNA_PLACE = 0.6
#: The gap a planted repeat keeps from a gene.
GENE_GAP = 300

#: The mitochondrion's copies per nuclear genome, in every library.
MITO_COPIES = 5
#: Bleed-through: the share of a species' nanopore reads that come from another
#: species of the pool (barcodes misassigned between runs). Spar's run holds
#: enough of Cgla's reads to assemble some: contigs the screen must move or
#: drop. The others, at 1 %, are too thin to assemble.
BLEED_THROUGH = {"Scer": ("Spar", 0.01), "Spar": ("Cgla", 0.15), "Cgla": ("Scer", 0.01)}
#: The errors every nanopore read of a place carries: ONT_HOMOPOLYMER_SHARE of
#: them a homopolymer (4 bases or more) a base too long or too short, the
#: others a substitution. Never in a planted repeat, whose copies are as
#: planted.
ONT_HOMOPOLYMER_SHARE = 0.7
#: Hi-C's read length; Arima's enzymes: the site, and where the filled-in end
#: stops (GATC|, GANT|).
READ_LENGTH = 150
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
#: The S. cerevisiae gene a note names (N. glabratus's: "similar to
#: uniprot|P53192 Saccharomyces cerevisiae YGL020c").
SCER_GENE = re.compile(
    r"Saccharomyces cerevisiae (Y[A-P][LR]\d{3}[WC](?:-[A-Z])?)", re.IGNORECASE
)


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


def assembly(species: Species, cache: Path) -> dict[str, str]:
    """Every sequence of the species' assembly, by its accession."""
    with gzip.open(download(species, cache, "genomic.fna.gz"), "rt") as handle:
        return read_fasta(handle.read())


def reference(species: Species, records: dict[str, str]) -> dict[str, str]:
    """The species' kept sequences, named ``<species>_<chromosome>``."""
    return {
        f"{species.name}_{new}": records[old] for old, new in species.sequences.items()
    }


@dataclass
class Gene:
    """A gene of a reference's annotation, its first transcript: 1-based,
    inclusive coordinates; its product and note, from its CDS."""

    gene_id: str
    chromosome: str
    strand: str
    biotype: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)
    product: str = ""
    note: str = ""

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
            if fields[2] == "CDS":
                gene.product = attributes.get("product", gene.product)
                gene.note = attributes.get("note", gene.note)
    for gene in genes.values():
        gene.biotype = biotypes.get(gene.gene_id, gene.biotype)
        gene.exons.sort()
        gene.cds = merge_intervals(gene.cds)
    return sorted(
        (g for g in genes.values() if g.exons), key=lambda g: (g.chromosome, g.start)
    )


def scer_characterized(cache: Path) -> dict[str, bool]:
    """Each S. cerevisiae protein-coding gene, genome-wide: characterized, or
    its product "uncharacterized protein"."""
    known: dict[str, bool] = {}
    with gzip.open(download(SPECIES[0], cache, "genomic.gtf.gz"), "rt") as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            known.setdefault(
                attributes["gene_id"],
                attributes.get("product", "") != "uncharacterized protein",
            )
    return known


def gene_class(species: str, gene: Gene, scer: dict[str, bool]) -> tuple[str, str]:
    """Whether a reference gene is characterized, and why. S. cerevisiae's and
    S. paradoxus's annotations say so in their products ("uncharacterized
    protein", "hypothetical protein"); N. glabratus's call every gene
    uncharacterized, and their notes name its S. cerevisiae homolog, itself
    characterized or not."""
    if species == "Cgla":
        match = SCER_GENE.search(gene.note)
        if not match:
            return "uncharacterized", "no S. cerevisiae homolog"
        homolog = match.group(1).upper()
        if homolog not in scer:
            return "uncharacterized", f"S. cerevisiae homolog {homolog}, not annotated"
        verdict = "characterized" if scer[homolog] else "uncharacterized"
        return verdict, f"S. cerevisiae homolog {homolog}, {verdict}"
    unknown = {"Scer": "uncharacterized protein", "Spar": "hypothetical protein"}
    verdict = "uncharacterized" if gene.product == unknown[species] else "characterized"
    return verdict, f"product: {gene.product}"


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


def diverged(sequence: str, rate: float, rng: random.Random) -> str:
    """A sequence with substitutions at ``rate`` per base, and 1-base
    insertions and deletions at a twentieth of it each."""
    if not rate:
        return sequence
    out: list[str] = []
    for base in sequence:
        roll = rng.random()
        if roll < rate:
            out.append(
                rng.choice("ACGT".replace(base, "") if base in "ACGT" else "ACGT")
            )
        elif roll < rate * 1.05:
            continue
        elif roll < rate * 1.1:
            out.append(base + rng.choice("ACGT"))
        else:
            out.append(base)
    return "".join(out)


def repeat_families(
    records: dict[str, str], rng: random.Random
) -> dict[str, dict[str, str]]:
    """Each species' repeat families: S288C's sequences (``records``, its
    assembly), diverged by FAMILY_DIVERGENCE."""
    sources = {
        family: records[accession][start - 1 : end]
        for family, (accession, start, end) in FAMILY_SOURCES.items()
    }
    return {
        species.name: {
            family: diverged(sequence, FAMILY_DIVERGENCE[species.name], rng)
            for family, sequence in sources.items()
        }
        for species in SPECIES
    }


def telomere_length(sequence: str, side: str) -> int:
    """The telomeric repeat at a chromosome's end: C1-3A at its left, TG1-3 at
    its right; 0 where the reference has none."""
    if side == "left":
        match = re.match(r"(?:C{1,3}A)+", sequence)
    else:
        match = re.search(r"(?:TG{1,3})+$", sequence)
    return len(match.group()) if match else 0


def its_tract(side: str, rng: random.Random) -> str:
    """The telomeric repeat between two tandem Y' elements: C1-3A on the left
    arms' strand, TG1-3 on the right arms'."""
    length = rng.randint(*ITS_LENGTH)
    units = ("CA", "CCA", "CCCA") if side == "left" else ("TG", "TGG", "TGGG")
    tract = ""
    while len(tract) < length:
        tract += rng.choice(units)
    return tract


@dataclass(frozen=True)
class Insertion:
    """A planted repeat: where it goes in the reference (0-based, before that
    base), its sequence, and what to record of it."""

    point: int
    sequence: str
    record: dict


def spread(sites: list[tuple[int, str]], n: int, length: int) -> list[tuple[int, str]]:
    """``n`` of the sites, each the nearest still free to its share of the
    chromosome."""
    free = list(sites)
    chosen = []
    for k in range(1, n + 1):
        if not free:
            break
        target = length * k // (n + 1)
        site = min(free, key=lambda s: abs(s[0] - target))
        free.remove(site)
        chosen.append(site)
    return sorted(chosen)


def plant_repeats(
    genome: dict[str, str],
    genes: list[Gene],
    level: Level,
    families: dict[str, str],
    rng: random.Random,
) -> tuple[dict[str, str], list[Gene], list[dict]]:
    """The genome with the level's repeats planted, the genes moved past each,
    and where each sits (1-based, inclusive)."""
    nuclear = [name for name in genome if not name.endswith("_chrM")]
    # Where no repeat lands: each gene, GENE_GAP around it; a tRNA gene's
    # body only, since Ty1 lands just upstream of it.
    spans = {
        name: merge_intervals(
            [
                (g.start - pad, g.end + pad)
                for g in genes
                if g.chromosome == name
                for pad in [0 if g.biotype == "tRNA" else GENE_GAP]
            ]
        )
        for name in nuclear
    }

    def between_genes(name: str, point: int) -> bool:
        return all(not (start <= point <= end) for start, end in spans[name])

    def gaps(name: str) -> list[tuple[int, str]]:
        """Midpoints of the gaps between genes: where a Ty1 lands when no tRNA
        gene is left."""
        return [
            ((left + right) // 2, "between genes")
            for (_, left), (right, _) in itertools.pairwise(spans[name])
            if right > left
        ]

    insertions: dict[str, list[Insertion]] = {name: [] for name in nuclear}
    ty1 = families["Ty1"]
    for name in nuclear:
        sequence = genome[name]
        length = len(sequence)
        inside = range(TY1_MARGIN, length - TY1_MARGIN)
        sites = []
        for gene in genes:
            if gene.chromosome != name or gene.biotype != "tRNA":
                continue
            distance = rng.randint(*TY1_UPSTREAM)
            point = (
                gene.start - 1 - distance if gene.strand == "+" else gene.end + distance
            )
            if point in inside and between_genes(name, point):
                sites.append((point, f"upstream of {gene.gene_id}"))
        wanted = 1 if level.ty1_spacing is None else max(1, length // level.ty1_spacing)
        wanted += 1 if level.ty1_cluster else 0
        chosen = spread(sites, wanted, length)
        if len(chosen) < wanted:
            spare = [g for g in gaps(name) if g[0] in inside and g not in chosen]
            chosen += spread(spare, wanted - len(chosen), length)
        chosen.sort()
        cluster = rng.randrange(len(chosen)) if level.ty1_cluster else None
        for index, (point, site) in enumerate(chosen):
            copies = level.ty1_cluster if index == cluster else 1
            strand = rng.choice("+-")
            elements = [
                diverged(ty1, COPY_DIVERGENCE["Ty1"], rng) for _ in range(copies)
            ]
            block = "".join(e if strand == "+" else revcomp(e) for e in elements)
            # The target site, duplicated: on both sides of the element.
            duplicate = sequence[point - TY1_TARGET_DUPLICATION : point]
            insertions[name].append(
                Insertion(
                    point,
                    block + duplicate,
                    {
                        "family": "Ty1",
                        "copies": copies,
                        "unit_length": len(ty1),
                        "tandem": copies > 1,
                        "strand": strand,
                        "site": site,
                    },
                )
            )
    ends = [(name, side) for name in nuclear for side in ("left", "right")]
    for name, side in sorted(rng.sample(ends, round(level.yprime_ends * len(ends)))):
        copies = rng.randint(*level.yprime_copies)
        units = []
        for _ in range(copies):
            unit = diverged(families["Y'"], COPY_DIVERGENCE["Y'"], rng)
            # As at chrXIV's left end: the Y', then the telomeric repeat on its
            # centromere's side; mirrored on a right arm.
            units.append(
                unit + its_tract(side, rng)
                if side == "left"
                else its_tract(side, rng) + revcomp(unit)
            )
        telomere = telomere_length(genome[name], side)
        point = telomere if side == "left" else len(genome[name]) - telomere
        insertions[name].append(
            Insertion(
                point,
                "".join(units),
                {
                    "family": "Y'",
                    "copies": copies,
                    "unit_length": len(families["Y'"]),
                    "tandem": copies > 1,
                    "strand": "+" if side == "left" else "-",
                    "site": f"{side} end",
                },
            )
        )
    if level.rdna_units:
        name = max(nuclear, key=lambda n: len(genome[n]))
        target = round(RDNA_PLACE * len(genome[name]))
        planted_points = [i.point for i in insertions[name]]
        candidates = [
            point
            for point, _ in gaps(name)
            if all(abs(point - p) > TY1_MARGIN for p in planted_points)
        ]
        point = min(candidates, key=lambda p: abs(p - target))
        insertions[name].append(
            Insertion(
                point,
                "".join(
                    diverged(families["rDNA"], COPY_DIVERGENCE["rDNA"], rng)
                    for _ in range(level.rdna_units)
                ),
                {
                    "family": "rDNA",
                    "copies": level.rdna_units,
                    "unit_length": len(families["rDNA"]),
                    "tandem": True,
                    "strand": "+",
                    "site": "between genes",
                },
            )
        )

    planted = dict(genome)
    repeats: list[dict] = []
    for name in nuclear:
        events = sorted(insertions[name], key=lambda i: i.point)
        points = [event.point for event in events]
        assert len(set(points)) == len(points), f"two repeats at one place on {name}"
        sequence = genome[name]
        # Planted from the chromosome's end, so that each point stands.
        for event in reversed(events):
            sequence = (
                sequence[: event.point] + event.sequence + sequence[event.point :]
            )
            genes = [
                g.shifted(event.point, len(event.sequence))
                if g.chromosome == name
                else g
                for g in genes
            ]
        planted[name] = sequence
        offset = 0
        for event in events:
            start = event.point + offset + 1
            end = start + len(event.sequence) - 1
            assert sequence[start - 1 : end] == event.sequence
            repeats.append(
                {"chromosome": name, "start": start, "end": end, **event.record}
            )
            offset += len(event.sequence)
    return planted, genes, repeats


def systematic_errors(
    genome: dict[str, str], repeats: list[dict], spacing: int, rng: random.Random
) -> tuple[dict[str, str], list[tuple[str, int, str]]]:
    """The genome as its nanopore reads all see it: one error every
    ``spacing`` bases, outside the planted repeats; and each error (the true
    genome's sequence, 1-based position, kind)."""
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
        n_errors = round(len(sequence) / spacing)
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


def write_gtf(
    genes: list[Gene], classes: dict[str, tuple[str, str]], path: Path
) -> None:
    """The genes, as GTF: a gene, a transcript, its exons and CDS; a
    protein-coding gene's class (characterized or not)."""
    with path.open("w", encoding="ascii") as sink:
        for g in genes:
            ids = f'gene_id "{g.gene_id}"; transcript_id "{g.gene_id}.t1";'
            attributes = f'{ids} gene_biotype "{g.biotype}";'
            if g.gene_id in classes:
                attributes += f' gene_class "{classes[g.gene_id][0]}";'
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


def nanopore_reads(
    species: Species, level: Level, ref: Path, seed: int, workdir: Path
) -> list[str]:
    """Badread's reads of a species, R10.4.1's error and quality models, the
    level's identity and lengths, as FASTQ records."""
    out = workdir / f"{species.name}.badread.fastq"
    with out.open("w") as sink:
        subprocess.run(
            [
                "badread", "simulate",
                "--reference", str(ref),
                "--quantity", f"{level.ont_depth:g}x",
                "--error_model", "nanopore2023",
                "--qscore_model", "nanopore2023",
                "--identity", ",".join(f"{x:g}" for x in level.ont_identity),
                "--length", ",".join(str(x) for x in level.ont_length),
                "--seed", str(seed),
            ],
            check=True, stdout=sink, stderr=subprocess.DEVNULL,
        )  # fmt: skip
    lines = out.read_text().splitlines()
    return ["\n".join(lines[i : i + 4]) + "\n" for i in range(0, len(lines), 4)]


def n50(lengths: list[int]) -> int:
    """The length of the read at which half the bases are in reads as long."""
    total = sum(lengths)
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if 2 * running >= total:
            return length
    return 0


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


def hic_pairs(
    genomes: dict[str, dict[str, str]], level: Level, seed: int, out: Path
) -> None:
    """The pool's Hi-C pairs, written to ``<out>_R1.fastq.gz`` and
    ``<out>_R2.fastq.gz`` as they are drawn."""
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

    raw_1, handle_1 = gzip_writer(Path(f"{out}_R1.fastq.gz"))
    raw_2, handle_2 = gzip_writer(Path(f"{out}_R2.fastq.gz"))
    with raw_1, handle_1, raw_2, handle_2:
        batch_1: list[str] = []
        batch_2: list[str] = []
        made = 0
        last = None
        while made < level.hic_pairs:
            if last is not None and rng.random() < level.hic_duplicates:
                read_1, read_2, label = last
            else:
                species = rng.choices([s.name for s in SPECIES], weights=weights)[0]
                name_a, index_a = random_site(species)
                insert = max(READ_LENGTH + 20, int(rng.gauss(350, 80)))
                roll = rng.random()
                if roll < level.hic_unligated:
                    start = sites[name_a].starts[index_a]
                    molecule = sequences[name_a][
                        max(0, start - insert // 2) : start + insert // 2
                    ]
                    label = f"{name_a}:{start}:unligated"
                else:
                    if roll < level.hic_unligated + level.hic_random_ligation:
                        # Random ligation: a partner from the whole pool, its
                        # own species included, by abundance. The pairs
                        # between species are then epsilon (1 - sum of p
                        # squared) of the ligated.
                        other = rng.choices([s.name for s in SPECIES], weights=weights)[
                            0
                        ]
                        name_b, index_b = random_site(other)
                    elif (
                        roll
                        < level.hic_unligated
                        + level.hic_random_ligation
                        + level.hic_trans
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
                    part_b = start_at(
                        sequences[name_b], site_b, enzyme_b, insert - split
                    )
                    # Either strand of the other fragment: both enzymes' sites
                    # are palindromes.
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
            batch_1.append(
                f"@hic_{made:07d} {label}\n{read_1}\n+\n{quality[: len(read_1)]}\n"
            )
            batch_2.append(
                f"@hic_{made:07d} {label}\n{read_2}\n+\n{quality[: len(read_2)]}\n"
            )
            made += 1
            if len(batch_1) == 20_000 or made == level.hic_pairs:
                handle_1.write("".join(batch_1).encode("ascii"))
                handle_2.write("".join(batch_2).encode("ascii"))
                batch_1.clear()
                batch_2.clear()


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


def make_set(
    level: Level,
    outdir: Path,
    threads: int,
    references: dict[str, dict[str, str]],
    annotations: dict[str, list[Gene]],
    classes: dict[str, dict[str, tuple[str, str]]],
    families: dict[str, dict[str, str]],
) -> None:
    """One set, in ``<outdir>/<level>/``."""
    root = outdir / level.name
    if root.exists():
        shutil.rmtree(root)
    for sub in ("ont", "pool", "references", "truth"):
        (root / sub).mkdir(parents=True, exist_ok=True)
    rng = random.Random(f"plexembler-benchmark-v{GENERATOR_VERSION}-{level.name}")

    genomes: dict[str, dict[str, str]] = {}
    genes: dict[str, list[Gene]] = {}
    repeats: dict[str, list[dict]] = {}
    for species in SPECIES:
        original = references[species.name]
        genome, genes[species.name], repeats[species.name] = plant_repeats(
            original, annotations[species.name], level, families[species.name], rng
        )
        # No repeat in a gene: each gene's sequence is the reference's.
        before = {g.gene_id: g for g in annotations[species.name]}
        for g in genes[species.name]:
            old = before[g.gene_id]
            assert spliced(genome[g.chromosome], g.exons, g.strand) == spliced(
                original[old.chromosome], old.exons, old.strand
            ), f"a repeat planted in {g.gene_id}"
        genomes[species.name] = genome
        write_reference(genome, root / "references" / f"{species.name}.fa")
        write_gtf(
            genes[species.name],
            classes[species.name],
            root / "references" / f"{species.name}.gtf",
        )

    # The nuclear protein-coding genes: their CDS (the M4 check finds them among
    # the pipeline's, the characterized ones counted) and their transcripts (the
    # RNA-seq is drawn from them).
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
        with (root / "references" / f"{species.name}.cds.fa").open("w") as sink:
            for g in coding:
                sink.write(
                    f">{g.gene_id}\n{spliced(genome[g.chromosome], g.cds, g.strand)}\n"
                )
        with (root / "truth" / f"{species.name}.genes.tsv").open("w") as sink:
            sink.write("gene\tchromosome\tstart\tend\tstrand\tcds_bp\tclass\tbasis\n")
            for g in coding:
                verdict, basis = classes[species.name][g.gene_id]
                cds_bp = sum(e - s + 1 for s, e in g.cds)
                sink.write(
                    f"{g.gene_id}\t{g.chromosome}\t{g.start}\t{g.end}\t{g.strand}\t"
                    f"{cds_bp}\t{verdict}\t{basis}\n"
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
                genomes[species.name],
                repeats[species.name],
                level.ont_error_spacing,
                rng,
            )
            write_reference(seen, work / f"{species.name}.fa", depth=True)
            with (root / "truth" / f"{species.name}.ont_errors.tsv").open("w") as sink:
                sink.write("sequence\tposition\tkind\n")
                for name, position, kind in errors[species.name]:
                    sink.write(f"{name}\t{position}\t{kind}\n")
        with ThreadPoolExecutor(max_workers=threads) as pool:
            jobs = {
                s.name: pool.submit(
                    nanopore_reads, s, level, work / f"{s.name}.fa", 1000 + i, work
                )
                for i, s in enumerate(SPECIES)
            }
            reads = {name: job.result() for name, job in jobs.items()}
        bleed = {}
        ont = {}
        for species in SPECIES:
            source, share = BLEED_THROUGH[species.name]
            extra = rng.sample(reads[source], round(share * len(reads[species.name])))
            bleed[species.name] = {"from": source, "share": share, "reads": len(extra)}
            lengths = [len(r.split("\n", 2)[1]) for r in reads[species.name]]
            size = sum(len(x) for x in genomes[species.name].values())
            ont[species.name] = {
                "depth": round(sum(lengths) / size, 1),
                "n50": n50(lengths),
            }
            write_gzip(
                root / "ont" / f"{species.name}.fastq.gz", reads[species.name] + extra
            )
        del reads

        r1, r2 = wgs_reads(genomes, threads, 2000, work)
        for read, name in ((r1, "R1"), (r2, "R2")):
            with read.open() as handle:
                write_gzip(root / "pool" / f"wgs_{name}.fastq.gz", handle)

    hic_pairs(genomes, level, 3000, root / "pool" / "hic")

    levels = {s.name: expression_levels(transcripts[s.name], rng) for s in SPECIES}
    for index, (condition, replicate) in enumerate(RNA_LIBRARIES):
        label = library_label(condition, replicate)
        rna_1, rna_2 = rnaseq_pairs(transcripts, levels, label, 4000 + index)
        stem = f"rna-{condition.lower()}-{replicate}"
        write_gzip(root / "pool" / f"{stem}_R1.fastq.gz", rna_1)
        write_gzip(root / "pool" / f"{stem}_R2.fastq.gz", rna_2)
    labels = [library_label(c, r) for c, r in RNA_LIBRARIES]
    for species in SPECIES:
        # The true TPM: a gene's share of the species' molecules, per million.
        totals = {label: sum(levels[species.name][label].values()) for label in labels}
        with (root / "truth" / f"{species.name}.tpm.tsv").open("w") as sink:
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

    (root / "species.csv").write_text(
        "species,pool,ont,chemistry,genome_size,ploidy,n_chromosomes,busco_lineage,genetic_code,telomere_repeat\n"
        + "".join(
            f"{s.name},bench,ont/{s.name}.fastq.gz,r10,"
            f"{sum(len(x) for x in genomes[s.name].values())},1,{s.n_chromosomes},"
            "saccharomycetes_odb12,1,NA\n"
            for s in SPECIES
        )
    )
    (root / "libraries.csv").write_text(
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
        "level": asdict(level),
        "repeat_families": {
            family: {
                "source": f"S288C {accession}:{start}-{end}",
                "length": end - start + 1,
                "divergence": FAMILY_DIVERGENCE,
                "copy_divergence": COPY_DIVERGENCE[family],
            }
            for family, (accession, start, end) in FAMILY_SOURCES.items()
        },
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
                "characterized": sum(
                    1
                    for g in genes[s.name]
                    if g.biotype == "protein_coding"
                    and g.cds
                    and not g.chromosome.endswith("_chrM")
                    and classes[s.name][g.gene_id][0] == "characterized"
                ),
                "expressed": len(transcripts[s.name]),
                "repeats": repeats[s.name],
                "ont": ont[s.name],
                "ont_errors": len(errors[s.name]),
                "wgs_depth": s.wgs_depth,
                "hic_share": s.hic_share,
                "rnaseq_share": s.rnaseq_share,
                "bleed_through": bleed[s.name],
            }
            for s in SPECIES
        },
        "mito_copies": MITO_COPIES,
        "ont": {
            "depth": level.ont_depth,
            "length_mean": level.ont_length[0],
            "length_sd": level.ont_length[1],
            "identity": list(level.ont_identity),
            "error_spacing": level.ont_error_spacing,
            "homopolymer_share": ONT_HOMOPOLYMER_SHARE,
        },
        "hic": {
            "pairs": level.hic_pairs,
            "enzymes": list(SITES),
            "trans": level.hic_trans,
            "random_ligation": level.hic_random_ligation,
            "between_species_expected": level.hic_random_ligation
            * (1 - sum(s.hic_share**2 for s in SPECIES)),
            "unligated": level.hic_unligated,
            "duplicates": level.hic_duplicates,
        },
        "rnaseq": {
            "libraries": labels,
            "pairs": RNA_PAIRS,
            "read_length": RNA_READ_LENGTH,
            "strandedness": "reverse",
            "changed_in_heat": RNA_CHANGED,
        },
    }
    (root / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")


def main(outdir: Path, threads: int, levels: list[str]) -> None:
    cache = outdir / ".cache"
    records = {s.name: assembly(s, cache) for s in SPECIES}
    references = {s.name: reference(s, records[s.name]) for s in SPECIES}
    annotations = {s.name: annotation(s, cache) for s in SPECIES}
    scer = scer_characterized(cache)
    classes = {
        s.name: {
            g.gene_id: gene_class(s.name, g, scer)
            for g in annotations[s.name]
            if g.biotype == "protein_coding"
        }
        for s in SPECIES
    }
    # Each species' repeat families, the same in every set.
    families = repeat_families(
        records["Scer"], random.Random(f"plexembler-benchmark-v{GENERATOR_VERSION}")
    )
    del records
    for level in LEVELS:
        if level.name in levels:
            make_set(level, outdir, threads, references, annotations, classes, families)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "-o", "--outdir", type=Path, default=Path(__file__).parent / "data"
    )
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument(
        "-l",
        "--level",
        action="append",
        choices=[level.name for level in LEVELS],
        help="a set to write (repeatable; all three by default)",
    )
    arguments = parser.parse_args()
    main(
        arguments.outdir,
        arguments.threads,
        arguments.level or [level.name for level in LEVELS],
    )
