"""The benchmark pool's small set: three yeasts, and the reads of their pool.

Two or three chromosomes and the mitochondrion of three yeasts with complete
references (docs/plan.md, T2): *Saccharomyces cerevisiae* S288C and its sister
species *S. paradoxus* CBS432, about 88 % identical, the pair of a genus whose
reads are hard to split; and *Nakaseomyces glabratus* (*Candida glabrata*)
CBS138, distant from both. For each species, its nanopore reads (Badread, R10's
error model), with a planted share of another species' reads (bleed-through);
for the pool, a WGS library (InSilicoSeq, NovaSeq's errors) and an Arima Hi-C
library (DpnII and HinfI) simulated here: each pair a ligation of two
restriction fragments' ends, its read through the junction when the molecule
is short, its contacts decaying with distance, some between chromosomes and a
few between species (random ligation), PCR duplicates and unligated molecules.

Every read's name says where it comes from, which the benchmark's checks
compare against. The set is deterministic: the same GENERATOR_VERSION writes the
same reads. It sits in ``<outdir>/small/``, with its samplesheets and the
references, and is published as the release ``benchmark-data-v<VERSION>``
(.github/workflows/benchmarks.yml).

    python benchmarks/make_benchmark_data.py --outdir benchmarks/data --threads 4

Needs badread (0.4.1) and iss (InSilicoSeq, 2.0.1) on PATH, and the network to
download the references from NCBI once.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import math
import random
import re
import shutil
import subprocess
import tempfile
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

#: Raised when what this script writes changes: it names the data's release.
GENERATOR_VERSION = 1
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
        wgs_depth=40,
        hic_share=0.4,
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
        wgs_depth=30,
        hic_share=0.35,
        n_chromosomes=3,
    ),
    Species(
        "Cgla",
        "Nakaseomyces glabratus CBS138",
        "GCF_000002545.3",
        "GCF/000/002/545/GCF_000002545.3_ASM254v2",
        {"NC_005967.2": "chrA", "NC_005968.1": "chrB", "NC_004691.1": "chrM"},
        ont_depth=40,
        wgs_depth=30,
        hic_share=0.25,
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
#: Hi-C: pairs, read length, and how the contacts fall.
HIC_PAIRS = 300_000
READ_LENGTH = 150
HIC_TRANS = 0.12  # between chromosomes of a species
HIC_RANDOM_LIGATION = 0.03  # a partner from anywhere in the pool, by abundance: epsilon
HIC_UNLIGATED = 0.05  # a molecule without a junction
HIC_DUPLICATES = 0.05  # PCR duplicates
#: Arima's enzymes: the site, and where the filled-in end stops (GATC|, GANT|).
SITES = {"DpnII": "GATC", "HinfI": "GA[ACGT]TC"}


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


def reference(species: Species, cache: Path) -> dict[str, str]:
    """The species' kept sequences, named ``<species>_<chromosome>``."""
    cache.mkdir(parents=True, exist_ok=True)
    name = species.ftp_path.rsplit("/", 1)[1]
    path = cache / f"{name}_genomic.fna.gz"
    if not path.exists():
        url = f"{NCBI}/{species.ftp_path}/{name}_genomic.fna.gz"
        with urllib.request.urlopen(url) as response, path.open("wb") as sink:
            shutil.copyfileobj(response, sink)
    with gzip.open(path, "rt") as handle:
        records = read_fasta(handle.read())
    return {
        f"{species.name}_{new}": records[old] for old, new in species.sequences.items()
    }


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


def nanopore_reads(species: Species, ref: Path, seed: int, workdir: Path) -> list[str]:
    """Badread's reads of a species, R10.4.1's error and quality models, as
    FASTQ records."""
    out = workdir / f"{species.name}.badread.fastq"
    with out.open("w") as sink:
        subprocess.run(
            [
                "badread", "simulate",
                "--reference", str(ref),
                "--quantity", f"{species.ont_depth}x",
                "--error_model", "nanopore2023",
                "--qscore_model", "nanopore2023",
                "--identity", "98,99.5,1",
                "--length", "9000,7000",
                "--seed", str(seed),
            ],
            check=True, stdout=sink, stderr=subprocess.DEVNULL,
        )  # fmt: skip
    lines = out.read_text().splitlines()
    return ["\n".join(lines[i : i + 4]) + "\n" for i in range(0, len(lines), 4)]


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


def mutate(read: str, rng: random.Random, rate: float = 0.002) -> str:
    return "".join(
        rng.choice("ACGT".replace(base, ""))
        if rng.random() < rate and base in "ACGT"
        else base
        for base in read
    )


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
        """A site at a distance drawn from P(s) ~ 1/s, from 1 kb to the chromosome."""
        here = sites[name].starts[index]
        length = len(sequences[name])
        distance = math.exp(rng.uniform(math.log(1_000), math.log(length)))
        target = here + rng.choice((-1, 1)) * distance
        target = min(max(target, 0), length - 1)
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


def main(outdir: Path, threads: int) -> None:
    small = outdir / "small"
    for sub in ("ont", "pool", "references"):
        (small / sub).mkdir(parents=True, exist_ok=True)
    genomes = {s.name: reference(s, outdir / ".cache") for s in SPECIES}
    for species in SPECIES:
        write_reference(
            genomes[species.name], small / "references" / f"{species.name}.fa"
        )

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for species in SPECIES:
            write_reference(
                genomes[species.name], work / f"{species.name}.fa", depth=True
            )
        with ThreadPoolExecutor(max_workers=threads) as pool:
            jobs = {
                s.name: pool.submit(
                    nanopore_reads, s, work / f"{s.name}.fa", 1000 + i, work
                )
                for i, s in enumerate(SPECIES)
            }
            reads = {name: job.result() for name, job in jobs.items()}
        rng = random.Random(f"plexembler-benchmark-v{GENERATOR_VERSION}")
        bleed = {}
        for species in SPECIES:
            source, share = BLEED_THROUGH[species.name]
            extra = rng.sample(reads[source], round(share * len(reads[species.name])))
            bleed[species.name] = {"from": source, "share": share, "reads": len(extra)}
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
                "ont_depth": s.ont_depth,
                "wgs_depth": s.wgs_depth,
                "hic_share": s.hic_share,
                "bleed_through": bleed[s.name],
            }
            for s in SPECIES
        },
        "mito_copies": MITO_COPIES,
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
