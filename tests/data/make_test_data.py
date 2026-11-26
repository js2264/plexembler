"""Make the test data: two small species, their nanopore reads, and a pool.

The pipeline's plumbing is tested on these files (conf/test.config, nf-test):
two species of a few dozen kilobases, each with telomeric repeats at its
chromosomes' ends; a few nanopore reads of each; and three paired-end
libraries made from the pool of both, Hi-C, WGS and RNA-seq, with a little
noise between the species, as a real pool has. Every read's name says where
it comes from (``species:chromosome:position:strand``), which tests of the
species' deconvolution compare against.

Run it from anywhere; it writes beside itself. It is deterministic: the same
version of this script writes the same bytes.

    python tests/data/make_test_data.py
"""

from __future__ import annotations

import argparse
import gzip
import math
import random
from dataclasses import dataclass
from pathlib import Path

#: Raised when what this script writes changes, so that tests can name it.
VERSION = 1
HERE = Path(__file__).resolve().parent
COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class Species:
    """A species of the pool: its chromosomes' lengths, GC and telomeric repeat."""

    name: str
    chromosomes: dict[str, int]
    gc: float
    telomere: str
    ont_depth: float


SPECIES = (
    Species(
        "sp1", {"chr1": 24_000, "chr2": 16_000}, gc=0.38, telomere="TGTGGG", ont_depth=8
    ),
    Species("sp2", {"chr1": 30_000}, gc=0.45, telomere="TTAGGG", ont_depth=8),
)


def revcomp(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def genome(species: Species, rng: random.Random) -> dict[str, str]:
    """Random chromosomes of the species' GC, capped by 300 bp of its repeat."""
    at, gc = (1 - species.gc) / 2, species.gc / 2
    chromosomes = {}
    for name, length in species.chromosomes.items():
        cap = (species.telomere * 60)[:300]
        body = "".join(rng.choices("ACGT", weights=(at, gc, gc, at), k=length - 600))
        chromosomes[name] = revcomp(cap) + body + cap
    return chromosomes


def mutate(sequence: str, rng: random.Random, rate: float) -> str:
    """Substitutions, insertions and deletions, each a third of ``rate``."""
    out = []
    for base in sequence:
        roll = rng.random()
        if roll < rate / 3:
            out.append(rng.choice("ACGT".replace(base, "")))
        elif roll < 2 * rate / 3:
            out.extend((base, rng.choice("ACGT")))
        elif roll >= rate:
            out.append(base)
    return "".join(out)


def quality(length: int, rng: random.Random, low: int, high: int) -> str:
    return "".join(chr(33 + rng.randint(low, high)) for _ in range(length))


def write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    """A gzipped FASTQ, without a timestamp in its header: the same bytes each time."""
    text = "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)
    with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(text.encode("ascii"))


def pick_position(genomes, rng):
    """A species, a chromosome and a position, each species by its genome's size."""
    flat = [
        (species, name, len(seq))
        for species, chromosomes in genomes.items()
        for name, seq in chromosomes.items()
    ]
    species, name, length = rng.choices(flat, weights=[n for *_, n in flat])[0]
    return species, name, length


def ont_reads(species: Species, chromosomes, rng):
    records = []
    total = sum(len(seq) for seq in chromosomes.values())
    target, made = species.ont_depth * total, 0
    while made < target:
        name = rng.choices(
            list(chromosomes), weights=[len(s) for s in chromosomes.values()]
        )[0]
        chromosome = chromosomes[name]
        length = min(len(chromosome), int(rng.lognormvariate(math.log(3_000), 0.5)))
        start = rng.randint(0, len(chromosome) - length)
        strand = rng.choice("+-")
        read = chromosome[start : start + length]
        read = mutate(read if strand == "+" else revcomp(read), rng, rate=0.03)
        records.append(
            (
                f"ont_{len(records):05d} {species.name}:{name}:{start}:{strand}",
                read,
                quality(len(read), rng, 10, 30),
            )
        )
        made += length
    return records


def pair(fragment: str, read_length: int) -> tuple[str, str]:
    """The two reads of a fragment: its start, and its end's reverse complement."""
    return fragment[:read_length], revcomp(fragment)[:read_length]


def wgs_pairs(genomes, rng, n=1_000, read_length=100):
    r1, r2 = [], []
    for index in range(n):
        species, name, length = pick_position(genomes, rng)
        size = int(rng.gauss(300, 30))
        start = rng.randint(0, length - size)
        fragment = genomes[species][name][start : start + size]
        a, b = pair(fragment, read_length)
        label = f"wgs_{index:05d} {species}:{name}:{start}:+"
        r1.append((label, mutate(a, rng, 0.005), quality(len(a), rng, 30, 40)))
        r2.append((label, mutate(b, rng, 0.005), quality(len(b), rng, 30, 40)))
    return r1, r2


def hic_pairs(genomes, rng, n=1_000, read_length=100, noise=0.05):
    """Contacts within a species, decaying with distance, and ``noise`` between species."""
    r1, r2 = [], []
    for index in range(n):
        species, name, length = pick_position(genomes, rng)
        pos1 = rng.randint(0, length - read_length)
        if rng.random() < noise:
            species2, name2, length2 = pick_position(genomes, rng)
            pos2 = rng.randint(0, length2 - read_length)
        elif rng.random() < 0.1 and len(genomes[species]) > 1:
            name2 = rng.choice([c for c in genomes[species] if c != name])
            species2, length2 = species, len(genomes[species][name2])
            pos2 = rng.randint(0, length2 - read_length)
        else:
            species2, name2 = species, name
            distance = int(math.exp(rng.uniform(math.log(200), math.log(length))))
            pos2 = min(
                max(0, pos1 + rng.choice((-1, 1)) * distance), length - read_length
            )
        a = genomes[species][name][pos1 : pos1 + read_length]
        b = revcomp(genomes[species2][name2][pos2 : pos2 + read_length])
        label = f"hic_{index:05d} {species}:{name}:{pos1}:+ {species2}:{name2}:{pos2}:-"
        r1.append((label, mutate(a, rng, 0.005), quality(len(a), rng, 30, 40)))
        r2.append((label, mutate(b, rng, 0.005), quality(len(b), rng, 30, 40)))
    return r1, r2


def rnaseq_pairs(genomes, rng, n=1_000, read_length=100):
    """Reads of a few genes per species, at different levels, from a stranded library."""
    genes = []
    for species, chromosomes in genomes.items():
        for name, seq in chromosomes.items():
            for start in range(2_000, len(seq) - 3_000, 6_000):
                genes.append(
                    (species, name, start, 1_500, rng.choice("+-"), rng.uniform(1, 20))
                )
    r1, r2 = [], []
    for index in range(n):
        species, name, start, length, strand, _ = rng.choices(
            genes, weights=[g[-1] for g in genes]
        )[0]
        transcript = genomes[species][name][start : start + length]
        if strand == "-":
            transcript = revcomp(transcript)
        size = int(rng.gauss(250, 25))
        offset = rng.randint(0, length - size)
        fragment = transcript[offset : offset + size]
        # dUTP: the first read is the transcript's reverse complement.
        b, a = pair(fragment, read_length)
        label = f"rna_{index:05d} {species}:{name}:{start}:{strand}"
        r1.append((label, mutate(a, rng, 0.005), quality(len(a), rng, 30, 40)))
        r2.append((label, mutate(b, rng, 0.005), quality(len(b), rng, 30, 40)))
    return r1, r2


def main(outdir: Path) -> None:
    rng = random.Random(f"plexembler-test-data-v{VERSION}")
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "truth").mkdir(exist_ok=True)
    genomes = {}
    for species in SPECIES:
        genomes[species.name] = genome(species, rng)
        with (outdir / "truth" / f"{species.name}.fa").open("w") as fasta:
            for name, seq in genomes[species.name].items():
                fasta.write(f">{name}\n")
                fasta.writelines(seq[i : i + 80] + "\n" for i in range(0, len(seq), 80))
        write_fastq(
            outdir / f"{species.name}.ont.fastq.gz",
            ont_reads(species, genomes[species.name], rng),
        )

    for library, make in (
        ("hic", hic_pairs),
        ("wgs", wgs_pairs),
        ("rnaseq", rnaseq_pairs),
    ):
        r1, r2 = make(genomes, rng)
        write_fastq(outdir / f"pool_{library}_R1.fastq.gz", r1)
        write_fastq(outdir / f"pool_{library}_R2.fastq.gz", r2)

    (outdir / "species.csv").write_text(
        "species,pool,ont,genome_size,ploidy,n_chromosomes,genetic_code,telomere_repeat\n"
        + "".join(
            f"{s.name},pool,{s.name}.ont.fastq.gz,{sum(s.chromosomes.values()) // 1000}k,"
            f"1,{len(s.chromosomes)},1,{s.telomere}\n"
            for s in SPECIES
        )
    )
    (outdir / "libraries.csv").write_text(
        "library,pool,type,fastq_1,fastq_2,condition,replicate,enzymes\n"
        "hic,pool,hic,pool_hic_R1.fastq.gz,pool_hic_R2.fastq.gz,,,DpnII\n"
        "wgs,pool,wgs,pool_wgs_R1.fastq.gz,pool_wgs_R2.fastq.gz,,,\n"
        "rna-ypd-1,pool,rnaseq,pool_rnaseq_R1.fastq.gz,pool_rnaseq_R2.fastq.gz,YPD,1,\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "-o", "--outdir", type=Path, default=HERE, help="where to write"
    )
    main(parser.parse_args().outdir)
