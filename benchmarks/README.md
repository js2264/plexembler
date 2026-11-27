# The benchmark pool

What each milestone's "Done when" ([docs/plan.md](../docs/plan.md)) is checked
on: a pool of yeasts whose truth is known, simulated from their references.

## The small set

`make_benchmark_data.py` writes it to `data/small/`, from three references at
NCBI:

| Species | Strain | Reference | Kept | Size |
|---|---|---|---|---|
| `Scer` | *Saccharomyces cerevisiae* S288C | GCF_000146045.2 | chrI, chrIII, chrVI, chrM | 0.90 Mb |
| `Spar` | *S. paradoxus* CBS432 | GCA_002079055.1 | chrI, chrIII, chrVI, chrM | 0.89 Mb |
| `Cgla` | *Nakaseomyces glabratus* CBS138 | GCF_000002545.3 | chrA, chrB, chrM | 1.01 Mb |

`Scer` and `Spar` are a pair of a genus, about 88 % identical: their pooled
short reads are the hard ones to split (docs/plan.md, Decisions). `Cgla` is far
from both.

| Reads | How | Truth in the read's name |
|---|---|---|
| Nanopore, per species (`ont/<species>.fastq.gz`) | Badread 0.4.1, R10.4.1's error and quality models (`nanopore2023`), 40×, reads of 98 % identity on average; bleed-through: 15 % of `Spar`'s reads are `Cgla`'s (enough to assemble some), 1 % of the others' another species' | the sequence, strand and position |
| WGS, pooled (`pool/wgs_R{1,2}.fastq.gz`) | InSilicoSeq 2.0.1, NovaSeq's errors, 40× of `Scer`, 30× of the others | the sequence |
| Hi-C, pooled (`pool/hic_R{1,2}.fastq.gz`) | here: Arima's DpnII and HinfI, 300,000 pairs; each pair two restriction fragments' ends ligated, a read through the junction when the molecule is short; contacts decaying as 1/s; 12 % between chromosomes, 3 % random ligations (a partner from anywhere in the pool, by abundance: ε), 5 % unligated, 5 % PCR duplicates | both ends, their sites and enzymes |

Each species' mitochondrion is 5 copies per nuclear genome in every library.
`truth.json` records every number above; `species.csv` and `libraries.csv`
are the run's samplesheets, and `references/` the truth.

The set is deterministic: the same `GENERATOR_VERSION` writes the same reads.
It is published, with that version, as the release `benchmark-data-v<N>`, which
[the benchmarks' workflow](../.github/workflows/benchmarks.yml) downloads, or
makes and publishes when it does not exist.

```bash
# The simulators: Badread and InSilicoSeq
micromamba create -p .bench -c conda-forge -c bioconda python=3.12 badread=0.4.1 insilicoseq=2.0.1
PATH=$PWD/.bench/bin:$PATH python benchmarks/make_benchmark_data.py --outdir benchmarks/data --threads 4
```

## A run, and its checks

```bash
nextflow run . -profile benchmark,docker --outdir results
python benchmarks/evaluate.py --results results --milestone 1
```

`evaluate.py` runs BUSCO on the references and QUAST against them (`--busco`,
`--quast`: how to run them, on PATH by default) and writes
`results/benchmark/m<N>.md`; it exits 1 when a check fails (a check with
nothing to judge reads "not judged").

## To come

The set grows with the milestones, each growth a new `GENERATOR_VERSION`:
RNA-seq of each species' genes (milestones 4 and 5), and misjoins planted in a
draft (milestone 3). `test_full` is real data: public reads of five to ten
yeasts with complete references, mixed in silico (docs/plan.md, T2).
