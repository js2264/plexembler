# plexembler

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![Pipeline](https://github.com/js2264/plexembler/actions/workflows/pipeline.yml/badge.svg)](https://github.com/js2264/plexembler/actions/workflows/pipeline.yml)
[![CI](https://github.com/js2264/plexembler/actions/workflows/ci.yml/badge.svg)](https://github.com/js2264/plexembler/actions/workflows/ci.yml)
[![Docs](https://github.com/js2264/plexembler/actions/workflows/docs.yml/badge.svg)](https://js2264.github.io/plexembler/)
[![Codecov](https://codecov.io/gh/js2264/plexembler/branch/main/graph/badge.svg)](https://codecov.io/gh/js2264/plexembler)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Reference genomes, gene annotations, Hi-C maps and expression tables for every
species of a pool, from per-species nanopore reads and pooled Hi-C, WGS and
RNA-seq.

Each species of a pool has its own nanopore reads; the Hi-C, WGS and RNA-seq
libraries are made once, from the cells of every species mixed in one tube.
plexembler assembles each species from its nanopore reads, splits the pooled
libraries between the species by competitive mapping, and gives each species
its polished, screened and scaffolded genome, its gene annotation, its Hi-C
map, its expression table and its report.

> [!NOTE]
> plexembler is being built, milestone by milestone: today it reads and checks
> its samplesheets and the reads, and makes each species' contigs, polished,
> checked, and chosen among its assemblers'. [The implementation plan](https://js2264.github.io/plexembler/plan/)
> says what comes next.

## Install

The pipeline needs [Nextflow](https://www.nextflow.io) (26.04.0 or newer) and
a container engine (Docker, Apptainer, Singularity) or conda:

```bash
nextflow pull js2264/plexembler
```

## Quick start

```bash
nextflow run js2264/plexembler -profile test,docker --outdir results
```

```bash
nextflow run js2264/plexembler -profile docker \
    --input species.csv --libraries libraries.csv --outdir results
```

## Documentation

<https://js2264.github.io/plexembler/>: [concepts](https://js2264.github.io/plexembler/concepts/),
[usage](https://js2264.github.io/plexembler/usage/),
[outputs](https://js2264.github.io/plexembler/output/),
[the implementation plan](https://js2264.github.io/plexembler/plan/) and
[how to cite plexembler](https://js2264.github.io/plexembler/citation/).
