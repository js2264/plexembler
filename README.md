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
> plexembler is being built, milestone by milestone: today it makes each
> species' contigs, splits the pooled libraries between the species, polishes
> and screens the contigs, scaffolds them into chromosomes with Hi-C, annotates
> their repeats and genes, and writes their expression tables, tracks, reports
> and a chromosomix store. Each milestone has steps still to come:
> [the implementation plan](https://js2264.github.io/plexembler/plan/) says
> which.

## Install

The pipeline needs [Nextflow](https://www.nextflow.io) (26.04.0 or newer) and
a container engine (Docker, Apptainer, Singularity) or conda:

```bash
nextflow pull js2264/plexembler
```

Or, on a Linux machine with nothing installed, let
[uv](https://docs.astral.sh/uv/) run it, with Nextflow and Java from PyPI and
the tools from conda-forge and bioconda (micromamba): `plexembler run` takes
`nextflow run`'s arguments.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh && . ~/.local/bin/env
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C ~/.local bin/micromamba
uvx --from 'plexembler[java]' plexembler run -profile test,micromamba --outdir results
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
