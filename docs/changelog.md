# Changelog

Every change a user can notice, by version, in
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/)'s format. Versions
follow [semantic versioning](https://semver.org).

## Unreleased

### Added

- **The repository**: the pipeline on nf-core's template (nf-core/tools 4.2.2)
  and the Python package on chromokit's (CONVENTIONS.md), with their docs, CI
  and tests.
- **Two samplesheets**: the species and their nanopore runs (`--input`), and the
  pooled Hi-C, WGS and RNA-seq libraries (`--libraries`), checked against each
  other. A relative path in them is read from the samplesheet's folder.
- **The species' nanopore chemistry**, a `chemistry` column: `r10`, `r9`, or
  `NA` to tell it from the reads. ([#1](https://github.com/js2264/plexembler/pull/1))
- **Unknown values**: `NA`, or an empty cell, in a species' `genome_size`,
  `ploidy`, `n_chromosomes`, `busco_lineage`, `genetic_code` and
  `telomere_repeat`; the step that needs one estimates it. A library's ID or
  condition that looks like a number (`37`) is read as text. ([#1](https://github.com/js2264/plexembler/pull/1))
- **The package's image**, `ghcr.io/js2264/plexembler`, published at each
  release, which the pipeline's own steps run in; PyPI's publishing is opt-in.
  ([#1](https://github.com/js2264/plexembler/pull/1))
- **The reads' QC**: NanoPlot for each species' nanopore reads, FastQC for each
  pooled library, and MultiQC.
- **`plexembler rename`**, which names each contig after its species
  (`<species>_<contig>`), and the command line's log, beside its outputs.
- **The test profile**: two species of a few dozen kilobases and their pool,
  written by `tests/data/make_test_data.py`.
- **Example samplesheets** of ten species in three pools, each pool sequenced
  differently (`assets/examples/`).
- **The implementation plan**, `docs/plan.md`.
- **Each species' contigs** (milestone 1): its nanopore reads as one FASTQ,
  filtered (chopper) and capped at `--ont_max_depth` (Filtlong); every
  assembler of `--assemblers` (Flye; hifiasm in its `--ont` mode) a candidate,
  named after the species, polished (`--ont_polisher`: Medaka, or `dorado
  polish`, experimental) and checked (gfastats, BUSCO, Merqury); one kept,
  `species/<species>/assembly/<species>.assembly.fa.gz`, with its table of
  candidates.
- **The species samplesheet's `ont` takes Dorado's unaligned BAMs**, their
  tags kept, and a new column, `assembler`, names the candidate to keep.
- **What the samplesheet leaves NA, estimated**: the chemistry from the reads'
  median quality, the genome size by LRGE; and GenomeScope2's fit of each
  species' nanopore k-mers, which says whether it looks heterozygous.
- **The QV of each candidate** against its pool's WGS k-mers, or its own
  nanopore k-mers when the pool has no WGS library; its completeness against
  its own nanopore k-mers.
- **`plexembler select`**, which keeps one candidate of a step by its checks
  and writes why; **`plexembler ani`**, the identity between the species of a
  pool (skani), with a warning above 95 %.
- **The test data**, version 2: 30× of nanopore reads per species and 10,000
  WGS pairs, so that the test profile assembles; its second species leaves its
  chemistry and genome size to the pipeline.
