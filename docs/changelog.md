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
