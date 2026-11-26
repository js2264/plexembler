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
  species' nanopore k-mers, a hint of its heterozygosity (nanopore errors
  inflate it).
- **The QV of each candidate** against its pool's WGS k-mers, or its own
  nanopore k-mers when the pool has no WGS library; its completeness against
  its own nanopore k-mers.
- **`plexembler select`**, which keeps one candidate of a step by its checks
  and writes why; **`plexembler ani`**, the identity between the species of a
  pool (skani), with a warning above 95 %.
- **The benchmark pool** (`-profile benchmark`): three yeasts, two of a genus,
  their nanopore reads and their pool's WGS and Arima Hi-C simulated from their
  references by `benchmarks/make_benchmark_data.py`, published as the release
  `benchmark-data-v1`; `benchmarks/evaluate.py` checks a milestone's "Done
  when" on its results.
- **The pool, split** (milestone 2): each pooled WGS and Hi-C library aligned
  to every species' contigs of its pool at once (bwa-mem2) and split between
  them by `plexembler split`: a pair to a species when both mates are unique
  (`--deconvolve_mapq`) and near-identical (`--deconvolve_min_identity`); the
  ambiguous pairs of each pair of species, the pairs between species and the
  Hi-C random-ligation rate (epsilon) counted, under
  `pools/<pool>/deconvolution/`.
- **Short-read polishing**: each species with its WGS pairs, or, in a pool
  without WGS, its Hi-C reads cut at their ligation junctions (`plexembler
  digest`); Polypolish then pypolca, kept only when the QV against pairs held
  out of the polishing rises (`plexembler select --step polish`);
  `--polish_reads`, `--sr_polisher`, `--polish_holdout`.
- **The contig screen**, `plexembler screen`: each contig kept, moved to
  another species of its pool, dropped, or named its species' mitochondrion
  (`<species>_chrM`), from each species' nanopore depth on it, its identity to
  the other species' contigs and its Hi-C contacts (which move a contig only
  when its own species' reads cover it thinly: `--hic-own-depth`); a moved
  contig is named after the species it goes to, its name kept; each species'
  contigs after it in `species/<species>/contigs/`.
- **The test data**, version 2: 30× of nanopore reads per species and 10,000
  WGS pairs, so that the test profile assembles; its second species leaves its
  chemistry and genome size to the pipeline.
- **Each species' genome** (milestone 3): the pool split again on the screened
  contigs; each species scaffolded with its own Hi-C pairs by YaHS, told its
  library's restriction motifs; `plexembler finalize`, which names the
  scaffolds, chromosomes first, longest first (`<species>_chr01`, …, the
  samplesheet's `n_chromosomes` of them or every scaffold of
  `--min_chromosome` bases or more), the mitochondrion (`<species>_chrM`) and
  the unplaced, and writes the genome, its sizes, the table of names and the
  AGP in the new names; the genome's gfastats, BUSCO and Merqury checks; all
  in `species/<species>/genome/`. A species of a pool without Hi-C is finished
  from its contigs.
- **`plexembler junctions`**, which scores each join of a species' scaffolds
  by the Hi-C pairs across it against the species' contact decay, from the
  pairs within its contigs: `supported`, `unsupported` (a misjoin, likely) or
  `thin` (`species/<species>/genome/qc/junctions/`).
- **Each species' Hi-C map** on its genome: its pairs aligned (bwa-mem2),
  parsed, sorted and deduplicated (pairtools), binned at `--map_resolution`
  (1 kb) and zoomed out, balanced (cooler): `species/<species>/hic/<species>.mcool`.
