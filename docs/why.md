# Why plexembler?

Building the genome of a species without a reference takes three kinds of
data: long reads to assemble it, Hi-C to put its contigs into chromosomes, and
short reads to polish it and to measure it. For ten species, that is thirty
libraries. plexembler takes the long reads of each species apart, and makes
each short-read library **once, for a pool of species**: one Hi-C library, one
WGS library and, if you want gene expression and better gene models, RNA-seq
libraries, each from the cells of all the species mixed in one tube.

From that, it makes for each species:

- a **genome reference**, assembled from the species' nanopore reads,
  polished, screened for contamination and scaffolded with its Hi-C pairs;
- its **gene annotation** (GTF and GFF3), with RNA-seq evidence when there is
  some, and its repeats, tRNAs, rDNA, telomeres, centromeres and mitochondrial
  genome;
- its **Hi-C map** (`.mcool`, `.hic`, `.pretext`);
- its **gene expression table**, one column per RNA-seq library;
- its **report**: what each step did, and how complete, correct and clean the
  genome is.

## What makes it different

- **The pool is split by competitive mapping, not by barcodes.** The contigs of
  every species of a pool are put together, and each pooled read goes to the
  species it maps best on ([concepts](concepts.md#deconvolution)).
- **The pool is its own control.** Hi-C pairs between species measure the
  library's random ligation; contigs whose reads and contacts point to another
  species are found and moved.
- **Tools compete.** Assemblers and scaffolders (YaHS and instaGRAAL among
  them) each make a candidate, and every candidate goes through the same
  checks before one is kept.
- **One pipeline for a workstation and a cluster**, built on nf-core's
  practices: containers, profiles, a checked samplesheet, tests.

## When to use something else

- **One species, with its own libraries:** a single-genome pipeline such as
  [nf-core/genomeassembler](https://nf-co.re/genomeassembler) or
  [sanger-tol/genomeassembly](https://github.com/sanger-tol/genomeassembly) is
  simpler.
- **A metagenome:** without a long-read run per species, the species have to be
  found from the contacts; [metaTOR](https://github.com/koszullab/metaTOR) does
  that.
- **Species closer than about 95 % identity:** the pool's short reads cannot
  tell them apart in their conserved regions. Pool species that are further
  apart, or sequence the close ones apart.
