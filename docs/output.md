# Outputs

Everything goes under `--outdir`, in three places: one folder per **species**,
one per **pool**, and the run's own files. The tree below is the one the
pipeline is built towards ([implementation plan](plan.md)); a path marked
*(now)* is written today.

```text
<outdir>/
├── chromosomix.store.json   every species' map, genome, tracks and workspace, for chromosomix
├── species/<species>/
│   ├── assembly/      (now) <species>.assembly.fa.gz: the kept contigs, polished;
│   │                  the table of candidates
│   ├── genome/        <species>.fa.gz (+ .fai, .gzi), <species>.softmasked.fa.gz,
│   │                  <species>.agp, <species>.chrom.sizes, <species>.chrM.fa.gz
│   ├── annotation/    <species>.gtf.gz, <species>.gff3.gz, <species>.proteins.faa.gz,
│   │                  <species>.cds.fna.gz, <species>.functions.tsv, repeats, ncRNA, tRNA,
│   │                  telomeres.bed, centromeres.bed, chrM.gff3
│   ├── hic/           <species>.mcool, <species>.hic, <species>.pretext, <species>.pairs.gz,
│   │                  maps (PNG), contact decay
│   ├── expression/    <species>.counts.tsv, <species>.tpm.tsv: genes × RNA-seq libraries
│   ├── variation/     heterozygous sites, ploidy and per-chromosome depth
│   ├── tracks/        coverage of the RNA-seq, WGS and nanopore reads (bigWig)
│   ├── contigs/       (now) <species>.contigs.fa.gz: polished with its short reads,
│   │                  screened, its organelle named <species>_chrM
│   ├── polish/        (now) the short-read polishing's choice; the Hi-C reads' digestion
│   ├── candidates/    assembly/<assembler>/ and polish/ (now): every candidate, its
│   │                  checks; later every scaffolder's, with its map
│   ├── estimates/     (now) what the samplesheet left NA: the genome size (LRGE);
│   │                  the k-mer spectrum's fit (GenomeScope2)
│   ├── reads/         (now, with --save_ont_reads) the nanopore reads as assembled
│   ├── qc/            nanoplot/ and filtlong/ (now), OMArk, contamination, …
│   ├── <species>.chromosomix.json   the species' view: map, genome, genes, tracks
│   └── <species>.report.html, .report.txt, .report_mqc.yaml
├── pools/<pool>/
│   ├── qc/            fastqc/ (now): the pooled libraries
│   ├── ani/           (now) the identity between the pool's species (skani)
│   ├── deconvolution/ (now) each library split between the species: their shares,
│   │                  ambiguous pairs, pairs between species, epsilon
│   ├── screen/        (now) each contig's verdict: kept, moved, dropped, organelle
│   ├── hic/           the pool's map, a .metacool whose genomes are its species
│   └── comparative/   orthogroups, species tree, synteny
├── multiqc/           multiqc_report.html (now): every step, every species and pool
└── pipeline_info/     (now) software versions, Nextflow's reports, parameters
```

The results open in [chromosomix](https://chromosomix.app) as they are: connect
the results folder (or open `chromosomix.store.json` once the results are on a
bucket), and each species is a card with its view.

## Now

### Nanopore reads

A species' runs become one FASTQ (Dorado's BAMs converted, their tags kept),
filtered by length and quality (chopper) and capped at `--ont_max_depth` times
the genome size (Filtlong, the longest and best reads kept).

- `species/<species>/qc/nanoplot/`: [NanoPlot](https://github.com/wdecoster/NanoPlot)'s
  report and statistics of the reads as they came, its runs together
  (`<species>_NanoPlot-report.html`, `<species>_NanoStats.txt`). Their median
  quality tells the chemistry when the samplesheet's is NA.
- `species/<species>/qc/filtlong/<species>.log`: what the cap kept.
- `species/<species>/estimates/`, for what the samplesheet left NA:
  `<species>.lrge.txt`, the genome size [LRGE](https://github.com/mbhall88/lrge)
  estimated from the reads, in bases, and its log with its 95 % interval; and,
  for every species, `genomescope2/`, [GenomeScope2](https://github.com/tbenavi1/genomescope2.0)'s
  fit of the nanopore reads' k-mer spectrum as a diploid's (`<species>_summary.txt`,
  plots). Read it with care: nanopore errors look like heterozygous k-mers, so
  a haploid species fits with a heterozygosity of about its reads' error rate,
  1 to 3 %; a heterozygous species, above. Milestone 2 fits each species' own
  short reads instead. A spectrum it cannot fit leaves no file, not a failed
  run.
- `species/<species>/reads/<species>.ont.fastq.gz`, with `--save_ont_reads`:
  the reads as they were assembled.

### Contigs

Every assembler of `--assemblers` makes a candidate; each is named after its
species ([`plexembler rename`](cli.md#rename)), polished with the species'
nanopore reads (Medaka, or `dorado polish`), and checked; one is kept by
[`plexembler select`](cli.md#select), or by the samplesheet's `assembler`.

- `species/<species>/assembly/`: the kept contigs, polished,
  `<species>.assembly.fa.gz`; the table of candidates,
  `<species>.assembly.selection.tsv` ([format](formats.md#selection-table)),
  its MultiQC table and the command's log.
- `species/<species>/candidates/assembly/<assembler>/`, for each candidate:
    - `<species>.<assembler>.contigs.fa`, its contigs named after the species,
      before polishing, and their names' table (`.contigs.names.tsv`,
      [format](formats.md#names-table));
    - `<species>.<assembler>.polished.fa.gz`, the polished contigs;
    - the assembler's own files: Flye's graph (`.assembly_graph.gfa.gz`),
      contig table (`.assembly_info.txt`), log and parameters; hifiasm's
      primary contigs' graph (`.bp.p_ctg.gfa`) and log;
    - `qc/gfastats/`: contiguity, [gfastats](https://github.com/vgl-hub/gfastats);
    - `qc/busco/`: [BUSCO](https://busco.ezlab.org)'s short summary
      (`short_summary.*.txt`, `.json`), batch summary and log;
    - `qc/merqury/`: [Merqury](https://github.com/marbl/merqury)'s QV
      (`<species>.<assembler>.wgs.qv` against the pool's WGS k-mers, or `.ont.qv`
      against the species' nanopore k-mers when the pool has no WGS library),
      its completeness against the species' nanopore k-mers
      (`.ont.completeness.stats`), and its spectra.

### Identity between the species of a pool

`pools/<pool>/ani/`, for a pool of two species or more:
[skani](https://github.com/bluenote-1577/skani)'s pairs of the kept contigs
(`<pool>.tsv`), and [`plexembler ani`](cli.md#ani)'s matrix (`<pool>.ani.tsv`),
pairs with their verdicts (`<pool>.ani_pairs.tsv`,
[format](formats.md#identity-pairs)), MultiQC table and log. A pair above
95 % is too close for the pool's short reads to be split well, and the run
says so.

### Pooled libraries

`pools/<pool>/qc/fastqc/`: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)'s
report of each read of each pooled library (`<library>_1_fastqc.html`,
`<library>_2_fastqc.html`), the library's runs together.

### The pool, split

`pools/<pool>/deconvolution/contigs/`: each WGS and Hi-C library of the pool,
aligned to every species' contigs at once (bwa-mem2) and split between the
species by [`plexembler split`](cli.md#split): `<pool>.<library>.split.tsv`,
each species' pairs, share, held-out pairs and duplicates;
`<pool>.<library>.split.json`, every number, with the ambiguous pairs of each
pair of species, the pairs between species and, for Hi-C, the random-ligation
rate epsilon; `<pool>.<library>.contacts.tsv` (Hi-C), the unique pairs between
contigs; and the command's log. The split reads are not kept: each step that
needs them makes them again from the alignments.

### Short-read polishing

Each species' contigs, polished with its own short reads: its WGS pairs, or,
when its pool has no WGS library, its Hi-C reads cut at their junctions
([`plexembler digest`](cli.md#digest)). Polypolish, then pypolca; each kept only
if its QV, against the k-mers of the pairs held out of the polishing, beats the
draft's.

- `species/<species>/polish/`: the choice, `<species>.polish.selection.tsv`
  ([format](formats.md#selection-table)), the kept sequence
  (`<species>.polish.fa.gz`), and, for Hi-C, the digestion's numbers
  (`<species>.digest.json`: the reads with 0, 1 and 2 or more junctions, a
  measure of the library's ligation).
- `species/<species>/candidates/polish/`: Polypolish's and pypolca's
  sequences, and their QVs and the draft's (`qc/merqury/*.holdout.qv`).

### The contig screen

[`plexembler screen`](cli.md#screen) weighs each contig's nanopore depth by
species, its identity to the other species' contigs and its Hi-C contacts.

- `pools/<pool>/screen/<pool>.screen.tsv`: each contig's signals, verdict
  (`keep`, `move`, `drop`, `organelle`), new name and reason.
- `species/<species>/contigs/<species>.contigs.fa.gz`: the species' contigs
  after the screen, those moved in from another species included, its
  mitochondrion named `<species>_chrM`: what milestone 3 scaffolds.

### MultiQC

`multiqc/multiqc_report.html`: the reads' QC, each candidate's BUSCO and
gfastats numbers, the table of candidates of every species, the identity
between the species of every pool, each library's split, each species'
polishing and its Hi-C reads' digestion, the contig screen, and the software
versions, in one report,
with the methods' description and its references.

### Pipeline information

`pipeline_info/`: Nextflow's execution report, timeline, trace and DAG; the
parameters of the run; the software versions
(`plexembler_software_mqc_versions.yml`).
