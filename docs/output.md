# Outputs

Everything goes under `--outdir`, in three places: one folder per **species**,
one per **pool**, and the run's own files. The tree below is the one the
pipeline is built towards ([implementation plan](plan.md)); a path marked
*(now)* is written today.

```text
<outdir>/
├── chromosomix.store.json   (now) every species' map, genome and tracks, for chromosomix
├── species/<species>/
│   ├── assembly/      (now) <species>.assembly.fa.gz: the kept contigs, polished;
│   │                  the table of candidates
│   ├── genome/        (now) <species>.fa.gz (+ .fai, .gzi), <species>.agp,
│   │                  <species>.chrom.sizes, <species>.names.tsv,
│   │                  <species>.softmasked.fa, qc/; later <species>.chrM.fa.gz
│   ├── annotation/    (now) <species>.gtf.gz, <species>.gff3.gz, <species>.proteins.faa.gz,
│   │                  <species>.cds.fna.gz, <species>.genes.tsv, repeats/, ncrna/, qc/;
│   │                  later <species>.functions.tsv, telomeres.bed, centromeres.bed, chrM.gff3
│   ├── rnaseq/        (now) <species>.<library>.bam (+ .bai): its pairs of each pooled
│   │                  RNA-seq library, on its genome
│   ├── hic/           (now) <species>.mcool, <species>.pairs.gz and their stats;
│   │                  later <species>.hic, <species>.pretext, maps (PNG), contact decay
│   ├── expression/    (now) <species>.counts.tsv, <species>.tpm.tsv: genes × RNA-seq libraries
│   ├── variation/     heterozygous sites, ploidy and per-chromosome depth
│   ├── tracks/        (now) coverage of the RNA-seq, by strand, and the nanopore reads'
│   │                  depth (bigWig); later the WGS'
│   ├── contigs/       (now) <species>.contigs.fa.gz: screened, its organelle named
│   │                  <species>_chrM; then polished with its short reads (polish/)
│   ├── polish/        (now) the short-read polishing's choice; the Hi-C reads' digestion
│   ├── candidates/    assembly/<assembler>/, polish/, scaffold/yahs/ and genes/ (now):
│   │                  every candidate, its checks; later every scaffolder's, with its map
│   ├── estimates/     (now) what the samplesheet left NA: the genome size (LRGE);
│   │                  the k-mer spectrum's fit (GenomeScope2)
│   ├── reads/         (now, with --save_ont_reads) the nanopore reads as assembled
│   ├── qc/            nanoplot/ and filtlong/ (now), OMArk, contamination, …
│   ├── <species>.chromosomix.json   the species' view: map, genome, genes, tracks (to come)
│   └── <species>.report.html, .report.txt, .report_mqc.yaml (now)
├── pools/<pool>/
│   ├── qc/            fastqc/ (now): the pooled libraries
│   ├── ani/           (now) the identity between the pool's species (skani)
│   ├── deconvolution/ (now) each library split between the species, on their
│   │                  contigs/ and on their screened/ contigs: their shares,
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
      before polishing, a circle the assembler wrote as copies of itself cut
      to one copy, and their names' table (`.contigs.names.tsv`: the copies
      cut, and which contigs are circular; [format](formats.md#names-table));
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

Each species' screened contigs (below), polished with its own short reads, of
the pool's second split, on the screened contigs: its WGS pairs, or, when its
pool has no WGS library, its Hi-C reads cut at their junctions
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
  mitochondrion named `<species>_chrM`: what is polished with its short reads
  (`species/<species>/polish/<species>.polish.fa.gz`), then scaffolded.

### The genome

The pool is split again, on the screened contigs
(`pools/<pool>/deconvolution/screened/`, as above), and each species is
scaffolded with its own Hi-C pairs: aligned to its contigs (bwa-mem2, each
mate on its own), then [YaHS](https://github.com/c-zhou/yahs), told the
library's restriction motifs. [`plexembler finalize`](cli.md#finalize) names
the scaffolds: the chromosomes, longest first, the samplesheet's
`n_chromosomes` of them, or every scaffold of `--min_chromosome` bases or more;
the mitochondrion; the unplaced. A species of a pool without Hi-C is finished
from its contigs, named the same way.

- `species/<species>/genome/`: the genome, `<species>.fa.gz` (bgzipped, with
  its `.fai` and `.gzi`); its sizes, `<species>.chrom.sizes`; the table of
  names, `<species>.names.tsv` ([format](formats.md#genome-names)); the AGP,
  `<species>.agp`, in the new names (scaffolded species); the command's log.
- `species/<species>/genome/qc/`: gfastats, BUSCO and Merqury, as for each
  candidate above; and `junctions/`, each join of the scaffolds scored by its
  Hi-C pairs against the species' contact decay
  ([`plexembler junctions`](cli.md#junctions): `<species>.junctions.tsv`, a
  row per join, `supported`, `unsupported` or `thin`).
- `species/<species>/candidates/scaffold/yahs/`: YaHS's scaffolds and AGP, in
  its names (`<species>.yahs_scaffolds_final.fa`, `.agp`), and its log.

The genome is what YaHS makes of the contigs, uncurated: a misjoin it makes
stays, the chromosomes are its longest scaffolds whatever they hold. The
joins' scores say where to look first, on the map.

### Hi-C maps

Each species' Hi-C pairs, aligned to its genome (bwa-mem2), parsed, sorted
and deduplicated ([pairtools](https://github.com/open2c/pairtools): each
mate's unique 5' alignment, MAPQ 30 or more), binned at `--map_resolution`
and zoomed out, balanced ([cooler](https://github.com/open2c/cooler)).

- `species/<species>/hic/<species>.mcool`: the map, every resolution from
  `--map_resolution` up, balanced; chromosomix and HiGlass open it.
- `species/<species>/hic/<species>.pairs.gz`: the unique pairs, in the
  [4DN pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md),
  and their numbers (`<species>.pairs.stat`); the parse's numbers, duplicates
  included (`<species>.pairsam.stat`).

### The pooled RNA-seq, split

Each pooled RNA-seq library, aligned by [STAR](https://github.com/alexdobin/STAR)
to its pool's genomes at once (introns capped at 5 kb), and split between the
species by [`plexembler split --type rnaseq`](cli.md#split): each species' pairs
are the alignments they were, on its genome already.

- `species/<species>/rnaseq/<species>.<library>.bam` (and `.bai`): the
  species' pairs of the library, sorted.
- `pools/<pool>/deconvolution/rnaseq/`: each library's split, as for the WGS
  and Hi-C libraries, and STAR's report of each (`star/<library>.Log.final.out`).

### Annotation

Each species' repeats, found de novo by
[RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler) and masked
soft by [RepeatMasker](https://www.repeatmasker.org) (`--repeats`); its genes,
predicted by [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) (`--fungus`)
with its RNA-seq alignments and `--orthodb_proteins`, or by
[AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) ab initio with
`--augustus_species` when it has neither; named and written by
[`plexembler genes`](cli.md#genes); its tRNAs and rRNAs.

- `species/<species>/annotation/`: the genes, `<species>.gff3.gz` and
  `<species>.gtf.gz`, named `<species>_000010`, … in the genome's order; the
  proteins, `<species>.proteins.faa.gz`, translated with the species'
  `genetic_code`; the CDS, `<species>.cds.fna.gz`; the table of genes,
  `<species>.genes.tsv`; its MultiQC table and the command's log.
- `species/<species>/annotation/repeats/`: RepeatModeler's families
  (`<species>.families.fa`, `.stk`, `.log`), RepeatMasker's table and its repeats (`<species>.out`,
  `.tbl`, `.gff`).
- `species/<species>/genome/<species>.softmasked.fa`: the genome, its repeats in
  lower case.
- `species/<species>/annotation/ncrna/`: [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE)'s
  tRNAs (`trnascan/<species>.tsv`, `.gff`, `.stats`) and
  [barrnap](https://github.com/tseemann/barrnap)'s rRNAs (`barrnap/<species>.gff`).
- `species/<species>/annotation/qc/busco/`: BUSCO on the proteins.
- `species/<species>/candidates/genes/`: the predictor's own files
  (`braker3/braker.gtf`, `.aa`, `.codingseq`, `.log`; `augustus/<species>.augustus.gff3`).

### Expression and tracks

Each species' read pairs of each pooled RNA-seq library, counted on its genes
by [featureCounts](https://subread.sourceforge.net) (exons, by gene,
`--rnaseq_strandedness`), and its tables written by
[`plexembler expression`](cli.md#expression).

- `species/<species>/expression/<species>.counts.tsv`: the pairs, genes by
  libraries, each library's column named `<condition>_<replicate>`, after each
  gene's length.
- `species/<species>/expression/<species>.tpm.tsv`: the TPM, within the
  species, from the counts and the genes' lengths.
- `species/<species>/expression/featurecounts/`: featureCounts' table and its
  summary.
- `species/<species>/tracks/`: [deepTools](https://deeptools.readthedocs.io)'
  bigWigs: each library's coverage of the genes on each strand, in counts per
  million (`<species>.<library>.forward.bigWig`, `.reverse.bigWig`), and the
  nanopore reads' depth on the genome (`<species>.ont.bigWig`).

### Reports

[`plexembler report`](cli.md#report) gathers what the pipeline made of each
species ([format](formats.md#the-report-of-a-species)).

- `species/<species>/<species>.report.html`: its page: the assembly's
  candidates and choice, the polishing, its share of each pooled library, the
  screen, its chromosomes, BUSCO and QV, its joins, genes and expression.
- `species/<species>/<species>.report.txt`: every number of it, for scripts.
- `species/<species>/<species>.report_mqc.yaml`: its row of MultiQC's General
  Statistics table.

### The browser

`chromosomix.store.json`, at the top of the results:
[chromosomix](https://chromosomix.app)'s store
([`plexembler browser`](cli.md#browser)), a card per species, its genome and
Hi-C map each with its genes and tracks placed, every file named by its path in
the results. Connect the results folder in chromosomix's app; or, once the
results are on a bucket, open chromosomix.app with `?store=<its URL>`, the
run's `--browser_base` its address.

### MultiQC

`multiqc/multiqc_report.html`: the reads' QC, each candidate's BUSCO and
gfastats numbers, the table of candidates of every species, the identity
between the species of every pool, each library's split (the first, on the
contigs), each species' polishing and its Hi-C reads' digestion, the contig
screen, each genome's BUSCO and gfastats numbers and its joins' verdicts, each
map's pairs, each RNA-seq library's split and STAR's report, each species'
genes and BUSCO on its proteins, featureCounts' assignments and each
library's pairs on genes, each species' row of the General Statistics (its
report), and the software versions, in one report,
with the methods' description and its references.

### Pipeline information

`pipeline_info/`: Nextflow's execution report, timeline, trace and DAG; the
parameters of the run; the software versions
(`plexembler_software_mqc_versions.yml`).
