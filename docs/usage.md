# Usage

```bash
nextflow run js2264/plexembler \
    -profile <docker/apptainer/singularity/conda/institute> \
    --input species.csv \
    --libraries libraries.csv \
    --outdir results
```

The pipeline takes two samplesheets: the **species** of the pools, with their
nanopore reads (`--input`), and the **pooled libraries**, Hi-C, WGS and
RNA-seq (`--libraries`). Every parameter is on the [parameters](parameters.md)
page; the results are on the [outputs](output.md) page. `plexembler run`
takes the same arguments, and runs the pipeline of its own version
([command line](cli.md#run)).

!!! note "Under development"

    The pipeline is being built milestone by milestone
    ([implementation plan](plan.md)). Today, it reads and checks both
    samplesheets, checks the reads (NanoPlot, FastQC, MultiQC), and makes each
    species' contigs (milestone 1): its nanopore reads filtered and capped,
    every assembler's contigs as candidates, polished and checked (BUSCO,
    Merqury, gfastats), one kept, and the identity between the species of each
    pool. The columns below are read and checked already; the steps that use the
    others come with their milestones.

## Paths in the samplesheets

A path in a samplesheet is a URL (`https://`, `s3://`, …), an absolute path,
or a path **relative to the samplesheet's folder**: a samplesheet beside the
data it lists names the files alone. The pipeline checks that each file exists
before anything runs.

!!! info "A deviation from nf-core, and why"

    nf-core's samplesheets read a relative path from the directory the pipeline
    is launched from. plexembler reads it from the samplesheet's folder, so that
    a samplesheet keeps working wherever it is launched from: its file paths are
    strings in the samplesheets' schemas, and the pipeline resolves and checks
    them itself (`subworkflows/local/utils_nfcore_plexembler_pipeline`).

## Species

One row per nanopore run of a species; the rows of a species are its runs, and
share its values. The [schema](https://github.com/js2264/plexembler/blob/main/assets/schema_input.json)
checks each column. What you do not know about a species' genome is `NA`, or an
empty cell: the step that needs it estimates it, and the report says so.
A species' runs may mix FASTQ files and Dorado's BAMs.

```csv title="species.csv"
species,pool,ont,genome_size,ploidy,n_chromosomes,genetic_code
Scer,poolA,ont/Scer.fastq.gz,12m,1,16,1
Calb,poolA,ont/Calb_run1.fastq.gz,14m,2,8,12
Calb,poolA,ont/Calb_run2.fastq.gz,14m,2,8,12
Huva,poolA,ont/Huva.fastq.gz,NA,NA,NA,NA
```

| Column | Required | Description | When it is `NA` |
|---|---|---|---|
| `species` | yes | The species' ID: a letter, then letters, digits and hyphens. It starts the name of every contig of the species, `<species>_<contig>`, and names its results. | |
| `pool` | no (`pool`) | The pool the species is in: the pooled libraries of the same pool are split between its species. | |
| `ont` | yes | The species' nanopore reads: FASTQ, gzipped or not, or Dorado's unaligned BAM (`.bam`), whose tags are kept. | |
| `chemistry` | no | The nanopore chemistry: `r10` (R10.4.1, Dorado SUP or HAC) or `r9` (R9.4.1 and older). It sets Flye's mode (`--nano-hq`, `--nano-raw`), Medaka's model and chopper's lowest quality. | told from the reads' median quality, NanoPlot's: Q15 and above is `r10` |
| `genome_size` | no | The expected genome size, in bases, with an optional `k`, `m` or `g` (`12m`). The nanopore reads are capped at `--ont_max_depth` times it. | estimated from the nanopore reads by LRGE (`species/<species>/estimates/`) |
| `ploidy` | no | The expected ploidy. A haploid species' assemblies keep no other haplotype's contigs (Flye's `--no-alt-contigs`, hifiasm's `-l0`). | haploid; GenomeScope2's fit of the nanopore k-mers hints at heterozygosity (nanopore errors inflate it), the species' own short reads tell it (milestone 2), and allele balance (milestone 6) |
| `n_chromosomes` | no | The expected number of nuclear chromosomes: the genome's longest scaffolds, that many, are named its chromosomes ([`plexembler finalize`](cli.md#finalize)); HapHiC's input (to come). | every scaffold of `--min_chromosome` bases (100 kb) or more is a chromosome |
| `busco_lineage` | no | A BUSCO lineage (`busco --list-datasets`), such as `saccharomycetes_odb12`. | `--busco_lineage`: `auto_euk` by default, BUSCO's choice |
| `genetic_code` | no | The NCBI translation table of the nuclear genes: `1`, `12` for the CUG-Ser clade, `26` for the CUG-Ala clade. The proteins are translated with it ([`plexembler genes`](cli.md#genes)). | `1`; to be inferred from the genome (Codetta, to come) |
| `telomere_repeat` | no | The telomeric repeat unit (`TGTGGG`). | found from the reads (milestone 3) |
| `assembler` | no | The assembler whose contigs to keep (`flye`, `hifiasm`), whatever the checks say: one of `--assemblers`. | the candidate the checks prefer ([`plexembler select`](cli.md#select)) |

## Pooled libraries

One row per pair of FASTQ files; the rows of a library are its runs (lanes),
and share its values. Every library belongs to a pool of the species
samplesheet. A pool without a Hi-C library, or without a WGS library, is
allowed: the steps that need one are skipped for its species, with a warning.
The [schema](https://github.com/js2264/plexembler/blob/main/assets/schema_libraries.json)
checks each column.

```csv title="libraries.csv"
library,pool,type,fastq_1,fastq_2,condition,replicate,enzymes
hic,poolA,hic,hic_R1.fastq.gz,hic_R2.fastq.gz,,,"DpnII,HinfI"
wgs,poolA,wgs,wgs_R1.fastq.gz,wgs_R2.fastq.gz,,,
rna-ypd-1,poolA,rnaseq,ypd1_R1.fastq.gz,ypd1_R2.fastq.gz,YPD,1,
rna-ypd-2,poolA,rnaseq,ypd2_R1.fastq.gz,ypd2_R2.fastq.gz,YPD,2,
```

| Column | Required | Description |
|---|---|---|
| `library` | yes | The library's ID: letters, digits, `_`, `.` and `-`. |
| `pool` | no (`pool`) | The pool the library was made from. |
| `type` | yes | `hic`, `wgs` or `rnaseq`. |
| `fastq_1`, `fastq_2` | yes | The pairs' two reads, as gzipped FASTQ. Pooled libraries are paired-end. |
| `condition`, `replicate` | no | An RNA-seq library's condition and replicate, which name the columns of the expression tables (the library's ID and `1` when absent). |
| `enzymes` | for `hic` | The Hi-C library's restriction enzymes, comma-separated, as Biopython names them (`DpnII,HinfI` for Arima's kit); `none` for a library made without one, such as Micro-C, which rules out instaGRAAL. Quote the value when it holds a comma. |

## An example: three pools

Ten species in three pools, each pool sequenced differently
([`assets/examples/`](https://github.com/js2264/plexembler/tree/main/assets/examples)):

```csv title="species.csv"
species,pool,ont,genome_size,ploidy,n_chromosomes,genetic_code
Scer,poolA,ont/Scer.fastq.gz,12m,1,16,1
Kmar,poolA,ont/Kmar_run1.fastq.gz,11m,1,8,1
Kmar,poolA,ont/Kmar_run2.fastq.gz,11m,1,8,1
Lthe,poolA,ont/Lthe.fastq.gz,10.4m,1,8,1
Huva,poolA,ont/Huva.fastq.gz,NA,NA,NA,NA
Calb,poolB,ont/Calb.fastq.gz,14m,2,8,12
Dhan,poolB,ont/Dhan.fastq.gz,12m,1,7,12
Ptan,poolB,ont/Ptan.fastq.gz,13m,1,NA,26
Spom,poolC,ont/Spom.fastq.gz,12.6m,1,3,1
Ylip,poolC,ont/Ylip.fastq.gz,20m,1,6,1
Pkud,poolC,ont/Pkud.fastq.gz,11m,2,5,1
```

```csv title="libraries.csv"
library,pool,type,fastq_1,fastq_2,condition,replicate,enzymes
A-hic,poolA,hic,poolA/hic_L001_R1.fastq.gz,poolA/hic_L001_R2.fastq.gz,,,"DpnII,HinfI"
A-hic,poolA,hic,poolA/hic_L002_R1.fastq.gz,poolA/hic_L002_R2.fastq.gz,,,"DpnII,HinfI"
A-wgs,poolA,wgs,poolA/wgs_R1.fastq.gz,poolA/wgs_R2.fastq.gz,,,
A-rna-ypd-1,poolA,rnaseq,poolA/rna_ypd1_R1.fastq.gz,poolA/rna_ypd1_R2.fastq.gz,YPD,1,
A-rna-ypd-2,poolA,rnaseq,poolA/rna_ypd2_R1.fastq.gz,poolA/rna_ypd2_R2.fastq.gz,YPD,2,
A-rna-heat-1,poolA,rnaseq,poolA/rna_heat1_R1.fastq.gz,poolA/rna_heat1_R2.fastq.gz,heat,1,
A-rna-heat-2,poolA,rnaseq,poolA/rna_heat2_R1.fastq.gz,poolA/rna_heat2_R2.fastq.gz,heat,2,
B-hic,poolB,hic,poolB/hic_R1.fastq.gz,poolB/hic_R2.fastq.gz,,,DpnII
B-wgs,poolB,wgs,poolB/wgs_R1.fastq.gz,poolB/wgs_R2.fastq.gz,,,
C-microc,poolC,hic,poolC/microc_R1.fastq.gz,poolC/microc_R2.fastq.gz,,,none
C-rna-ypd-1,poolC,rnaseq,poolC/rna_ypd1_R1.fastq.gz,poolC/rna_ypd1_R2.fastq.gz,YPD,1,
C-rna-ypd-2,poolC,rnaseq,poolC/rna_ypd2_R1.fastq.gz,poolC/rna_ypd2_R2.fastq.gz,YPD,2,
```

- `Kmar` has two nanopore runs, and `A-hic` two lanes: rows that share an ID
  are merged.
- `Calb` and `Dhan` (CUG-Ser clade) and `Ptan` (CUG-Ala) translate CUG
  differently: their `genetic_code` is 12 or 26.
- `NA` (or an empty cell) is what is not known: `Huva`'s genome size, ploidy,
  chromosome count and genetic code, `Ptan`'s chromosome count.

What each pool's species get, once the pipeline is complete
([plan](plan.md)):

| Pool | Hi-C | WGS | RNA-seq | Its species get |
|---|---|---|---|---|
| poolA | Arima (DpnII, HinfI), two lanes | yes | 2 conditions × 2 | everything: polished and scaffolded genomes, Hi-C maps, annotation with RNA evidence, expression tables |
| poolB | DpnII | yes | none | scaffolded genomes and maps; genes from protein evidence only; no expression table |
| poolC | Micro-C (`none`) | none | 1 condition × 2 | genomes polished with nanopore reads only, and a rougher QV, from nanopore k-mers; YaHS only (instaGRAAL needs an enzyme); expression tables |

A pool whose Hi-C library was made with a restriction enzyme can do without a
WGS library: its Hi-C reads, cut at their ligation junctions, polish its
species instead ([plan](plan.md), milestone 2). A Micro-C library has no
junction to cut at, which is why `poolC` stays with nanopore polishing.

## Contigs

Each species' contigs (milestone 1 of the [plan](plan.md)) are made from its
nanopore reads alone; [the parameters](parameters.md#contigs) set each step.

- **The reads:** filtered by length (`--ont_min_length`, 1 kb) and mean quality
  (`--ont_min_quality`: 10 for R10 reads, 7 for R9), then capped at
  `--ont_max_depth` (100×) times the genome size, the longest and best kept.
- **The candidates:** each assembler of `--assemblers` (`flye`; `flye,hifiasm`
  adds hifiasm in its `--ont` mode, for R10 reads) makes one; each is polished
  (`--ont_polisher medaka`; `dorado`, experimental, for the species whose
  reads are Dorado's BAMs, on a GPU; or `none`) and checked.
- **The one kept:** the most single-copy complete BUSCOs, then the highest QV,
  then the longest contig N50; the samplesheet's `assembler` names one instead.
- **The checks:** BUSCO with the species' lineage; Merqury's QV against the
  pool's WGS k-mers (an assembly k-mer missing from the whole pool is an
  error) or, without a WGS library, the species' nanopore k-mers, a rougher
  QV; Merqury's completeness against the species' nanopore k-mers; gfastats.
- **The pool:** skani's identity between each pair of its species. A pair
  above 95 % is too close for the pool's short reads to be split between them
  well, and the run warns about it; a pair from 85 % is of the same genus,
  which milestone 2 splits with an identity threshold as well.

BUSCO downloads its lineages when it runs; on a cluster whose nodes have no
network, download them beforehand (`busco --download saccharomycetes_odb12`)
and give their folder to `--busco_db`.

## The pool, split

The pooled libraries are split between the species of their pool (milestone 2
of the [plan](plan.md); [the parameters](parameters.md#the-pool-split)): every
species' contigs in one reference, each WGS and Hi-C library aligned to it at
once, each pair given to the species of its two mates when both are unique
(`--deconvolve_mapq`, 20) and near-identical (`--deconvolve_min_identity`,
98 %). The identity matters for a pair of species of a genus: a read of a
region missing from its own species' contigs aligns, unique, to the sister
species' copy, 85 to 95 % identical.

- **The screen:** each contig kept, moved to another species of its pool,
  dropped, or named its species' mitochondrion, from each species' nanopore
  depth on it, its identity to the other species' contigs and its Hi-C
  contacts ([`plexembler screen`](cli.md#screen)). The pool is then split
  again, on the screened contigs: a species whose DNA another species' contigs
  held (a barcode's bleed-through, assembled) gets back the pairs the first
  split found ambiguous between the two.
- **Polishing:** each species' screened contigs with its own short reads, of
  that second split (`--polish_reads auto`: its WGS pairs, or, in a pool
  without WGS, its Hi-C reads cut at their ligation junctions; `wgs` or `hic`
  to choose), Polypolish then pypolca, kept only when the QV against
  `--polish_holdout` (10 %) of the pairs, held out of the polishing, rises. A
  pool whose Hi-C library was made without an enzyme (Micro-C) and that has no
  WGS library is not polished with short reads. `--sr_polisher none` skips the
  step.

## Chromosomes

Each species' screened contigs become its genome (milestone 3 of the
[plan](plan.md); [the parameters](parameters.md#chromosomes)): the pool split
again on them, each species scaffolded with its own Hi-C pairs (YaHS, told the
restriction motifs of the library's `enzymes`), the scaffolds named, the
chromosomes first ([`plexembler finalize`](cli.md#finalize)), and its Hi-C map
made on the genome (pairtools, cooler; its finest bin `--map_resolution`,
1 kb). A Micro-C library (`enzymes` `none`) scaffolds too, without motifs. A
species whose pool has no Hi-C library is not scaffolded: its genome is its
contigs, named the same way.

The scaffolds are YaHS's, uncurated. The steps that check and correct them
(the junctions' support on the map, the telomeres, a second scaffolder) are to
come; until then, look at each species' map before trusting its chromosomes.

## Annotation

Each species' genome is annotated (milestone 4 of the [plan](plan.md);
[the parameters](parameters.md#annotation)). The pooled RNA-seq libraries are
split between the species first, on their genomes (STAR, then
[`plexembler split --type rnaseq`](cli.md#split)): each species' pairs are its
evidence. Its repeats are found and masked (`--repeats repeatmodeler`), then
its genes predicted (`--gene_predictor`):

- **`braker`** (the default): BRAKER3 `--fungus`, trained on the species
  itself, with its RNA-seq, `--orthodb_proteins` (OrthoDB's proteins of the
  fungi, `Fungi.fa.gz`, from
  [their partitions](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/)),
  or both. A species with neither is predicted by AUGUSTUS.
- **`augustus`**: AUGUSTUS ab initio, with `--augustus_species`, a species it
  was trained on (`augustus --species=help`), as close to yours as can be:
  `saccharomyces_cerevisiae_S288C` by default. Quick, and right for a species
  near one AUGUSTUS knows; BRAKER3 does better for the others.
- **`none`**: no genes.

The genes are named in the genome's order (`<species>_000010`, …) and their
proteins translated with the samplesheet's `genetic_code` (`1` when `NA`):
CUG is a sense codon in tables 1, 12 and 26 alike, so the gene structures do
not change, the proteins do. A wrong code shows as proteins with a stop codon
inside, which the run warns of.

## Expression

Each species' read pairs of each pooled RNA-seq library (milestone 5 of the
[plan](plan.md); [the parameters](parameters.md#expression-and-tracks)) are
counted on its genes, and its counts and TPM written, genes by libraries, each
library named `<condition>_<replicate>` from the libraries' samplesheet. The
libraries' strandedness is `--rnaseq_strandedness`: `reverse` (dUTP, as
Illumina's stranded kits make them) by default, `forward`, or `unstranded`.
The TPM are within each species, from its own pairs alone: a pair the split
could not give to one species is in neither's. There is no differential
expression: the tables are made for it.

Each species' tracks go beside them: each library's coverage on each strand,
and the nanopore reads' depth, as bigWig.

## Reports and the browser

Each species has its report (milestone 6): a page, its numbers as text, and its
row in MultiQC's General Statistics. The results open in
[chromosomix](https://chromosomix.app) through the store the run writes at
their top, `chromosomix.store.json`: connect the results folder in its app, or
put the results on a bucket and give the run their address
(`--browser_base https://…/results/`) to open them at chromosomix.app.

## Profiles

`-profile` takes a container engine (`docker`, `apptainer`, `singularity`,
`podman`, `charliecloud`), `conda`, an institution's profile from
[nf-core/configs](https://nf-co.re/configs/) (`maestro`, …), and `gpu` for the
steps that can use one. They combine: `-profile maestro,apptainer`. The
[guide on compute](guides/compute.md) says which to use on a workstation and on
a cluster.

`-profile test` runs on two species of a few dozen kilobases, from
`tests/data/`: the plumbing, in a few minutes.

## Releases

Pin a release, so that a run can be made again as it was:

```bash
nextflow run js2264/plexembler -r v0.1.0 -profile docker \
    --input species.csv --libraries libraries.csv --outdir results
```

`plexembler run` pins it already: plexembler 0.1.0 runs `-r v0.1.0`.

plexembler stays at 0.x ([plan](plan.md#versions-0x-for-good)): a patch
release (0.1.1) only fixes; a minor release (0.2.0) adds steps, and may change a
parameter, a column or an output, which its [changelog](changelog.md) says
under *Breaking changes*, with a page on moving a run to it.

## Running it again

`-resume` reuses every step whose inputs and code did not change. A pool's
deconvolution waits for every species of the pool: a species added to a pool
reruns the pool's deconvolution, and the steps after it, for every species of
the pool.
