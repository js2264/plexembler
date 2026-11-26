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
page; the results are on the [outputs](output.md) page.

!!! note "Under development"

    The pipeline is being built milestone by milestone
    ([implementation plan](plan.md)). Today, it reads and checks both
    samplesheets and checks the reads (NanoPlot, FastQC, MultiQC). The columns
    below are read and checked already; the steps that use them come with their
    milestones.

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
checks each column.

```csv title="species.csv"
species,pool,ont,genome_size,ploidy,n_chromosomes,genetic_code
Scer,poolA,ont/Scer.fastq.gz,12m,1,16,1
Calb,poolA,ont/Calb_run1.fastq.gz,14m,2,8,12
Calb,poolA,ont/Calb_run2.fastq.gz,14m,2,8,12
```

| Column | Required | Description | Used from |
|---|---|---|---|
| `species` | yes | The species' ID: a letter, then letters, digits and hyphens. It starts the name of every contig of the species, `<species>_<contig>`, and names its results. | now |
| `pool` | no (`pool`) | The pool the species is in: the pooled libraries of the same pool are split between its species. | now |
| `ont` | yes | The species' nanopore reads: FASTQ, gzipped or not. | now (QC); milestone 1 (assembly) |
| `genome_size` | no | The expected genome size, in bases, with an optional `k`, `m` or `g` (`12m`). | milestone 1 |
| `ploidy` | no | The expected ploidy. | milestones 1 and 6 |
| `n_chromosomes` | no | The expected number of nuclear chromosomes: a check of the scaffolds, and HapHiC's input. | milestone 3 |
| `busco_lineage` | no | A BUSCO lineage (`busco --list-datasets`); BUSCO picks one when absent. | milestone 1 |
| `genetic_code` | no (`1`) | The NCBI translation table of the nuclear genes: `12` for the CUG-Ser clade, `26` for the CUG-Ala clade. | milestone 4 |
| `telomere_repeat` | no | The telomeric repeat unit (`TGTGGG`); found from the reads when absent. | milestone 3 |

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
Huva,poolA,ont/Huva.fastq.gz,9m,,,1
Calb,poolB,ont/Calb.fastq.gz,14m,2,8,12
Dhan,poolB,ont/Dhan.fastq.gz,12m,1,7,12
Ptan,poolB,ont/Ptan.fastq.gz,13m,1,,26
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
- Empty cells are allowed where a column is optional: `Huva`'s ploidy and
  chromosome count are unknown.

What each pool's species get, once the pipeline is complete
([plan](plan.md)):

| Pool | Hi-C | WGS | RNA-seq | Its species get |
|---|---|---|---|---|
| poolA | Arima (DpnII, HinfI), two lanes | yes | 2 conditions × 2 | everything: polished and scaffolded genomes, Hi-C maps, annotation with RNA evidence, expression tables |
| poolB | DpnII | yes | none | scaffolded genomes and maps; genes from protein evidence only; no expression table |
| poolC | Micro-C (`none`) | none | 1 condition × 2 | genomes polished with nanopore reads only, and a rougher QV, from nanopore k-mers; YaHS only (instaGRAAL needs an enzyme); expression tables |

## Profiles

`-profile` takes a container engine (`docker`, `apptainer`, `singularity`,
`podman`, `charliecloud`), `conda`, an institution's profile from
[nf-core/configs](https://nf-co.re/configs/) (`maestro`, …), and `gpu` for the
steps that can use one. They combine: `-profile maestro,apptainer`. The
[guide on compute](guides/compute.md) says which to use on a workstation and on
a cluster.

`-profile test` runs on two species of a few dozen kilobases, from
`tests/data/`: the plumbing, in a few minutes.

## Running it again

`-resume` reuses every step whose inputs and code did not change. A pool's
deconvolution waits for every species of the pool: a species added to a pool
reruns the pool's deconvolution, and the steps after it, for every species of
the pool.
