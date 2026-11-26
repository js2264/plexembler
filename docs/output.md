# Outputs

Everything goes under `--outdir`, in three places: one folder per **species**,
one per **pool**, and the run's own files. The tree below is the one the
pipeline is built towards ([implementation plan](plan.md)); a path marked
*(now)* is written today.

```text
<outdir>/
├── species/<species>/
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
│   ├── candidates/    every assembler × scaffolder candidate, its checks and its map
│   ├── qc/            nanoplot/ (now), BUSCO, Merqury, OMArk, contamination, …
│   └── <species>.report.html, .report.txt, .report_mqc.yaml
├── pools/<pool>/
│   ├── qc/            fastqc/ (now): the pooled libraries
│   ├── deconvolution/ each species' share of each library; ambiguous reads; the
│   │                  nanopore cross-talk between species
│   ├── hic/           the pool's map: every species of the pool in one map
│   └── comparative/   identity between the species, orthogroups, species tree, synteny
├── multiqc/           multiqc_report.html (now): every step, every species and pool
└── pipeline_info/     (now) software versions, Nextflow's reports, parameters
```

## Now

### Nanopore reads

`species/<species>/qc/nanoplot/`: [NanoPlot](https://github.com/wdecoster/NanoPlot)'s
report and statistics of a species' nanopore reads, its runs together
(`<species>_NanoPlot-report.html`, `<species>_NanoStats.txt`).

### Pooled libraries

`pools/<pool>/qc/fastqc/`: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)'s
report of each read of each pooled library (`<library>_1_fastqc.html`,
`<library>_2_fastqc.html`), the library's runs together.

### MultiQC

`multiqc/multiqc_report.html`: the runs' QC and the software versions, in one
report, with the methods' description and its references.

### Pipeline information

`pipeline_info/`: Nextflow's execution report, timeline, trace and DAG; the
parameters of the run; the software versions
(`plexembler_software_mqc_versions.yml`).
