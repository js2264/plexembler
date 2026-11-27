# js2264/plexembler pipeline parameters

Reference genomes, gene annotations, Hi-C maps and expression tables for every species of a pool, from per-species nanopore reads and pooled Hi-C, WGS and RNA-seq

## Input/output options

The species and the pooled libraries, and where the results go.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | The species samplesheet: one row per nanopore run of a species. <details><summary>Help</summary><small>A CSV file with a header: species, pool, ont (FASTQ, or Dorado's unaligned BAM), and optionally chemistry, genome_size, ploidy, n_chromosomes, busco_lineage, genetic_code, telomere_repeat and assembler. A relative path in it is read from the samplesheet's folder. See the [usage page](https://js2264.github.io/plexembler/usage/#species).</small></details>| `string` |  | True |  |
| `libraries` | The pooled libraries samplesheet: one row per pair of FASTQ files of a Hi-C, WGS or RNA-seq library made from a pool of species. <details><summary>Help</summary><small>A CSV file with a header: library, pool, type (hic, wgs or rnaseq), fastq_1, fastq_2, and optionally condition and replicate (rnaseq), enzymes (hic). A relative path in it is read from the samplesheet's folder. See the [usage page](https://js2264.github.io/plexembler/usage/#pooled-libraries).</small></details>| `string` |  | True |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

## Contigs

Each species' nanopore reads, its assemblies as candidates, their polishing and checks, and the one kept (milestone 1).

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `ont_min_length` | The shortest nanopore read kept, in bases. | `integer` | 1000 |  |  |
| `ont_min_quality` | The lowest mean quality of a nanopore read kept; by default 10 for R10 reads, 7 for R9. | `number` |  |  |  |
| `ont_max_depth` | The depth the nanopore reads are capped at, from the genome size; the longest and best reads are kept. 0 for no cap. <details><summary>Help</summary><small>Filtlong keeps `ont_max_depth` times the genome size, in bases: the samplesheet's `genome_size`, or LRGE's estimate when it is NA.</small></details>| `integer` | 100 |  |  |
| `save_ont_reads` | Save each species' nanopore reads as they were assembled: filtered and capped. | `boolean` |  |  |  |
| `assemblers` | The assemblers that run, comma-separated: each one's contigs are a candidate (flye, hifiasm). <details><summary>Help</summary><small>Every candidate is renamed, polished and checked; `plexembler select` keeps one, by BUSCO's single-copy genes, then the QV, then the contig N50, unless the samplesheet's `assembler` names one. hifiasm runs in its `--ont` mode, for R10 reads.</small></details>| `string` | flye |  |  |
| `ont_polisher` | How the candidates are polished with their species' nanopore reads: medaka, dorado (experimental: a GPU, and Dorado's BAMs), or none. (accepted: `medaka`\|`dorado`\|`none`) <details><summary>Help</summary><small>`dorado` runs `dorado polish` for the species whose `ont` is Dorado's unaligned BAM, and Medaka for the others. It is experimental until it has been tested on a GPU (docs/plan.md).</small></details>| `string` | medaka |  |  |
| `medaka_model` | Medaka's model; by default r1041_e82_400bps_sup_v5.2.0 for R10 reads, r941_min_sup_g507 for R9. | `string` |  |  |  |
| `busco_lineage` | The BUSCO lineage of a species whose samplesheet's busco_lineage is NA: auto_euk to let BUSCO choose. | `string` | auto_euk |  |  |
| `busco_db` | A folder of BUSCO's lineages, downloaded beforehand (`busco --download`), for runs without network. | `string` |  |  |  |
| `kmer_size` | The k-mers' length, for Merqury's QV and completeness and GenomeScope2. | `integer` | 21 |  |  |

## The pool, split

Each pooled library split between its pool's species, and each species polished with its own short reads (milestone 2).

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `deconvolve_mapq` | A read's lowest mapping quality to be given to a species: below, it is ambiguous. | `integer` | 20 |  |  |
| `deconvolve_min_identity` | A read's lowest identity to its species' contig: below, it is a read of another species. <details><summary>Help</summary><small>A read of a region missing from its own species' contigs aligns, unique, to a sister species' copy, 85 to 95 % identical: the identity, not the mapping quality, catches it (docs/plan.md, Decisions).</small></details>| `number` | 0.98 |  |  |
| `polish_reads` | The short reads each species is polished with: auto (its WGS pairs, or its Hi-C reads cut at their junctions when its pool has no WGS library), wgs or hic. (accepted: `auto`\|`wgs`\|`hic`) | `string` | auto |  |  |
| `sr_polisher` | Short-read polishing: Polypolish, then pypolca, kept only if it raises the QV; or none. (accepted: `polypolish,pypolca`\|`none`) | `string` | polypolish,pypolca |  |  |
| `polish_holdout` | The share of each library's pairs kept apart from the polishing, whose k-mers the QV guard counts against. | `number` | 0.1 |  |  |
| `save_split_reads` | Save each species' split reads, of each pooled library. | `boolean` |  |  |  |

## Chromosomes

Each species scaffolded with its own Hi-C pairs, its scaffolds named, and its Hi-C map (milestone 3).

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `min_chromosome` | A species without n_chromosomes: the shortest scaffold named a chromosome, in bases. <details><summary>Help</summary><small>When the samplesheet gives a species' n_chromosomes, its longest scaffolds, that many, are its chromosomes; the others are unplaced (docs/output.md, The genome).</small></details>| `integer` | 100000 |  |  |
| `map_resolution` | The finest bin of each species' Hi-C map, in bases; the map is zoomed out from it. | `integer` | 1000 |  |  |

## Annotation

Each species' repeats, genes and non-coding RNAs, and each pooled RNA-seq library split between the species (milestone 4).

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `repeats` | Each species' repeats: found de novo (RepeatModeler) and masked soft (RepeatMasker), or none. (accepted: `repeatmodeler`\|`none`) | `string` | repeatmodeler |  |  |
| `gene_predictor` | Each species' genes: BRAKER3 (--fungus) with its RNA-seq alignments and --orthodb_proteins, AUGUSTUS ab initio for a species with neither; augustus for AUGUSTUS alone; or none. (accepted: `braker`\|`augustus`\|`none`) <details><summary>Help</summary><small>BRAKER3 trains GeneMark and AUGUSTUS on the species itself, from its RNA-seq, related proteins, or both. A species with neither is predicted ab initio by AUGUSTUS with --augustus_species, a species it was trained on.</small></details>| `string` | braker |  |  |
| `augustus_species` | The species AUGUSTUS predicts ab initio with: one it was trained on (augustus --species=help), as close as can be. | `string` | saccharomyces_cerevisiae_S288C |  |  |
| `orthodb_proteins` | Proteins of related species for BRAKER3, such as OrthoDB's partition of the fungi (Fungi.fa.gz). | `string` |  |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. (accepted: `symlink`\|`rellink`\|`link`\|`copy`\|`copyNoFollow`\|`move`) <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `trace_report_suffix` | Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss. | `string` |  |  | True |
| `help` | Display the help message. | `['boolean', 'string']` |  |  |  |
| `help_full` | Display the full detailed help message. | `boolean` |  |  |  |
| `show_hidden` | Display hidden parameters in the help message (only works when --help or --help_full are provided). | `boolean` |  |  |  |
