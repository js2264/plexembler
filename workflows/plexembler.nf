/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ONT_READS              } from '../subworkflows/local/ont_reads'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { KMERS                  } from '../subworkflows/local/kmers'
include { CONTIGS                } from '../subworkflows/local/contigs'
include { POOL_DECONVOLVE        } from '../subworkflows/local/pool_deconvolve'
include { POLISH_SR              } from '../subworkflows/local/polish_sr'
include { SCREEN                 } from '../subworkflows/local/screen'
include { POOL_DECONVOLVE as POOL_DECONVOLVE_SCREENED } from '../subworkflows/local/pool_deconvolve'
include { SCAFFOLD               } from '../subworkflows/local/scaffold'
include { GENOME_QC as GENOME_QC_FINAL } from '../subworkflows/local/genome_qc'
include { MAPS                   } from '../subworkflows/local/maps'
include { PLEXEMBLER_JUNCTIONS   } from '../modules/local/plexembler/junctions/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_plexembler_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PLEXEMBLER {

    take:
    ch_species   // channel: [ val(meta), [ path(ont) ] ], one item per species
    ch_libraries // channel: [ val(meta), [ path(fastq) ] ], one item per pooled library
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: each species' nanopore reads, as one FASTQ, checked, filtered
    // and capped; its chemistry and genome size when the samplesheet left them NA
    //
    ONT_READS(ch_species, params.ont_max_depth)
    ch_multiqc_files = ch_multiqc_files.mix(ONT_READS.out.multiqc_files)

    //
    // MODULE: the pooled libraries' reads as they come, a library's runs together
    //
    FASTQC(ch_libraries)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, zip -> zip }.flatten())

    //
    // SUBWORKFLOW: the k-mers of each species' nanopore reads, and of each pool's WGS reads
    //
    def ch_wgs_reads = ch_libraries
        .filter { meta, _reads -> meta.type == 'wgs' }
        .map { meta, reads -> [meta.pool, reads] }
        .groupTuple()
        .map { pool, reads -> [[id: pool], reads.flatten()] }
    KMERS(ONT_READS.out.reads, ch_wgs_reads, params.kmer_size)

    //
    // SUBWORKFLOW: each species' contigs, every assembler's a candidate, polished,
    // checked, one kept; and the identity between the species of each pool
    //
    def ch_bams = ch_species
        .map { meta, ont -> [meta, [ont].flatten().findAll { run -> run.name.endsWith('.bam') }] }
        .filter { _meta, bams -> bams }
    CONTIGS(
        ONT_READS.out.reads,
        ch_bams,
        KMERS.out.ont_db,
        KMERS.out.wgs_db,
        params.assemblers.tokenize(','),
        params.ont_polisher,
        params.busco_lineage,
        params.busco_db ? file(params.busco_db, checkIfExists: true) : [],
    )
    ch_multiqc_files = ch_multiqc_files.mix(CONTIGS.out.multiqc_files)

    //
    // SUBWORKFLOW: each pooled library split between its pool's species, on their contigs
    //
    POOL_DECONVOLVE(CONTIGS.out.contigs, ch_libraries, 'contigs')
    ch_multiqc_files = ch_multiqc_files.mix(POOL_DECONVOLVE.out.multiqc_files)

    //
    // SUBWORKFLOW: each species polished with its own short reads, when that raises its QV
    //
    def ch_polished = CONTIGS.out.contigs
    if (params.sr_polisher != 'none') {
        POLISH_SR(CONTIGS.out.contigs, POOL_DECONVOLVE.out.reads, params.polish_reads, params.kmer_size)
        ch_polished = POLISH_SR.out.contigs
        ch_multiqc_files = ch_multiqc_files.mix(POLISH_SR.out.multiqc_files)
    }

    //
    // SUBWORKFLOW: each contig kept, moved to another species of its pool, or
    // dropped; each species' organelle named
    //
    SCREEN(ch_polished, ONT_READS.out.reads, POOL_DECONVOLVE.out.contacts)
    ch_multiqc_files = ch_multiqc_files.mix(SCREEN.out.multiqc_files)

    //
    // SUBWORKFLOW: the pool split again, on the screened contigs: the pairs
    // each species is scaffolded with. Its tables are published, not given to
    // MultiQC, whose split is the first's (their files share their names).
    //
    POOL_DECONVOLVE_SCREENED(SCREEN.out.contigs, ch_libraries, 'screened')

    //
    // SUBWORKFLOW: each species scaffolded with its own Hi-C pairs, and finished
    //
    SCAFFOLD(SCREEN.out.contigs, POOL_DECONVOLVE_SCREENED.out.reads)

    //
    // SUBWORKFLOW: the final genome's checks
    //
    GENOME_QC_FINAL(
        SCAFFOLD.out.genome.map { meta, fasta -> [meta + [candidate: 'genome'], fasta] },
        KMERS.out.ont_db,
        KMERS.out.wgs_db,
        params.busco_lineage,
        params.busco_db ? file(params.busco_db, checkIfExists: true) : [],
    )
    ch_multiqc_files = ch_multiqc_files.mix(GENOME_QC_FINAL.out.multiqc_files)

    //
    // SUBWORKFLOW: each species' Hi-C map on its genome
    //
    MAPS(SCAFFOLD.out.genome, SCAFFOLD.out.sizes, SCAFFOLD.out.hic, params.map_resolution)
    ch_multiqc_files = ch_multiqc_files.mix(MAPS.out.multiqc_files)

    //
    // MODULE: each join of each species' scaffolds, scored by its Hi-C pairs
    //
    PLEXEMBLER_JUNCTIONS(
        MAPS.out.pairs
            .map { meta, pairs -> [meta.id, meta, pairs] }
            .join(SCAFFOLD.out.agp.map { meta, agp -> [meta.id, agp] })
            .map { _id, meta, pairs, agp -> [meta, pairs, agp] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(PLEXEMBLER_JUNCTIONS.out.multiqc.map { _meta, table -> table })

    //
    // The steps to come, milestone by milestone (docs/plan.md): annotation,
    // expression and the report of each species.
    //

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'plexembler_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'plexembler'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
