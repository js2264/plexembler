//
// Subworkflow with functionality specific to the js2264/plexembler pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to the species samplesheet (--input)
    libraries         //  string: Path to the pooled libraries samplesheet (--libraries)
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input species.csv --libraries libraries.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // The species and the pooled libraries, read and checked against each other
    //
    def species_rows = readSpecies(input)
    def library_rows = readLibraries(libraries)
    validatePools(species_rows, library_rows)

    emit:
    species   = channel.fromList(species_rows) // channel: [ val(meta), [ path(ont) ] ], one item per species
    libraries = channel.fromList(library_rows) // channel: [ val(meta), [ path(fastq) ] ], one item per library
    versions  = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:

    //
    // Completion email and summary
    //
    workflow.onComplete {

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// A path of a samplesheet: a URL or an absolute path as it is, a relative path
// read from the samplesheet's folder (docs/usage.md). It must exist.
//
def resolveSheetPath(value, sheet) {
    def path = (value ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/ || value.startsWith('/'))
        ? file(value)
        : sheet.parent.resolve(value)
    if (!path.exists()) {
        error("${sheet.name}: ${value} does not exist (a relative path is read from the samplesheet's folder, ${sheet.parent})")
    }
    return path
}

//
// The keys on which the rows of one species, or of one library, disagree
//
def disagreements(metas) {
    return metas[0].keySet().findAll { key -> metas.collect { meta -> meta[key] }.unique().size() > 1 }.join(', ')
}

//
// A species' meta map, every key present: what the samplesheet left empty is null
//
def speciesMeta(meta) {
    return [
        id              : meta.id,
        pool            : meta.pool ?: 'pool',
        genome_size     : meta.genome_size ?: null,
        ploidy          : meta.ploidy ?: null,
        n_chromosomes   : meta.n_chromosomes ?: null,
        busco_lineage   : meta.busco_lineage ?: null,
        genetic_code    : meta.genetic_code ?: 1,
        telomere_repeat : meta.telomere_repeat ? meta.telomere_repeat.toUpperCase() : null,
    ]
}

//
// A pooled library's meta map, every key present; single_end for nf-core's modules
//
def libraryMeta(meta) {
    def type = meta.type
    return [
        id         : meta.id,
        pool       : meta.pool ?: 'pool',
        type       : type,
        condition  : type == 'rnaseq' ? (meta.condition ?: meta.id) : null,
        replicate  : type == 'rnaseq' ? (meta.replicate ?: 1) : null,
        enzymes    : type == 'hic' ? meta.enzymes : null,
        single_end : false,
    ]
}

//
// The species samplesheet (--input): its rows grouped by species, a species'
// nanopore runs together. [ meta, [ ont, ... ] ]
//
def readSpecies(input) {
    def sheet = file(input, checkIfExists: true)
    def rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
        .collect { meta, ont -> [ speciesMeta(meta), resolveSheetPath(ont, sheet) ] }
    return rows
        .groupBy { meta, _ont -> meta.id }
        .collect { id, runs ->
            def metas = runs.collect { meta, _ont -> meta }.unique()
            if (metas.size() > 1) {
                error("${sheet.name}: the rows of species ${id} disagree on ${disagreements(metas)}: a species' rows are its runs, and share its values")
            }
            [ metas[0], runs.collect { _meta, ont -> ont } ]
        }
}

//
// The pooled libraries samplesheet (--libraries): its rows grouped by library, a
// library's runs together. [ meta, [ fastq_1, fastq_2, ... ] ]
//
def readLibraries(libraries) {
    def sheet = file(libraries, checkIfExists: true)
    def rows = samplesheetToList(libraries, "${projectDir}/assets/schema_libraries.json")
        .collect { meta, fastq_1, fastq_2 ->
            [ libraryMeta(meta), [ resolveSheetPath(fastq_1, sheet), resolveSheetPath(fastq_2, sheet) ] ]
        }
    return rows
        .groupBy { meta, _reads -> meta.id }
        .collect { id, runs ->
            def metas = runs.collect { meta, _reads -> meta }.unique()
            if (metas.size() > 1) {
                error("${sheet.name}: the rows of library ${id} disagree on ${disagreements(metas)}: a library's rows are its runs, and share its values")
            }
            [ metas[0], runs.collect { _meta, reads -> reads }.flatten() ]
        }
}

//
// Every library belongs to a pool of species; a pool without a Hi-C or a WGS
// library gets the steps that need one skipped, with a warning
//
def validatePools(species, libraries) {
    def pools = species.collect { meta, _ont -> meta.pool }.unique()
    def orphans = libraries.findAll { meta, _reads -> !(meta.pool in pools) }
    if (orphans) {
        def named = orphans.collect { meta, _reads -> "${meta.id} (pool ${meta.pool})" }.join(', ')
        error("--libraries: ${named} belong to no pool of --input, whose pools are ${pools.join(', ')}")
    }
    pools.each { pool ->
        def types = libraries.findAll { meta, _reads -> meta.pool == pool }.collect { meta, _reads -> meta.type }
        ['hic', 'wgs'].findAll { type -> !(type in types) }.each { type ->
            log.warn("Pool ${pool} has no ${type} library: the steps that need one will be skipped for its species")
        }
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // The tools of the run, as CITATIONS.md lists them: one line per tool, added with its step.
    def citation_text = [
            "Tools used in the workflow included:",
            "NanoPlot (De Coster et al. 2023),",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // The references of toolCitationText(), in its order.
    def reference_text = [
            "<li>De Coster, W., & Rademakers, R. (2023). NanoPack2: population-scale evaluation of long-read sequencing data. Bioinformatics, 39(5), btad311. doi: 10.1093/bioinformatics/btad311</li>",
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
