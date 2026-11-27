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
    assemblers        //  string: The assemblers that run, comma-separated (--assemblers)
    gene_predictor    //  string: What predicts the genes (--gene_predictor)
    ont_polisher      //  string: What polishes the contigs with nanopore reads (--ont_polisher)
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
    validateAssemblers(species_rows, assemblers.tokenize(','))
    // The steps' environments are conda's (micromamba's), not a container's.
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba', 'micromamba']).size() >= 1) {
        validateCondaSteps(gene_predictor, ont_polisher)
    }

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
// A value a samplesheet may leave unknown: NA, or an empty cell, is null
//
def known(value) {
    return (value == null || value == [] || value == '' || value == 'NA') ? null : value
}

//
// A species' meta map, every key present: what the samplesheet left unknown (NA,
// or empty) is null, and the steps that need it estimate it (docs/usage.md). The
// schema checks the values: a genetic code is an NCBI translation table.
//
def speciesMeta(meta) {
    def genetic_code = known(meta.genetic_code)
    def genome_size = known(meta.genome_size)
    def telomere_repeat = known(meta.telomere_repeat)
    return [
        id              : meta.id,
        pool            : known(meta.pool) ?: 'pool',
        chemistry       : known(meta.chemistry),
        genome_size     : genome_size?.toString(),
        ploidy          : known(meta.ploidy)?.toInteger(),
        n_chromosomes   : known(meta.n_chromosomes)?.toInteger(),
        busco_lineage   : known(meta.busco_lineage),
        genetic_code    : genetic_code?.toInteger(),
        telomere_repeat : telomere_repeat ? telomere_repeat.toUpperCase() : null,
        assembler       : known(meta.assembler),
        // Nanopore reads are single reads, as nf-core's modules need to know.
        single_end      : true,
    ]
}

//
// A species may name the assembler whose contigs it keeps: one that runs
//
def validateAssemblers(species, assemblers) {
    species.each { meta, _ont ->
        if (meta.assembler && !(meta.assembler in assemblers)) {
            error("--input: species ${meta.id} keeps the contigs of ${meta.assembler}, which --assemblers (${assemblers.join(',')}) does not run")
        }
    }
}

//
// A pooled library's meta map, every key present; single_end for nf-core's modules
//
//
// The steps that have no conda environment, refused with conda (-profile conda,
// mamba, micromamba): BRAKER3, whose authors support it in their container
// alone, and Dorado, which has no conda package.
//
def validateCondaSteps(gene_predictor, ont_polisher) {
    if (gene_predictor == 'braker') {
        error("--gene_predictor braker: BRAKER3 runs in a container alone, and conda's profiles run none: choose --gene_predictor augustus, or a container engine (docs/usage.md)")
    }
    if (ont_polisher == 'dorado') {
        error("--ont_polisher dorado: Dorado has no conda package: choose --ont_polisher medaka, or a container engine (docs/usage.md)")
    }
}

def libraryMeta(meta) {
    def type = meta.type
    return [
        id         : meta.id.toString(),
        pool       : meta.pool ?: 'pool',
        type       : type,
        condition  : type == 'rnaseq' ? (meta.condition ?: meta.id).toString() : null,
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
            "NanoPlot and chopper (De Coster et al. 2023),",
            "SAMtools (Danecek et al. 2021),",
            "Filtlong (Wick),",
            "LRGE (Hall and Coin 2024),",
            "Flye (Kolmogorov et al. 2019),",
            "hifiasm (Cheng et al. 2021),",
            "gfatools (Li),",
            "Medaka and Dorado (Oxford Nanopore Technologies),",
            "gfastats (Formenti et al. 2022),",
            "BUSCO (Manni et al. 2021),",
            "Merqury and meryl (Rhie et al. 2020),",
            "GenomeScope 2.0 (Ranallo-Benavidez et al. 2020),",
            "skani (Shaw and Yu 2023),",
            "bwa-mem2 (Vasimuddin et al. 2019),",
            "BWA (Li and Durbin 2009),",
            "Polypolish (Wick and Holt 2022),",
            "pypolca (Bouras et al. 2024),",
            "minimap2 (Li 2018),",
            "YaHS (Zhou et al. 2023),",
            "pairtools (Open2C et al. 2024),",
            "cooler (Abdennur and Mirny 2020),",
            "STAR (Dobin et al. 2013),",
            "RepeatModeler2 (Flynn et al. 2020),",
            "RepeatMasker (Smit et al. 2015),",
            "BRAKER3 (Gabriel et al. 2024),",
            "AUGUSTUS (Stanke et al. 2008),",
            "tRNAscan-SE (Chan et al. 2021),",
            "barrnap (Seemann),",
            "featureCounts (Liao et al. 2014),",
            "deepTools (Ramírez et al. 2016),",
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
            "<li>Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: 10.1093/gigascience/giab008</li>",
            "<li>Wick, R. R. Filtlong: quality filtering tool for long reads. URL: https://github.com/rrwick/Filtlong</li>",
            "<li>Hall, M. B., & Coin, L. J. M. (2024). Genome size estimation from long read overlaps. bioRxiv. doi: 10.1101/2024.11.27.625777</li>",
            "<li>Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37(5), 540-546. doi: 10.1038/s41587-019-0072-8</li>",
            "<li>Cheng, H., Concepcion, G. T., Feng, X., Zhang, H., & Li, H. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nature Methods, 18(2), 170-175. doi: 10.1038/s41592-020-01056-5</li>",
            "<li>Li, H. gfatools. URL: https://github.com/lh3/gfatools</li>",
            "<li>Oxford Nanopore Technologies. Medaka and Dorado. URL: https://github.com/nanoporetech</li>",
            "<li>Formenti, G., et al. (2022). Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs. Bioinformatics, 38(17), 4214-4216. doi: 10.1093/bioinformatics/btac460</li>",
            "<li>Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A., & Zdobnov, E. M. (2021). BUSCO Update. Molecular Biology and Evolution, 38(10), 4647-4654. doi: 10.1093/molbev/msab199</li>",
            "<li>Rhie, A., Walenz, B. P., Koren, S., & Phillippy, A. M. (2020). Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. Genome Biology, 21(1), 245. doi: 10.1186/s13059-020-02134-9</li>",
            "<li>Ranallo-Benavidez, T. R., Jaron, K. S., & Schatz, M. C. (2020). GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. Nature Communications, 11(1), 1432. doi: 10.1038/s41467-020-14998-3</li>",
            "<li>Shaw, J., & Yu, Y. W. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nature Methods, 20(11), 1661-1665. doi: 10.1038/s41592-023-02018-3</li>",
            "<li>Vasimuddin, M., Misra, S., Li, H., & Aluru, S. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE IPDPS, 314-324. doi: 10.1109/IPDPS.2019.00041</li>",
            "<li>Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25(14), 1754-1760. doi: 10.1093/bioinformatics/btp324</li>",
            "<li>Wick, R. R., & Holt, K. E. (2022). Polypolish: Short-read polishing of long-read bacterial genome assemblies. PLoS Computational Biology, 18(1), e1009802. doi: 10.1371/journal.pcbi.1009802</li>",
            "<li>Bouras, G., et al. (2024). How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. Microbial Genomics, 10(6), 001254. doi: 10.1099/mgen.0.001254</li>",
            "<li>Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100. doi: 10.1093/bioinformatics/bty191</li>",
            "<li>Zhou, C., McCarthy, S. A., & Durbin, R. (2023). YaHS: yet another Hi-C scaffolding tool. Bioinformatics, 39(1), btac808. doi: 10.1093/bioinformatics/btac808</li>",
            "<li>Open2C, Abdennur, N., Fudenberg, G., Flyamer, I. M., Galitsyna, A. A., Goloborodko, A., Imakaev, M., & Venev, S. V. (2024). Pairtools: from sequencing data to chromosome contacts. PLoS Computational Biology, 20(5), e1012164. doi: 10.1371/journal.pcbi.1012164</li>",
            "<li>Abdennur, N., & Mirny, L. A. (2020). Cooler: scalable storage for Hi-C data and other genomically labeled arrays. Bioinformatics, 36(1), 311-316. doi: 10.1093/bioinformatics/btz540</li>",
            "<li>Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15-21. doi: 10.1093/bioinformatics/bts635</li>",
            "<li>Flynn, J. M., et al. (2020). RepeatModeler2 for automated genomic discovery of transposable element families. PNAS, 117(17), 9451-9457. doi: 10.1073/pnas.1921046117</li>",
            "<li>Smit, A. F. A., Hubley, R., & Green, P. (2013-2015). RepeatMasker Open-4.0. URL: https://www.repeatmasker.org</li>",
            "<li>Gabriel, L., et al. (2024). BRAKER3: Fully automated genome annotation using RNA-seq and protein evidence with GeneMark-ETP, AUGUSTUS, and TSEBRA. Genome Research, 34(5), 769-777. doi: 10.1101/gr.278090.123</li>",
            "<li>Stanke, M., Diekhans, M., Baertsch, R., & Haussler, D. (2008). Using native and syntenically mapped cDNA alignments to improve de novo gene finding. Bioinformatics, 24(5), 637-644. doi: 10.1093/bioinformatics/btn013</li>",
            "<li>Chan, P. P., Lin, B. Y., Mak, A. J., & Lowe, T. M. (2021). tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. Nucleic Acids Research, 49(16), 9077-9096. doi: 10.1093/nar/gkab688</li>",
            "<li>Seemann, T. barrnap: BAsic Rapid Ribosomal RNA Predictor. URL: https://github.com/tseemann/barrnap</li>",
            "<li>Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930. doi: 10.1093/bioinformatics/btt656</li>",
            "<li>Ramírez, F., et al. (2016). deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Research, 44(W1), W160-W165. doi: 10.1093/nar/gkw257</li>",
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
