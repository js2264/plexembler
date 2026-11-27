process PLEXEMBLER_BROWSER {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // Each species and its pool, SPECIES=POOL; each of its files, by its path
    // in the results, SPECIES=PATH.
    val species
    val files

    output:
    path "chromosomix.store.json", emit: store
    path "*.log"                 , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def species_args = species.collect { pair -> "-s ${pair}" }.join(' ')
    def file_args = files.collect { pair -> "-f ${pair}" }.join(' ')
    """
    plexembler browser \\
        ${species_args} \\
        ${file_args} \\
        -o . \\
        ${args}
    """

    stub:
    """
    echo '{"kind": "chromosomix.store", "version": 1, "entries": []}' > chromosomix.store.json
    touch browser.plexembler_20260101000000.log
    """
}
