process PLEXEMBLER_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A species (meta.id, meta.pool), and the files the pipeline wrote for it
    // and its pool, each known by its name.
    tuple val(meta), path(files, stageAs: 'inputs/*')

    output:
    tuple val(meta), path("${meta.id}.report.txt")      , emit: txt
    tuple val(meta), path("${meta.id}.report.html")     , emit: html
    tuple val(meta), path("${meta.id}.report_mqc.yaml") , emit: multiqc
    tuple val(meta), path("*.log")                      , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pool = meta.pool ? "--pool ${meta.pool}" : ''
    """
    plexembler report \\
        inputs/* \\
        -p ${meta.id} \\
        ${pool} \\
        -o . \\
        ${args}
    """

    stub:
    """
    printf '[genome]\\n' > ${meta.id}.report.txt
    printf '<html></html>\\n' > ${meta.id}.report.html
    printf 'id: plexembler_report\\n' > ${meta.id}.report_mqc.yaml
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
