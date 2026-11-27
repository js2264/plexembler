process PLEXEMBLER_EXPRESSION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A species' featureCounts table; meta.labels, each library's LIBRARY=LABEL.
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${meta.id}.counts.tsv")        , emit: counts
    tuple val(meta), path("${meta.id}.tpm.tsv")           , emit: tpm
    tuple val(meta), path("${meta.id}.expression_mqc.tsv"), emit: multiqc
    tuple val(meta), path("*.log")                        , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def labels = (meta.labels ?: []).collect { label -> "-l ${label}" }.join(' ')
    """
    plexembler expression \\
        ${counts} \\
        -p ${meta.id} \\
        -o . \\
        ${labels} \\
        ${args}
    """

    stub:
    """
    printf 'gene\\tlength_bp\\n' > ${meta.id}.counts.tsv
    printf 'gene\\n' > ${meta.id}.tpm.tsv
    printf '# id: plexembler_expression\\nSample\\tPairs on genes\\n' > ${meta.id}.expression_mqc.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
