process LRGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lrge:1.0.0--hf1b6044_0' :
        'biocontainers/lrge:1.0.0--hf1b6044_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    // The estimate, in bases, alone on its line; and LRGE's log, with its interval.
    tuple val(meta), path("${prefix}.lrge.txt"), emit: estimate
    tuple val(meta), path("${prefix}.lrge.log"), emit: log
    tuple val("${task.process}"), val('lrge'), eval("lrge --version | sed 's/^lrge //'"), emit: versions_lrge, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    lrge \\
        -t ${task.cpus} \\
        -o ${prefix}.lrge.txt \\
        ${args} \\
        ${reads} \\
        2> >(tee ${prefix}.lrge.log >&2)
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo 12000000 > ${prefix}.lrge.txt
    touch ${prefix}.lrge.log
    """
}
