process PLEXEMBLER_JUNCTIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A species' pairs on its genome, and the genome's AGP.
    tuple val(meta), path(pairs), path(agp)

    output:
    tuple val(meta), path("${meta.id}.junctions.tsv")     , emit: table
    tuple val(meta), path("${meta.id}.junctions_mqc.tsv") , emit: multiqc
    tuple val(meta), path("*.log")                        , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    plexembler junctions \\
        ${pairs} \\
        -a ${agp} \\
        -p ${meta.id} \\
        -o . \\
        ${args}
    """

    stub:
    """
    printf 'scaffold\\tposition\\tleft_contig\\tright_contig\\tobserved\\texpected\\tratio\\tverdict\\n' > ${meta.id}.junctions.tsv
    printf '# id: plexembler_junctions\\nSample\\tJoins\\n${meta.id}\\t0\\n' > ${meta.id}.junctions_mqc.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
