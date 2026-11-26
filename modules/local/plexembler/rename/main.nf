process PLEXEMBLER_RENAME {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.fa")       , emit: fasta
    tuple val(meta), path("${prefix}.names.tsv"), emit: names
    tuple val(meta), path("*.log")              , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${fasta}" == "${prefix}.fa") {
        error("PLEXEMBLER_RENAME: the input and the output are both named ${fasta}: set ext.prefix")
    }
    """
    plexembler rename \\
        ${fasta} \\
        -p ${meta.id} \\
        -o ${prefix}.fa \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf ">${meta.id}_contig_1\\nACGT\\n" > ${prefix}.fa
    printf 'name\\toriginal_name\\tlength_bp\\n${meta.id}_contig_1\\tcontig_1\\t4\\n' > ${prefix}.names.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
