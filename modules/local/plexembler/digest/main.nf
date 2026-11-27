process PLEXEMBLER_DIGEST {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A species' Hi-C reads, and their library's enzymes (meta.enzymes).
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.digested.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.digest.json")      , emit: json
    tuple val(meta), path("${meta.id}.digest_mqc.tsv")   , emit: multiqc
    tuple val(meta), path("*.log")                       , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    plexembler digest \\
        ${reads} \\
        -e ${meta.enzymes} \\
        -p ${meta.id} \\
        -o . \\
        ${args}
    """

    stub:
    """
    echo "" | gzip -n > ${meta.id}.digested.fastq.gz
    echo '{"reads": 0}' > ${meta.id}.digest.json
    printf '# id: plexembler_digest\\nSample\\tReads\\n' > ${meta.id}.digest_mqc.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
