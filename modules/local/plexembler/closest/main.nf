process PLEXEMBLER_CLOSEST {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A species (meta.id), and sourmash gather's matches of its genome.
    tuple val(meta), path(gather)
    // The database's taxonomy (sourmash's lineages, CSV), or [].
    path lineages

    output:
    tuple val(meta), path("${meta.id}.closest.tsv")     , emit: table
    // The closest genome's sequence, from NCBI: none without network.
    tuple val(meta), path("${meta.id}.closest.fa.gz")   , emit: fasta, optional: true
    tuple val(meta), path("${meta.id}.closest_mqc.tsv") , emit: multiqc
    tuple val(meta), path("*.log")                      , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def taxonomy = lineages ? "--lineages ${lineages}" : ''
    """
    plexembler closest \\
        ${gather} \\
        -p ${meta.id} \\
        -o . \\
        ${taxonomy} \\
        ${args}
    """

    stub:
    """
    printf 'rank\\taccession\\tname\\tani\\tspecies_containment\\tmatch_containment\\tshared_bp\\tlineage\\tsequence\\n' > ${meta.id}.closest.tsv
    printf '1\\tGCA_000000000.1\\tGCA_000000000.1 stub\\t99.00\\t90.00\\t90.00\\t1000\\tNA\\tdownloaded\\n' >> ${meta.id}.closest.tsv
    printf '>stub\\nACGT\\n' | gzip > ${meta.id}.closest.fa.gz
    printf '# id: plexembler_closest\\nSample\\tClosest genome\\n' > ${meta.id}.closest_mqc.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
