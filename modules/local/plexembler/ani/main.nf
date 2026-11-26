process PLEXEMBLER_ANI {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A pool (meta.id) and its species (meta.species), and skani's pairs.
    tuple val(meta), path(skani)

    output:
    tuple val(meta), path("${meta.id}.ani.tsv")      , emit: matrix
    tuple val(meta), path("${meta.id}.ani_pairs.tsv"), emit: pairs
    tuple val(meta), path("${meta.id}.ani_mqc.tsv")  , emit: multiqc
    tuple val(meta), path("*.log")                   , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def species = (meta.species ?: []).collect { name -> "-s ${name}" }.join(' ')
    """
    plexembler ani \\
        ${skani} \\
        -p ${meta.id} \\
        -o . \\
        ${species} \\
        ${args}
    """

    stub:
    """
    printf 'species\\n' > ${meta.id}.ani.tsv
    printf 'species_a\\tspecies_b\\tani\\taligned_fraction_a\\taligned_fraction_b\\tverdict\\n' > ${meta.id}.ani_pairs.tsv
    printf '# id: plexembler_ani\\nSample\\tPool\\n' > ${meta.id}.ani_mqc.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
