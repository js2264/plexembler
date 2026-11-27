process PLEXEMBLER_DOTPLOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A species (meta.id): minimap2's alignments of its genome to its closest
    // known genome, the two genomes, and plexembler closest's table, which
    // names the second.
    tuple val(meta), path(paf), path(genome), path(reference), path(closest)

    output:
    tuple val(meta), path("${meta.id}.dotplot.svg")  , emit: svg
    tuple val(meta), path("${meta.id}.dotplot.json") , emit: figure
    tuple val(meta), path("${meta.id}.dotplot.tsv")  , emit: table
    tuple val(meta), path("*.log")                   , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    plexembler dotplot \\
        ${paf} \\
        -p ${meta.id} \\
        -o . \\
        -g ${genome} \\
        -r ${reference} \\
        -c ${closest} \\
        ${args}
    """

    stub:
    """
    printf '<svg xmlns="http://www.w3.org/2000/svg"></svg>\\n' > ${meta.id}.dotplot.svg
    printf '{"data": [], "layout": {}}\\n' > ${meta.id}.dotplot.json
    printf 'genome\\tsequence\\tlabel\\tlength_bp\\taligned_bp\\taligned_fraction\\tidentity\\tbest_match\\tbest_match_share\\tstrand\\n' > ${meta.id}.dotplot.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
