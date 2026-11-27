process AUGUSTUS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augustus:3.5.0--pl5321h5653ebf_10' :
        'biocontainers/augustus:3.5.0--pl5321h5653ebf_10' }"

    input:
    // A species' genome, soft-masked, gzipped or not.
    tuple val(meta), path(fasta)
    // One of AUGUSTUS's trained species (augustus --species=help).
    val species

    output:
    tuple val(meta), path("${prefix}.gff3"), emit: gff3
    tuple val("${task.process}"), val('augustus'), eval("augustus --version 2>&1 | sed -n 's/AUGUSTUS (\\(.*\\)) is a gene .*/\\1/p'"), emit: versions_augustus, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.augustus"
    def plain = fasta.name.endsWith('.gz') ? "gunzip -c ${fasta} > genome.fa" : "cp ${fasta} genome.fa"
    """
    ${plain}
    augustus \\
        --species=${species} \\
        --gff3=on \\
        --softmasking=1 \\
        ${args} \\
        genome.fa \\
        > ${prefix}.gff3
    rm genome.fa
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.augustus"
    """
    printf '##gff-version 3\\n' > ${prefix}.gff3
    """
}
