process PLEXEMBLER_FINALIZE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A species' scaffolds, and the scaffolder's AGP; meta.n_chromosomes, when known.
    tuple val(meta), path(fasta), path(agp)

    output:
    tuple val(meta), path("${meta.id}.fa.gz")          , emit: fasta
    tuple val(meta), path("${meta.id}.fa.gz.fai")      , emit: fai
    tuple val(meta), path("${meta.id}.fa.gz.gzi")      , emit: gzi
    tuple val(meta), path("${meta.id}.chrom.sizes")    , emit: sizes
    tuple val(meta), path("${meta.id}.names.tsv")      , emit: names
    tuple val(meta), path("${meta.id}.agp")            , emit: agp, optional: true
    tuple val(meta), path("*.log")                     , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def agp_arg = agp ? "-a ${agp}" : ''
    def chromosomes = meta.n_chromosomes ? "-n ${meta.n_chromosomes}" : ''
    """
    plexembler finalize \\
        ${fasta} \\
        -p ${meta.id} \\
        -o finalized \\
        ${agp_arg} \\
        ${chromosomes} \\
        ${args}
    mv finalized/* .
    # bgzip, so that the genome can be indexed and read where it lies.
    zcat ${meta.id}.fa.gz | bgzip -@ ${task.cpus} > genome.fa.gz
    mv genome.fa.gz ${meta.id}.fa.gz
    samtools faidx ${meta.id}.fa.gz
    """

    stub:
    """
    printf '>${meta.id}_chr01\\nACGT\\n' | gzip -n > ${meta.id}.fa.gz
    touch ${meta.id}.fa.gz.fai ${meta.id}.fa.gz.gzi
    printf '${meta.id}_chr01\\t4\\n' > ${meta.id}.chrom.sizes
    printf 'name\\toriginal_name\\tlength_bp\\trole\\n' > ${meta.id}.names.tsv
    touch ${meta.id}.agp ${meta.id}.plexembler_20260101000000.log
    """
}
