process PLEXEMBLER_GENES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A species' genome and its gene models (GTF or GFF3); meta.genetic_code.
    tuple val(meta), path(genome), path(models)
    // What predicted them: the GFF3's second column.
    val source

    output:
    tuple val(meta), path("${meta.id}.gff3.gz")         , emit: gff3
    tuple val(meta), path("${meta.id}.gtf.gz")          , emit: gtf
    tuple val(meta), path("${meta.id}.proteins.faa.gz") , emit: proteins
    tuple val(meta), path("${meta.id}.cds.fna.gz")      , emit: cds
    tuple val(meta), path("${meta.id}.genes.tsv")       , emit: table
    tuple val(meta), path("${meta.id}.genes_mqc.tsv")   , emit: multiqc
    tuple val(meta), path("*.log")                      , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genetic_code = meta.genetic_code ?: 1
    """
    plexembler genes \\
        ${genome} \\
        ${models} \\
        -p ${meta.id} \\
        -o . \\
        -g ${genetic_code} \\
        -s ${source} \\
        ${args}
    """

    stub:
    """
    for suffix in gff3 gtf proteins.faa cds.fna; do
        echo "" | gzip -n > ${meta.id}.\${suffix}.gz
    done
    # A gene, so that the steps after it run in a stub run too.
    printf 'gene\\ttranscripts\\n${meta.id}_000010\\t1\\n' > ${meta.id}.genes.tsv
    printf '# id: plexembler_genes\\nSample\\tGenes\\n${meta.id}\\t0\\n' > ${meta.id}.genes_mqc.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
