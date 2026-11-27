process PLEXEMBLER_SPLIT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A pooled library (meta.id, meta.pool, meta.type), and its alignments to
    // the pool's contigs as bwa-mem2 writes them: a pair's records together.
    tuple val(meta), path(alignments)

    output:
    tuple val(meta), path("*_R{1,2}.fastq.gz")              , emit: reads, optional: true
    tuple val(meta), path("${meta.pool}.${meta.id}.split.tsv")   , emit: table
    tuple val(meta), path("${meta.pool}.${meta.id}.split.json")  , emit: json
    tuple val(meta), path("${meta.pool}.${meta.id}.split_mqc.tsv"), emit: multiqc
    tuple val(meta), path("${meta.pool}.${meta.id}.contacts.tsv"), emit: contacts, optional: true
    tuple val(meta), path("*.log")                          , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools view -h ${alignments} \\
        | plexembler split - \\
            -p ${meta.pool} \\
            -o . \\
            -l ${meta.id} \\
            -t ${meta.type} \\
            ${args}
    """

    stub:
    """
    for species in ${(meta.species ?: ['sp']).join(' ')}; do
        echo "" | gzip -n > \${species}.${meta.id}_R1.fastq.gz
        echo "" | gzip -n > \${species}.${meta.id}_R2.fastq.gz
        echo "" | gzip -n > \${species}.${meta.id}.holdout_R1.fastq.gz
        echo "" | gzip -n > \${species}.${meta.id}.holdout_R2.fastq.gz
    done
    printf 'species\\tpairs\\tshare\\tholdout\\tduplicates\\n' > ${meta.pool}.${meta.id}.split.tsv
    echo '{"pairs": 0, "species": {}}' > ${meta.pool}.${meta.id}.split.json
    printf '# id: plexembler_split\\nSample\\tPairs\\n' > ${meta.pool}.${meta.id}.split_mqc.tsv
    if [ "${meta.type}" = "hic" ]; then printf 'contig_a\\tcontig_b\\tpairs\\n' > ${meta.pool}.${meta.id}.contacts.tsv; fi
    touch ${meta.pool}.plexembler_20260101000000.log
    """
}
