process PYPOLCA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypolca:0.5.0--pyhdfd78af_0' :
        'biocontainers/pypolca:0.5.0--pyhdfd78af_0' }"

    input:
    // A draft, and its short reads: a pair's two files, or single reads.
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("${prefix}.fa.gz") , emit: fasta
    tuple val(meta), path("${prefix}.report"), emit: report
    tuple val("${task.process}"), val('pypolca'), eval("pypolca --version 2>/dev/null | sed -n 's/^.*version //p'"), emit: versions_pypolca, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def files = [reads].flatten()
    def reads_args = files.size() == 2 ? "-1 ${files[0]} -2 ${files[1]}" : "-1 ${files[0]}"
    """
    pypolca run \\
        -a ${assembly} \\
        ${reads_args} \\
        -t ${task.cpus} \\
        -o pypolca \\
        -p ${prefix} \\
        ${args}
    gzip -cn pypolca/${prefix}_corrected.fasta > ${prefix}.fa.gz
    cp pypolca/${prefix}.report ${prefix}.report
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>${meta.id}_contig_1\\nACGT\\n' | gzip -n > ${prefix}.fa.gz
    touch ${prefix}.report
    """
}
