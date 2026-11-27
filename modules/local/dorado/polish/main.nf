process DORADO_POLISH {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    // Dorado has no conda package: a container engine is needed.
    container "docker.io/ontresearch/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f"

    input:
    // A species' unaligned reads (Dorado's BAMs, which name the basecalling
    // model) and its draft contigs.
    tuple val(meta), path(reads, stageAs: 'reads/*'), path(draft)

    output:
    tuple val(meta), path("${prefix}.fa.gz"), emit: fasta
    tuple val("${task.process}"), val('dorado'), eval("dorado --version 2>&1 | head -n 1 | sed 's/^.*dorado //; s/+.*//'"), emit: versions_dorado, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.containerEngine == null) {
        error("DORADO_POLISH: Dorado runs in its container only; use -profile docker, apptainer or singularity")
    }
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def device = task.accelerator ? '--device cuda:all' : '--device cpu'
    """
    gzip -cdf ${draft} > draft.fa
    dorado aligner \\
        --threads ${task.cpus} \\
        --output-dir aligned \\
        ${args} \\
        draft.fa \\
        reads/

    dorado polish \\
        --threads ${task.cpus} \\
        ${device} \\
        ${args2} \\
        aligned \\
        draft.fa \\
        > ${prefix}.fa
    gzip -n ${prefix}.fa
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>${meta.id}_contig_1\\nACGT\\n' | gzip -n > ${prefix}.fa.gz
    """
}
