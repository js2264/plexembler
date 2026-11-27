process PLEXEMBLER_SELECT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A species' candidates, each list in the order of `candidates`; a check a
    // candidate lacks is the placeholder NO_FILE (assets/), written as NA.
    tuple val(meta), val(candidates), path(fastas, stageAs: 'fasta?/*'), path(buscos, stageAs: 'busco?/*'), path(qvs, stageAs: 'qv?/*'), path(stats, stageAs: 'stats?/*')
    val step

    output:
    tuple val(meta), path("${meta.id}.${step}.fa.gz")          , emit: fasta
    tuple val(meta), path("${meta.id}.${step}.selection.tsv")  , emit: selection
    tuple val(meta), path("${meta.id}.${step}_selection_mqc.tsv"), emit: multiqc
    tuple val(meta), path("*.log")                             , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // The samplesheet's choice, when it named one (its `assembler` column).
    def keep = meta.keep ? "--keep ${meta.keep}" : ''
    def columns = [candidates, fastas, buscos, qvs, stats].collect { item -> [item].flatten() }
    if (columns*.size().unique().size() != 1) {
        error("PLEXEMBLER_SELECT: ${meta.id} has ${columns*.size()} candidates, sequences, BUSCO summaries, QVs and statistics")
    }
    def rows = columns.transpose().collect { name, fasta, busco, qv, stat ->
        "printf '%s\\t%s\\t%s\\t%s\\t%s\\n' '${name}' '${fasta}' '${checkCell(busco)}' '${checkCell(qv)}' '${checkCell(stat)}' >> candidates.tsv"
    }.join('\n    ')
    """
    printf 'candidate\\tfasta\\tbusco\\tqv\\tstats\\n' > candidates.tsv
    ${rows}

    plexembler select \\
        candidates.tsv \\
        --step ${step} \\
        -p ${meta.id} \\
        -o . \\
        ${keep} \\
        ${args}
    # bgzip, so that the kept sequence can be indexed where it lies.
    zcat ${meta.id}.${step}.fa.gz | bgzip -@ ${task.cpus} > bgzipped.fa.gz
    mv bgzipped.fa.gz ${meta.id}.${step}.fa.gz
    """

    stub:
    """
    printf '>${meta.id}_contig_1\\nACGT\\n' | gzip -n > ${meta.id}.${step}.fa.gz
    printf 'candidate\\tkept\\treason\\n' > ${meta.id}.${step}.selection.tsv
    printf '# id: plexembler_${step}_selection\\nSample\\tKept\\n' > ${meta.id}.${step}_selection_mqc.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}

//
// A check's cell in the table of candidates: NA for the placeholder NO_FILE
//
def checkCell(check) {
    def staged = check.toString()
    return staged.tokenize('/')[-1] == 'NO_FILE' ? 'NA' : staged
}
