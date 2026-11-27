process PLEXEMBLER_SCREEN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/js2264/plexembler:${workflow.manifest.version}"

    input:
    // A pool (meta.id): its contigs, each species' nanopore reads aligned to
    // them (PAF, in the order of meta.species), the contigs aligned to each
    // other, and its Hi-C libraries' contacts between contigs; meta.haploid,
    // its species of ploidy 1.
    tuple val(meta), path(contigs), path(reads, stageAs: 'reads/*'), path(identity), path(contacts, stageAs: 'contacts/*')

    output:
    tuple val(meta), path("*.screened.fa.gz")          , emit: fasta
    tuple val(meta), path("${meta.id}.screen.tsv")     , emit: table
    tuple val(meta), path("${meta.id}.screen_mqc.tsv") , emit: multiqc
    tuple val(meta), path("*.log")                     , emit: log
    tuple val("${task.process}"), val('plexembler'), eval("plexembler --version | sed 's/.* //'"), emit: versions_plexembler, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pafs = [reads].flatten()
    def species_reads = [meta.species, pafs].transpose().collect { species, paf -> "-r ${species}=${paf}" }.join(' ')
    def identity_arg = identity ? "-i ${identity}" : ''
    def contacts_args = [contacts].flatten().findAll { table -> table }.collect { table -> "-c ${table}" }.join(' ')
    def haploid_args = (meta.haploid ?: []).collect { species -> "--haploid ${species}" }.join(' ')
    """
    plexembler screen \\
        ${contigs} \\
        ${species_reads} \\
        ${identity_arg} \\
        ${contacts_args} \\
        ${haploid_args} \\
        -p ${meta.id} \\
        -o . \\
        ${args}
    # bgzip, so that each species' contigs can be indexed where they lie.
    for fasta in *.screened.fa.gz; do
        zcat \${fasta} | bgzip -@ ${task.cpus} > bgzipped.fa.gz
        mv bgzipped.fa.gz \${fasta}
    done
    """

    stub:
    """
    for species in ${meta.species.join(' ')}; do
        printf ">\${species}_contig_1\\nACGT\\n" | gzip -n > \${species}.screened.fa.gz
    done
    printf 'contig\\tspecies\\tverdict\\n' > ${meta.id}.screen.tsv
    printf '# id: plexembler_screen\\nSample\\tKept\\n' > ${meta.id}.screen_mqc.tsv
    touch ${meta.id}.plexembler_20260101000000.log
    """
}
