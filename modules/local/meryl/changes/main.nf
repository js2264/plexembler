process MERYL_CHANGES {
    tag "${meta.id}"
    label 'process_low'

    // meryl --version reports a git revision, not the release: the versions
    // output says 1.4.1, the container's, as nf-core's meryl modules do.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7df6a8db5e64dc18ecc1aa31b12308fc9d5e07ddefcc64f5ffffa1716c029da8/data'
        : 'community.wave.seqera.io/library/meryl:1.4.1--2481ef1dd13fafa3'}"

    input:
    // A polishing, the draft it started from, and the k-mers of the reads held
    // out of the polishing (a meryl database of k = kmer_size).
    tuple val(meta), path(draft, stageAs: 'draft/*'), path(polished, stageAs: 'polished/*'), path(reads_db)
    val kmer_size

    output:
    tuple val(meta), path("${prefix}.kmers.tsv"), emit: kmers
    tuple val("${task.process}"), val('meryl'), val('1.4.1'), emit: versions_meryl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "memory=${task.memory.toGiga()}" : ''
    """
    # The draft's k-mers and the polishing's, then those the held-out reads lack.
    meryl count k=${kmer_size} threads=${task.cpus} ${memory} ${args} ${draft} output draft.meryl
    meryl count k=${kmer_size} threads=${task.cpus} ${memory} ${args} ${polished} output polished.meryl
    meryl difference draft.meryl ${reads_db} output draft.absent.meryl
    meryl difference polished.meryl ${reads_db} output polished.absent.meryl
    # Where the two differ: the draft's absent k-mers that the polishing
    # removed (fixed), and the absent ones it made (introduced).
    meryl difference draft.absent.meryl polished.absent.meryl output fixed.meryl
    meryl difference polished.absent.meryl draft.absent.meryl output introduced.meryl
    # The first `distinct` line is the count (the histogram's header has one
    # too); awk reads on to the end, so that meryl never writes to a closed pipe.
    fixed=\$(meryl statistics fixed.meryl | awk '\$1 == "distinct" && !n++ { print \$2 }')
    introduced=\$(meryl statistics introduced.meryl | awk '\$1 == "distinct" && !n++ { print \$2 }')
    printf 'fixed\\tintroduced\\n%s\\t%s\\n' "\${fixed:-0}" "\${introduced:-0}" > ${prefix}.kmers.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'fixed\\tintroduced\\n0\\t0\\n' > ${prefix}.kmers.tsv
    """
}
