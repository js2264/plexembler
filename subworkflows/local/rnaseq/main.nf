//
// Each pooled RNA-seq library split between its pool's species (docs/plan.md,
// M4): the pool's genomes as one STAR index; each library aligned to it once
// (introns capped at 5 kb: compact genomes get false splices otherwise);
// `plexembler split --type rnaseq`, each pair to the species of its two mates'
// unique, near-identical alignments, written as that species' alignments,
// which lie on its genome already. No second alignment.
//

include { FIND_CONCATENATE as CAT_GENOMES        } from '../../../modules/nf-core/find/concatenate/main'
include { STAR_GENOMEGENERATE                    } from '../../../modules/nf-core/star/genomegenerate/main'
include { CAT_FASTQ as CAT_FASTQ_RNASEQ          } from '../../../modules/nf-core/cat/fastq/main'
include { STAR_ALIGN                             } from '../../../modules/nf-core/star/align/main'
include { PLEXEMBLER_SPLIT as PLEXEMBLER_SPLIT_RNASEQ } from '../../../modules/local/plexembler/split/main'

workflow RNASEQ {

    take:
    ch_genomes   // channel: [ val(meta), path(fasta) ], each species' genome (meta.id, meta.pool)
    ch_libraries // channel: [ val(meta), [ path(fastq) ] ], the pooled libraries

    main:

    //
    // MODULES: each pool's genomes, as one plain FASTA, and its STAR index
    //
    def ch_pools = ch_genomes
        .map { meta, fasta -> [meta.pool, meta.id, fasta] }
        .groupTuple()
        .map { pool, species, fastas ->
            def sorted = [species, fastas].transpose().sort { pair -> pair[0] }
            [[id: "${pool}.genomes".toString(), pool: pool, species: sorted.collect { pair -> pair[0] }], sorted.collect { pair -> pair[1] }]
        }
    CAT_GENOMES(ch_pools)
    STAR_GENOMEGENERATE(CAT_GENOMES.out.file_out, [[:], []])

    //
    // MODULES: each RNA-seq library, its runs as one pair of files, aligned
    // to its pool's genomes
    //
    def ch_rnaseq = ch_libraries
        .filter { meta, _reads -> meta.type == 'rnaseq' }
        .branch { _meta, reads ->
            as_is: reads.size() == 2
            merge: true
        }
    CAT_FASTQ_RNASEQ(ch_rnaseq.merge)
    def ch_align = ch_rnaseq.as_is
        .mix(CAT_FASTQ_RNASEQ.out.reads)
        .map { meta, reads -> [meta.pool, meta, reads] }
        .combine(STAR_GENOMEGENERATE.out.index.map { meta, index -> [meta.pool, meta, index] }, by: 0)
        .multiMap { _pool, meta, reads, reference, index ->
            reads: [meta + [species: reference.species], reads]
            index: [reference, index]
        }
    STAR_ALIGN(ch_align.reads, ch_align.index, [[:], []], true)

    //
    // MODULE: each library split between its pool's species, as BAMs
    //
    PLEXEMBLER_SPLIT_RNASEQ(STAR_ALIGN.out.bam)

    // Each species' alignments of each library: its meta, the library's.
    def ch_bams = PLEXEMBLER_SPLIT_RNASEQ.out.bam
        .join(PLEXEMBLER_SPLIT_RNASEQ.out.bai)
        .flatMap { meta, bams, bais ->
            def indexes = [bais].flatten().collectEntries { bai -> [bai.name - '.bai', bai] }
            [bams].flatten().collect { bam ->
                def species = bam.name.tokenize('.')[0]
                [
                    [id: species, pool: meta.pool, library: meta.id, condition: meta.condition, replicate: meta.replicate, single_end: false],
                    bam,
                    indexes[bam.name],
                ]
            }
        }

    def ch_multiqc_files = PLEXEMBLER_SPLIT_RNASEQ.out.multiqc.map { _meta, table -> table }
        .mix(STAR_ALIGN.out.log_final.map { _meta, log -> log })

    emit:
    bams          = ch_bams                             // channel: [ val(meta), path(bam), path(bai) ], each species' alignments of each library
    table         = PLEXEMBLER_SPLIT_RNASEQ.out.table  // channel: [ val(meta), path(tsv) ]
    multiqc_files = ch_multiqc_files                    // channel: path
}
