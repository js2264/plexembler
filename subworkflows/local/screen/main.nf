//
// The contig screen (docs/plan.md, M2): each species' nanopore reads aligned to
// every species' contigs of its pool, the pool's contigs aligned to each other,
// and the pool's Hi-C contacts between contigs, weighed by plexembler screen:
// each contig kept, moved to another species, dropped, or named an organelle.
//

include { FIND_CONCATENATE as CAT_POOL         } from '../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as CAT_POOL_TARGET  } from '../../../modules/nf-core/find/concatenate/main'
include { MINIMAP2_ALIGN as MINIMAP2_READS     } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_IDENTITY  } from '../../../modules/nf-core/minimap2/align/main'
include { PLEXEMBLER_SCREEN                    } from '../../../modules/local/plexembler/screen/main'

workflow SCREEN {

    take:
    ch_contigs  // channel: [ val(meta), path(fasta) ], each species' contigs
    ch_reads    // channel: [ val(meta), path(fastq) ], each species' nanopore reads
    ch_contacts // channel: [ val(meta), path(tsv) ], each pool's Hi-C contacts between contigs (meta.pool)

    main:

    //
    // MODULE: the pool's contigs together
    //
    def ch_pools = ch_contigs
        .map { meta, fasta -> [meta.pool, meta.id, fasta] }
        .groupTuple()
        .map { pool, species, fastas ->
            def sorted = [species, fastas].transpose().sort { pair -> pair[0] }
            [[id: pool, species: sorted.collect { pair -> pair[0] }], sorted.collect { pair -> pair[1] }]
        }
    CAT_POOL(ch_pools)
    // The same contigs again, under another name: the all-against-all's
    // target, which cannot be staged under its query's name.
    CAT_POOL_TARGET(ch_pools)

    //
    // MODULES: each species' reads on the pool's contigs; the contigs on each other
    //
    def ch_reads_on_pool = ch_reads
        .map { meta, reads -> [meta.pool, meta, reads] }
        .combine(CAT_POOL.out.file_out.map { meta, fasta -> [meta.id, meta, fasta] }, by: 0)
        .multiMap { _pool, meta, reads, pool_meta, fasta ->
            reads: [meta + [pool: pool_meta.id], reads]
            reference: [pool_meta, fasta]
        }
    MINIMAP2_READS(ch_reads_on_pool.reads, ch_reads_on_pool.reference, false, '', false, false)
    def ch_all_against_all = CAT_POOL.out.file_out
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(CAT_POOL_TARGET.out.file_out.map { meta, fasta -> [meta.id, fasta] })
        .multiMap { _pool, meta, fasta, target ->
            query: [meta, fasta]
            target: [meta, target]
        }
    MINIMAP2_IDENTITY(ch_all_against_all.query, ch_all_against_all.target, false, '', true, false)

    //
    // MODULE: the verdicts
    //
    def ch_screen = CAT_POOL.out.file_out
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(
            MINIMAP2_READS.out.paf
                .map { meta, paf -> [meta.pool, meta.id, paf] }
                .groupTuple()
        )
        .join(MINIMAP2_IDENTITY.out.paf.map { meta, paf -> [meta.id, paf] })
        .join(
            ch_contacts.map { meta, table -> [meta.pool, table] }.groupTuple(),
            remainder: true
        )
        // A pool without contacts joins with a null; a pool's contacts without
        // the rest (its contigs failed) come as a shorter row, left out.
        .filter { row -> row.size() == 7 && row[1] != null }
        .map { _pool, meta, fasta, species, pafs, identity, contacts ->
            // The PAFs in the order of meta.species.
            def ordered = [species, pafs].transpose().sort { pair -> meta.species.indexOf(pair[0]) }
            [meta, fasta, ordered.collect { pair -> pair[1] }, identity, contacts ?: []]
        }
    PLEXEMBLER_SCREEN(ch_screen)

    // Each species' screened contigs, its meta as it came.
    def ch_screened = ch_contigs
        .map { meta, _fasta -> [meta.id, meta] }
        .join(
            PLEXEMBLER_SCREEN.out.fasta
                .flatMap { _meta, fastas -> [fastas].flatten().collect { fasta -> [fasta.name.tokenize('.')[0], fasta] } }
        )
        .map { _id, meta, fasta -> [meta, fasta] }

    emit:
    contigs       = ch_screened                       // channel: [ val(meta), path(fasta) ]
    table         = PLEXEMBLER_SCREEN.out.table       // channel: [ val(meta), path(tsv) ], per pool
    multiqc_files = PLEXEMBLER_SCREEN.out.multiqc.map { _meta, table -> table } // channel: path
}
