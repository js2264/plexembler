//
// A pool's libraries split between its species (docs/concepts.md,
// deconvolution): every species' contigs in one reference; each WGS and Hi-C
// library aligned to it at once (bwa-mem2; -5SP for Hi-C); each pair given to
// the species of its mates, unique and near-identical (plexembler split).
// Reused at each stage whose contigs change (docs/plan.md, M2 and M3).
//

include { CAT_FASTQ        } from '../../../modules/nf-core/cat/fastq/main'
include { FIND_CONCATENATE as CAT_REFERENCE } from '../../../modules/nf-core/find/concatenate/main'
include { BWAMEM2_INDEX    } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM      } from '../../../modules/nf-core/bwamem2/mem/main'
include { PLEXEMBLER_SPLIT } from '../../../modules/local/plexembler/split/main'

workflow POOL_DECONVOLVE {

    take:
    ch_contigs   // channel: [ val(meta), path(fasta) ], each species' contigs, named <species>_<contig>
    ch_libraries // channel: [ val(meta), [ path(fastq) ] ], the pooled libraries
    stage        // string: the contigs' stage, which names the pool's reference (contigs, screened, genome)

    main:

    //
    // MODULE: the pool's reference, every species' contigs together
    //
    def ch_references = ch_contigs
        .map { meta, fasta -> [meta.pool, meta.id, fasta] }
        .groupTuple()
        .map { pool, species, fastas ->
            def sorted = [species, fastas].transpose().sort { pair -> pair[0] }
            [[id: "${pool}.${stage}".toString(), pool: pool, species: sorted.collect { pair -> pair[0] }], sorted.collect { pair -> pair[1] }]
        }
    CAT_REFERENCE(ch_references)
    BWAMEM2_INDEX(CAT_REFERENCE.out.file_out, 'bwamem2')

    //
    // MODULE: a library's runs as one pair of files
    //
    def ch_split_libraries = ch_libraries
        .filter { meta, _reads -> meta.type in ['wgs', 'hic'] }
        .branch { _meta, reads ->
            as_is: reads.size() == 2
            merge: true
        }
    CAT_FASTQ(ch_split_libraries.merge)
    def ch_paired = ch_split_libraries.as_is.mix(CAT_FASTQ.out.reads)

    //
    // MODULE: each library aligned to its pool's reference; Hi-C's mates on their own
    //
    def ch_to_align = ch_paired
        .map { meta, reads -> [meta.pool, meta, reads] }
        .combine(
            CAT_REFERENCE.out.file_out
                .join(BWAMEM2_INDEX.out.index)
                .map { meta, fasta, index -> [meta.pool, meta, fasta, index] },
            by: 0
        )
        .multiMap { _pool, meta, reads, reference, fasta, index ->
            reads: [meta + [species: reference.species, stage: stage], reads]
            index: [reference, index]
            fasta: [reference, fasta]
        }
    BWAMEM2_MEM(ch_to_align.reads, ch_to_align.index, ch_to_align.fasta, false)

    //
    // MODULE: each pair given to a species, or counted
    //
    PLEXEMBLER_SPLIT(BWAMEM2_MEM.out.bam)

    // Each species' pairs of each library: <species>.<library>[.holdout]_R{1,2}.fastq.gz
    def ch_species_reads = PLEXEMBLER_SPLIT.out.reads
        .flatMap { meta, files ->
            [files].flatten().collect { fastq ->
                def parts = fastq.name =~ /^([^.]+)\.(.+?)(\.holdout)?_R([12])\.fastq\.gz$/
                if (!parts) {
                    error("PLEXEMBLER_SPLIT wrote ${fastq.name}: not <species>.<library>[.holdout]_R<1|2>.fastq.gz")
                }
                def holdout = parts[0][3] != null
                ["${parts[0][1]}|${meta.id}|${holdout}".toString(), meta, parts[0][4], fastq]
            }
        }
        .groupTuple(size: 2)
        .map { key, metas, reads, files ->
            def (species, library, holdout) = key.tokenize('|')
            def library_meta = metas[0]
            def ordered = [reads, files].transpose().sort { pair -> pair[0] }.collect { pair -> pair[1] }
            [
                [
                    id      : species,
                    pool    : library_meta.pool,
                    library : library,
                    type    : library_meta.type,
                    enzymes : library_meta.enzymes,
                    holdout : holdout == 'true',
                    stage   : stage,
                ],
                ordered,
            ]
        }

    emit:
    reads      = ch_species_reads                // channel: [ val(meta), [ path(R1), path(R2) ] ], per species, library and part (meta.holdout)
    split      = PLEXEMBLER_SPLIT.out.json       // channel: [ val(meta), path(json) ], per library
    table      = PLEXEMBLER_SPLIT.out.table      // channel: [ val(meta), path(tsv) ]
    reference  = CAT_REFERENCE.out.file_out            // channel: [ val(meta), path(fasta) ], per pool
    contacts   = PLEXEMBLER_SPLIT.out.contacts   // channel: [ val(meta), path(tsv) ], per Hi-C library: its contacts between contigs
    multiqc_files = PLEXEMBLER_SPLIT.out.multiqc.map { _meta, table -> table } // channel: path
}
