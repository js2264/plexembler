//
// Each species' Hi-C map, on its final genome (docs/plan.md, M3): its Hi-C pairs
// aligned to the genome (bwa-mem2 -5SP), parsed, sorted and deduplicated
// (pairtools), binned and zoomed (cooler): <species>.mcool.
//

include { BWAMEM2_INDEX as BWAMEM2_INDEX_GENOME } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_GENOME     } from '../../../modules/nf-core/bwamem2/mem/main'
include { PAIRTOOLS_PARSE                       } from '../../../modules/nf-core/pairtools/parse/main'
include { PAIRTOOLS_SORT                        } from '../../../modules/nf-core/pairtools/sort/main'
include { PAIRTOOLS_DEDUP                       } from '../../../modules/nf-core/pairtools/dedup/main'
include { COOLER_CLOAD                          } from '../../../modules/nf-core/cooler/cload/main'
include { COOLER_ZOOMIFY                        } from '../../../modules/nf-core/cooler/zoomify/main'

workflow MAPS {

    take:
    ch_genome  // channel: [ val(meta), path(fasta) ], each species' genome
    ch_sizes   // channel: [ val(meta), path(sizes) ], its chromosome sizes
    ch_hic     // channel: [ val(meta), [ path(R1), path(R2) ] ], its Hi-C pairs
    resolution // integer: the finest bin, in bases

    main:

    //
    // MODULES: the Hi-C pairs on the genome
    //
    BWAMEM2_INDEX_GENOME(ch_genome, 'bwamem2')
    def ch_align = ch_genome
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(BWAMEM2_INDEX_GENOME.out.index.map { meta, index -> [meta.id, index] })
        .join(ch_hic.map { meta, reads -> [meta.id, reads] })
        .multiMap { _id, meta, fasta, index, reads ->
            reads: [meta, reads]
            index: [meta, index]
            fasta: [meta, fasta]
        }
    BWAMEM2_MEM_GENOME(ch_align.reads, ch_align.index, ch_align.fasta, false)

    //
    // MODULES: the pairs, parsed, sorted, their duplicates dropped
    //
    def ch_parse = BWAMEM2_MEM_GENOME.out.bam
        .map { meta, bam -> [meta.id, meta, bam] }
        .join(ch_sizes.map { meta, sizes -> [meta.id, sizes] })
        .multiMap { _id, meta, bam, sizes ->
            bam: [meta, bam]
            sizes: sizes
        }
    PAIRTOOLS_PARSE(ch_parse.bam, ch_parse.sizes)
    PAIRTOOLS_SORT(PAIRTOOLS_PARSE.out.pairsam)
    PAIRTOOLS_DEDUP(PAIRTOOLS_SORT.out.sorted)

    //
    // MODULES: the map, at every resolution from the finest
    //
    def ch_cload = PAIRTOOLS_DEDUP.out.pairs
        .map { meta, pairs -> [meta.id, meta, pairs] }
        .join(ch_sizes.map { meta, sizes -> [meta.id, meta, sizes] })
        .multiMap { _id, meta, pairs, sizes_meta, sizes ->
            pairs: [meta, pairs, []]
            sizes: [sizes_meta, sizes]
        }
    COOLER_CLOAD(ch_cload.pairs, ch_cload.sizes, 'pairs', resolution)
    COOLER_ZOOMIFY(COOLER_CLOAD.out.cool)

    def ch_multiqc_files = PAIRTOOLS_PARSE.out.stat.map { _meta, stat -> stat }
        .mix(PAIRTOOLS_DEDUP.out.stat.map { _meta, stat -> stat })

    emit:
    pairs         = PAIRTOOLS_DEDUP.out.pairs // channel: [ val(meta), path(pairs) ]
    mcool         = COOLER_ZOOMIFY.out.mcool  // channel: [ val(meta), path(mcool) ]
    multiqc_files = ch_multiqc_files          // channel: path
}
