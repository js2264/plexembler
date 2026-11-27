//
// Each species' genome, scaffolded with its own Hi-C pairs and finished
// (docs/plan.md, M3): its screened contigs; its Hi-C reads from the pool's
// split of the same contigs, aligned to them (bwa-mem2 -5SP); YaHS, its motifs
// from the library's enzymes (none for Micro-C); the scaffolds named,
// chromosomes first
// (plexembler finalize).
//

include { CAT_FASTQ as CAT_FASTQ_HIC              } from '../../../modules/nf-core/cat/fastq/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_SPECIES  } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_SPECIES      } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CONTIGS } from '../../../modules/nf-core/samtools/faidx/main'
include { YAHS                                    } from '../../../modules/nf-core/yahs/main'
include { PLEXEMBLER_FINALIZE                     } from '../../../modules/local/plexembler/finalize/main'

workflow SCAFFOLD {

    take:
    ch_contigs // channel: [ val(meta), path(fasta) ], each species' screened contigs
    ch_reads   // channel: [ val(meta), [ path(R1), path(R2) ] ], each species' split reads (POOL_DECONVOLVE.out.reads)

    main:

    //
    // MODULE: each species' Hi-C pairs, held-out ones included, as one pair of files
    //
    def ch_hic = ch_reads
        .filter { meta, _reads -> meta.type == 'hic' }
        .map { meta, reads -> [meta.id, meta.enzymes, reads] }
        .groupTuple()
        .map { species, enzymes, reads ->
            [[id: species, enzymes: enzymes.unique().join(','), single_end: false], reads.flatten()]
        }
        .branch { _meta, files ->
            merge: files.size() > 2
            as_is: true
        }
    CAT_FASTQ_HIC(ch_hic.merge)
    def ch_hic_pairs = ch_hic.as_is.mix(CAT_FASTQ_HIC.out.reads)

    //
    // MODULES: each species' Hi-C pairs on its own contigs; the contigs' index
    //
    def ch_to_scaffold = ch_contigs
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_hic_pairs.map { meta, reads -> [meta.id, meta.enzymes, reads] })
        .map { _id, meta, fasta, enzymes, reads -> [meta + [enzymes: enzymes], fasta, reads] }
    BWAMEM2_INDEX_SPECIES(ch_to_scaffold.map { meta, fasta, _reads -> [meta, fasta] }, 'bwamem2')
    def ch_aligning = ch_to_scaffold
        .map { meta, fasta, reads -> [meta.id, meta, fasta, reads] }
        .join(BWAMEM2_INDEX_SPECIES.out.index.map { meta, index -> [meta.id, index] })
        .multiMap { _id, meta, fasta, reads, index ->
            reads: [meta, reads]
            index: [meta, index]
            fasta: [meta, fasta]
        }
    BWAMEM2_MEM_SPECIES(ch_aligning.reads, ch_aligning.index, ch_aligning.fasta, false)
    SAMTOOLS_FAIDX_CONTIGS(ch_to_scaffold.map { meta, fasta, _reads -> [meta, fasta, []] }, false)

    //
    // MODULE: YaHS
    //
    def ch_yahs = ch_to_scaffold
        .map { meta, fasta, _reads -> [meta.id, meta, fasta] }
        .join(SAMTOOLS_FAIDX_CONTIGS.out.fai.map { meta, fai -> [meta.id, fai] })
        .join(BWAMEM2_MEM_SPECIES.out.bam.map { meta, bam -> [meta.id, bam] })
        .map { _id, meta, fasta, fai, bam -> [meta, fasta, fai, bam, []] }
    YAHS(ch_yahs)

    //
    // MODULE: the scaffolds named: chromosomes by size, the organelle, the
    // unplaced. A species without Hi-C: its contigs, named as a genome all the same.
    //
    def ch_without_hic = ch_contigs
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_hic_pairs.map { meta, _reads -> [meta.id, true] }, remainder: true)
        .filter { _id, meta, _fasta, has_hic -> meta != null && has_hic == null }
        .map { _id, meta, fasta, _has_hic -> [meta, fasta, []] }
    PLEXEMBLER_FINALIZE(
        YAHS.out.scaffolds_fasta
            .join(YAHS.out.scaffolds_agp)
            .mix(ch_without_hic)
    )

    emit:
    genome  = PLEXEMBLER_FINALIZE.out.fasta  // channel: [ val(meta), path(fasta) ], bgzipped, every species
    fai     = PLEXEMBLER_FINALIZE.out.fai    // channel: [ val(meta), path(fai) ]
    sizes   = PLEXEMBLER_FINALIZE.out.sizes  // channel: [ val(meta), path(sizes) ]
    agp     = PLEXEMBLER_FINALIZE.out.agp    // channel: [ val(meta), path(agp) ], the scaffolded species
    hic     = ch_hic_pairs                   // channel: [ val(meta), [ path(R1), path(R2) ] ], each species' Hi-C pairs
}
