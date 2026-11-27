//
// Each species' genome, scaffolded with its own Hi-C pairs and finished
// (docs/plan.md, M3): its screened contigs; its Hi-C reads from the pool's
// split of the same contigs, aligned to them (bwa-mem2 -5SP); YaHS, its motifs
// from the library's enzymes (none for Micro-C), then YaHS again on its own
// scaffolds (-a), which joins what its first pass left apart; YaHS's joins
// scored by the same pairs, parsed and deduplicated on the contigs and placed
// on the scaffolds through YaHS's AGP (plexembler junctions --contigs); the
// scaffolds named, chromosomes first, once the joins the pairs do not support
// are cut (plexembler finalize --cut).
//

include { CAT_FASTQ as CAT_FASTQ_HIC              } from '../../../modules/nf-core/cat/fastq/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_SPECIES  } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_SPECIES      } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CONTIGS } from '../../../modules/nf-core/samtools/faidx/main'
include { YAHS                                    } from '../../../modules/nf-core/yahs/main'
include { YAHS as YAHS_RESCAFFOLD                 } from '../../../modules/nf-core/yahs/main'
include { PAIRTOOLS_PARSE as PAIRTOOLS_PARSE_CONTIGS } from '../../../modules/nf-core/pairtools/parse/main'
include { PAIRTOOLS_SORT as PAIRTOOLS_SORT_CONTIGS   } from '../../../modules/nf-core/pairtools/sort/main'
include { PAIRTOOLS_DEDUP as PAIRTOOLS_DEDUP_CONTIGS } from '../../../modules/nf-core/pairtools/dedup/main'
include { PLEXEMBLER_JUNCTIONS as JUNCTIONS_SCAFFOLDS } from '../../../modules/local/plexembler/junctions/main'
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
    // MODULES: YaHS, then YaHS on its own scaffolds (-a): a scaffold of the
    // first pass is one piece of the second, whose ends gather the links that
    // its contigs' ends split (conf/modules.config, its resolutions)
    //
    def ch_yahs = ch_to_scaffold
        .map { meta, fasta, _reads -> [meta.id, meta, fasta] }
        .join(SAMTOOLS_FAIDX_CONTIGS.out.fai.map { meta, fai -> [meta.id, fai] })
        .join(BWAMEM2_MEM_SPECIES.out.bam.map { meta, bam -> [meta.id, bam] })
        .map { _id, meta, fasta, fai, bam -> [meta, fasta, fai, bam, []] }
    YAHS(ch_yahs)
    def ch_rescaffold = ch_yahs
        .map { meta, fasta, fai, bam, _agp -> [meta.id, meta, fasta, fai, bam] }
        .join(YAHS.out.scaffolds_agp.map { meta, agp -> [meta.id, agp] })
        .map { _id, meta, fasta, fai, bam, agp -> [meta, fasta, fai, bam, agp] }
    YAHS_RESCAFFOLD(ch_rescaffold)

    //
    // MODULES: YaHS's joins, each scored by the pairs across it against the
    // species' contact decay: the pairs on the contigs, parsed and deduplicated
    // as the genome's are (MAPS), each mate placed on its scaffold by YaHS's AGP
    //
    def ch_parse = BWAMEM2_MEM_SPECIES.out.bam
        .map { meta, bam -> [meta.id, meta, bam] }
        .join(SAMTOOLS_FAIDX_CONTIGS.out.fai.map { meta, fai -> [meta.id, fai] })
        .multiMap { _id, meta, bam, fai ->
            bam: [meta, bam]
            sizes: fai
        }
    PAIRTOOLS_PARSE_CONTIGS(ch_parse.bam, ch_parse.sizes)
    PAIRTOOLS_SORT_CONTIGS(PAIRTOOLS_PARSE_CONTIGS.out.pairsam)
    PAIRTOOLS_DEDUP_CONTIGS(PAIRTOOLS_SORT_CONTIGS.out.sorted)
    JUNCTIONS_SCAFFOLDS(
        PAIRTOOLS_DEDUP_CONTIGS.out.pairs
            .join(YAHS_RESCAFFOLD.out.scaffolds_agp)
    )

    //
    // MODULE: the scaffolds named: chromosomes by size, the organelle, the
    // unplaced; first, each join the pairs do not support cut. A species
    // without Hi-C: its contigs, named as a genome all the same.
    //
    def ch_without_hic = ch_contigs
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_hic_pairs.map { meta, _reads -> [meta.id, true] }, remainder: true)
        .filter { _id, meta, _fasta, has_hic -> meta != null && has_hic == null }
        .map { _id, meta, fasta, _has_hic -> [meta, fasta, [], []] }
    PLEXEMBLER_FINALIZE(
        YAHS_RESCAFFOLD.out.scaffolds_fasta
            .join(YAHS_RESCAFFOLD.out.scaffolds_agp)
            .join(JUNCTIONS_SCAFFOLDS.out.table)
            .mix(ch_without_hic)
    )

    emit:
    genome  = PLEXEMBLER_FINALIZE.out.fasta  // channel: [ val(meta), path(fasta) ], bgzipped, every species
    fai     = PLEXEMBLER_FINALIZE.out.fai    // channel: [ val(meta), path(fai) ]
    sizes   = PLEXEMBLER_FINALIZE.out.sizes  // channel: [ val(meta), path(sizes) ]
    names   = PLEXEMBLER_FINALIZE.out.names  // channel: [ val(meta), path(tsv) ], each genome's names, lengths and roles
    agp     = PLEXEMBLER_FINALIZE.out.agp    // channel: [ val(meta), path(agp) ], the scaffolded species
    joins   = JUNCTIONS_SCAFFOLDS.out.table  // channel: [ val(meta), path(tsv) ], YaHS's joins scored, before the cut
    hic     = ch_hic_pairs                   // channel: [ val(meta), [ path(R1), path(R2) ] ], each species' Hi-C pairs
}
