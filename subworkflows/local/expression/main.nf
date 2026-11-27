//
// Each species' expression and tracks (docs/plan.md, M5): featureCounts on its
// alignments of each pooled RNA-seq library, its genes' counts and TPM, genes
// by libraries (plexembler expression); the libraries' coverage, by strand,
// and the nanopore reads' depth on its genome, as bigWig (deepTools).
//

include { SUBREAD_FEATURECOUNTS                         } from '../../../modules/nf-core/subread/featurecounts/main'
include { PLEXEMBLER_EXPRESSION                         } from '../../../modules/local/plexembler/expression/main'
include { DEEPTOOLS_BAMCOVERAGE as BAMCOVERAGE_FORWARD  } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { DEEPTOOLS_BAMCOVERAGE as BAMCOVERAGE_REVERSE  } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { MINIMAP2_ALIGN as MINIMAP2_ONT_GENOME         } from '../../../modules/nf-core/minimap2/align/main'
include { DEEPTOOLS_BAMCOVERAGE as BAMCOVERAGE_ONT      } from '../../../modules/nf-core/deeptools/bamcoverage/main'

workflow EXPRESSION {

    take:
    ch_bams      // channel: [ val(meta), path(bam), path(bai) ], each species' alignments of each RNA-seq library (meta.id: the species)
    ch_gtf       // channel: [ val(meta), path(gtf) ], each species' genes, a species without any left out
    ch_genomes   // channel: [ val(meta), path(fasta) ], each species' genome
    ch_ont       // channel: [ val(meta), path(fastq) ], each species' nanopore reads
    strandedness // string: the RNA-seq libraries' strandedness: reverse, forward or unstranded

    main:

    //
    // MODULES: each species' genes counted in each library, then its tables
    //
    def ch_count = ch_bams
        .map { meta, bam, _bai -> [meta.id, bam, "${meta.library}=${meta.condition}_${meta.replicate}".toString()] }
        .groupTuple()
        .join(ch_gtf.map { meta, gtf -> [meta.id, meta, gtf] })
        .map { _id, bams, labels, meta, gtf ->
            [meta + [strandedness: strandedness, labels: labels, single_end: false], bams, gtf]
        }
    SUBREAD_FEATURECOUNTS(ch_count)
    PLEXEMBLER_EXPRESSION(SUBREAD_FEATURECOUNTS.out.counts)

    //
    // MODULES: the tracks: each library's coverage by strand; the nanopore
    // reads' depth
    //
    BAMCOVERAGE_FORWARD(ch_bams, [], [], [[:], []])
    BAMCOVERAGE_REVERSE(ch_bams, [], [], [[:], []])
    def ch_ont_align = ch_ont
        .map { meta, reads -> [meta.id, meta, reads] }
        .join(ch_genomes.map { meta, fasta -> [meta.id, fasta] })
        .multiMap { _id, meta, reads, fasta ->
            reads: [meta, reads]
            fasta: [meta, fasta]
        }
    MINIMAP2_ONT_GENOME(ch_ont_align.reads, ch_ont_align.fasta, true, 'bai', false, false)
    BAMCOVERAGE_ONT(MINIMAP2_ONT_GENOME.out.bam.join(MINIMAP2_ONT_GENOME.out.index), [], [], [[:], []])

    def ch_multiqc_files = SUBREAD_FEATURECOUNTS.out.summary.map { _meta, summary -> summary }
        .mix(PLEXEMBLER_EXPRESSION.out.multiqc.map { _meta, table -> table })

    emit:
    counts        = PLEXEMBLER_EXPRESSION.out.counts   // channel: [ val(meta), path(tsv) ]
    tpm           = PLEXEMBLER_EXPRESSION.out.tpm      // channel: [ val(meta), path(tsv) ]
    rnaseq_tracks = BAMCOVERAGE_FORWARD.out.bigwig.mix(BAMCOVERAGE_REVERSE.out.bigwig) // channel: [ val(meta), path(bigwig) ]
    ont_tracks    = BAMCOVERAGE_ONT.out.bigwig         // channel: [ val(meta), path(bigwig) ]
    multiqc_files = ch_multiqc_files                   // channel: path
}
