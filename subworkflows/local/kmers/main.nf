//
// The k-mers the genome checks count against (docs/plan.md, M1): each species'
// own nanopore reads, and each pool's WGS reads; and the spectrum of each
// species' nanopore k-mers, which GenomeScope2 fits for its heterozygosity.
//

include { MERYL_COUNT as MERYL_COUNT_ONT } from '../../../modules/nf-core/meryl/count/main'
include { MERYL_COUNT as MERYL_COUNT_WGS } from '../../../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM                 } from '../../../modules/nf-core/meryl/unionsum/main'
include { MERYL_HISTOGRAM                } from '../../../modules/nf-core/meryl/histogram/main'
include { GENOMESCOPE2                   } from '../../../modules/nf-core/genomescope2/main'

workflow KMERS {

    take:
    ch_ont_reads // channel: [ val(meta), path(fastq) ], one item per species
    ch_wgs_reads // channel: [ val(meta), [ path(fastq) ] ], one item per pool with WGS reads, meta.id its pool
    kmer_size    // integer: k

    main:

    //
    // MODULE: each species' nanopore k-mers
    //
    MERYL_COUNT_ONT(ch_ont_reads, kmer_size)

    //
    // MODULE: each pool's WGS k-mers, every file of its WGS libraries together
    //
    MERYL_COUNT_WGS(ch_wgs_reads, kmer_size)
    MERYL_UNIONSUM(MERYL_COUNT_WGS.out.meryl_db, kmer_size)

    //
    // MODULES: the spectrum of each species' nanopore k-mers, fitted by GenomeScope2
    //
    MERYL_HISTOGRAM(MERYL_COUNT_ONT.out.meryl_db, kmer_size)
    GENOMESCOPE2(MERYL_HISTOGRAM.out.hist)

    emit:
    ont_db       = MERYL_COUNT_ONT.out.meryl_db // channel: [ val(meta), path(meryl) ], one per species
    wgs_db       = MERYL_UNIONSUM.out.meryl_db  // channel: [ val(meta), path(meryl) ], one per pool with WGS reads
    histogram    = MERYL_HISTOGRAM.out.hist     // channel: [ val(meta), path(hist) ]
    genomescope  = GENOMESCOPE2.out.summary     // channel: [ val(meta), path(txt) ]
}
