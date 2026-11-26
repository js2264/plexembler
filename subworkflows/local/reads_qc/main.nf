//
// The reads as they come: the nanopore runs of each species, and each pooled library
//

include { NANOPLOT } from '../../../modules/nf-core/nanoplot/main'
include { FASTQC   } from '../../../modules/nf-core/fastqc/main'

workflow READS_QC {

    take:
    ch_species   // channel: [ val(meta), [ path(ont) ] ], one item per species
    ch_libraries // channel: [ val(meta), [ path(fastq) ] ], one item per pooled library

    main:

    //
    // MODULE: the nanopore reads of a species, its runs together
    //
    NANOPLOT(ch_species)

    //
    // MODULE: the reads of a pooled library, its runs together
    //
    FASTQC(ch_libraries)

    // One file per item: FastQC gives a library's two archives together.
    ch_multiqc_files = NANOPLOT.out.txt.map { _meta, txt -> txt }
        .mix(FASTQC.out.zip.map { _meta, zip -> zip })
        .flatten()

    emit:
    nanoplot      = NANOPLOT.out.txt  // channel: [ val(meta), path(txt) ]
    fastqc_zip    = FASTQC.out.zip    // channel: [ val(meta), [ path(zip) ] ]
    multiqc_files = ch_multiqc_files  // channel: path
}
