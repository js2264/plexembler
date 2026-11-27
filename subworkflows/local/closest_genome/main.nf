//
// Each species' closest known genome (docs/plan.md, M6): its genome's k-mers
// (sourmash sketch) matched against a database of genomes (sourmash gather),
// the best match named and, the network allowing, its sequence downloaded from
// NCBI (plexembler closest); then the species' genome aligned to it (minimap2)
// and drawn against it (plexembler dotplot).
//
// The download is the only step of the pipeline that needs the network. A
// species whose closest genome's sequence cannot be had (no network, an
// accession NCBI does not know, no match at all) has no dotplot: the run warns,
// and goes on.
//

include { SOURMASH_SKETCH                    } from '../../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_GATHER                    } from '../../../modules/nf-core/sourmash/gather/main'
include { PLEXEMBLER_CLOSEST                 } from '../../../modules/local/plexembler/closest/main'
include { MINIMAP2_ALIGN as MINIMAP2_CLOSEST } from '../../../modules/nf-core/minimap2/align/main'
include { PLEXEMBLER_DOTPLOT                 } from '../../../modules/local/plexembler/dotplot/main'

workflow CLOSEST_GENOME {

    take:
    ch_genomes  // channel: [ val(meta), path(fasta) ], each species' genome
    sourmash_db // path: a sourmash database of genomes, at the k-mer size gather is given
    lineages    // path: the database's taxonomy, sourmash's lineages (CSV); [] without one

    main:

    //
    // MODULES: each genome's k-mers, matched against the database's genomes
    //
    SOURMASH_SKETCH(ch_genomes, true)
    SOURMASH_GATHER(SOURMASH_SKETCH.out.signatures, sourmash_db, false, false, false, false)

    //
    // MODULE: each species' closest genome, named, and its sequence downloaded
    //
    PLEXEMBLER_CLOSEST(SOURMASH_GATHER.out.result, lineages)

    // Its warnings, in the run's log: a species without its closest genome's
    // sequence has no dotplot, and the run says why.
    PLEXEMBLER_CLOSEST.out.log
        .flatMap { _meta, file -> file.readLines().findAll { line -> line.contains(' :: WARNING :: ') } }
        .subscribe { line -> log.warn(line.split(' :: WARNING :: ', 2)[1]) }

    //
    // MODULES: each genome aligned to its closest, when that one's sequence
    // came, and drawn against it
    //
    def ch_references = PLEXEMBLER_CLOSEST.out.fasta.map { meta, fasta -> [meta.id, fasta] }
    def ch_pairs = ch_genomes
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_references)
        .multiMap { _id, meta, genome, reference ->
            query: [meta, genome]
            target: [meta, reference]
        }
    MINIMAP2_CLOSEST(ch_pairs.query, ch_pairs.target, false, '', true, false)
    PLEXEMBLER_DOTPLOT(
        MINIMAP2_CLOSEST.out.paf
            .map { meta, paf -> [meta.id, meta, paf] }
            .join(ch_genomes.map { meta, fasta -> [meta.id, fasta] })
            .join(ch_references)
            .join(PLEXEMBLER_CLOSEST.out.table.map { meta, table -> [meta.id, table] })
            .map { _id, meta, paf, genome, reference, table -> [meta, paf, genome, reference, table] }
    )

    emit:
    table         = PLEXEMBLER_CLOSEST.out.table   // channel: [ val(meta), path(tsv) ], each species' closest genomes
    dotplot       = PLEXEMBLER_DOTPLOT.out.svg     // channel: [ val(meta), path(svg) ], the species with their closest genome's sequence
    figure        = PLEXEMBLER_DOTPLOT.out.figure  // channel: [ val(meta), path(json) ], the dotplot, for the report
    alignment     = PLEXEMBLER_DOTPLOT.out.table   // channel: [ val(meta), path(tsv) ], how the two genomes align
    multiqc_files = PLEXEMBLER_CLOSEST.out.multiqc.map { _meta, table -> table } // channel: path
}
