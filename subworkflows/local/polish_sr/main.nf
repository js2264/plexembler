//
// Each species' contigs polished with its own short reads (docs/plan.md, M2):
// its WGS pairs; or, in a pool without WGS (or with --polish_reads hic), its
// Hi-C reads cut at their ligation junctions (plexembler digest), as single
// reads. Polypolish on every alignment of each read (bwa mem -a), then pypolca.
// The guard: the pairs plexembler split held out, which the polishing never
// saw. Where a polishing changed the draft, the draft's k-mers they lack that
// it removed, and those it made that they lack (MERYL_CHANGES); the polishing
// kept is the one that removed the most more than it made, if one did
// (plexembler select). Merqury's QV of each, against the same pairs, reported.
//

include { PLEXEMBLER_DIGEST                        } from '../../../modules/local/plexembler/digest/main'
include { CAT_FASTQ as CAT_FASTQ_POLISH            } from '../../../modules/nf-core/cat/fastq/main'
include { BWA_INDEX                                } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM                                  } from '../../../modules/nf-core/bwa/mem/main'
include { FIND_CONCATENATE as CAT_DRAFT_PLAIN      } from '../../../modules/nf-core/find/concatenate/main'
include { POLYPOLISH_POLISH                        } from '../../../modules/nf-core/polypolish/polish/main'
include { PYPOLCA                                  } from '../../../modules/local/pypolca/main'
include { MERYL_COUNT as MERYL_COUNT_HOLDOUT       } from '../../../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM as MERYL_UNIONSUM_HOLDOUT } from '../../../modules/nf-core/meryl/unionsum/main'
include { MERQURY_MERQURY as MERQURY_POLISH        } from '../../../modules/nf-core/merqury/merqury/main'
include { MERYL_CHANGES                            } from '../../../modules/local/meryl/changes/main'
include { PLEXEMBLER_SELECT as SELECT_POLISH       } from '../../../modules/local/plexembler/select/main'

workflow POLISH_SR {

    take:
    ch_contigs    // channel: [ val(meta), path(fasta) ], each species' draft
    ch_reads      // channel: [ val(meta), [ path(R1), path(R2) ] ], each species' split reads (POOL_DECONVOLVE.out.reads)
    polish_reads  // string: auto, wgs or hic
    kmer_size     // integer: k, for the QV guard

    main:

    //
    // The reads each species is polished with, and the pairs held out from them
    //
    def ch_by_type = ch_reads
        .map { meta, reads -> [meta.id, meta, reads] }
        .groupTuple()
        .flatMap { species, metas, reads ->
            def items = [metas, reads].transpose()
            def has_wgs = items.any { meta, _r -> meta.type == 'wgs' }
            def has_hic = items.any { meta, _r -> meta.type == 'hic' && meta.enzymes && meta.enzymes != 'none' }
            def chosen = polish_reads == 'wgs' || (polish_reads == 'auto' && has_wgs) ? 'wgs'
                : (has_hic ? 'hic' : null)
            if (chosen == null) {
                return []
            }
            items
                .findAll { meta, _r -> meta.type == chosen }
                .collect { meta, files -> [species, meta, files] }
        }
        .branch { _species, meta, _files ->
            holdout: meta.holdout
            wgs: meta.type == 'wgs'
            hic: true
        }

    //
    // MODULE: Hi-C reads cut at their junctions, every piece a single read
    //
    PLEXEMBLER_DIGEST(ch_by_type.hic.map { _species, meta, files -> [meta, files] })

    // One item per species: its polishing reads, as a pair (WGS) or single reads (Hi-C)
    def ch_polish_reads = ch_by_type.wgs
        .map { _species, meta, files -> [meta, files] }
        .mix(PLEXEMBLER_DIGEST.out.reads.map { meta, pieces -> [meta + [single_end: true], [pieces]] })
        .map { meta, files -> [meta.id, meta, files] }
        .groupTuple()
        .map { species, metas, files ->
            def pieces = metas[0].type == 'hic'
            [[id: species, pool: metas[0].pool, polish_reads: metas[0].type, single_end: pieces], files.flatten()]
        }
        .branch { _meta, files ->
            merge: files.size() > (_meta.single_end ? 1 : 2)
            as_is: true
        }
    CAT_FASTQ_POLISH(ch_polish_reads.merge)
    def ch_polishing = ch_polish_reads.as_is.mix(CAT_FASTQ_POLISH.out.reads.map { meta, reads -> [meta, [reads].flatten()] })

    //
    // MODULES: Polypolish on every alignment of each read, then pypolca
    //
    def ch_drafts = ch_contigs
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_polishing.map { meta, reads -> [meta.id, meta.polish_reads, reads] })
        .map { _id, meta, fasta, kind, reads -> [meta + [polish_reads: kind], fasta, reads] }
    BWA_INDEX(ch_drafts.map { meta, fasta, _reads -> [meta, fasta] }, 'bwa')
    def ch_alignments = ch_drafts
        .map { meta, fasta, reads -> [meta.id, meta, fasta, reads] }
        .join(BWA_INDEX.out.index.map { meta, index -> [meta.id, index] })
        .flatMap { _id, meta, fasta, reads, index ->
            reads.withIndex().collect { file, n -> [meta + [read: n + 1, single_end: true], file, index, fasta] }
        }
        .multiMap { meta, file, index, fasta ->
            reads: [meta, file]
            index: [meta, index]
            fasta: [meta, fasta]
        }
    BWA_MEM(ch_alignments.reads, ch_alignments.index, ch_alignments.fasta, false)
    def ch_sams = BWA_MEM.out.sam
        .map { meta, sam -> [meta.id, sam] }
        .groupTuple()
    // Polypolish reads a gzipped FASTA's first member alone: a bgzipped draft
    // (the screen's, a member per 64 kb) goes to it plain.
    CAT_DRAFT_PLAIN(ch_drafts.map { meta, fasta, _reads -> [meta, [fasta]] })
    def ch_polypolish = CAT_DRAFT_PLAIN.out.file_out
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_sams)
        .multiMap { _id, meta, fasta, sams ->
            fasta: [meta, fasta]
            sams: [meta, sams]
        }
    POLYPOLISH_POLISH(ch_polypolish.fasta, ch_polypolish.sams, false)
    PYPOLCA(
        POLYPOLISH_POLISH.out.fasta
            .map { meta, fasta -> [meta.id, meta, fasta] }
            .join(ch_drafts.map { meta, _fasta, reads -> [meta.id, reads] })
            .map { _id, meta, fasta, reads -> [meta, fasta, reads] }
    )

    //
    // MODULES: the guard, the QV of the draft and of each polishing against the held-out pairs
    //
    MERYL_COUNT_HOLDOUT(
        ch_by_type.holdout
            .map { species, meta, files -> [species, meta, files] }
            .groupTuple()
            .map { species, metas, files -> [[id: species, pool: metas[0].pool], files.flatten()] },
        kmer_size,
    )
    MERYL_UNIONSUM_HOLDOUT(MERYL_COUNT_HOLDOUT.out.meryl_db, kmer_size)
    def ch_candidates = ch_drafts.map { meta, fasta, _reads -> [meta + [candidate: 'draft'], fasta] }
        .mix(POLYPOLISH_POLISH.out.fasta.map { meta, fasta -> [meta + [candidate: 'polypolish'], fasta] })
        .mix(PYPOLCA.out.fasta.map { meta, fasta -> [meta + [candidate: 'pypolca'], fasta] })
    MERQURY_POLISH(
        ch_candidates
            .map { meta, fasta -> [meta.id, meta, fasta] }
            .combine(MERYL_UNIONSUM_HOLDOUT.out.meryl_db.map { meta, db -> [meta.id, db] }, by: 0)
            .map { _id, meta, fasta, db -> [meta, db, fasta] }
    )

    //
    // MODULE: the k-mers each polishing changed, that the held-out pairs lack
    //
    MERYL_CHANGES(
        ch_candidates
            .filter { meta, _fasta -> meta.candidate != 'draft' }
            .map { meta, fasta -> [meta.id, meta, fasta] }
            .combine(CAT_DRAFT_PLAIN.out.file_out.map { meta, draft -> [meta.id, draft] }, by: 0)
            .combine(MERYL_UNIONSUM_HOLDOUT.out.meryl_db.map { meta, db -> [meta.id, db] }, by: 0)
            .map { _id, meta, fasta, draft, db -> [meta, draft, fasta, db] },
        kmer_size,
    )

    //
    // MODULE: the polishing kept, if the held-out pairs favour one
    //
    def no_file = file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    def order = ['draft', 'polypolish', 'pypolca']
    def ch_selection = ch_candidates
        .map { meta, fasta -> [[meta.id, meta.candidate], meta, fasta] }
        .join(MERQURY_POLISH.out.assembly_qv.map { meta, qv -> [[meta.id, meta.candidate], qv] }, remainder: true)
        .join(MERYL_CHANGES.out.kmers.map { meta, kmers -> [[meta.id, meta.candidate], kmers] }, remainder: true)
        .filter { _key, meta, _fasta, _qv, _kmers -> meta != null }
        .map { _key, meta, fasta, qv, kmers ->
            def species = meta.findAll { key, _value -> key != 'candidate' }
            [groupKey(species, order.size()), meta.candidate, fasta, qv ?: no_file, kmers ?: no_file]
        }
        .groupTuple()
        .map { species, names, fastas, qvs, kmers ->
            def rows = [names, fastas, qvs, kmers].transpose().sort { row -> order.indexOf(row[0]) }
            [
                species.getGroupTarget(),
                rows.collect { row -> row[0] },
                rows.collect { row -> row[1] },
                rows.collect { _row -> no_file },
                rows.collect { row -> row[2] },
                rows.collect { _row -> no_file },
                rows.collect { row -> row[3] },
            ]
        }
    SELECT_POLISH(ch_selection, 'polish')

    // Every species: polished when a polishing raised its QV; as it came when
    // its pool has no reads to polish it with.
    def ch_polished = ch_contigs
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(SELECT_POLISH.out.fasta.map { meta, fasta -> [meta.id, fasta] }, remainder: true)
        .filter { _id, meta, _draft, _polished -> meta != null }
        .map { _id, meta, draft, polished -> [meta, polished ?: draft] }

    emit:
    contigs       = ch_polished                 // channel: [ val(meta), path(fasta) ], each species', meta as it came
    selection     = SELECT_POLISH.out.selection // channel: [ val(meta), path(tsv) ]
    digest        = PLEXEMBLER_DIGEST.out.json  // channel: [ val(meta), path(json) ]
    multiqc_files = SELECT_POLISH.out.multiqc.map { _meta, table -> table }
        .mix(PLEXEMBLER_DIGEST.out.multiqc.map { _meta, table -> table }) // channel: path
}
