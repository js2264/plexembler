//
// The checks of a species' sequence, at every stage (docs/plan.md, M1): its
// contiguity (gfastats), its gene completeness (BUSCO), and its consensus
// quality and k-mer completeness (Merqury).
//
// The QV is counted against the pool's WGS k-mers: an assembly k-mer missing
// from the whole pool is an error of this species' sequence. Without a WGS
// library, against the species' own nanopore k-mers: a rougher QV. The
// completeness is always the species' own nanopore k-mers': the pool's short
// reads hold every species' k-mers.
//

include { GFASTATS                                } from '../../../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO                             } from '../../../modules/nf-core/busco/busco/main'
include { MERQURY_MERQURY as MERQURY_QV           } from '../../../modules/nf-core/merqury/merqury/main'
include { MERQURY_MERQURY as MERQURY_COMPLETENESS } from '../../../modules/nf-core/merqury/merqury/main'

workflow GENOME_QC {

    take:
    ch_sequences  // channel: [ val(meta), path(fasta) ], meta.id the species, meta.pool its pool
    ch_ont_db     // channel: [ val(meta), path(meryl) ], one per species
    ch_wgs_db     // channel: [ val(meta), path(meryl) ], one per pool with WGS reads, meta.id the pool
    busco_lineage // string: the BUSCO lineage of a species without one (auto_euk: BUSCO chooses)
    busco_db      // path: BUSCO's lineages, downloaded beforehand; [] to download them

    main:

    //
    // MODULE: contiguity
    //
    GFASTATS(ch_sequences, '', '', '', [[:], []], [[:], []], [[:], []], [[:], []])

    //
    // MODULE: gene completeness, each species with its own lineage
    //
    def ch_busco = ch_sequences.multiMap { meta, fasta ->
        fasta: [meta, fasta]
        lineage: meta.busco_lineage ?: busco_lineage
    }
    BUSCO_BUSCO(ch_busco.fasta, 'genome', ch_busco.lineage, busco_db, [], true)

    //
    // MODULES: Merqury, the QV against the pool's WGS k-mers when it has some;
    // the completeness against the species' nanopore k-mers
    //
    // Every pool's WGS k-mers, by pool: a pool has several species, and a
    // species several sequences, each counted against them.
    def ch_wgs_by_pool = ch_wgs_db
        .map { meta, db -> [meta.id, db] }
        .toList()
        .map { pairs -> pairs.collectEntries() }
    def ch_with_kmers = ch_sequences
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .combine(ch_ont_db.map { meta, db -> [meta.id, db] }, by: 0)
        .combine(ch_wgs_by_pool)
        .map { _id, meta, fasta, ont_db, wgs_by_pool -> [meta.pool, meta, fasta, ont_db, wgs_by_pool[meta.pool]] }
        .branch { _pool, _meta, _fasta, _ont_db, wgs_db ->
            wgs: wgs_db != null
            ont_only: true
        }
    MERQURY_QV(
        ch_with_kmers.wgs
            .map { _pool, meta, fasta, _ont_db, wgs_db -> [meta + [qv_kmers: 'wgs'], wgs_db, fasta] }
            .mix(ch_with_kmers.ont_only.map { _pool, meta, fasta, ont_db, _wgs_db -> [meta + [qv_kmers: 'ont'], ont_db, fasta] })
    )
    MERQURY_COMPLETENESS(
        ch_with_kmers.wgs.map { _pool, meta, fasta, ont_db, _wgs_db -> [meta + [qv_kmers: 'ont'], ont_db, fasta] }
    )
    // The completeness against the species' nanopore k-mers: from the second
    // run when the QV's was against the pool's WGS k-mers, else from the first.
    def ch_completeness = MERQURY_COMPLETENESS.out.stats
        .mix(MERQURY_QV.out.stats.filter { meta, _stats -> meta.qv_kmers == 'ont' })
        .map { meta, stats -> [meta.findAll { key, _value -> key != 'qv_kmers' }, stats] }
    def ch_qv = MERQURY_QV.out.assembly_qv
        .map { meta, qv -> [meta.findAll { key, _value -> key != 'qv_kmers' }, qv, meta.qv_kmers] }

    def ch_multiqc_files = BUSCO_BUSCO.out.short_summaries_txt.map { _meta, txt -> txt }
        .mix(GFASTATS.out.assembly_summary.map { _meta, summary -> summary })
        .flatten()

    emit:
    stats         = GFASTATS.out.assembly_summary      // channel: [ val(meta), path(summary) ]
    busco         = BUSCO_BUSCO.out.short_summaries_txt // channel: [ val(meta), path(txt) ]
    busco_json    = BUSCO_BUSCO.out.short_summaries_json // channel: [ val(meta), path(json) ]
    qv            = ch_qv                              // channel: [ val(meta), path(qv), val(kmers: wgs or ont) ]
    completeness  = ch_completeness                    // channel: [ val(meta), path(stats) ]
    multiqc_files = ch_multiqc_files                   // channel: path
}
