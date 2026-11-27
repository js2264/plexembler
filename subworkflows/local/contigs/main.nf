//
// Each species' contigs (docs/plan.md, M1): every assembler's, named after the
// species and polished with its nanopore reads, as candidates; each one checked
// (GENOME_QC); one kept, by the checks or by the samplesheet's `assembler`; and
// the identity between the species of each pool.
//

include { FLYE              } from '../../../modules/nf-core/flye/main'
include { HIFIASM           } from '../../../modules/nf-core/hifiasm/main'
include { GFATOOLS_GFA2FA   } from '../../../modules/nf-core/gfatools/gfa2fa/main'
include { PLEXEMBLER_RENAME } from '../../../modules/local/plexembler/rename/main'
include { MEDAKA            } from '../../../modules/nf-core/medaka/main'
include { DORADO_POLISH     } from '../../../modules/local/dorado/polish/main'
include { GENOME_QC         } from '../genome_qc'
include { PLEXEMBLER_SELECT } from '../../../modules/local/plexembler/select/main'
include { SKANI_TRIANGLE    } from '../../../modules/nf-core/skani/triangle/main'
include { PLEXEMBLER_ANI    } from '../../../modules/local/plexembler/ani/main'

workflow CONTIGS {

    take:
    ch_reads      // channel: [ val(meta), path(fastq) ], one item per species (ONT_READS.out.reads)
    ch_bams       // channel: [ val(meta), [ path(bam) ] ], the species whose reads are Dorado's BAMs
    ch_ont_db     // channel: [ val(meta), path(meryl) ], one per species (KMERS.out.ont_db)
    ch_wgs_db     // channel: [ val(meta), path(meryl) ], one per pool with WGS reads (KMERS.out.wgs_db)
    assemblers    // list: the assemblers that run, each a candidate (flye, hifiasm)
    polisher      // string: medaka, dorado or none
    busco_lineage // string: the BUSCO lineage of a species without one
    busco_db      // path: BUSCO's lineages, downloaded beforehand; [] to download them

    main:

    //
    // MODULES: each assembler's contigs, a candidate; Flye's mode by chemistry
    //
    def ch_flye = ch_reads
        .filter { _meta, _reads -> 'flye' in assemblers }
        .multiMap { meta, reads ->
            reads: [meta + [candidate: 'flye'], reads]
            mode: meta.chemistry == 'r9' ? '--nano-raw' : '--nano-hq'
        }
    FLYE(ch_flye.reads, ch_flye.mode)

    HIFIASM(
        ch_reads
            .filter { _meta, _reads -> 'hifiasm' in assemblers }
            .map { meta, reads -> [meta + [candidate: 'hifiasm'], reads, []] },
        [[:], [], []],
        [[:], [], []],
        [[:], []],
    )
    GFATOOLS_GFA2FA(HIFIASM.out.primary_contigs)

    //
    // MODULE: every contig named after its species, <species>_<contig>; a
    // circular contig that the assembler wrote as copies of itself (Flye's
    // doubled mitochondrion) cut to one copy, by its sequence. Flye's table of
    // contigs says which are linear, and stay whole; hifiasm writes none.
    //
    PLEXEMBLER_RENAME(
        FLYE.out.fasta
            .join(FLYE.out.txt)
            .mix(GFATOOLS_GFA2FA.out.fasta.map { meta, fasta -> [meta, fasta, []] })
    )

    //
    // MODULES: the candidates polished with their species' nanopore reads:
    // Medaka; dorado polish for a species whose reads are Dorado's BAMs, when asked
    //
    def ch_polished = PLEXEMBLER_RENAME.out.fasta
    if (polisher != 'none') {
        def ch_bams_by_species = ch_bams
            .map { meta, bams -> [meta.id, bams] }
            .toList()
            .map { pairs -> pairs.collectEntries() }
        def ch_to_polish = PLEXEMBLER_RENAME.out.fasta
            .map { meta, contigs -> [meta.id, meta, contigs] }
            .combine(ch_reads.map { meta, reads -> [meta.id, reads] }, by: 0)
            .combine(ch_bams_by_species)
            .branch { _id, meta, _contigs, _reads, bams_by_species ->
                dorado: polisher == 'dorado' && bams_by_species[meta.id]
                medaka: true
            }
        DORADO_POLISH(
            ch_to_polish.dorado.map { _id, meta, contigs, _reads, bams_by_species -> [meta, bams_by_species[meta.id], contigs] }
        )
        MEDAKA(
            ch_to_polish.medaka.map { _id, meta, contigs, reads, _bams_by_species -> [meta, reads, contigs] }
        )
        ch_polished = MEDAKA.out.assembly.mix(DORADO_POLISH.out.fasta)
    }

    //
    // SUBWORKFLOW: each candidate's checks
    //
    GENOME_QC(ch_polished, ch_ont_db, ch_wgs_db, busco_lineage, busco_db)

    //
    // MODULE: one candidate kept per species, by its checks or by the samplesheet
    //
    // A check a candidate lacks (BUSCO finds no gene in a small test genome) is
    // a placeholder, which the table of candidates writes as NA.
    def no_file = file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    def ch_candidates = ch_polished
        .map { meta, fasta -> [candidateKey(meta), meta, fasta] }
        .join(GENOME_QC.out.busco.map { meta, txt -> [candidateKey(meta), [txt].flatten()[0]] }, remainder: true)
        .join(GENOME_QC.out.qv.map { meta, qv, _kmers -> [candidateKey(meta), qv] }, remainder: true)
        .join(GENOME_QC.out.stats.map { meta, stats -> [candidateKey(meta), stats] }, remainder: true)
        .filter { _key, meta, _fasta, _busco, _qv, _stats -> meta != null }
        .map { _key, meta, fasta, busco, qv, stats ->
            def species = meta.findAll { key, _value -> key != 'candidate' }
            [groupKey(species, assemblers.size()), meta.candidate, fasta, busco ?: no_file, qv ?: no_file, stats ?: no_file]
        }
        .groupTuple()
        .map { species, candidates, fastas, buscos, qvs, stats ->
            def meta = species.getGroupTarget()
            // No k-mer changes: those are the polishing's.
            [meta + [keep: meta.assembler], candidates, fastas, buscos, qvs, stats, candidates.collect { _name -> no_file }]
        }
    PLEXEMBLER_SELECT(ch_candidates, 'assembly')

    // The species' meta, as it came: the selection's key `keep` left out.
    def ch_contigs = PLEXEMBLER_SELECT.out.fasta
        .map { meta, fasta -> [meta.findAll { key, _value -> key != 'keep' }, fasta] }

    //
    // MODULES: the identity between the species of each pool of two or more
    //
    def ch_pools = ch_contigs
        .map { meta, fasta -> [meta.pool, meta.id, fasta] }
        .groupTuple()
        .filter { _pool, species, _fastas -> species.size() > 1 }
        .map { pool, species, fastas -> [[id: pool, species: species.sort()], fastas] }
    SKANI_TRIANGLE(ch_pools)
    PLEXEMBLER_ANI(SKANI_TRIANGLE.out.triangle)

    def ch_multiqc_files = GENOME_QC.out.multiqc_files
        .mix(PLEXEMBLER_SELECT.out.multiqc.map { _meta, table -> table })
        .mix(PLEXEMBLER_ANI.out.multiqc.map { _meta, table -> table })

    emit:
    contigs       = ch_contigs                       // channel: [ val(meta), path(fasta) ], the kept contigs of each species
    candidates    = ch_polished                      // channel: [ val(meta), path(fasta) ], every candidate, meta.candidate its name
    graphs        = FLYE.out.gfa.mix(HIFIASM.out.primary_contigs) // channel: [ val(meta), path(gfa) ]
    selection     = PLEXEMBLER_SELECT.out.selection  // channel: [ val(meta), path(tsv) ]
    busco         = GENOME_QC.out.busco              // channel: [ val(meta), path(txt) ]
    qv            = GENOME_QC.out.qv                 // channel: [ val(meta), path(qv), val(kmers) ]
    ani           = PLEXEMBLER_ANI.out.pairs         // channel: [ val(meta), path(tsv) ], one per pool
    multiqc_files = ch_multiqc_files                 // channel: path
}

//
// A candidate's key: its species and its name
//
def candidateKey(meta) {
    return [meta.id, meta.candidate]
}
