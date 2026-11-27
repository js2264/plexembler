//
// Each species' annotation (docs/plan.md, M4): its repeats, de novo
// (RepeatModeler) and masked soft (RepeatMasker); its genes, by BRAKER3 with
// its RNA-seq alignments and the proteins given, or AUGUSTUS ab initio with
// a trained species when it has neither; named and written by
// `plexembler genes`, translated with its genetic code; BUSCO on its
// proteins; its tRNAs (tRNAscan-SE) and rRNAs (barrnap).
//

include { FIND_CONCATENATE as CAT_GENOME_PLAIN } from '../../../modules/nf-core/find/concatenate/main'
include { REPEATMODELER_BUILDDATABASE        } from '../../../modules/nf-core/repeatmodeler/builddatabase/main'
include { REPEATMODELER_REPEATMODELER        } from '../../../modules/nf-core/repeatmodeler/repeatmodeler/main'
include { REPEATMASKER_REPEATMASKER          } from '../../../modules/nf-core/repeatmasker/repeatmasker/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_RNASEQ } from '../../../modules/nf-core/samtools/merge/main'
include { BRAKER3                            } from '../../../modules/nf-core/braker3/main'
include { AUGUSTUS                           } from '../../../modules/local/augustus/main'
include { PLEXEMBLER_GENES                   } from '../../../modules/local/plexembler/genes/main'
include { BUSCO_BUSCO as BUSCO_PROTEINS      } from '../../../modules/nf-core/busco/busco/main'
include { TRNASCANSE                         } from '../../../modules/nf-core/trnascanse/main'
include { BARRNAP                            } from '../../../modules/nf-core/barrnap/main'

workflow ANNOTATION {

    take:
    ch_genomes       // channel: [ val(meta), path(fasta) ], each species' genome
    ch_bams          // channel: [ val(meta), path(bam), path(bai) ], each species' RNA-seq alignments, a library each (meta.id: the species)
    repeats          // string: repeatmodeler, or none
    predictor        // string: braker, augustus, or none
    augustus_species // string: AUGUSTUS's trained species, for ab initio predictions
    proteins         // path: proteins of related species for BRAKER (OrthoDB's), or []
    busco_lineage    // string: the lineage of a species whose samplesheet gives none
    busco_db         // path: BUSCO's lineages, downloaded, or []

    main:

    //
    // MODULES: the repeats, found de novo and masked soft
    //
    def ch_masked = channel.empty()
    def ch_multiqc_files = channel.empty()
    if (repeats == 'repeatmodeler') {
        REPEATMODELER_BUILDDATABASE(ch_genomes)
        REPEATMODELER_REPEATMODELER(REPEATMODELER_BUILDDATABASE.out.db)
        def ch_mask = ch_genomes
            .map { meta, fasta -> [meta.id, meta, fasta] }
            .join(REPEATMODELER_REPEATMODELER.out.fasta.map { meta, library -> [meta.id, library] })
            .multiMap { _id, meta, fasta, library ->
                fasta: [meta, fasta]
                library: library
            }
        REPEATMASKER_REPEATMASKER(ch_mask.fasta, ch_mask.library)
        ch_masked = REPEATMASKER_REPEATMASKER.out.masked
    } else {
        // Unmasked, plain: BRAKER reads its genome as text.
        CAT_GENOME_PLAIN(ch_genomes.map { meta, fasta -> [meta, [fasta]] })
        ch_masked = CAT_GENOME_PLAIN.out.file_out
    }

    //
    // MODULE: each species' RNA-seq alignments, its libraries merged
    //
    def ch_species_bams = ch_bams
        .map { meta, bam, bai -> [meta.id, bam, bai] }
        .groupTuple()
        .branch { _id, bams, _bais ->
            merge: bams.size() > 1
            as_is: true
        }
    SAMTOOLS_MERGE_RNASEQ(
        ch_species_bams.merge.map { id, bams, bais -> [[id: id], bams, bais] },
        [[:], [], [], []],
        'bai',
    )
    def ch_rnaseq = ch_species_bams.as_is
        .map { id, bams, _bais -> [id, bams[0]] }
        .mix(SAMTOOLS_MERGE_RNASEQ.out.bam.map { meta, bam -> [meta.id, bam] })

    //
    // MODULES: the genes. BRAKER3 for a species with RNA-seq or proteins to
    // go by; AUGUSTUS ab initio, with a trained species, for the others.
    //
    def ch_models = channel.empty()
    if (predictor != 'none') {
        def ch_predict = ch_masked
            .map { meta, fasta -> [meta.id, meta, fasta] }
            .join(ch_rnaseq, remainder: true)
            .filter { _id, meta, _fasta, _bam -> meta != null }
            .branch { _id, _meta, _fasta, bam ->
                braker: predictor == 'braker' && (bam != null || proteins)
                augustus: true
            }
        def ch_braker = ch_predict.braker
            .multiMap { _id, meta, fasta, bam ->
                fasta: [meta, fasta]
                bam: bam ?: []
            }
        BRAKER3(ch_braker.fasta, ch_braker.bam, [], [], proteins, [])
        AUGUSTUS(ch_predict.augustus.map { _id, meta, fasta, _bam -> [meta, fasta] }, augustus_species)
        ch_models = BRAKER3.out.gtf.map { meta, gtf -> [meta.id, gtf, 'BRAKER3'] }
            .mix(AUGUSTUS.out.gff3.map { meta, gff3 -> [meta.id, gff3, 'AUGUSTUS'] })
    }

    //
    // MODULE: the genes named, translated and written
    //
    def ch_annotate = ch_genomes
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_models)
        .multiMap { _id, meta, fasta, models, source ->
            models: [meta, fasta, models]
            source: source
        }
    PLEXEMBLER_GENES(ch_annotate.models, ch_annotate.source)
    ch_multiqc_files = ch_multiqc_files.mix(PLEXEMBLER_GENES.out.multiqc.map { _meta, table -> table })

    //
    // MODULE: BUSCO on each species' proteins
    //
    // A species without genes (a predictor found none) has no proteins to check.
    def ch_with_genes = PLEXEMBLER_GENES.out.table
        .filter { _meta, table -> table.readLines().size() > 1 }
        .map { meta, _table -> [meta.id, true] }
    def ch_busco = PLEXEMBLER_GENES.out.proteins
        .map { meta, faa -> [meta.id, meta, faa] }
        .join(ch_with_genes)
        .multiMap { _id, meta, faa, _genes ->
            fasta: [meta + [candidate: 'proteins'], faa]
            lineage: meta.busco_lineage ?: busco_lineage
        }
    BUSCO_PROTEINS(ch_busco.fasta, 'proteins', ch_busco.lineage, busco_db, [], true)
    ch_multiqc_files = ch_multiqc_files.mix(BUSCO_PROTEINS.out.short_summaries_txt.map { _meta, summary -> summary })

    //
    // MODULES: the non-coding RNAs: tRNAs, the rDNA locus
    //
    TRNASCANSE(ch_genomes, true, false, true, false)
    BARRNAP(ch_genomes.map { meta, fasta -> [meta, fasta, 'euk'] })

    emit:
    masked        = ch_masked                         // channel: [ val(meta), path(fasta) ], plain, soft-masked when repeats are masked
    gff3          = PLEXEMBLER_GENES.out.gff3         // channel: [ val(meta), path(gff3) ]
    gtf           = PLEXEMBLER_GENES.out.gtf          // channel: [ val(meta), path(gtf) ], every species, with genes or not
    proteins      = PLEXEMBLER_GENES.out.proteins     // channel: [ val(meta), path(faa) ]
    proteins_busco = BUSCO_PROTEINS.out.short_summaries_txt // channel: [ val(meta), path(txt) ]
    genes         = PLEXEMBLER_GENES.out.table        // channel: [ val(meta), path(tsv) ]
    with_genes    = ch_with_genes                     // channel: [ val(id), true ], the species with genes
    trnas         = TRNASCANSE.out.gff                // channel: [ val(meta), path(gff) ]
    rrnas         = BARRNAP.out.gff                   // channel: [ val(meta), path(gff) ]
    multiqc_files = ch_multiqc_files                  // channel: path
}
