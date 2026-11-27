//
// A species' nanopore reads, ready to assemble (docs/plan.md, M1): its runs as
// one FASTQ, Dorado's BAMs converted; their QC; the chemistry and the genome
// size the samplesheet left unknown, estimated; filtered by length and quality,
// and capped at a depth, the longest and best reads kept.
//

include { SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main'
include { CAT_FASTQ      } from '../../../modules/nf-core/cat/fastq/main'
include { NANOPLOT       } from '../../../modules/nf-core/nanoplot/main'
include { CHOPPER        } from '../../../modules/nf-core/chopper/main'
include { LRGE           } from '../../../modules/local/lrge/main'
include { FILTLONG       } from '../../../modules/nf-core/filtlong/main'

workflow ONT_READS {

    take:
    ch_species // channel: [ val(meta), [ path(ont) ] ], one item per species, its runs in the samplesheet's order
    max_depth  // integer: the depth the reads are capped at, from the genome size; 0 for no cap

    main:

    //
    // Each run on its own: Dorado's unaligned BAMs to FASTQ, their tags kept
    //
    def ch_runs = ch_species
        .flatMap { meta, files ->
            def runs = [files].flatten()
            runs.withIndex().collect { run, index -> [meta, index, runs.size(), run] }
        }
        .branch { _meta, _index, _n, run ->
            bam: run.name.endsWith('.bam')
            fastq: true
        }

    SAMTOOLS_FASTQ(
        ch_runs.bam.map { meta, index, _n, bam -> [meta + [id: "${meta.id}_run${index + 1}".toString(), species: meta.id], bam] },
        false,
    )
    // A read of an unaligned BAM is neither the first nor the second of a pair:
    // samtools writes it to its "other" file.
    def ch_converted = ch_runs.bam
        .map { meta, index, n, _bam -> ["${meta.id}_run${index + 1}".toString(), meta, index, n] }
        .join(SAMTOOLS_FASTQ.out.other.map { meta, fastq -> [meta.id, fastq] })
        .map { _run, meta, index, n, fastq -> [meta, index, n, fastq] }

    //
    // A species' runs as one FASTQ, gzipped: a single gzipped FASTQ as it is
    //
    def ch_grouped = ch_runs.fastq
        .mix(ch_converted)
        .map { meta, index, n, fastq -> [groupKey(meta.id, n), meta, index, fastq] }
        .groupTuple()
        .map { _key, metas, indexes, fastqs ->
            [metas[0], [indexes, fastqs].transpose().sort { pair -> pair[0] }.collect { pair -> pair[1] }]
        }
        .branch { _meta, fastqs ->
            as_is: fastqs.size() == 1 && fastqs[0].name.endsWith('.gz')
            merge: true
        }
    CAT_FASTQ(ch_grouped.merge)
    def ch_raw = ch_grouped.as_is
        .map { meta, fastqs -> [meta, fastqs[0]] }
        .mix(CAT_FASTQ.out.reads)

    //
    // MODULE: the reads as they come, and their median quality
    //
    NANOPLOT(ch_raw)

    // The chemistry: the samplesheet's, else told from the reads' median quality.
    def ch_profiled = ch_raw
        .join(NANOPLOT.out.txt.map { meta, txt -> [meta, [txt].flatten().find { stats -> stats.name.endsWith('NanoStats.txt') }] })
        .map { meta, reads, nanostats ->
            def median = medianReadQuality(nanostats)
            def profiled = meta + [
                chemistry           : meta.chemistry ?: chemistryFromQuality(median),
                chemistry_estimated : !meta.chemistry,
                median_read_quality : median,
            ]
            [profiled, reads]
        }

    //
    // MODULE: the reads filtered by length and quality
    //
    CHOPPER(ch_profiled, [])

    //
    // MODULE: the genome size the samplesheet left unknown, from the reads
    //
    def ch_filtered = CHOPPER.out.fastq.branch { meta, _reads ->
        known: meta.genome_size
        unknown: true
    }
    LRGE(ch_filtered.unknown)
    def ch_sized = ch_filtered.known
        .map { meta, reads -> [meta + [genome_size_bp: genomeSizeBases(meta.genome_size), genome_size_estimated: false], reads] }
        .mix(
            ch_filtered.unknown
                .join(LRGE.out.estimate)
                .map { meta, reads, estimate ->
                    [meta + [genome_size_bp: estimate.text.trim().toBigDecimal().toLong(), genome_size_estimated: true], reads]
                }
        )

    //
    // MODULE: the reads capped at max_depth times the genome size
    //
    def ch_to_cap = ch_sized.branch { _meta, _reads ->
        cap: max_depth > 0
        as_is: true
    }
    FILTLONG(ch_to_cap.cap.map { meta, reads -> [meta + [target_bases: meta.genome_size_bp * max_depth], [], reads] })
    def ch_reads = ch_to_cap.as_is.mix(FILTLONG.out.reads)

    def ch_multiqc_files = NANOPLOT.out.txt.map { _meta, txt -> txt }.flatten()

    emit:
    reads         = ch_reads              // channel: [ val(meta), path(fastq) ], meta with chemistry and genome_size_bp
    raw           = ch_raw                // channel: [ val(meta), path(fastq) ], the runs as one FASTQ
    nanostats     = NANOPLOT.out.txt      // channel: [ val(meta), path(txt) ]
    genome_size   = LRGE.out.log          // channel: [ val(meta), path(log) ], LRGE's estimate and interval
    filtlong_log  = FILTLONG.out.log      // channel: [ val(meta), path(log) ]
    multiqc_files = ch_multiqc_files      // channel: path
}

//
// A genome size as the samplesheet gives it (12m, 10.4m, 40k, 12000000), in bases
//
def genomeSizeBases(size) {
    def text = size.toString().trim().toLowerCase()
    def unit = ['k': 1e3, 'm': 1e6, 'g': 1e9][text[-1]] ?: 1
    def number = unit == 1 ? text : text[0..-2]
    return Math.round(number.toBigDecimal() * unit) as Long
}

//
// The median read quality of NanoPlot's NanoStats.txt, or null when it has none
//
def medianReadQuality(nanostats) {
    def line = nanostats.readLines().find { text -> text.startsWith('Median read quality:') }
    return line ? line.tokenize(':')[1].trim().replace(',', '').toBigDecimal() : null
}

//
// The nanopore chemistry told from the reads' median quality: R10.4.1, basecalled
// with Dorado's HAC or SUP models, reads at Q15 and above; R9.4.1 below. Unknown
// quality is taken as R10, the chemistry sold today.
//
def chemistryFromQuality(median) {
    return (median == null || median >= 15) ? 'r10' : 'r9'
}
