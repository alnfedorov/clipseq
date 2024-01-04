include { STAR_ALIGN as STAR_ALIGN_RNALIB } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_GENOME } from '../../modules/nf-core/star/align/main'
include { UMITOOLS_DEDUP } from '../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

include { FIX_BAM_FLAGS_TAGS } from '../../modules/local/fix_bam_flags_tags'
include { CROSS_LINKS } from '../../modules/local/xlink'
include { NMS } from '../../modules/local/nms'

workflow ALIGN_NMS_DEDUP {
    take:
    reads                 // channel: [ val(meta), path(reads)  ]
    rnalib_index           // value: path (STAR index)
    genome_index          // value: paht(STAR index)

    main:
    versions = Channel.empty()

    // Align to the RNA collection
    STAR_ALIGN_RNALIB(
        reads,
        rnalib_index.map { [ [:], it ] },
        Channel.value([ [:], file('NA')]),
        true,
        false,
        false
    )
    versions = versions.mix(STAR_ALIGN_RNALIB.out.versions.first())

    // Align to the genome
    STAR_ALIGN_GENOME(
        reads,
        genome_index.map { [ [:], it ] },
        Channel.value([ [:], file('NA')]),
        true,
        false,
        false
    )
    versions = versions.mix(STAR_ALIGN_GENOME.out.versions.first())

    // Select the most suitable alignment in each case
    nms = reads \
        .join(STAR_ALIGN_RNALIB.out.bam, failOnMismatch: true) \
        .join(STAR_ALIGN_GENOME.out.bam, failOnMismatch: true)
    NMS(nms)
    versions = versions.mix(NMS.out.versions.first())

    // Deduplicate libraries
    UMITOOLS_DEDUP(NMS.out.alignment, true)
    versions = versions.mix(UMITOOLS_DEDUP.out.versions.first())

    SAMTOOLS_INDEX ( UMITOOLS_DEDUP.out.bam )
    versions = versions.mix(SAMTOOLS_INDEX.out.versions.first())

    alignments = UMITOOLS_DEDUP.out.bam \
        .join(SAMTOOLS_INDEX.out.bai)
    
    // Fix BAM mapq/NH tags to properly reflect the real number of alignments (they are broken by NMS)
    FIX_BAM_FLAGS_TAGS(alignments)
    versions = versions.mix(FIX_BAM_FLAGS_TAGS.out.versions.first())

    // Generate crosslinks
    CROSS_LINKS(FIX_BAM_FLAGS_TAGS.out.alignment)
    versions = versions.mix(CROSS_LINKS.out.versions.first())

    emit:
    // Alignment
    rnalib_align_log = STAR_ALIGN_RNALIB.out.log_final
    genome_align_log = STAR_ALIGN_GENOME.out.log_final
    nms_log = NMS.out.summary

    // Deduplication
    dedup_log = UMITOOLS_DEDUP.out.log
    dedup_edit_distance = UMITOOLS_DEDUP.out.tsv_edit_distance
    dedup_per_umi = UMITOOLS_DEDUP.out.tsv_per_umi
    dedup_umi_per_position = UMITOOLS_DEDUP.out.tsv_umi_per_position

    // Final alignment
    alignments = FIX_BAM_FLAGS_TAGS.out.alignment

    // Cross-link sites
    xlink_bed = CROSS_LINKS.out.bed
    xlink_bedgraph = CROSS_LINKS.out.bedgraph

    versions = versions.ifEmpty(null) // channel: [ path(versions.yml) ]
}
