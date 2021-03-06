//
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
// Alignment with STAR
//
include { STAR_ALIGN        } from '../../modules/local/star_align'
include { BAM_SORT_SAMTOOLS } from '../nf-core/bam_sort_samtools'
workflow ALIGN_STAR {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/star/index/
    gtf   // channel: /path/to/genome.gtf
    main:
    ch_versions = Channel.empty()
    //
    // Map reads with STAR
    //
params.outdir_custom = "${params.outdir}/align/star"
params.outdir_custom = "${params.outdir}/align/star"
params.outdir_custom = "${params.outdir}/align/star"
params.outdir_custom = "${params.outdir}/align/star"
    STAR_ALIGN ( reads, index, gtf )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
params.outdir_custom = "${params.outdir}/samtools/sort/bam"
params.outdir_custom = "${params.outdir}/samtools/sort/bam"
params.outdir_custom = "${params.outdir}/samtools/sort/bam"
params.outdir_custom = "${params.outdir}/samtools/sort/bam"
    BAM_SORT_SAMTOOLS ( STAR_ALIGN.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)
    emit:
    orig_bam       = STAR_ALIGN.out.bam             // channel: [ val(meta), bam            ]
    log_final      = STAR_ALIGN.out.log_final       // channel: [ val(meta), log_final      ]
    log_out        = STAR_ALIGN.out.log_out         // channel: [ val(meta), log_out        ]
    log_progress   = STAR_ALIGN.out.log_progress    // channel: [ val(meta), log_progress   ]
    bam_sorted     = STAR_ALIGN.out.bam_sorted      // channel: [ val(meta), bam_sorted     ]
    bam_transcript = STAR_ALIGN.out.bam_transcript  // channel: [ val(meta), bam_transcript ]
    fastq          = STAR_ALIGN.out.fastq           // channel: [ val(meta), fastq          ]
    tab            = STAR_ALIGN.out.tab             // channel: [ val(meta), tab            ]
    bam            = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai            = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi            = BAM_SORT_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats          = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    versions       = ch_versions                    // channel: [ versions.yml ]
}
