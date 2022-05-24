//
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
// Run SAMtools stats, flagstat and idxstats
//
include { SAMTOOLS_STATS    } from '../../modules/nf-core/modules/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/modules/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/modules/samtools/flagstat/main'
workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai // channel: [ val(meta), [ bam ], [bai/csi] ]
    main:
    ch_versions = Channel.empty()
params.outdir_custom = "${params.outdir}/stats/samtools"
params.outdir_custom = "${params.outdir}/stats/samtools"
params.outdir_custom = "${params.outdir}/stats/samtools"
params.outdir_custom = "${params.outdir}/stats/samtools"
    SAMTOOLS_STATS ( ch_bam_bai, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
params.outdir_custom = "${params.outdir}/flagstat/samtools"
params.outdir_custom = "${params.outdir}/flagstat/samtools"
params.outdir_custom = "${params.outdir}/flagstat/samtools"
params.outdir_custom = "${params.outdir}/flagstat/samtools"
    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())
params.outdir_custom = "${params.outdir}/idxstats/samtools"
params.outdir_custom = "${params.outdir}/idxstats/samtools"
params.outdir_custom = "${params.outdir}/idxstats/samtools"
params.outdir_custom = "${params.outdir}/idxstats/samtools"
    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())
    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    versions = ch_versions                    // channel: [ versions.yml ]
}
