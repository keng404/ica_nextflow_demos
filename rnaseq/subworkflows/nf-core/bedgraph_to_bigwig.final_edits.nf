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
// Clip over-running ends from bedGraph file and convert to bigWig
//
include { UCSC_BEDCLIP          } from '../../modules/nf-core/modules/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'
workflow BEDGRAPH_TO_BIGWIG {
    take:
    bedgraph // channel: [ val(meta), [ bedgraph ] ]
    sizes    //    path: chrom.sizes
    main:
    ch_versions = Channel.empty()
    //
    // Clip bedGraph file
    //
params.outdir_custom = "${params.outdir}/bedclip/ucsc"
params.outdir_custom = "${params.outdir}/bedclip/ucsc"
params.outdir_custom = "${params.outdir}/bedclip/ucsc"
params.outdir_custom = "${params.outdir}/bedclip/ucsc"
    UCSC_BEDCLIP ( bedgraph, sizes )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions.first())
    //
    // Convert bedGraph to bigWig
    //
params.outdir_custom = "${params.outdir}/bedgraphtobigwig/ucsc"
params.outdir_custom = "${params.outdir}/bedgraphtobigwig/ucsc"
params.outdir_custom = "${params.outdir}/bedgraphtobigwig/ucsc"
params.outdir_custom = "${params.outdir}/bedgraphtobigwig/ucsc"
    UCSC_BEDGRAPHTOBIGWIG ( UCSC_BEDCLIP.out.bedgraph, sizes )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())
    emit:
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    bedgraph = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), [ bedgraph ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
