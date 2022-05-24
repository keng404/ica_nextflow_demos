//
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
// Read QC, UMI extraction and trimming
//
include { FASTQC           } from '../../modules/nf-core/modules/fastqc/main'
include { UMITOOLS_EXTRACT } from '../../modules/nf-core/modules/umitools/extract/main'
include { TRIMGALORE       } from '../../modules/nf-core/modules/trimgalore/main'
workflow FASTQC_UMITOOLS_TRIMGALORE {
    take:
    reads            // channel: [ val(meta), [ reads ] ]
    skip_fastqc      // boolean: true/false
    with_umi         // boolean: true/false
    skip_trimming    // boolean: true/false
    umi_discard_read // integer: 0, 1 or 2
    main:
    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()
    if (!skip_fastqc) {
params.outdir_custom = "${params.outdir}/fastqc"
params.outdir_custom = "${params.outdir}/fastqc"
params.outdir_custom = "${params.outdir}/fastqc"
params.outdir_custom = "${params.outdir}/fastqc"
        FASTQC ( reads ).html.set { fastqc_html }
        fastqc_zip  = FASTQC.out.zip
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }
    umi_reads = reads
    umi_log   = Channel.empty()
    if (with_umi) {
        UMITOOLS_EXTRACT ( reads ).reads.set { umi_reads }
        umi_log     = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())
        // Discard R1 / R2 if required
        if (umi_discard_read in [1,2]) {
            UMITOOLS_EXTRACT
                .out
                .reads
                .map { meta, reads ->
                    if (!meta.single_end) {
                        meta['single_end'] = true
                        reads = reads[umi_discard_read % 2]
                    }
                    return [ meta, reads ]
                }
                .set { umi_reads }
        }
    }
    trim_reads = umi_reads
    trim_html  = Channel.empty()
    trim_zip   = Channel.empty()
    trim_log   = Channel.empty()
    if (!skip_trimming) {
params.outdir_custom = "${params.outdir}/trimgalore"
params.outdir_custom = "${params.outdir}/trimgalore"
params.outdir_custom = "${params.outdir}/trimgalore"
params.outdir_custom = "${params.outdir}/trimgalore"
        TRIMGALORE ( umi_reads ).reads.set { trim_reads }
        trim_html   = TRIMGALORE.out.html
        trim_zip    = TRIMGALORE.out.zip
        trim_log    = TRIMGALORE.out.log
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    }
    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]
    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]
    umi_log            // channel: [ val(meta), [ log ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]
    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
