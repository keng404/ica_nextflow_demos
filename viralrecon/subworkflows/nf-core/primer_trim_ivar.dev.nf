params = [:]
process = [:]
process.ext = [:]
params.protocol = null
params.ivar_trim_noprimer = false
params.ivar_trim_offset = null
//
// iVar trim, sort, index BAM file and run samtools stats, flagstat and idxstats
//
include { IVAR_TRIM         } from '../../modules/nf-core/modules/ivar/trim/main'
include { BAM_SORT_SAMTOOLS } from './bam_sort_samtools'
workflow PRIMER_TRIM_IVAR {
    take:
    bam // channel: [ val(meta), [ bam ], [bai] ]
    bed // path   : bed
    main:
    ch_versions = Channel.empty()
    //
    // iVar trim primers
    //
if (!params.skip_ivar_trim && params.protocol == 'amplicon') {
    process.ext.args = [ '-m 30 -q 20', params.ivar_trim_noprimer ? '' : '-e', params.ivar_trim_offset ? "-x ${params.ivar_trim_offset}" : '' ].join(' ').trim()
    process.ext.prefix =  "$meta.id.ivar_trim" 
}
    IVAR_TRIM (
        bam,
        bed
    )
    ch_versions = ch_versions.mix(IVAR_TRIM.out.versions.first())
    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_SAMTOOLS (
        IVAR_TRIM.out.bam
    )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)
    emit:
    bam_orig = IVAR_TRIM.out.bam              // channel: [ val(meta), bam   ]
    log_out  = IVAR_TRIM.out.log              // channel: [ val(meta), log   ]
    bam      = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    versions = ch_versions                    // channel: [ versions.yml ]
}
