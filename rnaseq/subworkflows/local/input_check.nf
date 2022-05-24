//
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
// Check input samplesheet and get read channels
//
include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    main:
params.outdir_custom = "${params.outdir}/check/samplesheet"
params.outdir_custom = "${params.outdir}/check/samplesheet"
publishDir =  path: { "${params.outdir}/pipeline_info" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
params.outdir_custom = "${params.outdir}/check/samplesheet"
publishDir =  path: { "${params.outdir}/pipeline_info" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
params.outdir_custom = "${params.outdir}/check/samplesheet"
publishDir =  path: { "${params.outdir}/pipeline_info" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }
    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    meta.strandedness = row.strandedness
    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}
