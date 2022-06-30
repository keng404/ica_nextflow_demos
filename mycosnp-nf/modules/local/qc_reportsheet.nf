process QC_REPORTSHEET {
    container 'ubuntu:20.04'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-small'
    errorStrategy 'ignore'
    publishDir  enabled: true,mode: "${params.publish_dir_mode}",path: { "${params.outdir}/stats/qc_report" },pattern: "*",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    time '1day'
    input:
    path(qc_lines)
    output:
    path("qc_report.txt"), emit: qc_reportsheet
    script:
    """
    printf \"Sample Name\\t# Reads Before Trimming\\tGC Before Trimming\\tAverage Phred Before Trimming\\tCoverage Before Trimming\\t# Reads After Trimming\\t# Paired Reads After Trimming\\t# Unpaired Reads After Trimming\\tGC After Trimming\\tAverage Phred After Trimming\\tCoverage After Trimming\\n\" > qc_report.txt
    sort ${qc_lines} > sorted.txt
    cat sorted.txt >> qc_report.txt
    """
}
