process QC_REPORTSHEET {
	publishDir  path: { "${params.outdir}/reportsheetqc"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    container 'library/ubuntu:20.04'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-small'
    
cpus 1
    
memory '6 GB'
    errorStrategy 'ignore'
    time '1day'
    maxForks 10
    input:
    path(qc_lines)
    output:
    path("qc_report.txt"), emit: qc_reportsheet
    script:
    """
    printf \"Sample Name\\tReads Before Trimming\\tGC Before Trimming\\tAverage Q Score Before Trimming\\tReference Length Coverage Before Trimming\\tReads After Trimming\\tPaired Reads After Trimming\\tUnpaired Reads After Trimming\\tGC After Trimming\\tAverage Q Score After Trimming\\tReference Length Coverage After Trimming\\tMean Coverage Depth\\tReads Mapped\\n\" > qc_report.txt
    sort ${qc_lines} > sorted.txt
    cat sorted.txt >> qc_report.txt
    """
}
