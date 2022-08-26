process VCF_QC {
	publishDir  path: { "${params.outdir}/qcvcf"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "vcf-qc"
    container 'library/ubuntu:20.04'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-small'
    
cpus 1
    
memory '6 GB'
    errorStrategy 'ignore'
    time '1day'
    maxForks 10
    input:
    path(vcffasta)
    output:
    path("vcf-qc-report.txt"),   emit: vcf_qc_report
script:
"""
    printf "Sample Name\\tLength\\tNumber-N\\n" > vcf-qc-report.txt
    awk '\$0 ~ ">" {if (NR > 1) {print c "\\t" d;} c=0;d=0;printf substr(\$0,2,200) "\\t"; } \$0 !~ ">" {c+=length(\$0);d+=gsub(/N/, "");d+=gsub(/n/, "")} END { print c "\\t" d; }' $vcffasta >> vcf-qc-report.txt
"""
}
