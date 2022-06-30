process VCF_QC {
    tag "vcf-qc"
    container 'ubuntu:20.04'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: true,mode: "${params.publish_dir_mode}",path: { "${params.outdir}/combined/vcf-qc-report" },pattern: "*.{txt}" 
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
