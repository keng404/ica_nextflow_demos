process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir enabled: true,mode: "${params.publish_dir_mode}",path: { "${params.outdir}/samples/${meta.id}/variant_calling/haplotypecaller"},pattern: "*{vcf.gz,vcf.gz.tbi}"
    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path fasta
    path fai
    path dict
    path dbsnp
    path dbsnp_tbi
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path "versions.yml"              , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = "-ERC GVCF --sample-ploidy \"${params.sample_ploidy}"
    def prefix = "${meta.id}.g"
    def interval_option = intervals ? "-L ${intervals}" : ""
    def dbsnp_option    = dbsnp ? "-D ${dbsnp}" : ""
    def avail_mem       = 16
    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        HaplotypeCaller \\
        -R $fasta \\
        -I $input \\
        ${dbsnp_option} \\
        ${interval_option} \\
        -O ${prefix}.vcf.gz \\
        $args \\
        --tmp-dir .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
