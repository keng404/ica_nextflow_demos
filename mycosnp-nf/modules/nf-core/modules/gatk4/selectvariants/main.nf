process GATK4_SELECTVARIANTS {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: true,mode: "${params.publish_dir_mode}",path: { "${params.outdir}/combined/selectedsnps" },pattern: "*{vcf.gz,vcf.gz.tbi}"
    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    output:
    tuple val(meta), path("*.selectvariants.vcf.gz")       , emit: vcf
    tuple val(meta), path("*.selectvariants.vcf.gz.tbi")   , emit: tbi
    path "versions.yml" , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = '--select-type-to-include "SNP"'
    def prefix = "combined_genotype_filtered_snps"
    def avail_mem = 16
    """
    gatk --java-options "-Xmx${avail_mem}G" SelectVariants \\
        -V $vcf \\
        -O ${prefix}.selectvariants.vcf.gz \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
