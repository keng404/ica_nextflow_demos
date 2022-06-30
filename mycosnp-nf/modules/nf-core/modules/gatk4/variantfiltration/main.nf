process GATK4_VARIANTFILTRATION {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: true,mode: "${params.publish_dir_mode}",path: { "${params.outdir}/combined/filteredgvcfs" },pattern: "*{vcf.gz,vcf.gz.tbi}"
    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    path  fasta
    path  fai
    path  dict
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path "versions.yml", emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = "--filter-expression '${params.gvcfs_filter}' --filter-name filter"
    def prefix = "combined_genotype_filtered"
    def avail_mem = 16
    """
    gatk --java-options "-Xmx${avail_mem}G" VariantFiltration \\
        -R $fasta \\
        -V $vcf \\
        -O ${prefix}.vcf.gz \\
        $args
    cat <<-END_VERSIONS > versions.yml

    """
}
