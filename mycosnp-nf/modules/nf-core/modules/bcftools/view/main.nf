process BCFTOOLS_VIEW {
    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::bcftools=1.14' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.14--h88f3f91_0' :
        'quay.io/biocontainers/bcftools:1.14--h88f3f91_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: true,mode: "${params.publish_dir_mode}",path: { "${params.outdir}/combined/finalfiltered" },pattern: "*{vcf.gz,vcf.gz.tbi}"
    input:
    tuple val(meta), path(vcf), path(index)
    path(regions)
    path(targets)
    path(samples)
    output:
    tuple val(meta), path("*.gz") , emit: vcf
    path "versions.yml"           , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = "-Oz"
    def prefix = "finalfiltered"
    def regions_file  = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file =  samples ? "--samples-file ${samples}" : ""
    """
    bcftools view \\
        --output ${prefix}.vcf.gz \\
        ${regions_file} \\
        ${targets_file} \\
        ${samples_file} \\
        $args \\
        --threads 12 \\
        ${vcf}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
