process SPLIT_VCF {
    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::bcftools=1.14' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.14--h88f3f91_0' :
        'quay.io/biocontainers/bcftools:1.14--h88f3f91_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: true,mode: "${params.publish_dir_mode}",path: { "${params.outdir}/combined/splitvcf" },pattern: "*{txt,vcf.gz}"
    input:
    tuple val(meta), path(vcf), path(tbi)
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcfs
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # First get a list of samples
    bcftools query \\
        -l \\
        $vcf > samplelist.txt
    for SAMPLE in \$(cat samplelist.txt); do
        bcftools view -Oz -s \$SAMPLE -o \$SAMPLE.vcf.gz $vcf
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}