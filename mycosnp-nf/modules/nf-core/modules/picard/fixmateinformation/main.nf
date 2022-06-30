process PICARD_FIXMATEINFORMATION {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::picard=2.26.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.9--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.9--hdfd78af_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: "${params.save_debug}",mode: "${params.publish_dir_mode}",path: { "${params.outdir}/samples/${meta.id}/picard_fixmate" },pattern: "*.bam"
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_fixmate"
    def STRINGENCY = "LENIENT"
    def avail_mem = 16
    """
    picard \\
        FixMateInformation \\
        -Xmx${avail_mem}g \\
        -I ${bam} \\
        -O ${prefix}.bam \\
        --VALIDATION_STRINGENCY ${STRINGENCY}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard FixMateInformation --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
