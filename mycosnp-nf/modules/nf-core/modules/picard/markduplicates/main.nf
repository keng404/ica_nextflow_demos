process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: "${params.save_debug}",mode: "${params.publish_dir_mode}",path: { "${params.outdir}/samples/${meta.id}/picard_markduplicates" },pattern: "*.bam"
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.bai")        , optional:true, emit: bai
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    path  "versions.yml"                  , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = "REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT" 
    def prefix = "${meta.id}_markdups"
    def avail_mem = 16
    """
    picard \\
        -Xmx${avail_mem}g \\
        MarkDuplicates \\
        $args \\
        I=$bam \\
        O=${prefix}.bam \\
        M=${prefix}.MarkDuplicates.metrics.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
