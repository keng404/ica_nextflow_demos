process TABIX_BGZIP {
    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::tabix=1.11' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    input:
    tuple val(meta), path(input)
    output:
    tuple val(meta), path("*.gz"), emit: gz
    path  "versions.yml"         , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bgzip -c $args $input > ${prefix}.${input.getExtension()}.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
