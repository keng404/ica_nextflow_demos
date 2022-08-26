process UNTAR {
    tag "$archive"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    cpus 6    
    memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"
    input:
    tuple val(meta), path(archive)
    output:
    tuple val(meta), path("$untar"), emit: untar
    path "versions.yml"            , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    untar     = archive.toString().tokenize('.')[0]
    """
    mkdir output
    tar \\
        -C output --strip-components 1 \\
        -xzvf \\
        $args \\
        $archive \\
        $args2
    mv output ${untar}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
