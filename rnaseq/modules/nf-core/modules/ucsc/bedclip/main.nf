def VERSION = '377' // Version information not provided by tool on CLI
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
process UCSC_BEDCLIP {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::ucsc-bedclip=377" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedclip:377--h0b8a92a_2' :
        'quay.io/biocontainers/ucsc-bedclip:377--h0b8a92a_2' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    cpus 16
    memory '30 GB'
    input:
    tuple val(meta), path(bedgraph)
    path  sizes
    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path "versions.yml"                , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedClip \\
        $bedgraph \\
        $sizes \\
        ${prefix}.bedGraph
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
