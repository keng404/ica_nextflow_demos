def VERSION = '2.2.0' // Version information not provided by tool on CLI
process HISAT2_EXTRACTSPLICESITES {
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$gtf"
    conda (params.enable_conda ? 'bioconda::hisat2=2.2.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h1b792b2_3' :
        'quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    cpus 16
    memory '30 GB'
    input:
    path gtf
    output:
    path "*.splice_sites.txt", emit: txt
    path "versions.yml"      , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.splice_sites.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: $VERSION
    END_VERSIONS
    """
}
