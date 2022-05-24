process BEDTOOLS_GENOMECOV {
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path("*.forward.bedGraph"), emit: bedgraph_forward
    tuple val(meta), path("*.reverse.bedGraph"), emit: bedgraph_reverse
    path "versions.yml"                        , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def prefix_forward = "${prefix}.forward"
    def prefix_reverse = "${prefix}.reverse"
    if (meta.strandedness == 'reverse') {
        prefix_forward = "${prefix}.reverse"
        prefix_reverse = "${prefix}.forward"
    }
    """
    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand + \\
        $args \\
        | bedtools sort > ${prefix_forward}.bedGraph
    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand - \\
        $args \\
        | bedtools sort > ${prefix_reverse}.bedGraph
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
