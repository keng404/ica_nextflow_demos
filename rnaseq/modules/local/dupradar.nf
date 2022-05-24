process DUPRADAR {
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::bioconductor-dupradar=1.18.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dupradar:1.18.0--r40_1' :
        'quay.io/biocontainers/bioconductor-dupradar:1.18.0--r40_1' }"
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
    path  gtf
    output:
    tuple val(meta), path("*.pdf")    , emit: pdf
    tuple val(meta), path("*.txt")    , emit: txt
    tuple val(meta), path("*_mqc.txt"), emit: multiqc
    path "versions.yml"               , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    def paired_end = meta.single_end ? 'single' :  'paired'
    """
    dupradar.r \\
        $bam \\
        $prefix \\
        $gtf \\
        $strandedness \\
        $paired_end \\
        $task.cpus
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dupradar: \$(Rscript -e "library(dupRadar); cat(as.character(packageVersion('dupRadar')))")
    END_VERSIONS
    """
}
