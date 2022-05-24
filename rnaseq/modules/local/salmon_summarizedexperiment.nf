process SALMON_SUMMARIZEDEXPERIMENT {
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$tx2gene"
    conda (params.enable_conda ? "bioconda::bioconductor-summarizedexperiment=1.20.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.20.0--r40_0' :
        'quay.io/biocontainers/bioconductor-summarizedexperiment:1.20.0--r40_0' }"
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
    path counts
    path tpm
    path tx2gene
    output:
    path "*.rds"       , emit: rds
    path "versions.yml", emit: versions
    when:
    task.ext.when == null || task.ext.when
    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_summarizedexperiment.r \\
        NULL \\
        $counts \\
        $tpm
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-summarizedexperiment: \$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}
