process QC_REPORT {
	publishDir  path: { "${params.outdir}/reportqc"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' : 
        'quay.io/biocontainers/pandas:1.1.5' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    
cpus 6
    
memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    maxForks 10
    input:
    tuple val(meta), path(txt), path(results) //input values are from channel that joins FAQCS("txt") and QUALIMAP("results") outputs
    path reference
    output:
    path("*_output.txt"), emit: qc_line
    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python $projectDir/bin/qc_report_stats.py \\
        --sample ${meta.id} \\
        --stats ${meta.id}.stats.txt \\
        --base_content_before_trim qa.${meta.id}.base_content.txt \\
        --base_content_after_trim ${meta.id}.base_content.txt \\
        --qual_scores_before_trim qa.${meta.id}.for_qual_histogram.txt \\
        --qual_scores_after_trim ${meta.id}.for_qual_histogram.txt \\
        --reference ${reference} \\
        --bam_coverage ${meta.id}/genome_results.txt > ${meta.id}_output.txt
    """
}
