process GENERATE_PIPELINE_STATS_FAILURE {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'
    input:
    tuple val(meta), path(trimmed_reads), \
    path(fastp_raw_qc), \
    path(fastp_total_qc), \
    path(kraken2_trimd_report), \
    path(krona_trimd), \
    path(kraken2_trimd_summary), \
    path(taxID), \
    val(spades_outcome)
    output:
    tuple val(meta), path('*.synopsis'), emit: pipeline_stats
    when:
    "${spades_outcome[0]}" == "run_failure" || "${spades_outcome[1]}" == "no_scaffolds" || "${spades_outcome[2]}" == "no_contigs"
    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
bash ${workflow.launchDir}/bin/pipeline_stats_writer.sh \\
        -a $fastp_raw_qc \\
        -b $fastp_total_qc \\
        -c ${trimmed_reads[0]} \\
        -d ${trimmed_reads[1]} \\
        -e $kraken2_trimd_report \\
        -f $kraken2_trimd_summary \\
        -g $krona_trimd
    """
}
