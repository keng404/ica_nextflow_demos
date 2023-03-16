process GENERATE_PIPELINE_STATS {
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
    path(assembly_scaffolds), \
    path(filtered_assembly), \
    path(mlst_file), \
    path(gamma_HV), \
    path(gamma_AR), \
    path(gamma_replicon), \
    path(quast_report), \
    path(krona_weighted), \
    path(kraken2_weighted_report), \
    path(kraken2_weighted_summary), \
    path(taxID), \
    path(fastANI_formatted_file), \
    path(assembly_ratio_file), \
    path(amr_file), \
    path(gc_content)
    output:
    tuple val(meta), path('*.synopsis'), emit: pipeline_stats
    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) {
        terra = ""
    } else if (params.terra==true) {
        terra = "-5 terra"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
bash ${workflow.launchDir}/bin/pipeline_stats_writer.sh \\
        -a $fastp_raw_qc \\
        -b $fastp_total_qc \\
        -c ${trimmed_reads[0]} \\
        -d ${trimmed_reads[1]} \\
        -e $kraken2_trimd_report \\
        -f $kraken2_trimd_summary \\
        -g $krona_trimd \\
        -h $assembly_scaffolds \\
        -i $filtered_assembly \\
        -m $kraken2_weighted_report \\
        -n $kraken2_weighted_summary \\
        -o $krona_weighted \\
        -p $quast_report \\
        -q $taxID \\
        -r $assembly_ratio_file \\
        -t $fastANI_formatted_file \\
        -u $gamma_AR \\
        -v $gamma_replicon \\
        -w $gamma_HV \\
        -y $mlst_file \\
        -2 $amr_file \\
        -4 $gc_content \\
        $terra
    """
}
