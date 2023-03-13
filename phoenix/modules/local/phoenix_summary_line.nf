process CREATE_SUMMARY_LINE {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(trimmed_qc_data_file), \
    path(mlst_file), \
    path(hypervirulence_gamma_file), \
    path(ar_gamma_file), \
    path(pf_gamma_file), \
    path(quast_report), \
    path(ratio_file), \
    path(synopsis), \
    path(taxonomy_file), \
    path(trimd_ksummary), \
    path(amr_report)

    output:
    path '*_summaryline.tsv'           , emit: line_summary
    path "versions.yml"                , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Phoenix_summary_line.py \\
        -q $quast_report \\
        -t $trimmed_qc_data_file \\
        -a $ar_gamma_file \\
        -v $hypervirulence_gamma_file \\
        -p $pf_gamma_file \\
        -r $ratio_file \\
        -m $mlst_file \\
        -u $amr_report \\
        -n ${prefix} \\
        -s $synopsis \\
        -x $taxonomy_file \\
        -k $trimd_ksummary \\
        -o ${prefix}_summaryline.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
