process GATHER_SUMMARY_LINES {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'
    input:
    path(summary_line_files)
    path(outdir_path)
    val(busco_val)
    output:
    path 'Phoenix_Output_Report.tsv'  , emit: summary_report
    path "versions.yml"               , emit: versions
    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def busco_parameter = busco_val ? "--busco" : ""
    """
    #check if the empty summaryline file exists - comes from fetch_failed_summaries.nf module
    if [ -f "empty_summaryline.tsv" ]; then
        rm empty_summaryline.tsv
        new_summary_line_files=\$(echo $summary_line_files | sed 's/empty_summaryline.tsv //')
        if [ -f "placeholder.txt" ]; then
            rm placeholder.txt
            new_summary_line_files=\$(echo \$new_summary_line_files | sed 's/placeholder.txt //')
        fi
python ${workflow.launchDir}/bin/Create_phoenix_summary_tsv.py --out Phoenix_Output_Report.tsv $busco_parameter \$new_summary_line_files
    elif [ -f "placeholder.txt" ]; then
        rm placeholder.txt
        new_summary_line_files=\$(echo $summary_line_files | sed 's/placeholder.txt //')
        if [ -f "empty_summaryline.tsv" ]; then
            rm empty_summaryline.tsv
            new_summary_line_files=\$(echo \$new_summary_line_files | sed 's/empty_summaryline.tsv //')
        fi
python ${workflow.launchDir}/bin/Create_phoenix_summary_tsv.py --out Phoenix_Output_Report.tsv $busco_parameter \$new_summary_line_files
    else
python ${workflow.launchDir}/bin/Create_phoenix_summary_tsv.py \\
            --out Phoenix_Output_Report.tsv \\
            $busco_parameter \\
            $summary_line_files
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
