process FETCH_FAILED_SUMMARIES {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(directory)
    path(failed_summaries)
    path(summaries)

    output:
    path('*_summaryline.tsv'), emit: spades_failure_summary_line

    script:
    """
    #for each summaryline_failure.tsv file check to see if 'SPAdes_Failure' is in the file.
    if [ -f ${directory}/*/*_summaryline_failure.tsv ]; then
        for file in ${directory}/*/*_summaryline_failure.tsv; do 
            if grep -q SPAdes_Failure "\$file"; then
                # if so then add the sample name to the front of the file and move it to the correct place.
                fname=\$(basename \$file _summaryline_failure.tsv)
                cp \$file \${fname}_summaryline.tsv
                mv \$file ${directory}/\${fname}/\${fname}_summaryline.tsv
            fi
        done
    # If the summarylines file doesn't exist then just create an empty file. 
    else
        touch empty_summaryline.tsv
    fi
    """
}