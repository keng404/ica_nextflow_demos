process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'
    input:
    path samplesheet
    output:
    path '*.valid.csv' , emit: csv
    path "versions.yml", emit: versions
    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    """
python ${workflow.launchDir}/bin/check_samplesheet.py \\
    $samplesheet \\
    samplesheet.valid.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
