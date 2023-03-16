process DETERMINE_TAXA_ID {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'
    input:
    tuple val(meta), path(kraken_weighted), path(formatted_ani_file), path(k2_bh_summary)
    path(taxa_file)
    output:
    tuple val(meta), path('*.tax'), emit: taxonomy
    path "versions.yml"           , emit: versions
    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
bash ${workflow.launchDir}/bin/determine_taxID.sh -k $kraken_weighted -r $k2_bh_summary -s $meta.id -f $formatted_ani_file -d $taxa_file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NCBI Taxonomy Reference File: $taxa_file
    END_VERSIONS
    """
}
