process DETERMINE_TOP_TAXA {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.0.0'
    input:
    tuple val(meta), path(mash_dists), path(assembly_scaffolds)
    output:
    tuple val(meta), path('*_best_MASH_hits.txt'), emit: top_taxa_list
    tuple val(meta), path('*_genomic.fna.gz'),     emit: reference_files
    path("versions.yml"),                          emit: versions
    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) {
        terra = ""
    } else if (params.terra==true) {
        terra = "-t terra"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
bash ${workflow.launchDir}/bin/sort_and_prep_dist.sh -a $assembly_scaffolds -x $mash_dists $terra
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Date of RefSeq Pull: \$(date +"%d-%m-%y")
    END_VERSIONS
    """
}
