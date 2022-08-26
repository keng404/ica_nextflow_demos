process BCFTOOLS_QUERY {
    tag "$meta.id"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    publishDir  path: { "${params.outdir}/bcftools_query"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    cpus 6    
    memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    conda (params.enable_conda ? "bioconda::bcftools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.15.1--h0ea216a_0':
        'quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0' }"
    input:
    tuple val(meta), path(vcf), path(tbi)
    path regions
    path targets
    path samples
    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?:  "-H -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t[%DP\\t]\\t[%REF_DP\\t]\\t[%ALT_DP\\t]\\n'"
    def prefix = task.ext.prefix ?: "${meta.id}.bcftools_query"
    def regions_file = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file =  samples ? "--samples-file ${samples}" : ""
    """
    bcftools query \\
        --output ${prefix}.txt \\
        $regions_file \\
        $targets_file \\
        $samples_file \\
        $args \\
        $vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
