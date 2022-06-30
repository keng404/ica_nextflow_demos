process FILTER_GATK_GENOTYPES {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::scipy=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scipy%3A1.1.0' :
        'quay.io/biocontainers/scipy:1.1.0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: "${params.save_debug}", mode: "${params.publish_dir_mode}",path: { "${params.outdir}/combined/selectedsnpsfiltered" },pattern: "*{vcf.gz,vcf.gz.tbi}"
    input:
    tuple val(meta), path(vcf)
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    // path "versions.yml", emit: versions
    script:
    def args = "${params.vcftools_filter}"
    def prefix = "combined_genotype_filtered_snps_filtered"
    def is_compressed_vcf = vcf.getName().endsWith(".gz") ? true : false
    def vcf_name = vcf.getName().replace(".gz", "")
    """
    if [ "$is_compressed_vcf" == "true" ]; then
        gzip -c -d $vcf > $vcf_name
    fi
    python $projectDir/bin/broad-vcf-filter/filterGatkGenotypes.py  $vcf_name \\
                            $args \\
                           > ${prefix}.vcf
    gzip ${prefix}.vcf
    """
}
