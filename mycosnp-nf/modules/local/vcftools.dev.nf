process FILTER_GATK_GENOTYPES {
	publishDir  path: { "${params.outdir}/genotypesgatkfilter"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::scipy=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scipy%3A1.1.0' :
        'quay.io/biocontainers/scipy:1.1.0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    
cpus 6
    
memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    maxForks 10
    input:
    tuple val(meta), path(vcf)
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    // path "versions.yml", emit: versions
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
