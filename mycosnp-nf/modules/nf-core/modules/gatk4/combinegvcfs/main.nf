process GATK4_COMBINEGVCFS {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    path (fasta)
    path (fasta_fai)
    path (fasta_dict)
    output:
    tuple val(meta), path("*.combined.g.vcf.gz"), emit: combined_gvcf
    path "versions.yml"                         , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem       = 16
    def input_files = vcf.collect{"-V ${it}"}.join(' ') // add '-V' to each vcf file
    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        CombineGVCFs \\
        -R ${fasta} \\
        -O ${prefix}.combined.g.vcf.gz \\
        ${args} \\
        ${input_files}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
