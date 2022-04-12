process SEQKIT_REPLACE {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::seqkit=2.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.1.0--h9ee0642_0':
        'quay.io/biocontainers/seqkit:2.1.0--h9ee0642_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: true,mode: "${params.publish_dir_mode}",path: { "${params.outdir}/combined/vcf-to-fasta" },pattern: "vcf-to-fasta.fasta"
    input:
    tuple val(meta), path(fastx)
    output:
    tuple val(meta), path("*.fast*"), emit: fastx
    path "versions.yml"             , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = "-s -p '\\*' -r '-'" 
    def prefix = "vcf-to-fasta"
    def extension = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz/) {
        extension = "fasta"
    }
    def endswith = "fasta" ?: "${extension}.gz"
    """
    seqkit \\
        replace \\
        ${args} \\
        --threads 12 \\
        -i ${fastx} \\
        -o ${prefix}.${endswith}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """
}
