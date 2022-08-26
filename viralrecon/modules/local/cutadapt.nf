process CUTADAPT {
	publishDir  path: { "${params.outdir}/cutadapt"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::cutadapt=3.5' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0' :
        'quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    
cpus 6
    
memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    input:
    tuple val(meta), path(reads)
    path adapters
    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: '--overlap 5 --minimum-length 30 --error-rate 0.1'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "-a file:adapters.sub.fa" : "-a file:adapters.sub.fa -A file:adapters.sub.fa"
    def trimmed = meta.single_end ? "-o ${prefix}.fastq.gz" : "-o ${prefix}_1.fastq.gz -p ${prefix}_2.fastq.gz"
    """
    sed -r '/^[ACTGactg]+\$/ s/\$/X/g' $adapters > adapters.sub.fa
    cutadapt \\
        --cores $task.cpus \\
        $args \\
        $paired \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
