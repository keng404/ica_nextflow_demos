process KRAKEN2_KRAKEN2 {
	publishDir  path: { "${params.outdir}/kraken2kraken2"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::kraken2=2.1.2 conda-forge::pigz=2.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    
cpus 6
    
memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    input:
    tuple val(meta), path(reads)
    path  db
    val save_output_fastqs
    val save_reads_assignment
    output:
    tuple val(meta), path('*classified*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*unclassified*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classifiedreads*'), optional:true, emit: classified_reads_assignment
    tuple val(meta), path('*report.txt')                     , emit: report
    path "versions.yml"                                      , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def classified_command = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_command = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_command = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : ""
    def compress_reads_command = save_output_fastqs ? "pigz -p $task.cpus *.fastq" : ""
    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        $unclassified_command \\
        $classified_command \\
        $readclassification_command \\
        $paired \\
        $args \\
        $reads
    $compress_reads_command
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}