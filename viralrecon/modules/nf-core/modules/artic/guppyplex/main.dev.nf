process ARTIC_GUPPYPLEX {
	publishDir  path: { "${params.outdir}/guppyplexartic"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::artic=1.2.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.2--pyhdfd78af_0' :
        'quay.io/biocontainers/artic:1.2.2--pyhdfd78af_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    
cpus 6
    
memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    input:
    tuple val(meta), path(fastq_dir)
    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "versions.yml"               , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.2' // WARN: Version information provided by tool on CLI is incorrect. Please update this string when bumping container versions.
    """
    artic \\
        guppyplex \\
        $args \\
        --directory $fastq_dir \\
        --output ${prefix}.fastq
    pigz -p $task.cpus *.fastq
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: $VERSION
    END_VERSIONS
    """
}
