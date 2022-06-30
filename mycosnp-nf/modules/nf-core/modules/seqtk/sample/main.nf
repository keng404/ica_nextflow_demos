process SEQTK_SAMPLE {
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: "${params.save_debug}",mode: "${params.publish_dir_mode}",path: { "${params.outdir}/samples/${meta.id}/seqkit_sample" },pattern: "*.fastq.gz"
    input:
    tuple val(meta), path(reads)
    val sample_size
    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = "-s12237"
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        seqtk \\
            sample \\
            $args \\
            $reads \\
            $sample_size \\
            | gzip --no-name > ${prefix}.fastq.gz \\
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        if (!(args ==~ /.*-s[0-9]+.*/)) {
            args += " -s100"
        }
        """
        seqtk \\
            sample \\
            $args \\
            ${reads[0]} \\
            $sample_size \\
            | gzip --no-name > ${prefix}_1.fastq.gz \\
        seqtk \\
            sample \\
            $args \\
            ${reads[1]} \\
            $sample_size \\
            | gzip --no-name > ${prefix}_2.fastq.gz \\
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }
}
