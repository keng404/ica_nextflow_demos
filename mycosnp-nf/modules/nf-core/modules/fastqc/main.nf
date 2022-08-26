process FASTQC {
	publishDir  path: { "${params.outdir}/fastqc"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "keng404/fastqc:0.11.9--0"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    
cpus 6
    
memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    maxForks 10
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions
    script:
    def args = task.ext.args ?: '--quiet'
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${meta.id}.fastq.gz
        fastqc $args --threads $task.cpus ${meta.id}.fastq.gz
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
        END_VERSIONS
        // mv ${meta.id}_fastqc.zip ${reads[0].baseName}.fastqc.zip
        // mv ${meta.id}_fastqc.html ${reads[0].baseName}_fastqc.html
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${meta.id}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${meta.id}_2.fastq.gz
        fastqc $args --threads $task.cpus ${meta.id}_1.fastq.gz ${meta.id}_2.fastq.gz
        mv ${meta.id}_1_fastqc.zip ${reads[0].baseName}.1.fastqc.zip
        mv ${meta.id}_2_fastqc.zip ${reads[1].baseName}.2.fastqc.zip
        mv ${meta.id}_1_fastqc.html ${reads[0].baseName}.1.fastqc.html
        mv ${meta.id}_2_fastqc.html ${reads[1].baseName}.2.fastqc.html
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
        END_VERSIONS
        """
    }
}