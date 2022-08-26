process TABIX_TABIX {
	publishDir  path: { "${params.outdir}/tabixtabix"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    conda (params.enable_conda ? 'bioconda::tabix=1.11' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    
cpus 6
    
memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    input:
    tuple val(meta), path(tab)
    output:
    tuple val(meta), path("*.tbi"), optional:true, emit: tbi
    tuple val(meta), path("*.csi"), optional:true, emit: csi
    path  "versions.yml"          , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: '-p vcf -f'
    """
    tabix $args $tab
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${tab}.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
