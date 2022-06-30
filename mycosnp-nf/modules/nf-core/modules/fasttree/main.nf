process FASTTREE {
    conda (params.enable_conda ? "bioconda::fasttree=2.1.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fasttree:2.1.10--h516909a_4' :
        'quay.io/biocontainers/fasttree:2.1.10--h516909a_4' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: true,mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.endsWith(".tre") ? "fasttree_phylogeny.nh" : filename  },path: { "${params.outdir}/combined/phylogeny/fasttree" },pattern: "*"
    input:
    path alignment
    output:
    path "*.tre",         emit: phylogeny
    path "versions.yml" , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = "-gtr -gamma -fastest"
    """
    fasttree \\
        $args \\
        -log fasttree_phylogeny.tre.log \\
        -nt $alignment \\
        > fasttree_phylogeny.tre
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$(fasttree -help 2>&1 | head -1  | sed 's/^FastTree \\([0-9\\.]*\\) .*\$/\\1/')
    END_VERSIONS
    """
}
