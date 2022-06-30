process RAXMLNG {
    conda (params.enable_conda ? 'bioconda::raxml-ng=1.0.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:1.0.3--h32fcf60_0' :
        'quay.io/biocontainers/raxml-ng:1.0.3--h32fcf60_0' }"
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    publishDir  enabled: true,mode: "${params.publish_dir_mode}",saveAs: { filename -> if( filename.endsWith(".bestTree")) { return "raxmlng_bestTree.nh" } else if ( filename.endsWith(".support") ) { return "raxmlng_support.nh" } else { return filename }  },path: { "${params.outdir}/combined/phylogeny/raxmlng" },pattern: "*"
    input:
    path alignment
    output:
    path "*.raxml.bestTree", emit: phylogeny
    path "*.raxml.support" , optional:true, emit: phylogeny_bootstrapped
    path "versions.yml"    , emit: versions
    when:
    task.ext.when == null || task.ext.when
    script:
    def args = "--all --model GTR+G --bs-trees 1000"
    """
    raxml-ng \\
        $args \\
        --msa $alignment \\
        --threads 8 \\
        --prefix output
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxmlng: \$(echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//')
    END_VERSIONS
    """
}
