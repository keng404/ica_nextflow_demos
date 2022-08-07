#!/usr/bin/env nextflow


println("Currently using the Grandeur workflow for use with microbial sequencing. The view is great from 8299 feet (2530 meters) above sea level.\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v2.0.20220610")
println("")

nextflow.enable.dsl               = 2
// gunzip process

params.prokka        = true
params.blobtools     = true



params.outdir                     = workflow.launchDir + '/out/grandeur'
params.maxcpus                    = 16
params.medcpus                    = 8
params.gff_files = null

// core workflow of fastq to contig
params.fastp_options              = "--detect_adapter_for_pe"
params.bbduk_options              = "k=31 hdist=1"
params.spades_options             = '--isolate'

// fastq information
params.fastq_processes            = ['fastp', 'bbduk', 'spades', 'fastqc', 'cg_pipeline', 'mash', 'kraken2', 'summary', 'multiqc', 'shigatyper']
params.fastqc_options             = ''
params.cg_pipeline_options        = '--qual_offset 33 --minLength 1'
params.shigatyper_options         = ''
params.plasmidfinder_options      = ''
params.kraken2_db_dir = null
if(params.kraken2_db_dir != null){
  params.kraken2_db                 = "${workflow.launchDir/params.kraken2_db_dir}"
} else{
  params.kraken2_db = null
}
params.kraken2_options            = ''
// WARNING : DO NOT CHANGE params.mash_reference UNLESS YOU ARE USING A DIFFERENT MASH CONTAINER
params.mash_reference             = '/db/RefSeqSketchesDefaults.msh'
params.mash_options               = '-v 0 -d 0.5'
// duplicates as seqsero2 and serotypefinder defaults are for contigs
// params.seqsero2_options        = '-t 2 -m a -b mem'
// params.serotypefinder_options  = ''

// contig information
params.contig_processes           = ['amrfinderplus', 'kleborate', 'fastani', 'mlst', 'quast', 'serotypefinder', 'blobtools', 'summary', 'multiqc', 'plasmidfinder', 'seqsero2', 'kraken2', 'mash']
params.amrfinderplus_options      = ''
params.fastani_options            = ''
params.kleborate_options          = '-all'
params.mlst_options               = ''
params.quast_options              = ''
params.serotypefinder_options     = ''
params.seqsero2_options           = '-m a -b mem'
// // duplicate as kraken2 with fastq reads is default
// // params.kraken2_options         = ''

// for blobtools
params.blast_db_dir = null
params.blast_db                   = null
if(params.blast_db_dir != null){
  params.blast_db = "${workflow.launchDir/params.blast_db_dir}"
}
params.local_db_type              = 'ref_prok_rep_genomes'
params.blobtools_create_options   = ''
params.blobtools_view_options     = ''
params.blobtools_plot_options     = '--format png -r species'
params.bwa_options                = ''
params.samtools_sort_options      = ''

// summary
params.multiqc_options            = ''

// phylogenetic analysis
params.phylogenetic_processes     = []
params.iqtree2_options            = '-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000'
params.outgroup                   = ''
params.prokka_options             = '--mincontiglen 500 --compliant --locustag locus_tag --centre STAPHB'
params.roary_options              = ''
params.snp_dists_options          = '-c'

include { de_novo_alignment }     from './subworkflows/de_novo_alignment.nf'     addParams( outdir:                     params.outdir,
                                                                                            fastq_processes:            params.fastq_processes,
                                                                                            spades_options:             params.spades_options,
                                                                                            bbduk_optons:               params.bbduk_options,
                                                                                            fastp_options:              params.fastp_options)
include { fastq_information }     from './subworkflows/fastq_information.nf'     addParams( outdir:                     params.outdir,
                                                                                            fastq_processes:            params.fastq_processes,
                                                                                            fastqc_options:             params.fastqc_options,
                                                                                            cg_pipeline_options:        params.cg_pipeline_options,
                                                                                            shigatyper_options:         params.shigatyper_options,
                                                                                            kraken2_options:            params.kraken2_options,
                                                                                            mash_reference:             params.mash_reference,
                                                                                            mash_options:               params.mash_options,
                                                                                            serotypefinder_options:     params.serotypefinder_options,
                                                                                            plasmidfinder_options:      params.plasmidfinder_options,
                                                                                            seqsero2_options:           params.seqsero2_options)
include { contig_information }    from './subworkflows/contig_information.nf'    addParams( outdir:                     params.outdir,
                                                                                            contig_processes:           params.contig_processes,
                                                                                            amrfinderplus_options:      params.amrfinderplus_options,
                                                                                            fastani_options:            params.fastani_options,
                                                                                            kleborate_options:          params.kleborate_options,
                                                                                            mlst_options:               params.mlst_options,
                                                                                            quast_options:              params.quast_options,
                                                                                            serotypefinder_options:     params.serotypefinder_options,
                                                                                            seqsero2_options:           params.seqsero2_options,
                                                                                            plasmidfinder_options:      params.plasmidfinder_options,
                                                                                            kraken2_options:            params.kraken2_options)
include { blobtools }             from './subworkflows/blobtools.nf'             addParams( local_db_type:              params.local_db_type,
                                                                                            blobtools_create_options:   params.blobtools_create_options,
                                                                                            blobtools_view_options:     params.blobtools_view_options,
                                                                                            blobtools_plot_options:     params.blobtools_plot_options,
                                                                                            bwa_options:                params.bwa_options,
                                                                                            samtools_sort_options:      params.samtools_sort_options)
include { phylogenetic_analysis } from './subworkflows/phylogenetic_analysis.nf' addParams( outdir:                     params.outdir,
                                                                                            phylogenetic_processes:     params.phylogenetic_processes,
                                                                                            iqtree2_options:            params.iqtree2_options,
                                                                                            outgroup:                   params.outgroup,
                                                                                            prokka_options:             params.prokka_options,
                                                                                            roary_options:              params.roary_options,
                                                                                            snp_dists_options:          params.snp_dists_options)
include { mash_dist as mash }     from './modules/mash'                          addParams( fastq_processes:            params.fastq_processes,
                                                                                            mash_reference:             params.mash_reference,
                                                                                            mash_options:               params.mash_options)
include { summary }               from './modules/summary'                       addParams( fastq_processes:            params.fastq_processes,
                                                                                            contig_processes:           params.contig_processes,
                                                                                            phylogenetic_processes:     params.phylogenetic_processes)
include { multiqc }               from './modules/multiqc'                       addParams( multiqc_options:            params.multiqc_options,
                                                                                            fastq_processes:            params.fastq_processes,
                                                                                            contig_processes:           params.contig_processes,
                                                                                            phylogenetic_processes:     params.phylogenetic_processes)

// TODO : frp_plasmid
// TODO : pointfinder
// TODO : mubsuite
// TODO : sistr
// TODO : plasmidseeker
// TODO : socru?
// TODO : kaptive
// TODO : ngmaster

println("The files and directory for results is " + params.outdir)
println("The maximum number of CPUS for any one process is ${params.maxcpus}")
params.genome_sizes = workflow.launchDir + "/configs/genome_sizes.json"
params.fastani_refs = workflow.launchDir + "/configs/fastani_ref.tar.gz"

include { GRANDEUR } from './workflows/grandeur'
workflow {
  GRANDEUR()
}


workflow.onError {
    // copy intermediate files + directories
    println("Getting intermediate files from ICA")
    ['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
    // return trace files
    println("Returning workflow run-metric reports from ICA")
    ['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null', '|', 'xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()

}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("MultiQC report can be found at ${params.outdir}/multiqc/multiqc_report.html")
    println("Summary can be found at ${params.outdir}/grandeur_results.tsv")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
