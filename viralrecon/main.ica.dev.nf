#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/viralrecon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/viralrecon
    Website: https://nf-co.re/viralrecon
    Slack  : https://nfcore.slack.com/channels/viralrecon
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def primer_set         = 'artic'
def primer_set_version = 1
params.input = null
params.platform = null
params.protocol = null
params.genome = null
params.primer_set_version = primer_set_version
params.primer_set = primer_set
params.primer_fasta = null
params.primer_left_suffix = '_LEFT'
params.primer_right_suffix = '_RIGHT'
params.save_reference = false
params.fastq_dir = null
params.fast5_dir = null
params.sequencing_summary = null
params.min_barcode_reads = 100
params.min_guppyplex_reads = 10
params.artic_minion_caller = 'nanopolish'
params.artic_minion_aligner = 'minimap2'
params.artic_minion_medaka_model = null
params.skip_pycoqc = false
params.skip_nanoplot = false
params.asciigenome_read_depth = 50
params.asciigenome_window_size = 50
params.multiqc_title = null
params.multiqc_config = null
params.max_multiqc_email_size = '25.MB'
params.skip_mosdepth = false
params.skip_pangolin = false
params.skip_nextclade = false
params.skip_variants_quast = false
params.skip_snpeff = false
params.skip_asciigenome = false
params.skip_variants_long_table = false
params.skip_multiqc = false
// params.kraken2_db = 'ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz'
params.kraken2_db_name = 'human'
params.kraken2_variants_host_filter = false
params.kraken2_assembly_host_filter = true
params.save_trimmed_fail = false
params.skip_fastqc = false
params.skip_kraken2 = false
params.skip_fastp = false
params.skip_cutadapt = false
params.variant_caller = null
params.consensus_caller = 'bcftools'
params.min_mapped_reads = 1000
params.ivar_trim_noprimer = false
params.ivar_trim_offset = null
params.filter_duplicates = false
params.save_unaligned = false
params.save_mpileup = false
params.skip_ivar_trim = false
params.skip_markduplicates = true
params.skip_picard_metrics = false
params.skip_consensus_plots = false
params.skip_consensus = false
params.skip_variants = false
params.assemblers = 'spades'
params.spades_mode = 'rnaviral'
params.spades_hmm = null
params.blast_db = null
params.skip_bandage = false
params.skip_blast = false
params.skip_abacas = false
params.skip_plasmidid = false
params.skip_assembly_quast = false
params.skip_assembly = false
params.outdir = null
params.tracedir = "${params.outdir}/pipeline_info"
params.publish_dir_mode = 'copy'
params.email = null
params.email_on_fail = null
params.plaintext_email = false
params.monochrome_logs = false
params.help = false
params.validate_params = true
params.show_hidden_params = false
params.schema_ignore_params = 'genomes,igenomes_base'
params.enable_conda = false
params.custom_config_version = 'master'
params.custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
params.config_profile_description = null
params.config_profile_contact = null
params.config_profile_url = null
params.config_profile_name = null
params.max_memory = '128.GB'
params.max_cpus = 16
params.max_time = '240.h'
if (params.platform == 'illumina' && params.protocol == 'amplicon') {
    primer_set         = params.primer_set
    primer_set_version = params.primer_set_version
} else if (params.platform == 'nanopore') {
    primer_set          = 'artic'
    primer_set_version  = params.primer_set_version
    params.artic_scheme = WorkflowMain.getGenomeAttribute(params, 'scheme', log, primer_set, primer_set_version)
}

params.genomes = [:]

params.genomes['NC_045512.2'] = [:]
params.genomes['NC_045512.2'].fasta            = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/NC_045512.2/GCF_009858895.2_ASM985889v3_genomic.200409.fna.gz'
params.genomes['NC_045512.2'].gff              = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/NC_045512.2/GCF_009858895.2_ASM985889v3_genomic.200409.gff.gz'
params.genomes['NC_045512.2'].nextclade_dataset           = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/nextclade_sars-cov-2_MN908947_2022-01-18T12_00_00Z.tar.gz'
params.genomes['NC_045512.2'].nextclade_dataset_name      = 'sars-cov-2'
params.genomes['NC_045512.2'].nextclade_dataset_reference = 'MN908947'
params.genomes['NC_045512.2'].nextclade_dataset_tag       = '2022-01-18T12:00:00Z'

params.genomes['MN908947.3'] = [:]
params.genomes['MN908947.3'].fasta            = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.fna.gz'
params.genomes['MN908947.3'].gff              = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz'
params.genomes['MN908947.3'].nextclade_dataset           = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/nextclade_sars-cov-2_MN908947_2022-01-18T12_00_00Z.tar.gz'
params.genomes['MN908947.3'].nextclade_dataset_name      = 'sars-cov-2'
params.genomes['MN908947.3'].nextclade_dataset_reference = 'MN908947'
params.genomes['MN908947.3'].nextclade_dataset_tag       = '2022-01-18T12:00:00Z'
params.genomes['MN908947.3'].primer_sets = [:]
params.genomes['MN908947.3'].primer_sets.artic = [:]

params.genomes['MN908947.3'].primer_sets.artic['1'] = [:]
params.genomes['MN908947.3'].primer_sets.artic['1'].fasta      = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V1/nCoV-2019.reference.fasta'
params.genomes['MN908947.3'].primer_sets.artic['1'].gff        = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz'
params.genomes['MN908947.3'].primer_sets.artic['1'].primer_bed = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V1/nCoV-2019.primer.bed'
params.genomes['MN908947.3'].primer_sets.artic['1'].scheme     = 'nCoV-2019'
params.genomes['MN908947.3'].primer_sets.artic['2'] = [:] 
params.genomes['MN908947.3'].primer_sets.artic['2'].fasta      = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V2/nCoV-2019.reference.fasta'
params.genomes['MN908947.3'].primer_sets.artic['2'].gff        = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz'
params.genomes['MN908947.3'].primer_sets.artic['2'].primer_bed = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V2/nCoV-2019.primer.bed'
params.genomes['MN908947.3'].primer_sets.artic['2'].scheme     = 'nCoV-2019'

params.genomes['MN908947.3'].primer_sets.artic['3'] = [:]
params.genomes['MN908947.3'].primer_sets.artic['3'].fasta      = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta'
params.genomes['MN908947.3'].primer_sets.artic['3'].gff        = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz'
params.genomes['MN908947.3'].primer_sets.artic['3'].primer_bed = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V3/nCoV-2019.primer.bed'
params.genomes['MN908947.3'].primer_sets.artic['3'].scheme     = 'nCoV-2019'

params.genomes['MN908947.3'].primer_sets.artic['4'] = [:]
params.genomes['MN908947.3'].primer_sets.artic['4'].fasta      = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V4/SARS-CoV-2.reference.fasta'
params.genomes['MN908947.3'].primer_sets.artic['4'].gff        = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz'
params.genomes['MN908947.3'].primer_sets.artic['4'].primer_bed = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V4/SARS-CoV-2.scheme.bed'
params.genomes['MN908947.3'].primer_sets.artic['4'].scheme     = 'SARS-CoV-2'

params.genomes['MN908947.3'].primer_sets.artic['4.1'] = [:]
params.genomes['MN908947.3'].primer_sets.artic['4.1'].fasta      = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.reference.fasta'
params.genomes['MN908947.3'].primer_sets.artic['4.1'].gff        = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz'
params.genomes['MN908947.3'].primer_sets.artic['4.1'].primer_bed = 'https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.scheme.bed'
params.genomes['MN908947.3'].primer_sets.artic['4.1'].scheme     = 'SARS-CoV-2'
  
params.genomes['MN908947.3'].primer_sets.artic['1200'] = [:]
params.genomes['MN908947.3'].primer_sets.artic['1200'].fasta      = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/primer_schemes/artic/nCoV-2019/V1200/nCoV-2019.reference.fasta'
params.genomes['MN908947.3'].primer_sets.artic['1200'].gff        = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz'
params.genomes['MN908947.3'].primer_sets.artic['1200'].primer_bed = 'https://github.com/nf-core/test-datasets/raw/viralrecon/genome/MN908947.3/primer_schemes/artic/nCoV-2019/V1200/nCoV-2019.bed'
params.genomes['MN908947.3'].primer_sets.artic['1200'].scheme     = 'nCoV-2019'

params.fasta         = WorkflowMain.getGenomeAttribute(params, 'fasta'     , log, primer_set, primer_set_version)
params.gff           = WorkflowMain.getGenomeAttribute(params, 'gff'       , log, primer_set, primer_set_version)
params.bowtie2_index = WorkflowMain.getGenomeAttribute(params, 'bowtie2'   , log, primer_set, primer_set_version)
params.primer_bed    = WorkflowMain.getGenomeAttribute(params, 'primer_bed', log, primer_set, primer_set_version)
params.nextclade_dataset           = WorkflowMain.getGenomeAttribute(params, 'nextclade_dataset'          , log, primer_set, primer_set_version)
params.nextclade_dataset_name      = WorkflowMain.getGenomeAttribute(params, 'nextclade_dataset_name'     , log, primer_set, primer_set_version)
params.nextclade_dataset_reference = WorkflowMain.getGenomeAttribute(params, 'nextclade_dataset_reference', log, primer_set, primer_set_version)
params.nextclade_dataset_tag       = WorkflowMain.getGenomeAttribute(params, 'nextclade_dataset_tag'      , log, primer_set, primer_set_version)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
WorkflowMain.initialise(workflow, params, log)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
if (params.platform == 'illumina') {
    include { ILLUMINA } from './workflows/illumina'
} else if (params.platform == 'nanopore') {
    include { NANOPORE } from './workflows/nanopore'
}
workflow NFCORE_VIRALRECON {
    //
    // WORKFLOW: Variant and de novo assembly analysis for Illumina data
    //
    if (params.platform == 'illumina') {
        if(params.kraken2_db_file){
            params.kraken2_db = "${workflow.launchDir}/" + params.kraken2_db_file
        }
        println("my KRAKEN2 DB is here: " +  params.kraken2_db)
        ILLUMINA ()
    //
    // WORKFLOW: Variant analysis for Nanopore data
    //
    } else if (params.platform == 'nanopore') {
                if(params.kraken2_db_file){
            params.kraken2_db = "${workflow.launchDir}/" + params.kraken2_db_file
        }
        println("my KRAKEN2 DB is here: " +  params.kraken2_db)
        NANOPORE ()
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_VIRALRECON ()
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onError {
// copy intermediate files + directories
println("Getting intermediate files from ICA")
['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
// return trace files
println("Returning workflow run-metric reports from ICA")
['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null','|','xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()
}
