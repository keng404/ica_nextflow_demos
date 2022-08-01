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
def primer_set         = ''
def primer_set_version = 0
if (params.platform == 'illumina' && params.protocol == 'amplicon') {
    primer_set         = params.primer_set
    primer_set_version = params.primer_set_version
params.input = null
params.platform = null
params.protocol = null
params.genome = null
params.primer_set = null
params.primer_set_version = null
params.primer_fasta = null
params.primer_left_suffix = '_LEFT'
params.primer_right_suffix = '_RIGHT'
params.save_reference = false
params.fastq_dir = null
params.fast5_dir = null
params.sequencing_summary = null
params.min_barcode_reads = "100"
params.min_guppyplex_reads = "10"
params.artic_minion_caller = 'nanopolish'
params.artic_minion_aligner = 'minimap2'
params.artic_minion_medaka_model = null
params.skip_pycoqc = false
params.skip_nanoplot = false
params.asciigenome_read_depth = "50"
params.asciigenome_window_size = "50"
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
params.kraken2_db = 's3://nf-core-awsmegatests/viralrecon/input_data/kraken2_human.tar.gz'
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
params.min_mapped_reads = "1000"
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
params.max_cpus = "16"
params.max_time = '240.h'
params.params.enable_conda = true
} else if (params.platform == 'nanopore') {
    primer_set          = 'artic'
    primer_set_version  = params.primer_set_version
    params.artic_scheme = WorkflowMain.getGenomeAttribute(params, 'scheme', log, primer_set, primer_set_version)
}
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
        ILLUMINA ()
    //
    // WORKFLOW: Variant analysis for Nanopore data
    //
    } else if (params.platform == 'nanopore') {
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
