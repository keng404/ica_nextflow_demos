/*
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def valid_params = [
    aligners       : ['star_salmon', 'star_rsem', 'hisat2'],
    pseudoaligners : ['salmon'],
    rseqc_modules  : ['bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation', 'junction_saturation', 'read_distribution', 'read_duplication', 'tin']
]
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
// Validate input parameters
WorkflowRnaseq.initialise(params, log, valid_params)
// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.transcript_fasta, params.additional_fasta,
    params.gtf, params.gff, params.gene_bed,
    params.ribo_database_manifest, params.splicesites,
    params.star_index, params.hisat2_index, params.rsem_index, params.salmon_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
// Check rRNA databases for sortmerna
if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}
// Check if file with list of fastas is provided when running BBSplit
if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
    ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list, checkIfExists: true)
    if (ch_bbsplit_fasta_list.isEmpty()) {exit 1, "File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!"}
}
// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_bbsplit)   { prepareToolIndices << 'bbsplit'             }
if (!params.skip_alignment) { prepareToolIndices << params.aligner        }
if (params.pseudo_aligner)  { prepareToolIndices << params.pseudo_aligner }
// Get RSeqC modules to run
def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
if (params.bam_csi_index) {
    for (rseqc_module in ['read_distribution', 'inner_distance', 'tin']) {
        if (rseqc_modules.contains(rseqc_module)) {
            rseqc_modules.remove(rseqc_module)
        }
    }
}
// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}
// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
// Header files for MultiQC
ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Loaded from modules/local/
//
include { BEDTOOLS_GENOMECOV                 } from '../modules/local/bedtools_genomecov'
include { DESEQ2_QC as DESEQ2_QC_STAR_SALMON } from '../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_RSEM        } from '../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_SALMON      } from '../modules/local/deseq2_qc'
include { DUPRADAR                           } from '../modules/local/dupradar'
include { MULTIQC                            } from '../modules/local/multiqc'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../modules/local/multiqc_custom_biotype'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_FAIL_MAPPED  } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_STRAND_CHECK } from '../modules/local/multiqc_tsv_from_list'
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ALIGN_STAR     } from '../subworkflows/local/align_star'
include { QUANTIFY_RSEM  } from '../subworkflows/local/quantify_rsem'
include { QUANTIFY_SALMON as QUANTIFY_STAR_SALMON } from '../subworkflows/local/quantify_salmon'
include { QUANTIFY_SALMON as QUANTIFY_SALMON      } from '../subworkflows/local/quantify_salmon'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { BBMAP_BBSPLIT               } from '../modules/nf-core/modules/bbmap/bbsplit/main'
include { SAMTOOLS_SORT               } from '../modules/nf-core/modules/samtools/sort/main'
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/modules/preseq/lcextrap/main'
include { QUALIMAP_RNASEQ             } from '../modules/nf-core/modules/qualimap/rnaseq/main'
include { SORTMERNA                   } from '../modules/nf-core/modules/sortmerna/main'
include { STRINGTIE                   } from '../modules/nf-core/modules/stringtie/stringtie/main'
include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/modules/subread/featurecounts/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastqc_umitools_trimgalore'
include { ALIGN_HISAT2               } from '../subworkflows/nf-core/align_hisat2'
include { BAM_SORT_SAMTOOLS          } from '../subworkflows/nf-core/bam_sort_samtools'
include { MARK_DUPLICATES_PICARD     } from '../subworkflows/nf-core/mark_duplicates_picard'
include { RSEQC                      } from '../subworkflows/nf-core/rseqc'
include { DEDUP_UMI_UMITOOLS as DEDUP_UMI_UMITOOLS_GENOME        } from '../subworkflows/nf-core/dedup_umi_umitools'
include { DEDUP_UMI_UMITOOLS as DEDUP_UMI_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/dedup_umi_umitools'
include { BEDGRAPH_TO_BIGWIG as BEDGRAPH_TO_BIGWIG_FORWARD       } from '../subworkflows/nf-core/bedgraph_to_bigwig'
include { BEDGRAPH_TO_BIGWIG as BEDGRAPH_TO_BIGWIG_REVERSE       } from '../subworkflows/nf-core/bedgraph_to_bigwig'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Info required for completion email and summary
def multiqc_report      = []
def pass_percent_mapped = [:]
def fail_percent_mapped = [:]
workflow RNASEQ {
    ch_versions = Channel.empty()
    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    PREPARE_GENOME (
        prepareToolIndices,
        biotype
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    // Check if contigs in genome fasta file > 512 Mbp
    PREPARE_GENOME
        .out
        .fai
        .map { WorkflowRnaseq.checkMaxContigSize(it, log) }
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
params.outdir_custom = "${params.outdir}/check/input"
params.outdir_custom = "${params.outdir}/check/input"
params.outdir_custom = "${params.outdir}/check/input"
params.outdir_custom = "${params.outdir}/check/input"
params.outdir_custom = "${params.outdir}/check/input"
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
params.outdir_custom = "${params.outdir}/fastq/cat"
params.outdir_custom = "${params.outdir}/fastq/cat"
publishDir =  path: { "${params.outdir}/fastq" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }, enabled: params.save_merged_fastq
params.outdir_custom = "${params.outdir}/fastq/cat"
publishDir =  path: { "${params.outdir}/fastq" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }, enabled: params.save_merged_fastq
params.outdir_custom = "${params.outdir}/fastq/cat"
publishDir =  path: { "${params.outdir}/fastq" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }, enabled: params.save_merged_fastq
params.outdir_custom = "${params.outdir}/fastq/cat"
publishDir =  path: { "${params.outdir}/fastq" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }, enabled: params.save_merged_fastq
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))
    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //
params.outdir_custom = "${params.outdir}/trimgalore/umitools/fastqc"
params.outdir_custom = "${params.outdir}/trimgalore/umitools/fastqc"
params.outdir_custom = "${params.outdir}/trimgalore/umitools/fastqc"
params.outdir_custom = "${params.outdir}/trimgalore/umitools/fastqc"
params.outdir_custom = "${params.outdir}/trimgalore/umitools/fastqc"
    FASTQC_UMITOOLS_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_trimming,
        params.umi_discard_read
    )
    ch_versions = ch_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.versions)
    //
    // MODULE: Remove genome contaminant reads
    //
    ch_filtered_reads = FASTQC_UMITOOLS_TRIMGALORE.out.reads
    if (!params.skip_bbsplit) {
params.outdir_custom = "${params.outdir}/bbsplit/bbmap"
if (!params.skip_bbsplit) {
    process.ext.args = 'build=1ambiguous2=allmaxindel=150000'
}
params.outdir_custom = "${params.outdir}/bbsplit/bbmap"
if (!params.skip_bbsplit) {
    process.ext.args = 'build=1ambiguous2=allmaxindel=150000'
    publishDir =   path: { "${params.outdir}/bbsplit" }, mode: params.publish_dir_mode, pattern: '*.txt' ],  path: { "${params.outdir}/bbsplit" }, mode: params.publish_dir_mode, pattern: '*.fastq.gz', enabled: params.save_bbsplit_reads
}
params.outdir_custom = "${params.outdir}/bbsplit/bbmap"
if (!params.skip_bbsplit) {
    process.ext.args = 'build=1ambiguous2=allmaxindel=150000'
    publishDir =   path: { "${params.outdir}/bbsplit" }, mode: params.publish_dir_mode, pattern: '*.txt' ],  path: { "${params.outdir}/bbsplit" }, mode: params.publish_dir_mode, pattern: '*.fastq.gz', enabled: params.save_bbsplit_reads
}
params.outdir_custom = "${params.outdir}/bbsplit/bbmap"
if (!params.skip_bbsplit) {
    process.ext.args = 'build=1ambiguous2=allmaxindel=150000'
    publishDir =   path: { "${params.outdir}/bbsplit" }, mode: params.publish_dir_mode, pattern: '*.txt' ],  path: { "${params.outdir}/bbsplit" }, mode: params.publish_dir_mode, pattern: '*.fastq.gz', enabled: params.save_bbsplit_reads
}
params.outdir_custom = "${params.outdir}/bbsplit/bbmap"
if (!params.skip_bbsplit) {
    process.ext.args = 'build=1ambiguous2=allmaxindel=150000'
    publishDir =   path: { "${params.outdir}/bbsplit" }, mode: params.publish_dir_mode, pattern: '*.txt' ],  path: { "${params.outdir}/bbsplit" }, mode: params.publish_dir_mode, pattern: '*.fastq.gz', enabled: params.save_bbsplit_reads
}
        BBMAP_BBSPLIT (
            ch_filtered_reads,
            PREPARE_GENOME.out.bbsplit_index,
            [],
            [ [], [] ],
            false
        )
        .primary_fastq
        .set { ch_filtered_reads }
        ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())
    }
    //
    // MODULE: Remove ribosomal RNA reads
    //
    ch_sortmerna_multiqc = Channel.empty()
    if (params.remove_ribo_rna) {
        ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()
params.outdir_custom = "${params.outdir}/sortmerna"
if (params.remove_ribo_rna) {
    process.ext.args = '--num_alignments1--fastx-v'
}
params.outdir_custom = "${params.outdir}/sortmerna"
if (params.remove_ribo_rna) {
    process.ext.args = '--num_alignments1--fastx-v'
    publishDir =   path: { "${params.outdir}/sortmerna" }, mode: params.publish_dir_mode, pattern: "*.log" ],  path: { "${params.outdir}/sortmerna" }, mode: params.publish_dir_mode, pattern: "*.fastq.gz", enabled: params.save_non_ribo_reads
}
params.outdir_custom = "${params.outdir}/sortmerna"
if (params.remove_ribo_rna) {
    process.ext.args = '--num_alignments1--fastx-v'
    publishDir =   path: { "${params.outdir}/sortmerna" }, mode: params.publish_dir_mode, pattern: "*.log" ],  path: { "${params.outdir}/sortmerna" }, mode: params.publish_dir_mode, pattern: "*.fastq.gz", enabled: params.save_non_ribo_reads
}
params.outdir_custom = "${params.outdir}/sortmerna"
if (params.remove_ribo_rna) {
    process.ext.args = '--num_alignments1--fastx-v'
    publishDir =   path: { "${params.outdir}/sortmerna" }, mode: params.publish_dir_mode, pattern: "*.log" ],  path: { "${params.outdir}/sortmerna" }, mode: params.publish_dir_mode, pattern: "*.fastq.gz", enabled: params.save_non_ribo_reads
}
params.outdir_custom = "${params.outdir}/sortmerna"
if (params.remove_ribo_rna) {
    process.ext.args = '--num_alignments1--fastx-v'
    publishDir =   path: { "${params.outdir}/sortmerna" }, mode: params.publish_dir_mode, pattern: "*.log" ],  path: { "${params.outdir}/sortmerna" }, mode: params.publish_dir_mode, pattern: "*.fastq.gz", enabled: params.save_non_ribo_reads
}
        SORTMERNA (
            ch_filtered_reads,
            ch_sortmerna_fastas
        )
        .reads
        .set { ch_filtered_reads }
        ch_sortmerna_multiqc = SORTMERNA.out.log
        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
    }
    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_salmon') {
params.outdir_custom = "${params.outdir}/star/align"
params.outdir_custom = "${params.outdir}/star/align"
params.outdir_custom = "${params.outdir}/star/align"
params.outdir_custom = "${params.outdir}/star/align"
params.outdir_custom = "${params.outdir}/star/align"
        ALIGN_STAR (
            ch_filtered_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_STAR.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {
            // Deduplicate genome BAM file before downstream analysis
            DEDUP_UMI_UMITOOLS_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0])
            )
            ch_genome_bam        = DEDUP_UMI_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index  = DEDUP_UMI_UMITOOLS_GENOME.out.bai
            ch_samtools_stats    = DEDUP_UMI_UMITOOLS_GENOME.out.stats
            ch_samtools_flagstat = DEDUP_UMI_UMITOOLS_GENOME.out.flagstat
            ch_samtools_idxstats = DEDUP_UMI_UMITOOLS_GENOME.out.idxstats
            if (params.bam_csi_index) {
                ch_genome_bam_index  = DEDUP_UMI_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(DEDUP_UMI_UMITOOLS_GENOME.out.versions)
            // Co-ordinate sort, index and run stats on transcriptome BAM
params.outdir_custom = "${params.outdir}/samtools/sort/bam"
params.outdir_custom = "${params.outdir}/samtools/sort/bam"
params.outdir_custom = "${params.outdir}/samtools/sort/bam"
params.outdir_custom = "${params.outdir}/samtools/sort/bam"
params.outdir_custom = "${params.outdir}/samtools/sort/bam"
            BAM_SORT_SAMTOOLS (
                ch_transcriptome_bam
            )
            ch_transcriptome_sorted_bam = BAM_SORT_SAMTOOLS.out.bam
            ch_transcriptome_sorted_bai = BAM_SORT_SAMTOOLS.out.bai
            // Deduplicate transcriptome BAM file before read counting with Salmon
params.outdir_custom = "${params.outdir}/transcriptome/umitools/umi/dedup"
params.outdir_custom = "${params.outdir}/transcriptome/umitools/umi/dedup"
params.outdir_custom = "${params.outdir}/transcriptome/umitools/umi/dedup"
params.outdir_custom = "${params.outdir}/transcriptome/umitools/umi/dedup"
params.outdir_custom = "${params.outdir}/transcriptome/umitools/umi/dedup"
            DEDUP_UMI_UMITOOLS_TRANSCRIPTOME (
                ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0])
            )
            // Name sort BAM before passing to Salmon
params.outdir_custom = "${params.outdir}/sort/samtools"
params.outdir_custom = "${params.outdir}/sort/samtools"
params.outdir_custom = "${params.outdir}/sort/samtools"
params.outdir_custom = "${params.outdir}/sort/samtools"
params.outdir_custom = "${params.outdir}/sort/samtools"
            SAMTOOLS_SORT (
                DEDUP_UMI_UMITOOLS_TRANSCRIPTOME.out.bam
            )
            ch_transcriptome_bam = SAMTOOLS_SORT.out.bam
        }
        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
params.outdir_custom = "${params.outdir}/salmon/star/quantify"
params.outdir_custom = "${params.outdir}/salmon/star/quantify"
params.outdir_custom = "${params.outdir}/salmon/star/quantify"
params.outdir_custom = "${params.outdir}/salmon/star/quantify"
params.outdir_custom = "${params.outdir}/salmon/star/quantify"
        QUANTIFY_STAR_SALMON (
            ch_transcriptome_bam,
            ch_dummy_file,
            PREPARE_GENOME.out.transcript_fasta,
            PREPARE_GENOME.out.gtf,
            true,
            params.salmon_quant_libtype ?: ''
        )
        ch_versions = ch_versions.mix(QUANTIFY_STAR_SALMON.out.versions)
        if (!params.skip_qc & !params.skip_deseq2_qc) {
params.outdir_custom = "${params.outdir}/salmon/star/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_salmon'
}
params.outdir_custom = "${params.outdir}/salmon/star/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_salmon'
    publishDir =  path: { "${params.outdir}/${params.aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
params.outdir_custom = "${params.outdir}/salmon/star/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_salmon'
    publishDir =  path: { "${params.outdir}/${params.aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
params.outdir_custom = "${params.outdir}/salmon/star/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_salmon'
    publishDir =  path: { "${params.outdir}/${params.aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
params.outdir_custom = "${params.outdir}/salmon/star/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_salmon'
    publishDir =  path: { "${params.outdir}/${params.aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
            DESEQ2_QC_STAR_SALMON (
                QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_aligner_pca_multiqc        = DESEQ2_QC_STAR_SALMON.out.pca_multiqc
            ch_aligner_clustering_multiqc = DESEQ2_QC_STAR_SALMON.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_STAR_SALMON.out.versions)
        }
    }
    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with RSEM
    //
    ch_rsem_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_rsem') {
params.outdir_custom = "${params.outdir}/rsem/quantify"
params.outdir_custom = "${params.outdir}/rsem/quantify"
params.outdir_custom = "${params.outdir}/rsem/quantify"
params.outdir_custom = "${params.outdir}/rsem/quantify"
params.outdir_custom = "${params.outdir}/rsem/quantify"
        QUANTIFY_RSEM (
            ch_filtered_reads,
            PREPARE_GENOME.out.rsem_index
        )
        ch_genome_bam        = QUANTIFY_RSEM.out.bam
        ch_genome_bam_index  = QUANTIFY_RSEM.out.bai
        ch_samtools_stats    = QUANTIFY_RSEM.out.stats
        ch_samtools_flagstat = QUANTIFY_RSEM.out.flagstat
        ch_samtools_idxstats = QUANTIFY_RSEM.out.idxstats
        ch_star_multiqc      = QUANTIFY_RSEM.out.logs
        ch_rsem_multiqc      = QUANTIFY_RSEM.out.stat
        if (params.bam_csi_index) {
            ch_genome_bam_index = QUANTIFY_RSEM.out.csi
        }
        ch_versions = ch_versions.mix(QUANTIFY_RSEM.out.versions)
        if (!params.skip_qc & !params.skip_deseq2_qc) {
params.outdir_custom = "${params.outdir}/rsem/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_rsem'
}
params.outdir_custom = "${params.outdir}/rsem/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_rsem'
    publishDir =  path: { "${params.outdir}/${params.aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
params.outdir_custom = "${params.outdir}/rsem/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_rsem'
    publishDir =  path: { "${params.outdir}/${params.aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
params.outdir_custom = "${params.outdir}/rsem/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_rsem'
    publishDir =  path: { "${params.outdir}/${params.aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
params.outdir_custom = "${params.outdir}/rsem/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'star_rsem'
    publishDir =  path: { "${params.outdir}/${params.aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
            DESEQ2_QC_RSEM (
                QUANTIFY_RSEM.out.merged_counts_gene,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_aligner_pca_multiqc        = DESEQ2_QC_RSEM.out.pca_multiqc
            ch_aligner_clustering_multiqc = DESEQ2_QC_RSEM.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_RSEM.out.versions)
        }
    }
    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    ch_hisat2_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'hisat2') {
params.outdir_custom = "${params.outdir}/hisat2/align"
params.outdir_custom = "${params.outdir}/hisat2/align"
params.outdir_custom = "${params.outdir}/hisat2/align"
params.outdir_custom = "${params.outdir}/hisat2/align"
params.outdir_custom = "${params.outdir}/hisat2/align"
        ALIGN_HISAT2 (
            ch_filtered_reads,
            PREPARE_GENOME.out.hisat2_index,
            PREPARE_GENOME.out.splicesites
        )
        ch_genome_bam        = ALIGN_HISAT2.out.bam
        ch_genome_bam_index  = ALIGN_HISAT2.out.bai
        ch_samtools_stats    = ALIGN_HISAT2.out.stats
        ch_samtools_flagstat = ALIGN_HISAT2.out.flagstat
        ch_samtools_idxstats = ALIGN_HISAT2.out.idxstats
        ch_hisat2_multiqc    = ALIGN_HISAT2.out.summary
        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_HISAT2.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_HISAT2.out.versions)
        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {
            DEDUP_UMI_UMITOOLS_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0])
            )
            ch_genome_bam        = DEDUP_UMI_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index  = DEDUP_UMI_UMITOOLS_GENOME.out.bai
            ch_samtools_stats    = DEDUP_UMI_UMITOOLS_GENOME.out.stats
            ch_samtools_flagstat = DEDUP_UMI_UMITOOLS_GENOME.out.flagstat
            ch_samtools_idxstats = DEDUP_UMI_UMITOOLS_GENOME.out.idxstats
            if (params.bam_csi_index) {
                ch_genome_bam_index = DEDUP_UMI_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(DEDUP_UMI_UMITOOLS_GENOME.out.versions)
        }
    }
    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    ch_fail_mapping_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_multiqc
            .map { meta, align_log -> [ meta ] + WorkflowRnaseq.getStarPercentMapped(params, align_log) }
            .set { ch_percent_mapped }
        ch_genome_bam
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam }
        ch_genome_bam_index
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam_index }
        ch_percent_mapped
            .branch { meta, mapped, pass ->
                pass: pass
                    pass_percent_mapped[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    fail_percent_mapped[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }
        def header = [
            "Sample",
            "STAR uniquely mapped reads (%)"
        ]
params.outdir_custom = "${params.outdir}/mapped/fail/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
}
params.outdir_custom = "${params.outdir}/mapped/fail/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/mapped/fail/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/mapped/fail/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/mapped/fail/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
        MULTIQC_TSV_FAIL_MAPPED (
            ch_pass_fail_mapped.fail.collect(),
            header,
            'fail_mapped_samples'
        )
        .set { ch_fail_mapping_multiqc }
    }
    //
    // MODULE: Run Preseq
    //
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc && !params.skip_preseq) {
params.outdir_custom = "${params.outdir}/lcextrap/preseq"
if (!params.skip_preseq) {
    process.ext.args = '-verbose-bam-seed1-seg_len100000000'
}
params.outdir_custom = "${params.outdir}/lcextrap/preseq"
if (!params.skip_preseq) {
    process.ext.args = '-verbose-bam-seed1-seg_len100000000'
    publishDir =   path: { "${params.outdir}/${params.aligner}/preseq" }, mode: params.publish_dir_mode, pattern: "*.txt" ],  path: { "${params.outdir}/${params.aligner}/preseq/log" }, mode: params.publish_dir_mode, pattern: "*.log"
}
params.outdir_custom = "${params.outdir}/lcextrap/preseq"
if (!params.skip_preseq) {
    process.ext.args = '-verbose-bam-seed1-seg_len100000000'
    publishDir =   path: { "${params.outdir}/${params.aligner}/preseq" }, mode: params.publish_dir_mode, pattern: "*.txt" ],  path: { "${params.outdir}/${params.aligner}/preseq/log" }, mode: params.publish_dir_mode, pattern: "*.log"
}
params.outdir_custom = "${params.outdir}/lcextrap/preseq"
if (!params.skip_preseq) {
    process.ext.args = '-verbose-bam-seed1-seg_len100000000'
    publishDir =   path: { "${params.outdir}/${params.aligner}/preseq" }, mode: params.publish_dir_mode, pattern: "*.txt" ],  path: { "${params.outdir}/${params.aligner}/preseq/log" }, mode: params.publish_dir_mode, pattern: "*.log"
}
params.outdir_custom = "${params.outdir}/lcextrap/preseq"
if (!params.skip_preseq) {
    process.ext.args = '-verbose-bam-seed1-seg_len100000000'
    publishDir =   path: { "${params.outdir}/${params.aligner}/preseq" }, mode: params.publish_dir_mode, pattern: "*.txt" ],  path: { "${params.outdir}/${params.aligner}/preseq/log" }, mode: params.publish_dir_mode, pattern: "*.log"
}
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_preseq_multiqc = PRESEQ_LCEXTRAP.out.ccurve
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }
    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    ch_markduplicates_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_markduplicates) {
params.outdir_custom = "${params.outdir}/picard/duplicates/mark"
params.outdir_custom = "${params.outdir}/picard/duplicates/mark"
params.outdir_custom = "${params.outdir}/picard/duplicates/mark"
params.outdir_custom = "${params.outdir}/picard/duplicates/mark"
params.outdir_custom = "${params.outdir}/picard/duplicates/mark"
        MARK_DUPLICATES_PICARD (
            ch_genome_bam
        )
        ch_genome_bam             = MARK_DUPLICATES_PICARD.out.bam
        ch_genome_bam_index       = MARK_DUPLICATES_PICARD.out.bai
        ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = MARK_DUPLICATES_PICARD.out.metrics
        if (params.bam_csi_index) {
            ch_genome_bam_index = MARK_DUPLICATES_PICARD.out.csi
        }
        ch_versions = ch_versions.mix(MARK_DUPLICATES_PICARD.out.versions)
    }
    //
    // MODULE: STRINGTIE
    //
    if (!params.skip_alignment && !params.skip_stringtie) {
params.outdir_custom = "${params.outdir}/stringtie"
if (!params.skip_stringtie) {
    process.ext.args = [ '-v', params.stringtie_ignore_gtf ? '' : '-e' ].join(' ').trim()
}
params.outdir_custom = "${params.outdir}/stringtie"
if (!params.skip_stringtie) {
    process.ext.args = [ '-v', params.stringtie_ignore_gtf ? '' : '-e' ].join(' ').trim()
    publishDir =  path: { "${params.outdir}/${params.aligner}/stringtie" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/stringtie"
if (!params.skip_stringtie) {
    process.ext.args = [ '-v', params.stringtie_ignore_gtf ? '' : '-e' ].join(' ').trim()
    publishDir =  path: { "${params.outdir}/${params.aligner}/stringtie" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/stringtie"
if (!params.skip_stringtie) {
    process.ext.args = [ '-v', params.stringtie_ignore_gtf ? '' : '-e' ].join(' ').trim()
    publishDir =  path: { "${params.outdir}/${params.aligner}/stringtie" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/stringtie"
if (!params.skip_stringtie) {
    process.ext.args = [ '-v', params.stringtie_ignore_gtf ? '' : '-e' ].join(' ').trim()
    publishDir =  path: { "${params.outdir}/${params.aligner}/stringtie" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
        STRINGTIE (
            ch_genome_bam,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(STRINGTIE.out.versions.first())
    }
    //
    // MODULE: Feature biotype QC using featureCounts
    //
    ch_featurecounts_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc && !params.skip_biotype_qc && biotype) {
        PREPARE_GENOME
            .out
            .gtf
            .map { WorkflowRnaseq.biotypeInGtf(it, biotype, log) }
            .set { biotype_in_gtf }
        // Prevent any samples from running if GTF file doesn't have a valid biotype
        ch_genome_bam
            .combine(PREPARE_GENOME.out.gtf)
            .combine(biotype_in_gtf)
            .filter { it[-1] }
            .map { it[0..<it.size()-1] }
            .set { ch_featurecounts }
params.outdir_custom = "${params.outdir}/featurecounts/subread"
if (!params.skip_biotype_qc && params.featurecounts_group_type) {
    process.ext.args = [ '-B -C', params.gencode ? "-g gene_type" : "-g $params.featurecounts_group_type", "-t $params.featurecounts_feature_type" ].join(' ').trim()
}
params.outdir_custom = "${params.outdir}/featurecounts/subread"
if (!params.skip_biotype_qc && params.featurecounts_group_type) {
    process.ext.args = [ '-B -C', params.gencode ? "-g gene_type" : "-g $params.featurecounts_group_type", "-t $params.featurecounts_feature_type" ].join(' ').trim()
    publishDir =  path: { "${params.outdir}/${params.aligner}/featurecounts" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/featurecounts/subread"
if (!params.skip_biotype_qc && params.featurecounts_group_type) {
    process.ext.args = [ '-B -C', params.gencode ? "-g gene_type" : "-g $params.featurecounts_group_type", "-t $params.featurecounts_feature_type" ].join(' ').trim()
    publishDir =  path: { "${params.outdir}/${params.aligner}/featurecounts" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/featurecounts/subread"
if (!params.skip_biotype_qc && params.featurecounts_group_type) {
    process.ext.args = [ '-B -C', params.gencode ? "-g gene_type" : "-g $params.featurecounts_group_type", "-t $params.featurecounts_feature_type" ].join(' ').trim()
    publishDir =  path: { "${params.outdir}/${params.aligner}/featurecounts" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/featurecounts/subread"
if (!params.skip_biotype_qc && params.featurecounts_group_type) {
    process.ext.args = [ '-B -C', params.gencode ? "-g gene_type" : "-g $params.featurecounts_group_type", "-t $params.featurecounts_feature_type" ].join(' ').trim()
    publishDir =  path: { "${params.outdir}/${params.aligner}/featurecounts" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
        SUBREAD_FEATURECOUNTS (
            ch_featurecounts
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())
params.outdir_custom = "${params.outdir}/biotype/custom/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
}
params.outdir_custom = "${params.outdir}/biotype/custom/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/biotype/custom/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/biotype/custom/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/biotype/custom/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
        MULTIQC_CUSTOM_BIOTYPE (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_biotypes_header_multiqc
        )
        ch_featurecounts_multiqc = MULTIQC_CUSTOM_BIOTYPE.out.tsv
        ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())
    }
    //
    // MODULE: Genome-wide coverage with BEDTools
    //
    if (!params.skip_alignment && !params.skip_bigwig) {
params.outdir_custom = "${params.outdir}/genomecov/bedtools"
if (!params.skip_bigwig) {
    process.ext.args = '-split-du'
}
params.outdir_custom = "${params.outdir}/genomecov/bedtools"
if (!params.skip_bigwig) {
    process.ext.args = '-split-du'
    publishDir =  path: { "${params.outdir}/bedtools/${meta.id}" }, enabled: false
}
params.outdir_custom = "${params.outdir}/genomecov/bedtools"
if (!params.skip_bigwig) {
    process.ext.args = '-split-du'
    publishDir =  path: { "${params.outdir}/bedtools/${meta.id}" }, enabled: false
}
params.outdir_custom = "${params.outdir}/genomecov/bedtools"
if (!params.skip_bigwig) {
    process.ext.args = '-split-du'
    publishDir =  path: { "${params.outdir}/bedtools/${meta.id}" }, enabled: false
}
params.outdir_custom = "${params.outdir}/genomecov/bedtools"
if (!params.skip_bigwig) {
    process.ext.args = '-split-du'
    publishDir =  path: { "${params.outdir}/bedtools/${meta.id}" }, enabled: false
}
        BEDTOOLS_GENOMECOV (
            ch_genome_bam
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())
        //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
params.outdir_custom = "${params.outdir}/forward/bigwig/to/bedgraph"
params.outdir_custom = "${params.outdir}/forward/bigwig/to/bedgraph"
params.outdir_custom = "${params.outdir}/forward/bigwig/to/bedgraph"
params.outdir_custom = "${params.outdir}/forward/bigwig/to/bedgraph"
params.outdir_custom = "${params.outdir}/forward/bigwig/to/bedgraph"
        BEDGRAPH_TO_BIGWIG_FORWARD (
            BEDTOOLS_GENOMECOV.out.bedgraph_forward,
            PREPARE_GENOME.out.chrom_sizes
        )
        ch_versions = ch_versions.mix(BEDGRAPH_TO_BIGWIG_FORWARD.out.versions)
params.outdir_custom = "${params.outdir}/reverse/bigwig/to/bedgraph"
params.outdir_custom = "${params.outdir}/reverse/bigwig/to/bedgraph"
params.outdir_custom = "${params.outdir}/reverse/bigwig/to/bedgraph"
params.outdir_custom = "${params.outdir}/reverse/bigwig/to/bedgraph"
params.outdir_custom = "${params.outdir}/reverse/bigwig/to/bedgraph"
        BEDGRAPH_TO_BIGWIG_REVERSE (
            BEDTOOLS_GENOMECOV.out.bedgraph_reverse,
            PREPARE_GENOME.out.chrom_sizes
        )
    }
    //
    // MODULE: Downstream QC steps
    //
    ch_qualimap_multiqc           = Channel.empty()
    ch_dupradar_multiqc           = Channel.empty()
    ch_bamstat_multiqc            = Channel.empty()
    ch_inferexperiment_multiqc    = Channel.empty()
    ch_innerdistance_multiqc      = Channel.empty()
    ch_junctionannotation_multiqc = Channel.empty()
    ch_junctionsaturation_multiqc = Channel.empty()
    ch_readdistribution_multiqc   = Channel.empty()
    ch_readduplication_multiqc    = Channel.empty()
    ch_fail_strand_multiqc        = Channel.empty()
    ch_tin_multiqc                = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc) {
        if (!params.skip_qualimap) {
params.outdir_custom = "${params.outdir}/rnaseq/qualimap"
params.outdir_custom = "${params.outdir}/rnaseq/qualimap"
if (!params.skip_qualimap) {
    publishDir =  path: { "${params.outdir}/${params.aligner}/qualimap" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/rnaseq/qualimap"
if (!params.skip_qualimap) {
    publishDir =  path: { "${params.outdir}/${params.aligner}/qualimap" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/rnaseq/qualimap"
if (!params.skip_qualimap) {
    publishDir =  path: { "${params.outdir}/${params.aligner}/qualimap" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/rnaseq/qualimap"
if (!params.skip_qualimap) {
    publishDir =  path: { "${params.outdir}/${params.aligner}/qualimap" }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
            QUALIMAP_RNASEQ (
                ch_genome_bam,
                PREPARE_GENOME.out.gtf
            )
            ch_qualimap_multiqc = QUALIMAP_RNASEQ.out.results
            ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions.first())
        }
        if (!params.skip_dupradar) {
params.outdir_custom = "${params.outdir}/dupradar"
params.outdir_custom = "${params.outdir}/dupradar"
if (!params.skip_dupradar) {
    publishDir =   path: { "${params.outdir}/${params.aligner}/dupradar/scatter_plot" }, mode: params.publish_dir_mode, pattern: "*Dens.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/box_plot" }, mode: params.publish_dir_mode, pattern: "*Boxplot.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/histogram" }, mode: params.publish_dir_mode, pattern: "*Hist.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/gene_data" }, mode: params.publish_dir_mode, pattern: "*Matrix.txt" ],  path: { "${params.outdir}/${params.aligner}/dupradar/intercepts_slope" }, mode: params.publish_dir_mode, pattern: "*slope.txt"
}
params.outdir_custom = "${params.outdir}/dupradar"
if (!params.skip_dupradar) {
    publishDir =   path: { "${params.outdir}/${params.aligner}/dupradar/scatter_plot" }, mode: params.publish_dir_mode, pattern: "*Dens.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/box_plot" }, mode: params.publish_dir_mode, pattern: "*Boxplot.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/histogram" }, mode: params.publish_dir_mode, pattern: "*Hist.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/gene_data" }, mode: params.publish_dir_mode, pattern: "*Matrix.txt" ],  path: { "${params.outdir}/${params.aligner}/dupradar/intercepts_slope" }, mode: params.publish_dir_mode, pattern: "*slope.txt"
}
params.outdir_custom = "${params.outdir}/dupradar"
if (!params.skip_dupradar) {
    publishDir =   path: { "${params.outdir}/${params.aligner}/dupradar/scatter_plot" }, mode: params.publish_dir_mode, pattern: "*Dens.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/box_plot" }, mode: params.publish_dir_mode, pattern: "*Boxplot.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/histogram" }, mode: params.publish_dir_mode, pattern: "*Hist.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/gene_data" }, mode: params.publish_dir_mode, pattern: "*Matrix.txt" ],  path: { "${params.outdir}/${params.aligner}/dupradar/intercepts_slope" }, mode: params.publish_dir_mode, pattern: "*slope.txt"
}
params.outdir_custom = "${params.outdir}/dupradar"
if (!params.skip_dupradar) {
    publishDir =   path: { "${params.outdir}/${params.aligner}/dupradar/scatter_plot" }, mode: params.publish_dir_mode, pattern: "*Dens.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/box_plot" }, mode: params.publish_dir_mode, pattern: "*Boxplot.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/histogram" }, mode: params.publish_dir_mode, pattern: "*Hist.pdf" ],  path: { "${params.outdir}/${params.aligner}/dupradar/gene_data" }, mode: params.publish_dir_mode, pattern: "*Matrix.txt" ],  path: { "${params.outdir}/${params.aligner}/dupradar/intercepts_slope" }, mode: params.publish_dir_mode, pattern: "*slope.txt"
}
            DUPRADAR (
                ch_genome_bam,
                PREPARE_GENOME.out.gtf
            )
            ch_dupradar_multiqc = DUPRADAR.out.multiqc
            ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())
        }
        if (!params.skip_rseqc && rseqc_modules.size() > 0) {
params.outdir_custom = "${params.outdir}/rseqc"
params.outdir_custom = "${params.outdir}/rseqc"
params.outdir_custom = "${params.outdir}/rseqc"
params.outdir_custom = "${params.outdir}/rseqc"
params.outdir_custom = "${params.outdir}/rseqc"
            RSEQC (
                ch_genome_bam,
                ch_genome_bam_index,
                PREPARE_GENOME.out.gene_bed,
                rseqc_modules
            )
            ch_bamstat_multiqc            = RSEQC.out.bamstat_txt
            ch_inferexperiment_multiqc    = RSEQC.out.inferexperiment_txt
            ch_innerdistance_multiqc      = RSEQC.out.innerdistance_freq
            ch_junctionannotation_multiqc = RSEQC.out.junctionannotation_log
            ch_junctionsaturation_multiqc = RSEQC.out.junctionsaturation_rscript
            ch_readdistribution_multiqc   = RSEQC.out.readdistribution_txt
            ch_readduplication_multiqc    = RSEQC.out.readduplication_pos_xls
            ch_tin_multiqc                = RSEQC.out.tin_txt
            ch_versions = ch_versions.mix(RSEQC.out.versions)
            ch_inferexperiment_multiqc
                .map { meta, strand_log -> [ meta ] + WorkflowRnaseq.getInferexperimentStrandedness(strand_log, 30) }
                .filter { it[0].strandedness != it[1] }
                .map { meta, strandedness, sense, antisense, undetermined ->
                    [ "$meta.id\t$meta.strandedness\t$strandedness\t$sense\t$antisense\t$undetermined" ]
                }
                .set { ch_fail_strand }
            def header = [
                "Sample",
                "Provided strandedness",
                "Inferred strandedness",
                "Sense (%)",
                "Antisense (%)",
                "Undetermined (%)"
            ]
params.outdir_custom = "${params.outdir}/check/strand/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
}
params.outdir_custom = "${params.outdir}/check/strand/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/check/strand/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/check/strand/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/check/strand/tsv/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
            MULTIQC_TSV_STRAND_CHECK (
                ch_fail_strand.collect(),
                header,
                'fail_strand_check'
            )
            .set { ch_fail_strand_multiqc }
        }
    }
    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //
    ch_salmon_multiqc                   = Channel.empty()
    ch_pseudoaligner_pca_multiqc        = Channel.empty()
    ch_pseudoaligner_clustering_multiqc = Channel.empty()
    if (params.pseudo_aligner == 'salmon') {
params.outdir_custom = "${params.outdir}/salmon/quantify"
params.outdir_custom = "${params.outdir}/salmon/quantify"
params.outdir_custom = "${params.outdir}/salmon/quantify"
params.outdir_custom = "${params.outdir}/salmon/quantify"
params.outdir_custom = "${params.outdir}/salmon/quantify"
        QUANTIFY_SALMON (
            ch_filtered_reads,
            PREPARE_GENOME.out.salmon_index,
            ch_dummy_file,
            PREPARE_GENOME.out.gtf,
            false,
            params.salmon_quant_libtype ?: ''
        )
        ch_salmon_multiqc = QUANTIFY_SALMON.out.results
        ch_versions = ch_versions.mix(QUANTIFY_SALMON.out.versions)
        if (!params.skip_qc & !params.skip_deseq2_qc) {
params.outdir_custom = "${params.outdir}/salmon/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'salmon'
}
params.outdir_custom = "${params.outdir}/salmon/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'salmon'
    publishDir =  path: { "${params.outdir}/${params.pseudo_aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
params.outdir_custom = "${params.outdir}/salmon/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'salmon'
    publishDir =  path: { "${params.outdir}/${params.pseudo_aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
params.outdir_custom = "${params.outdir}/salmon/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'salmon'
    publishDir =  path: { "${params.outdir}/${params.pseudo_aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
params.outdir_custom = "${params.outdir}/salmon/qc/deseq2"
if (!params.skip_qc & !params.skip_deseq2_qc) {
    process.ext.args = [ "--id_col 1", "--sample_suffix ''", "--outprefix deseq2", "--count_col 3", params.deseq2_vst ? '--vst TRUE' : '' ].join(' ').trim()
    process.ext.args2 = 'salmon'
    publishDir =  path: { "${params.outdir}/${params.pseudo_aligner}/deseq2_qc" }, mode: params.publish_dir_mode, pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log}"
}
            DESEQ2_QC_SALMON (
                QUANTIFY_SALMON.out.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_pseudoaligner_pca_multiqc        = DESEQ2_QC_SALMON.out.pca_multiqc
            ch_pseudoaligner_clustering_multiqc = DESEQ2_QC_SALMON.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_SALMON.out.versions)
        }
    }
    //
    // MODULE: Pipeline reporting
    //
params.outdir_custom = "${params.outdir}/dumpsoftwareversions/custom"
params.outdir_custom = "${params.outdir}/dumpsoftwareversions/custom"
publishDir =  path: { "${params.outdir}/pipeline_info" }, mode: params.publish_dir_mode, pattern: '*_versions.yml'
params.outdir_custom = "${params.outdir}/dumpsoftwareversions/custom"
publishDir =  path: { "${params.outdir}/pipeline_info" }, mode: params.publish_dir_mode, pattern: '*_versions.yml'
params.outdir_custom = "${params.outdir}/dumpsoftwareversions/custom"
publishDir =  path: { "${params.outdir}/pipeline_info" }, mode: params.publish_dir_mode, pattern: '*_versions.yml'
params.outdir_custom = "${params.outdir}/dumpsoftwareversions/custom"
publishDir =  path: { "${params.outdir}/pipeline_info" }, mode: params.publish_dir_mode, pattern: '*_versions.yml'
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowRnaseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)
params.outdir_custom = "${params.outdir}/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
}
params.outdir_custom = "${params.outdir}/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
params.outdir_custom = "${params.outdir}/multiqc"
if (!params.skip_multiqc) {
    process.ext.args = params.multiqc_title?"--title\"$params.multiqc_title\"":''
    publishDir =  path: {  params.skip_alignment? '' : "/${params.aligner}" ].join('') }, mode: params.publish_dir_mode, saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
}
        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_fail_mapping_multiqc.ifEmpty([]),
            ch_fail_strand_multiqc.ifEmpty([]),
            FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
            ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
            ch_star_multiqc.collect{it[1]}.ifEmpty([]),
            ch_hisat2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_rsem_multiqc.collect{it[1]}.ifEmpty([]),
            ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
            ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]),
            ch_aligner_pca_multiqc.collect().ifEmpty([]),
            ch_aligner_clustering_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_pca_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_clustering_multiqc.collect().ifEmpty([]),
            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
            ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
            ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
            ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readduplication_multiqc.collect{it[1]}.ifEmpty([]),
            ch_tin_multiqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, fail_percent_mapped)
    }
    NfcoreTemplate.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
