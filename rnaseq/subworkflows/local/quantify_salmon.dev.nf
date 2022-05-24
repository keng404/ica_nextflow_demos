//
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    publishDir path: { "${params.outdir_custom}" },mode: "${params.publish_dir_mode}",saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
// Pseudo-alignment and quantification with Salmon
//
include { SALMON_QUANT    } from '../../modules/nf-core/modules/salmon/quant/main'
include { SALMON_TX2GENE  } from '../../modules/local/salmon_tx2gene'
include { SALMON_TXIMPORT } from '../../modules/local/salmon_tximport'
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE               } from '../../modules/local/salmon_summarizedexperiment'
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE_LENGTH_SCALED } from '../../modules/local/salmon_summarizedexperiment'
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_GENE_SCALED        } from '../../modules/local/salmon_summarizedexperiment'
include { SALMON_SUMMARIZEDEXPERIMENT as SALMON_SE_TRANSCRIPT         } from '../../modules/local/salmon_summarizedexperiment'
workflow QUANTIFY_SALMON {
    take:
    reads            // channel: [ val(meta), [ reads ] ]
    index            // channel: /path/to/salmon/index/
    transcript_fasta // channel: /path/to/transcript.fasta
    gtf              // channel: /path/to/genome.gtf
    alignment_mode   //    bool: Run Salmon in alignment mode
    lib_type         //     val: String to override salmon library type
    main:
    ch_versions = Channel.empty()
    //
    // Quantify and merge counts across samples
    //
params.outdir_custom = "${params.outdir}/quant/salmon"
params.outdir_custom = "${params.outdir}/quant/salmon"
params.outdir_custom = "${params.outdir}/quant/salmon"
params.outdir_custom = "${params.outdir}/quant/salmon"
    SALMON_QUANT ( reads, index, gtf, transcript_fasta, alignment_mode, lib_type )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
params.outdir_custom = "${params.outdir}/tx2gene/salmon"
params.outdir_custom = "${params.outdir}/tx2gene/salmon"
params.outdir_custom = "${params.outdir}/tx2gene/salmon"
params.outdir_custom = "${params.outdir}/tx2gene/salmon"
    SALMON_TX2GENE ( SALMON_QUANT.out.results.collect{it[1]}, gtf )
    ch_versions = ch_versions.mix(SALMON_TX2GENE.out.versions)
params.outdir_custom = "${params.outdir}/tximport/salmon"
params.outdir_custom = "${params.outdir}/tximport/salmon"
params.outdir_custom = "${params.outdir}/tximport/salmon"
params.outdir_custom = "${params.outdir}/tximport/salmon"
    SALMON_TXIMPORT ( SALMON_QUANT.out.results.collect{it[1]}, SALMON_TX2GENE.out.tsv.collect() )
    ch_versions = ch_versions.mix(SALMON_TXIMPORT.out.versions)
params.outdir_custom = "${params.outdir}/gene/se/salmon"
params.outdir_custom = "${params.outdir}/gene/se/salmon"
params.outdir_custom = "${params.outdir}/gene/se/salmon"
params.outdir_custom = "${params.outdir}/gene/se/salmon"
    SALMON_SE_GENE (
        SALMON_TXIMPORT.out.counts_gene,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.tsv.collect()
    )
    ch_versions = ch_versions.mix(SALMON_SE_GENE.out.versions)
params.outdir_custom = "${params.outdir}/scaled/length/gene/se/salmon"
params.outdir_custom = "${params.outdir}/scaled/length/gene/se/salmon"
params.outdir_custom = "${params.outdir}/scaled/length/gene/se/salmon"
params.outdir_custom = "${params.outdir}/scaled/length/gene/se/salmon"
    SALMON_SE_GENE_LENGTH_SCALED (
        SALMON_TXIMPORT.out.counts_gene_length_scaled,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.tsv.collect()
    )
params.outdir_custom = "${params.outdir}/scaled/gene/se/salmon"
params.outdir_custom = "${params.outdir}/scaled/gene/se/salmon"
params.outdir_custom = "${params.outdir}/scaled/gene/se/salmon"
params.outdir_custom = "${params.outdir}/scaled/gene/se/salmon"
    SALMON_SE_GENE_SCALED (
        SALMON_TXIMPORT.out.counts_gene_scaled,
        SALMON_TXIMPORT.out.tpm_gene,
        SALMON_TX2GENE.out.tsv.collect()
    )
params.outdir_custom = "${params.outdir}/transcript/se/salmon"
params.outdir_custom = "${params.outdir}/transcript/se/salmon"
params.outdir_custom = "${params.outdir}/transcript/se/salmon"
params.outdir_custom = "${params.outdir}/transcript/se/salmon"
    SALMON_SE_TRANSCRIPT (
        SALMON_TXIMPORT.out.counts_transcript,
        SALMON_TXIMPORT.out.tpm_transcript,
        SALMON_TX2GENE.out.tsv.collect()
    )
    emit:
    results                       = SALMON_QUANT.out.results                      // channel: [ val(meta), results_dir ]
    tpm_gene                      = SALMON_TXIMPORT.out.tpm_gene                  // channel: [ val(meta), counts ]
    counts_gene                   = SALMON_TXIMPORT.out.counts_gene               // channel: [ val(meta), counts ]
    counts_gene_length_scaled     = SALMON_TXIMPORT.out.counts_gene_length_scaled // channel: [ val(meta), counts ]
    counts_gene_scaled            = SALMON_TXIMPORT.out.counts_gene_scaled        // channel: [ val(meta), counts ]
    tpm_transcript                = SALMON_TXIMPORT.out.tpm_transcript            // channel: [ val(meta), counts ]
    counts_transcript             = SALMON_TXIMPORT.out.counts_transcript         // channel: [ val(meta), counts ]
    merged_gene_rds               = SALMON_SE_GENE.out.rds                        //    path: *.rds
    merged_gene_rds_length_scaled = SALMON_SE_GENE_LENGTH_SCALED.out.rds          //    path: *.rds
    merged_gene_rds_scaled        = SALMON_SE_GENE_SCALED.out.rds                 //    path: *.rds
    merged_counts_transcript      = SALMON_TXIMPORT.out.counts_transcript         //    path: *.transcript_counts.tsv
    merged_tpm_transcript         = SALMON_TXIMPORT.out.tpm_transcript            //    path: *.transcript_tpm.tsv
    merged_transcript_rds         = SALMON_SE_TRANSCRIPT.out.rds                  //    path: *.rds
    versions                      = ch_versions                                   // channel: [ versions.yml ]
}
