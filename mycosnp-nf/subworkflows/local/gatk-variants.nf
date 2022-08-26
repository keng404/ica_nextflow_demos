process = [:]
process.ext = [:]
params.max_amb_samples = "10000000"
params.max_perc_amb_samples = "10"
params.min_depth = "50"
/*
========================================================================================
    GATK Variants Sub-Workflow
========================================================================================
*/
include { GATK4_GENOTYPEGVCFS                     } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'
include { GATK4_VARIANTFILTRATION                 } from '../../modules/nf-core/modules/gatk4/variantfiltration/main'
include { GATK4_SELECTVARIANTS                    } from '../../modules/nf-core/modules/gatk4/selectvariants/main'
include { FILTER_GATK_GENOTYPES                   } from '../../modules/local/vcftools.nf'
include { BCFTOOLS_INDEX                          } from '../../modules/nf-core/modules/bcftools/index/main'
include { BCFTOOLS_VIEW  as BCFTOOLS_VIEW_CONVERT } from '../../modules/nf-core/modules/bcftools/view/main'
include { SPLIT_VCF                               } from '../../modules/local/splitvcf.nf'
include { VCF_TO_FASTA                            } from '../../modules/local/vcftofasta.nf'
include { VCF_QC                                  } from '../../modules/local/vcfqc.nf'
include { BCFTOOLS_QUERY                          } from '../../modules/nf-core/modules/bcftools/query/main'
include { VCF_CONSENSUS                           } from '../../modules/local/vcfconsensus.nf'
combined_gvcf = Channel.empty()
workflow GATK_VARIANTS {
    take:
    fasta
    fai
    bai
    dict
    thismeta
    vcffile
    vcfidx
    main:
    ch_versions = Channel.empty()
    combined_gvcf = thismeta.combine(vcffile).combine(vcfidx)
    process.ext.prefix = "combined_genotype"
    GATK4_GENOTYPEGVCFS(
        combined_gvcf.map{meta, vcf, idx -> [ meta, vcf, idx, [] ]},
        fasta, 
        fai, 
        dict, [], []
        )
    process.ext.args =  "--filter-expression '${params.gvcfs_filter}' --filter-name filter"
    process.ext.prefix = "combined_genotype_filtered"
    GATK4_VARIANTFILTRATION(
                            GATK4_GENOTYPEGVCFS.out.vcf.combine(GATK4_GENOTYPEGVCFS.out.tbi).map{ meta1, vcf, meta2, tbi->[meta1, vcf, tbi]},
                            fasta, 
                            fai, 
                            dict
                            )
    process.ext.args =  '--select-type-to-include "SNP"'
    process.ext.prefix = "combined_genotype_filtered_snps"
    GATK4_SELECTVARIANTS(
                         GATK4_VARIANTFILTRATION.out.vcf.combine(GATK4_VARIANTFILTRATION.out.tbi).map{ meta1, vcf, meta2, tbi->[meta1, vcf, tbi]}
                        )
    process.ext.args = params.gatkgenotypes_filter
    process.ext.prefix = "combined_genotype_filtered_snps_filtered"
    FILTER_GATK_GENOTYPES(GATK4_SELECTVARIANTS.out.vcf)
    fin_comb_vcf = FILTER_GATK_GENOTYPES.out.vcf.first()
    process.ext.args = "-Oz"
    process.ext.prefix = "finalfiltered"
    BCFTOOLS_VIEW_CONVERT(fin_comb_vcf.map{meta, vcf->[ meta, vcf, [] ] }, [], [], []  )   // Convert to bgzip
    BCFTOOLS_INDEX(BCFTOOLS_VIEW_CONVERT.out.vcf)
    SPLIT_VCF(
              BCFTOOLS_VIEW_CONVERT.out.vcf.combine(BCFTOOLS_INDEX.out.csi).map{meta1, vcf, meta2, csi-> [meta1, vcf, csi] }
             )
    final_vcf_txt = Channel.empty()
    fin_comb_vcf.combine(SPLIT_VCF.out.txt).map{meta1, vcf, meta2, txt -> 
                    [ meta1, vcf, txt, params.max_amb_samples, params.max_perc_amb_samples, params.min_depth ]}.set{final_vcf_txt}
    VCF_CONSENSUS(
        BCFTOOLS_VIEW_CONVERT.out.vcf.combine(BCFTOOLS_INDEX.out.csi).map{meta1, vcf, meta2, csi-> [meta1, vcf, csi] },
        fasta
    )
    VCF_TO_FASTA(final_vcf_txt, fasta)
    VCF_QC(VCF_TO_FASTA.out.fasta.map{meta, fa-> fa })
    // Collect versions information
    ch_versions = ch_versions.mix(  GATK4_GENOTYPEGVCFS.out.versions, 
                                    GATK4_VARIANTFILTRATION.out.versions, 
                                    GATK4_SELECTVARIANTS.out.versions,
                                    BCFTOOLS_VIEW_CONVERT.out.versions,
                                    BCFTOOLS_INDEX.out.versions,
                                    SPLIT_VCF.out.versions,
                                    VCF_CONSENSUS.out.versions
                                )
    emit:
    snps_fasta = VCF_TO_FASTA.out.fasta  // channel: [ val(meta), fasta ]
    versions = ch_versions               // channel: [ versions.yml ]
    // filtered_vcf    = BROAD_VCFFILTER.out
    // split_vcf_broad = SPLITVCF.out        --> the broad vcf file
    // variants        = GATK4_SELECTVARIANTS.out
    // split_vcf_gatk4 = SPLITVCF.out        --> the gatk4 vcf file
    // consensus_fasta = BCFTOOLS_CONSENSUS.out
    // qc_report       = VCF_QCREPORT.out
    
}
