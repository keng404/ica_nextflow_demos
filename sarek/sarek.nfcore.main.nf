#!/usr/bin/env nextflow
// Global default params, used in configs

// Workflow flags:

// Mandatory arguments
params.input = null // No default input TSV file that defines input files and associated metadata
params.step = 'mapping' // Starts with mapping
reads = Channel.fromPath(params.reads) // dummy variable

println params.reads

// Genome and references options
params.genome = 'GRCh38'
params.igenomes_base = 's3://ngi-igenomes/igenomes/'
params.igenomes_ignore = false
params.genomes_base = null // Disabled by default
params.save_reference = null // Built references not saved

// Main options
params.help = false
params.no_intervals = null // Intervals will be built from the fasta file
params.nucleotides_per_second = 1000 // Default interval size
params.sentieon = null // Not using Sentieon by default
params.skip_qc = null // All QC tools are used
params.target_bed = false // No default TargetBED file for targeted sequencing
params.tools = null // No default Variant_Calling or Annotation tools

// Modify fastqs (trim/split)
params.trim_fastq = false // No trimming
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0
params.trim_nextseq = 0
params.save_trimmed = true
params.split_fastq = null // Fastq files will not be split by default

// Preprocessing
params.aligner = 'bwa-mem'
params.markdup_java_options = '"-Xms4000m -Xmx7g"' // Established values for markDuplicates memory consumption, see https://github.com/SciLifeLab/Sarek/pull/689 for details
params.use_gatk_spark = false // GATK Spark implementation of their tools in local mode not used by default
params.save_bam_mapped = null // Mapped BAMs not saved
params.skip_markduplicates = null // Do not skip markDuplicates by default

// Variant Calling
params.ascat_ploidy = null // Use default value
params.ascat_purity = null // Use default value
params.cf_coeff = 0.05                    // default value for Control-FREEC
params.cf_contamination = null            // by default not specified in Control-FREEC
params.cf_contamination_adjustment = null // by default we are not using this in Control-FREEC
params.cf_ploidy = 2                      // you can use 2,3,4
params.cf_window = null                   // by default we are not using this in Control-FREEC
params.generate_gvcf = null // g.vcf are not produced by HaplotypeCaller by default
params.no_strelka_bp = null // Strelka will use Manta candidateSmallIndels if available
params.pon = false // No default PON (Panel of Normals) file for GATK Mutect2 / Sentieon TNscope
params.pon_index = false // No default PON index for GATK Mutect2 / Sentieon TNscope
params.ignore_soft_clipped_bases = null // no --dont-use-soft-clipped-bases for GATK Mutect2
params.umi = null // no umi
params.read_structure1 = null // no umi
params.read_structure2 = null // no umi

// Annotation
params.annotate_tools = null // Only with --step annotate
params.annotation_cache = null // Annotation cache disabled
params.cadd_cache = null // CADD cache disabled
params.cadd_indels = false // No CADD InDels file
params.cadd_indels_tbi = false // No CADD InDels index
params.cadd_wg_snvs = false // No CADD SNVs file
params.cadd_wg_snvs_tbi = false // No CADD SNVs index
params.genesplicer = null // genesplicer disabled within VEP
params.snpeff_cache = null // No directory for snpEff cache
params.vep_cache = null // No directory for VEP cache

// Custom config
params.config_profile_contact = false
params.config_profile_description = false
params.config_profile_url = false
params.custom_config_version = 'master'
params.custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"

// Other options
params.outdir = "${workflow.launchDir}/out"
params.publish_dir_mode = 'copy' // Default PublishDirMode (same as other nf-core pipelines)
params.sequencing_center = null // No sequencing center to be written in BAM header in MapReads process
params.multiqc_config = false
params.monochrome_logs = true // Monochrome logs disabled
params.email = "" // No default email
params.email_on_fail = false
params.plaintext_email = false // Plaintext email disabled
params.max_multiqc_email_size = 25.MB

params.hostnames = false
params.config_profile_name = null
params.config_profile_description = false
params.config_profile_contact = false
params.config_profile_url = false
params.validate_params = true
params.show_hidden_params = false
params.schema_ignore_params = 'genomes,input_paths'
params.tracedir = "${params.outdir}/pipeline_info"

// Base specifications
// Defaults only, expecting to be overwritten
params.cpus = 24 
params.max_cpus = 48
params.max_memory = 512.GB
params.max_time = 240.h
params.single_cpu_mem = 10.GB


// Container slugs
// Stable releases should specify release tag (ie: `2.5.2`)
// Developmental code should specify dev
pipeline_container = 'nfcore/sarek:2.7.1'
instance_size = 'himem-medium'
snpeff_container = {(params.annotation_cache && params.snpeff_cache) ? "${pipeline_container}" : "nfcore/sareksnpeff:2.7.1.${params.genome}"}
vep_container = {(params.annotation_cache && params.vep_cache) ? "${pipeline_container}" : "nfcore/sarekvep:2.7.1.${params.genome}"}
msisensor_container ="quay.io/biocontainers/msisensor-pro:1.1.a--hb3646a4_0" 
//////////


// illumina iGenomes reference file paths
params.genomes = [:]

// GRCh37
params.genomes['GRCh37']  = [:]
params.genomes['GRCh37'].ac_loci                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/ASCAT/1000G_phase3_20130502_SNP_maf0.3.loci"
params.genomes['GRCh37'].ac_loci_gc              = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/ASCAT/1000G_phase3_20130502_SNP_maf0.3.loci.gc"
params.genomes['GRCh37'].bwa                     = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Sequence/BWAIndex/human_g1k_v37_decoy.fasta.{amb,ann,bwt,pac,sa}"
params.genomes['GRCh37'].chr_dir                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Sequence/Chromosomes"
params.genomes['GRCh37'].chr_length              = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Sequence/Length/human_g1k_v37_decoy.len"
params.genomes['GRCh37'].dbsnp                   = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf"
params.genomes['GRCh37'].dbsnp_index             = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf.idx"
params.genomes['GRCh37'].dict                    = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict"
params.genomes['GRCh37'].fasta                   = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"
params.genomes['GRCh37'].fasta_fai               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai"
params.genomes['GRCh37'].germline_resource       = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh37.PASS.AC.AF.only.vcf.gz"
params.genomes['GRCh37'].germline_resource_index = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh37.PASS.AC.AF.only.vcf.gz.tbi"
params.genomes['GRCh37'].intervals               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/intervals/wgs_calling_regions_Sarek.list"
params.genomes['GRCh37'].known_indels            = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/{1000G_phase1,Mills_and_1000G_gold_standard}.indels.b37.vcf"
params.genomes['GRCh37'].known_indels_index      = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/{1000G_phase1,Mills_and_1000G_gold_standard}.indels.b37.vcf.idx"
params.genomes['GRCh37'].mappability             = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh37/Annotation/Control-FREEC/out100m2_hg19.gem"
params.genomes['GRCh37'].snpeff_db               = 'GRCh37.75'
params.genomes['GRCh37'].species                 = 'homo_sapiens'
params.genomes['GRCh37'].vep_cache_version       = '99'

// Ensembl.GRCh37
params.genomes['Ensembl.GRCh37']  = [:]
params.genomes['Ensembl.GRCh37'] .fasta                   = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Ensembl.GRCh37'] .bwa                     = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"

// GRCh38
params.genomes['GRCh38']  = [:]
params.genomes['GRCh38'].ac_loci                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/1000G_phase3_GRCh38_maf0.3.loci"
params.genomes['GRCh38'].ac_loci_gc              = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/1000G_phase3_GRCh38_maf0.3.loci.gc"
params.genomes['GRCh38'].bwa                     = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.{alt,amb,ann,bwt,pac,sa}"
params.genomes['GRCh38'].chr_dir                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes"
params.genomes['GRCh38'].chr_length              = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/Length/Homo_sapiens_assembly38.len"
params.genomes['GRCh38'].bsnp                   = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz"
params.genomes['GRCh38'].dbsnp_index             = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz.tbi"
params.genomes['GRCh38'].dict                    = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict"
params.genomes['GRCh38'].fasta                   = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.genomes['GRCh38'].fasta_fai               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai"
params.genomes['GRCh38'].germline_resource       = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz"
params.genomes['GRCh38'].germline_resource_index = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz.tbi"
params.genomes['GRCh38'].intervals               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/intervals/wgs_calling_regions.hg38.bed"
params.genomes['GRCh38'].known_indels            = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz"
params.genomes['GRCh38'].known_indels_index      = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz.tbi"
params.genomes['GRCh38'].mappability             = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/Control-FREEC/out100m2_hg38.gem"
params.genomes['GRCh38'].snpeff_db               = 'GRCh38.86'
params.genomes['GRCh38'].species                 = 'homo_sapiens'
params.genomes['GRCh38'].vep_cache_version       = '99'

// NCBI.GRCh38
params.genomes['NCBI.GRCh38']  = [:]
params.genomes['NCBI.GRCh38'].fasta                   = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['NCBI.GRCh38'].bwa                     = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"

// GRCm38
params.genomes['GRCm38']  = [:]
params.genomes['GRCm38'].bwa                     = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['GRCm38'].chr_dir                 = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/Chromosomes"
params.genomes['GRCm38'].chr_length              = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/Length/GRCm38.len"
params.genomes['GRCm38'].dbsnp                   = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
params.genomes['GRCm38'].dbsnp_index             = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.snps_all.dbSNP142.vcf.gz.tbi"
params.genomes['GRCm38'].dict                    = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.dict"
params.genomes['GRCm38'].fasta                   = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['GRCm38'].fasta_fai               = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa.fai"
params.genomes['GRCm38'].intervals               = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Annotation/intervals/GRCm38_calling_list.bed"
params.genomes['GRCm38'].known_indels            = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz"
params.genomes['GRCm38'].known_indels_index      = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz.tbi"
params.genomes['GRCm38'].mappability             = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Annotation/Control-FREEC/GRCm38_68_mm10.gem"
params.genomes['GRCm38'].snpeff_db               = 'GRCm38.86'
params.genomes['GRCm38'].species                 = 'mus_musculus'
params.genomes['GRCm38'].vep_cache_version       = '99'

// TAIR10
params.genomes['TAIR10']  = [:]
params.genomes['TAIR10'].bwa                     = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['TAIR10'].fasta                   = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa"

// EB2
params.genomes['EB2']  = [:]
params.genomes['EB2'].bwa                     = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['EB2'].fasta                   = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/WholeGenomeFasta/genome.fa"

// UMD3.1
params.genomes['UMD3.1']  = [:]
params.genomes['UMD3.1'].bwa                     = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['UMD3.1'].fasta                   = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/WholeGenomeFasta/genome.fa"

// bosTau8
params.genomes['bosTau8']  = [:]
params.genomes['bosTau8'].bwa                     = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['bosTau8'].fasta                   = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/WholeGenomeFasta/genome.fa"

// WBcel235
params.genomes['WBcel235']  = [:]
params.genomes['WBcel235'].bwa                     = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['WBcel235'].fasta                   = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['WBcel235'].snpeff_db               = 'WBcel235.86'
params.genomes['WBcel235'].species                 = 'caenorhabditis_elegans'
params.genomes['WBcel235'].vep_cache_version       = '99'

// ce10
params.genomes['ce10']  = [:]
params.genomes['ce10'].bwa                     = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['ce10'].fasta                   = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/WholeGenomeFasta/genome.fa"

// CanFam3.1
params.genomes['CanFam3.1']  = [:]
params.genomes['CanFam3.1'].bwa                     = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['CanFam3.1'].fasta                   = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/WholeGenomeFasta/genome.fa"

// canFam3
params.genomes['canFam3']  = [:]
params.genomes['canFam3'].bwa                     = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['canFam3'].fasta                   = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/WholeGenomeFasta/genome.fa"

// GRCz10
params.genomes['GRCz10']  = [:]
params.genomes['GRCz10'].bwa                     = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['GRCz10'].fasta                   = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/WholeGenomeFasta/genome.fa"

// danRer10
params.genomes['danRer10']  = [:]
params.genomes['danRer10'].bwa                     = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['danRer10'].fasta                   = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/WholeGenomeFasta/genome.fa"

// BDGP6
params.genomes['BDGP6']  = [:]
params.genomes['BDGP6'].bwa                     = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['BDGP6'].fasta                   = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/WholeGenomeFasta/genome.fa"

// dm6
params.genomes['dm6']  = [:]
params.genomes['dm6'].bwa                     = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['dm6'].fasta                   = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa"

//  EquCab2
params.genomes['EquCab2']  = [:]
params.genomes['EquCab2'].bwa                     = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['EquCab2'].fasta                   = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/WholeGenomeFasta/genome.fa"

// equCab2
params.genomes['equCab2']  = [:]
params.genomes['equCab2'].bwa                     = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['equCab2'].fasta                   = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/WholeGenomeFasta/genome.fa"

// EB1
params.genomes['EB1']  = [:]
params.genomes['EB1'].bwa                     = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['EB1'].fasta                   = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/WholeGenomeFasta/genome.fa"

// Galgal4
params.genomes['Galgal4']  = [:]
params.genomes['Galgal4'].bwa                     = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['Galgal4'].fasta                   = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/WholeGenomeFasta/genome.fa"

// galGal4
params.genomes['galGal4']  = [:]
params.genomes['galGal4'].bwa                     = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['galGal4'].fasta                   = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/WholeGenomeFasta/genome.fa"

// Gm01
params.genomes['Gm01']  = [:]
params.genomes['Gm01'].bwa                     = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['Gm01'].fasta                   = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/WholeGenomeFasta/genome.fa"

// hg38
params.genomes['hg38']  = [:]
params.genomes['hg38'].bwa                     = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['hg38'].fasta                   = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"

// hg19
params.genomes['hg19']  = [:]
params.genomes['hg19'].bwa                     = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['hg19'].fasta                   = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"

// Mmul_1
params.genomes['Mmul_1']  = [:]
params.genomes['Mmul_1'].bwa                     = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['Mmul_1'].fasta                   = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/WholeGenomeFasta/genome.fa"

// mm10
params.genomes['mm10']  = [:]
params.genomes['mm10'].bwa                     = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['mm10'].fasta                   = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"

// IRGSP-1.0
params.genomes['IRGSP-1.0']  = [:]
params.genomes['IRGSP-1.0'].bwa                     = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['IRGSP-1.0'].fasta                   = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/WholeGenomeFasta/genome.fa"

// CHIMP2.1.4
params.genomes['CHIMP2.1.4']  = [:]
params.genomes['CHIMP2.1.4'].bwa                     = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['CHIMP2.1.4'].fasta                   = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/WholeGenomeFasta/genome.fa"

// panTro4
params.genomes['panTro4']  = [:]
params.genomes['panTro4'].bwa                     = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['panTro4'].fasta                   = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/WholeGenomeFasta/genome.fa"

// Rnor_6.0
params.genomes['Rnor_6.0']  = [:]
params.genomes['Rnor_6.0'].bwa                     = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['Rnor_6.0'].fasta                   = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/WholeGenomeFasta/genome.fa"

// rn6
params.genomes['rn6']  = [:]
params.genomes['rn6'].bwa                     = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['rn6'].fasta                   = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa"

// R64-1-1
params.genomes['R64-1-1']  = [:]
params.genomes['R64-1-1'].bwa                     = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['R64-1-1'].fasta                   = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/WholeGenomeFasta/genome.fa"

// sacCer3
params.genomes['sacCer3']  = [:]
params.genomes['sacCer3'].bwa                     = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['sacCer3'].fasta                   = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa"

// EF2
params.genomes['EF2']  = [:]
params.genomes['EF2'].bwa                     = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['EF2'].fasta                   = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/WholeGenomeFasta/genome.fa"

// Sbi1
params.genomes['Sbi1']  = [:]
params.genomes['Sbi1'].bwa                     = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['Sbi1'].fasta                   = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/WholeGenomeFasta/genome.fa"

// Sscrofa10.2
params.genomes['Sscrofa10.2']  = [:]
params.genomes['Sscrofa10.2'].bwa                     = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['Sscrofa10.2'].fasta                   = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/WholeGenomeFasta/genome.fa"

// susScr3
params.genomes['susScr3']  = [:]
params.genomes['susScr3'].bwa                     = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['susScr3'].fasta                   = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/WholeGenomeFasta/genome.fa"

// AGPv3     
params.genomes['AGPv3']  = [:]
params.genomes['AGPv3'].bwa                     = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/BWAIndex/genome.fa.{amb,ann,bwt,pac,sa}"
params.genomes['AGPv3'].fasta                   = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/WholeGenomeFasta/genome.fa"

// Capture exit codes from upstream processes when piping
// process.shell = ['/bin/bash', '-euo', 'pipefail']

timeline = {}
timeline.enabled = true
timeline.file = "${params.tracedir}/execution_timeline.html"

report = {}
report.enabled = true
report.file = "${params.tracedir}/execution_report.html"

trace = {}
trace.enabled = true
trace.file = "${params.tracedir}/execution_trace.txt"

dag = {}
dag.enabled = true
dag.file = "${params.tracedir}/pipeline_dag.svg"


// Check if genome exists in the config
try {
    params.genomes[params.genome]
} catch(MissingPropertyException){
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. "
}

/*
================================================================================
                                  nf-core/sarek
================================================================================
Started March 2016.
Ported to nf-core May 2019.
--------------------------------------------------------------------------------
nf-core/sarek:
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://nf-co.re/sarek
--------------------------------------------------------------------------------
 @Documentation
 https://nf-co.re/sarek/docs
--------------------------------------------------------------------------------
*/

// log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////
//def json_schema = "$projectDir/nextflow_schema.json"
//if (params.help) {
//    def command = "nextflow run nf-core/sarek --input sample.tsv -profile docker"
//    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
//////    exit 0
// }

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////
//if (params.validate_params) {
//    NfcoreSchema.validateParameters(params, json_schema, log)
//}

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

/*
================================================================================
                         SET UP CONFIGURATION VARIABLES
================================================================================
*/

stepList = defineStepList()
step = params.step ? params.step.toLowerCase().replaceAll('-', '').replaceAll('_', '') : ''

// Handle deprecation
if (step == 'preprocessing') step = 'mapping'

if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (!checkParameterExistence(step, stepList)) exit 1, "Unknown step ${step}, see --help for more information"

toolList = defineToolList()
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
if (step == 'controlfreec') tools = ['controlfreec']
if (!checkParameterList(tools, toolList)) exit 1, 'Unknown tool(s), see --help for more information'

skipQClist = defineSkipQClist()
skipQC = params.skip_qc ? params.skip_qc == 'all' ? skipQClist : params.skip_qc.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
if (!checkParameterList(skipQC, skipQClist)) exit 1, 'Unknown QC tool(s), see --help for more information'

annoList = defineAnnoList()
annotate_tools = params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '')} : []
if (!checkParameterList(annotate_tools,annoList)) exit 1, 'Unknown tool(s) to annotate, see --help for more information'

// Check parameters
if ((params.ascat_ploidy && !params.ascat_purity) || (!params.ascat_ploidy && params.ascat_purity)) exit 1, 'Please specify both --ascat_purity and --ascat_ploidy, or none of them'
if (params.umi && !(params.read_structure1 && params.read_structure2)) exit 1, 'Please specify both --read_structure1 and --read_structure2, when using --umi'

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// MultiQC
// Stage config files
ch_multiqc_config = file("${workflow.launchDir/params.config_files}/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("${workflow.launchDir/params.config_files}/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("${workflow.launchDir/params.config_files}/docs/images/", checkIfExists: true)

// Handle input
tsvPath = null
if (params.input && (hasExtension(params.input, "tsv") || hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) tsvPath = params.input
if (params.input && (hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) step = "annotate"

save_bam_mapped = params.skip_markduplicates ? true : params.save_bam_mapped ? true : false

// If no input file specified, trying to get TSV files corresponding to step in the TSV directory
// only for steps preparerecalibration, recalibrate, variantcalling and controlfreec
if (!params.input && params.sentieon) {
    switch (step) {
        case 'mapping': break
        case 'recalibrate': tsvPath = "${params.outdir}/Preprocessing/TSV/sentieon_deduped.tsv"; break
        case 'variantcalling': tsvPath = "${params.outdir}/Preprocessing/TSV/sentieon_recalibrated.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (!params.input && !params.sentieon && !params.skip_markduplicates) {
    switch (step) {
        case 'mapping': break
        case 'preparerecalibration': tsvPath = "${params.outdir}/Preprocessing/TSV/duplicates_marked_no_table.tsv"; break
        case 'recalibrate': tsvPath = "${params.outdir}/Preprocessing/TSV/duplicates_marked.tsv"; break
        case 'variantcalling': tsvPath = "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"; break
        case 'controlfreec': tsvPath = "${params.outdir}/VariantCalling/TSV/control-freec_mpileup.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (!params.input && !params.sentieon && params.skip_markduplicates) {
    switch (step) {
        case 'mapping': break
        case 'preparerecalibration': tsvPath = "${params.outdir}/Preprocessing/TSV/mapped.tsv"; break
        case 'recalibrate': tsvPath = "${params.outdir}/Preprocessing/TSV/mapped_no_duplicates_marked.tsv"; break
        case 'variantcalling': tsvPath = "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"; break
        case 'controlfreec': tsvPath = "${params.outdir}/VariantCalling/TSV/control-freec_mpileup.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
}

inputSample = Channel.empty()
if (tsvPath) {
    tsvFile = file(tsvPath)
    switch (step) {
        case 'mapping': inputSample = extractFastq(tsvFile); break
        case 'preparerecalibration': inputSample = extractBam(tsvFile); break
        case 'recalibrate': inputSample = extractRecal(tsvFile); break
        case 'variantcalling': inputSample = extractBam(tsvFile); break
        case 'controlfreec': inputSample = extractPileup(tsvFile); break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (params.input && !hasExtension(params.input, "tsv")) {
    log.info "No TSV file"
    if (step != 'mapping') exit 1, 'No step other than "mapping" supports a directory as an input'
    log.info "Reading ${params.input} directory"
    log.warn "[nf-core/sarek] in ${params.input} directory, all fastqs are assuming to be from the same sample, which is assumed to be a germline one"
    inputSample = extractFastqFromDir(params.input)
    (inputSample, fastqTMP) = inputSample.into(2)
    fastqTMP.toList().subscribe onNext: {
        if (it.size() == 0) exit 1, "No FASTQ files found in --input directory '${params.input}'"
    }
    tsvFile = params.input  // used in the reports
} else if (tsvPath && step == 'annotate') {
    log.info "Annotating ${tsvPath}"
} else if (step == 'annotate') {
    log.info "Trying automatic annotation on files in the VariantCalling/ directory"
} else exit 1, 'No sample were defined, see --help'

(genderMap, statusMap, inputSample) = extractInfos(inputSample)

/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
// params.fasta has to be the first one
params.fasta = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta ?: null : null

// The rest can be sorted
params.ac_loci = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci ?: null : null
params.ac_loci_gc = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci_gc ?: null : null
params.bwa = params.genome && params.fasta && 'mapping' in step ? params.genomes[params.genome].bwa ?: null : null
params.chr_dir = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_dir ?: null : null
params.chr_length = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_length ?: null : null
params.dbsnp = params.genome && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || params.sentieon) ? params.genomes[params.genome].dbsnp ?: null : null
params.dbsnp_index = params.genome && params.dbsnp ? params.genomes[params.genome].dbsnp_index ?: null : null
params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
params.fasta_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_fai ?: null : null
params.germline_resource = params.genome && 'mutect2' in tools ? params.genomes[params.genome].germline_resource ?: null : null
params.germline_resource_index = params.genome && params.germline_resource ? params.genomes[params.genome].germline_resource_index ?: null : null
params.intervals = params.genome && !('annotate' in step) ? params.genomes[params.genome].intervals ?: null : null
params.known_indels = params.genome && ('mapping' in step || 'preparerecalibration' in step) ? params.genomes[params.genome].known_indels ?: null : null
params.known_indels_index = params.genome && params.known_indels ? params.genomes[params.genome].known_indels_index ?: null : null
params.mappability = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].mappability ?: null : null
params.snpeff_db = params.genome && ('snpeff' in tools || 'merge' in tools) ? params.genomes[params.genome].snpeff_db ?: null : null
params.species = params.genome && ('vep' in tools || 'merge' in tools) ? params.genomes[params.genome].species ?: null : null
params.vep_cache_version = params.genome && ('vep' in tools || 'merge' in tools) ? params.genomes[params.genome].vep_cache_version ?: null : null

// Initialize channels with files based on params
ch_ac_loci = params.ac_loci && 'ascat' in tools ? Channel.value(file(params.ac_loci)) : "null"
ch_ac_loci_gc = params.ac_loci_gc && 'ascat' in tools ? Channel.value(file(params.ac_loci_gc)) : "null"
ch_chr_dir = params.chr_dir && 'controlfreec' in tools ? Channel.value(file(params.chr_dir)) : "null"
ch_chr_length = params.chr_length && 'controlfreec' in tools ? Channel.value(file(params.chr_length)) : "null"
ch_dbsnp = params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || params.sentieon) ? Channel.value(file(params.dbsnp)) : "null"
ch_fasta = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
ch_fai = params.fasta_fai && !('annotate' in step) ? Channel.value(file(params.fasta_fai)) : "null"
ch_germline_resource = params.germline_resource && 'mutect2' in tools ? Channel.value(file(params.germline_resource)) : "null"
ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"
ch_known_indels = params.known_indels && ('mapping' in step || 'preparerecalibration' in step) ? Channel.value(file(params.known_indels)) : "null"
ch_mappability = params.mappability && 'controlfreec' in tools ? Channel.value(file(params.mappability)) : "null"

// Initialize channels with values based on params
ch_snpeff_cache = params.snpeff_cache ? Channel.value(file(params.snpeff_cache)) : "null"
ch_snpeff_db = params.snpeff_db ? Channel.value(params.snpeff_db) : "null"
ch_vep_cache_version = params.vep_cache_version ? Channel.value(params.vep_cache_version) : "null"
ch_vep_cache = params.vep_cache ? Channel.value(file(params.vep_cache)) : "null"

// Optional files, not defined within the params.genomes[params.genome] scope
ch_cadd_indels = params.cadd_indels ? Channel.value(file(params.cadd_indels)) : "null"
ch_cadd_indels_tbi = params.cadd_indels_tbi ? Channel.value(file(params.cadd_indels_tbi)) : "null"
ch_cadd_wg_snvs = params.cadd_wg_snvs ? Channel.value(file(params.cadd_wg_snvs)) : "null"
ch_cadd_wg_snvs_tbi = params.cadd_wg_snvs_tbi ? Channel.value(file(params.cadd_wg_snvs_tbi)) : "null"
ch_pon = params.pon ? Channel.value(file(params.pon)) : "null"
ch_target_bed = params.target_bed ? Channel.value(file(params.target_bed)) : "null"

// Optional values, not defined within the params.genomes[params.genome] scope
ch_read_structure1 = params.read_structure1 ? Channel.value(params.read_structure1) : "null"
ch_read_structure2 = params.read_structure2 ? Channel.value(params.read_structure2) : "null"

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
//log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}

summary['Input']             = params.input
summary['Step']              = step
summary['Genome']            = params.genome

if (params.no_intervals && step != 'annotate')  summary['Intervals']         = 'Do not use'
summary['Nucleotides/s']     = params.nucleotides_per_second
if (params.sentieon)            summary['Sention']                           = "Using Sentieon for Preprocessing and/or Variant Calling"
if (params.skip_qc)             summary['QC tools skipped']                  = skipQC.join(', ')
if (params.target_bed)          summary['Target BED']                        = params.target_bed
if (params.tools)               summary['Tools']                             = tools.join(', ')

if (params.trim_fastq || params.split_fastq) summary['Modify fastqs (trim/split)'] = ""

if (params.trim_fastq) {
    summary['Fastq trim']         = "Fastq trim selected"
    summary['Trim R1']            = "${params.clip_r1} bp"
    summary['Trim R2']            = "${params.clip_r2} bp"
    summary["Trim 3' R1"]         = "${params.three_prime_clip_r1} bp"
    summary["Trim 3' R2"]         = "${params.three_prime_clip_r2} bp"
    summary['NextSeq Trim']       = "${params.trim_nextseq} bp"
    summary['Saved Trimmed Fastq'] = params.save_trimmed ? 'Yes' : 'No'
}
if (params.split_fastq)          summary['Reads in fastq']                   = params.split_fastq

summary['MarkDuplicates'] = "Options"
summary['Java options'] = params.markdup_java_options
summary['GATK Spark']   = params.use_gatk_spark ? 'Yes' : 'No'

summary['Save BAMs mapped']   = params.save_bam_mapped ? 'Yes' : 'No'
summary['Skip MarkDuplicates']   = params.skip_markduplicates ? 'Yes' : 'No'

if ('ascat' in tools) {
    summary['ASCAT'] = "Options"
    if (params.ascat_purity) summary['purity'] = params.ascat_purity
    if (params.ascat_ploidy) summary['ploidy'] = params.ascat_ploidy
}

if ('controlfreec' in tools) {
    summary['Control-FREEC'] = "Options"
    if (params.cf_coeff)                    summary['coefficientOfVariation']  = params.cf_coeff
    if (params.cf_contamination)            summary['contamination']           = params.cf_contamination
    if (params.cf_contamination_adjustment) summary['contaminationAdjustment'] = params.cf_contamination_adjustment
    if (params.cf_ploidy)                   summary['ploidy']                  = params.cf_ploidy
    if (params.cf_window)                   summary['window']                  = params.cf_window
}

if ('haplotypecaller' in tools)             summary['GVCF']       = params.generate_gvcf ? 'Yes' : 'No'
if ('strelka' in tools && 'manta' in tools) summary['Strelka BP'] = params.no_strelka_bp ? 'No' : 'Yes'
if (params.pon && ('mutect2' in tools || (params.sentieon && 'tnscope' in tools))) summary['Panel of normals'] = params.pon

if (params.annotate_tools) summary['Tools to annotate'] = annotate_tools.join(', ')

if (params.annotation_cache) {
    summary['Annotation cache'] = "Enabled"
    if (params.snpeff_cache) summary['snpEff cache'] = params.snpeff_cache
    if (params.vep_cache)    summary['VEP cache']    = params.vep_cache
}

if (params.cadd_cache) {
    summary['CADD cache'] = "Enabled"
    if (params.cadd_indels)  summary['CADD indels']  = params.cadd_indels
    if (params.cadd_wg_snvs) summary['CADD wg snvs'] = params.cadd_wg_snvs
}

if (params.genesplicer) summary['genesplicer'] = "Enabled"

if (params.igenomes_base && !params.igenomes_ignore) summary['AWS iGenomes base'] = params.igenomes_base
if (params.igenomes_ignore)                          summary['AWS iGenomes']      = "Do not use"
if (params.genomes_base && !params.igenomes_ignore)  summary['Genomes base']      = params.genomes_base

summary['Save Reference']    = params.save_reference ? 'Yes' : 'No'

if (params.ac_loci)                 summary['Loci']                    = params.ac_loci
if (params.ac_loci_gc)              summary['Loci GC']                 = params.ac_loci_gc
if (params.bwa)                     summary['BWA indexes']             = params.bwa
if (params.chr_dir)                 summary['Chromosomes']             = params.chr_dir
if (params.chr_length)              summary['Chromosomes length']      = params.chr_length
if (params.dbsnp)                   summary['dbsnp']                   = params.dbsnp
if (params.dbsnp_index)             summary['dbsnpIndex']              = params.dbsnp_index
if (params.dict)                    summary['dict']                    = params.dict
if (params.fasta)                   summary['fasta reference']         = params.fasta
if (params.fasta_fai)               summary['fasta index']             = params.fasta_fai
if (params.germline_resource)       summary['germline resource']       = params.germline_resource
if (params.germline_resource_index) summary['germline resource index'] = params.germline_resource_index
if (params.intervals)               summary['intervals']               = params.intervals
if (params.known_indels)            summary['known indels']            = params.known_indels
if (params.known_indels_index)      summary['known indels index']      = params.known_indels_index
if (params.mappability)             summary['Mappability']             = params.mappability
if (params.snpeff_cache)            summary['snpEff cache']            = params.snpeff_cache
if (params.snpeff_db)               summary['snpEff DB']               = params.snpeff_db
if (params.species)                 summary['species']                 = params.species
if (params.vep_cache)               summary['VEP cache']               = params.vep_cache
if (params.vep_cache_version)       summary['VEP cache version']       = params.vep_cache_version

summary['Output dir']        = params.outdir
summary['Publish dir mode']  = params.publish_dir_mode
if (params.sequencing_center) summary['Sequenced by'] = params.sequencing_center

summary['Launch dir']  = workflow.launchDir
summary['Working dir'] = workflow.workDir
summary['Script dir']  = workflow.projectDir
summary['User']        = workflow.userName

if (params.multiqc_config) summary['MultiQC config'] = params.multiqc_config

summary['Config Profile'] = workflow.profile

if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url

summary['Config Files']   = workflow.configFiles.join(', ')


if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region'] = params.awsregion
    summary['AWS Queue']  = params.awsqueue
    summary['AWS CLI']    = params.awscli
}

if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

if ('mutect2' in tools && !(params.pon)) log.warn "[nf-core/sarek] Mutect2 was requested, but as no panel of normals were given, results will not be optimal"
if (params.sentieon) log.warn "[nf-core/sarek] Sentieon will be used, only works if Sentieon is available where nf-core/sarek is run"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'sarek-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/sarek Workflow Summary'
    section_href: 'https://github.com/nf-core/sarek'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

// Parse software version numbers

process get_software_versions {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    when: !('versions' in skipQC)

    script:
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    """
    alleleCounter --version &> v_allelecount.txt 2>&1 || true
    bcftools --version &> v_bcftools.txt 2>&1 || true
    ${aligner} version &> v_bwa.txt 2>&1 || true
    cnvkit.py version &> v_cnvkit.txt 2>&1 || true
    configManta.py --version &> v_manta.txt 2>&1 || true
    configureStrelkaGermlineWorkflow.py --version &> v_strelka.txt 2>&1 || true
    echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
    echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
    snpEff -version &> v_snpeff.txt 2>&1 || true
    fastqc --version &> v_fastqc.txt 2>&1 || true
    freebayes --version &> v_freebayes.txt 2>&1 || true
    freec &> v_controlfreec.txt 2>&1 || true
    gatk ApplyBQSR --help &> v_gatk.txt 2>&1 || true
    msisensor &> v_msisensor.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    qualimap --version &> v_qualimap.txt 2>&1 || true
    R --version &> v_r.txt 2>&1 || true
    R -e "library(ASCAT); help(package='ASCAT')" &> v_ascat.txt 2>&1 || true
    samtools --version &> v_samtools.txt 2>&1 || true
    tiddit &> v_tiddit.txt 2>&1 || true
    trim_galore -v &> v_trim_galore.txt 2>&1 || true
    vcftools --version &> v_vcftools.txt 2>&1 || true
    vep --help &> v_vep.txt 2>&1 || true

    python ${workflow.launchDir/params.config_files}/bin/scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

ch_software_versions_yaml = ch_software_versions_yaml.dump(tag:'SOFTWARE VERSIONS')

/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built
process BuildBWAindexes {

    tag "${fasta}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/BWAIndex/${it}" : null }

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta}.*") into bwa_built

    when: !(params.bwa) && params.fasta && 'mapping' in step

    script:
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    """
    ${aligner} index ${fasta}
    """
}


ch_bwa = params.bwa ? Channel.value(file(params.bwa)) : bwa_built
process BuildDict {

    tag "${fasta}"
    memory '16G'
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta.baseName}.dict") into dictBuilt

    when: !(params.dict) && params.fasta && !('annotate' in step) && !('controlfreec' in step)

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}

ch_dict = params.dict ? Channel.value(file(params.dict)) : dictBuilt
process BuildFastaFai {

    tag "${fasta}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta}.fai") into fai_built

    when: !(params.fasta_fai) && params.fasta && !('annotate' in step)

    script:
    """
    samtools faidx ${fasta}
    """
}

ch_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : fai_built

process BuildDbsnpIndex {

    tag "${dbsnp}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(dbsnp) from ch_dbsnp

    output:
        file("${dbsnp}.tbi") into dbsnp_tbi

    when: !(params.dbsnp_index) && params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || 'tnscope' in tools)

    script:
    """
    tabix -p vcf ${dbsnp}
    """
}

ch_dbsnp_tbi = params.dbsnp ? params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : dbsnp_tbi : "null"

process BuildGermlineResourceIndex {

    tag "${germlineResource}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(germlineResource) from ch_germline_resource

    output:
        file("${germlineResource}.tbi") into germline_resource_tbi

    when: !(params.germline_resource_index) && params.germline_resource && 'mutect2' in tools

    script:
    """
    tabix -p vcf ${germlineResource}
    """
}

ch_germline_resource_tbi = params.germline_resource ? params.germline_resource_index ? Channel.value(file(params.germline_resource_index)) : germline_resource_tbi : "null"

process BuildKnownIndelsIndex {

    tag "${knownIndels}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        each file(knownIndels) from ch_known_indels

    output:
        file("${knownIndels}.tbi") into known_indels_tbi

    when: !(params.known_indels_index) && params.known_indels && ('mapping' in step || 'preparerecalibration' in step)

    script:
    """
    tabix -p vcf ${knownIndels}
    """
}

ch_known_indels_tbi = params.known_indels ? params.known_indels_index ? Channel.value(file(params.known_indels_index)) : known_indels_tbi.collect() : "null"

process BuildPonIndex {

    tag "${pon}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(pon) from ch_pon

    output:
        file("${pon}.tbi") into pon_tbi

    when: !(params.pon_index) && params.pon && ('tnscope' in tools || 'mutect2' in tools)

    script:
    """
    tabix -p vcf ${pon}
    """
}

ch_pon_tbi = params.pon ? params.pon_index ? Channel.value(file(params.pon_index)) : pon_tbi : "null"

process BuildIntervals {

    tag "${fastaFai}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir params.outdir, mode: params.publish_dir_mode,
    saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(fastaFai) from ch_fai

    output:
        file("${fastaFai.baseName}.bed") into intervalBuilt

    when: !(params.intervals) && !('annotate' in step) && !('controlfreec' in step) 

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fastaFai} > ${fastaFai.baseName}.bed
    """
}

ch_intervals = params.no_intervals ? "null" : params.intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : intervalBuilt

/*
================================================================================
                                  PREPROCESSING
================================================================================
*/

// STEP 0: CREATING INTERVALS FOR PARALLELIZATION (PREPROCESSING AND VARIANT CALLING)

process CreateIntervalBeds {

    tag "${intervals}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    input:
        file(intervals) from ch_intervals

    output:
        file '*.bed' into bedIntervals mode flatten

    when: (!params.no_intervals) && step != 'annotate'

    script:
    // If the interval file is BED format, the fifth column is interpreted to
    // contain runtime estimates, which is then used to combine short-running jobs
    if (hasExtension(intervals, "bed"))
        """
        awk -vFS="\t" '{
          t = \$5  # runtime estimate
          if (t == "") {
            # no runtime estimate in this row, assume default value
            t = (\$3 - \$2) / ${params.nucleotides_per_second}
          }
          if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
            # start a new chunk
            name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
            chunk = 0
            longest = 0
          }
          if (t > longest)
            longest = t
          chunk += t
          print \$0 > name
        }' ${intervals}
        """
    else if (hasExtension(intervals, "interval_list"))
        """
        grep -v '^@' ${intervals} | awk -vFS="\t" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }'
        """
    else
        """
        awk -vFS="[:-]" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }' ${intervals}
        """
}

bedIntervals = bedIntervals
    .map { intervalFile ->
        def duration = 0.0
        for (line in intervalFile.readLines()) {
            final fields = line.split('\t')
            if (fields.size() >= 5) duration += fields[4].toFloat()
            else {
                start = fields[1].toInteger()
                end = fields[2].toInteger()
                duration += (end - start) / params.nucleotides_per_second
            }
        }
        [duration, intervalFile]
        }.toSortedList({ a, b -> b[0] <=> a[0] })
    .flatten().collate(2)
    .map{duration, intervalFile -> intervalFile}

bedIntervals = bedIntervals.dump(tag:'bedintervals')

if (params.no_intervals && step != 'annotate') {
    file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
    bedIntervals = Channel.from(file("${params.outdir}/no_intervals.bed"))
}

(intBaseRecalibrator, intApplyBQSR, intHaplotypeCaller, intFreebayesSingle, intMpileup, bedIntervalsSingle, bedIntervalsPair) = bedIntervals.into(7)

// PREPARING CHANNELS FOR PREPROCESSING AND QC

inputBam = Channel.create()
inputPairReads = Channel.create()

if (step in ['preparerecalibration', 'recalibrate', 'variantcalling', 'controlfreec', 'annotate']) {
    inputBam.close()
    inputPairReads.close()
} else inputSample.choice(inputPairReads, inputBam) {hasExtension(it[3], "bam") ? 1 : 0}

(inputBam, inputBamFastQC) = inputBam.into(2)

// Removing inputFile2 which is null in case of uBAM
inputBamFastQC = inputBamFastQC.map {
    idPatient, idSample, idRun, inputFile1, inputFile2 ->
    [idPatient, idSample, idRun, inputFile1]
}

if (params.split_fastq){
    inputPairReads = inputPairReads
        // newly splitfastq are named based on split, so the name is easier to catch
        .splitFastq(by: params.split_fastq, compress:true, file:"split", pe:true)
        .map {idPatient, idSample, idRun, reads1, reads2 ->
            // The split fastq read1 is the 4th element (indexed 3) its name is split_3
            // The split fastq read2's name is split_4
            // It's followed by which split it's acutally based on the mother fastq file
            // Index start at 1
            // Extracting the index to get a new IdRun
            splitIndex = reads1.fileName.toString().minus("split_3.").minus(".gz")
            newIdRun = idRun + "_" + splitIndex
            // Giving the files a new nice name
            newReads1 = file("${idSample}_${newIdRun}_R1.fastq.gz")
            newReads2 = file("${idSample}_${newIdRun}_R2.fastq.gz")
            [idPatient, idSample, newIdRun, reads1, reads2]}
}

inputPairReads = inputPairReads.dump(tag:'INPUT')

(inputPairReads, inputPairReadsTrimGalore, inputPairReadsFastQC, inputPairReadsUMI) = inputPairReads.into(4)

if (params.umi) inputPairReads.close()
else inputPairReadsUMI.close()

if (params.trim_fastq) inputPairReads.close()
else inputPairReadsTrimGalore.close()

// STEP 0.5: QC ON READS

// TODO: Use only one process for FastQC for FASTQ files and uBAM files
// FASTQ and uBAM files are renamed based on the sample name

process FastQCFQ {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'

    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz") from inputPairReadsFastQC

    output:
        file("*.{html,zip}") into fastQCFQReport

    when: !('fastqc' in skipQC)

    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz
    """
}

process FastQCBAM {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") from inputBamFastQC

    output:
        file("*.{html,zip}") into fastQCBAMReport

    when: !('fastqc' in skipQC)

    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}.bam
    """
}

fastQCReport = fastQCFQReport.mix(fastQCBAMReport)

fastQCReport = fastQCReport.dump(tag:'FastQC')

process TrimGalore {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-small'
    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/TrimGalore/${idSample}_${idRun}", mode: params.publish_dir_mode,
      saveAs: {filename ->
        if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
        else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
        else if (params.save_trimmed) filename
        else null
      }

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz") from inputPairReadsTrimGalore

    output:
        file("*.{html,zip,txt}") into trimGaloreReport
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1_val_1.fq.gz"), file("${idSample}_${idRun}_R2_val_2.fq.gz") into outputPairReadsTrimGalore

    when: params.trim_fastq

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
      cores = (task.cpus as int) - 4
      if (cores < 1) cores = 1
      if (cores > 4) cores = 4
      }
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
    """
    trim_galore \
         --cores ${cores} \
        --paired \
        --fastqc \
        --gzip \
        ${c_r1} ${c_r2} \
        ${tpc_r1} ${tpc_r2} \
        ${nextseq} \
        ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz

    mv *val_1_fastqc.html "${idSample}_${idRun}_R1.trimmed_fastqc.html"
    mv *val_2_fastqc.html "${idSample}_${idRun}_R2.trimmed_fastqc.html"
    mv *val_1_fastqc.zip "${idSample}_${idRun}_R1.trimmed_fastqc.zip"
    mv *val_2_fastqc.zip "${idSample}_${idRun}_R2.trimmed_fastqc.zip"
    """
}

/*
================================================================================
                            UMIs PROCESSING
================================================================================
*/

// UMI - STEP 1 - ANNOTATE
// the process needs to convert fastq to unmapped bam
// and while doing the conversion, tag the bam field RX with the UMI sequence

process UMIFastqToBAM {

    publishDir "${params.outdir}/Reports/${idSample}/UMI/${idSample}_${idRun}", mode: params.publish_dir_mode
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz") from inputPairReadsUMI
        val read_structure1 from ch_read_structure1
        val read_structure2 from ch_read_structure2

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi_converted.bam") into umi_converted_bams_ch

    when: params.umi

    // tmp folder for fgbio might be solved more elengantly?

    script:
    """
    mkdir tmp

    fgbio --tmp-dir=${PWD}/tmp \
    FastqToBam \
    -i "${idSample}_${idRun}_R1.fastq.gz" "${idSample}_${idRun}_R2.fastq.gz" \
    -o "${idSample}_umi_converted.bam" \
    --read-structures ${read_structure1} ${read_structure2} \
    --sample ${idSample} \
    --library ${idSample}
    """
}

// UMI - STEP 2 - MAP THE BAM FILE
// this is necessary because the UMI groups are created based on
// mapping position + same UMI tag

process UMIMapBamFile {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    input:
        set idPatient, idSample, idRun, file(convertedBam) from umi_converted_bams_ch
        file(bwaIndex) from ch_bwa
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi_unsorted.bam") into umi_aligned_bams_ch

    when: params.umi

    script:
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    """
    samtools bam2fq -T RX ${convertedBam} | \
    ${aligner} mem -p -t ${task.cpus} -C -M -R \"@RG\\tID:${idSample}\\tSM:${idSample}\\tPL:Illumina\" \
    ${fasta} - | \
    samtools view -bS - > ${idSample}_umi_unsorted.bam
    """
}

// UMI - STEP 3 - GROUP READS BY UMIs
// We have chose the Adjacency method, following the nice paper and blog explanation integrated in both
// UMItools and FGBIO
// https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/
// alternatively we can define this as input for the user to choose from

process GroupReadsByUmi {

    publishDir "${params.outdir}/Reports/${idSample}/UMI/${idSample}_${idRun}", mode: params.publish_dir_mode
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    input:
        set idPatient, idSample, idRun, file(alignedBam) from umi_aligned_bams_ch

    output:
        file("${idSample}_umi_histogram.txt") into umi_histogram_ch
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi-grouped.bam") into umi_grouped_bams_ch

    when: params.umi

    script:
    """
    mkdir tmp

    samtools view -h ${alignedBam} | \
    samblaster -M --addMateTags | \
    samtools view -Sb - >${idSample}_unsorted_tagged.bam

    fgbio --tmp-dir=${PWD}/tmp \
    GroupReadsByUmi \
    -s Adjacency \
    -i ${idSample}_unsorted_tagged.bam \
    -o ${idSample}_umi-grouped.bam \
    -f ${idSample}_umi_histogram.txt
    """
}

// UMI - STEP 4 - CALL MOLECULAR CONSENSUS
// Now that the reads are organised by UMI groups a molecular consensus will be created
// the resulting bam file will be again unmapped and therefore can be fed into the
// existing workflow from the step mapping
process CallMolecularConsensusReads {
    publishDir "${params.outdir}/Reports/${idSample}/UMI/${idSample}_${idRun}", mode: params.publish_dir_mode
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    input:
        set idPatient, idSample, idRun, file(groupedBamFile) from umi_grouped_bams_ch

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi-consensus.bam"), val("null") into consensus_bam_ch

    when: params.umi

    script:
    """
    mkdir tmp

    fgbio --tmp-dir=${PWD}/tmp \
    CallMolecularConsensusReads \
    -i $groupedBamFile \
    -o ${idSample}_umi-consensus.bam \
    -M 1 -S Coordinate
    """
}

// ################# END OF UMI READS PRE-PROCESSING
// from this moment on the generated uBam files can feed into the existing tools

// STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA MEM

input_pair_reads_sentieon = Channel.empty()

if (params.umi) {
    inputPairReads = inputPairReads.dump(tag:'INPUT before mapping')
    if (params.sentieon) input_pair_reads_sentieon = consensus_bam_ch
    else inputPairReads = consensus_bam_ch
}
else {
    if (params.trim_fastq) inputPairReads = outputPairReadsTrimGalore
    else inputPairReads = inputPairReads.mix(inputBam)
    inputPairReads = inputPairReads.dump(tag:'INPUT before mapping')

    (inputPairReads, input_pair_reads_sentieon) = inputPairReads.into(2)
    if (params.sentieon) inputPairReads.close()
    else input_pair_reads_sentieon.close()
}

process MapReads {

    container pipeline_container
    memory '70G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'
    tag "${idPatient}-${idRun}"

    input:
        set idPatient, idSample, idRun, file(inputFile1), file(inputFile2) from inputPairReads
        file(bwaIndex) from ch_bwa
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") into bamMapped
        set idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam") into bamMappedBamQC

    when: !(params.sentieon)

    script:
    // -K is an hidden option, used to fix the number of reads processed by bwa mem
    // Chunk size can affect bwa results, if not specified,
    // the number of threads can change which can give not deterministic result.
    // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
    // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    // adjust mismatch penalty for tumor samples
    status = statusMap[idPatient, idSample]
    extra = status == 1 ? "-B 3" : ""
    convertToFastq = hasExtension(inputFile1, "bam") ? "gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile1} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true | \\" : ""
    input = hasExtension(inputFile1, "bam") ? "-p /dev/stdin - 2> >(tee ${inputFile1}.bwa.stderr.log >&2)" : "${inputFile1} ${inputFile2}"
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    """
    ${convertToFastq}
    ${aligner} mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
    ${input} | \
    samtools sort --threads ${task.cpus} -m 2G - > ${idSample}_${idRun}.bam
    """
}

bamMapped = bamMapped.dump(tag:'Mapped BAM')
// Sort BAM whether they are standalone or should be merged

singleBam = Channel.create()
multipleBam = Channel.create()
bamMapped.groupTuple(by:[0, 1])
    .choice(singleBam, multipleBam) {it[2].size() > 1 ? 1 : 0}
singleBam = singleBam.map {
    idPatient, idSample, idRun, bam ->
    [idPatient, idSample, bam]
}
singleBam = singleBam.dump(tag:'Single BAM')

// STEP 1': MAPPING READS TO REFERENCE GENOME WITH SENTIEON BWA MEM

process Sentieon_MapReads {
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'

    tag "${idPatient}-${idRun}"

    input:
        set idPatient, idSample, idRun, file(inputFile1), file(inputFile2) from input_pair_reads_sentieon
        file(bwaIndex) from ch_bwa
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") into bam_sentieon_mapped

    when: params.sentieon

    script:
    // -K is an hidden option, used to fix the number of reads processed by bwa mem
    // Chunk size can affect bwa results, if not specified,
    // the number of threads can change which can give not deterministic result.
    // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
    // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    // adjust mismatch penalty for tumor samples
    status = statusMap[idPatient, idSample]
    extra = status == 1 ? "-B 3" : ""
    """
    sentieon bwa mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
    ${inputFile1} ${inputFile2} | \
    sentieon util sort -r ${fasta} -o ${idSample}_${idRun}.bam -t ${task.cpus} --sam2bam -i -
    """
}

bam_sentieon_mapped = bam_sentieon_mapped.dump(tag:'Sentieon Mapped BAM')
// Sort BAM whether they are standalone or should be merged

singleBamSentieon = Channel.create()
multipleBamSentieon = Channel.create()
bam_sentieon_mapped.groupTuple(by:[0, 1])
    .choice(singleBamSentieon, multipleBamSentieon) {it[2].size() > 1 ? 1 : 0}
singleBamSentieon = singleBamSentieon.map {
    idPatient, idSample, idRun, bam ->
    [idPatient, idSample, bam]
}
singleBamSentieon = singleBamSentieon.dump(tag:'Single BAM')

// STEP 1.5: MERGING BAM FROM MULTIPLE LANES

multipleBam = multipleBam.mix(multipleBamSentieon)

process MergeBamMapped {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'

    tag "${idPatient}-${idSample}"

    input:
        set idPatient, idSample, idRun, file(bam) from multipleBam

    output:
        set idPatient, idSample, file("${idSample}.bam") into bam_mapped_merged

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
    """
}

bam_mapped_merged = bam_mapped_merged.dump(tag:'Merged BAM')

bam_mapped_merged = bam_mapped_merged.mix(singleBam,singleBamSentieon)

(bam_mapped_merged, bam_sentieon_mapped_merged) = bam_mapped_merged.into(2)

if (!params.sentieon) bam_sentieon_mapped_merged.close()
else bam_mapped_merged.close()

bam_mapped_merged = bam_mapped_merged.dump(tag:'BAMs for MD')
bam_sentieon_mapped_merged = bam_sentieon_mapped_merged.dump(tag:'Sentieon BAMs to Index')

process IndexBamMergedForSentieon {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'

    tag "${idPatient}-${idSample}"

    input:
        set idPatient, idSample, file("${idSample}.bam") from bam_sentieon_mapped_merged

    output:
        set idPatient, idSample, file("${idSample}.bam"), file("${idSample}.bam.bai") into bam_sentieon_mapped_merged_indexed

    script:
    """
    samtools index ${idSample}.bam
    """
}

(bam_mapped_merged, bam_mapped_merged_to_index) = bam_mapped_merged.into(2)

process IndexBamFile {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (save_bam_mapped) "Preprocessing/${idSample}/Mapped/${it}"
            else null
        }

    input:
        set idPatient, idSample, file("${idSample}.bam") from bam_mapped_merged_to_index

    output:
        set idPatient, idSample, file("${idSample}.bam"), file("${idSample}.bam.bai") into bam_mapped_merged_indexed
        set idPatient, idSample into tsv_bam_indexed

    when: save_bam_mapped || !(params.known_indels)

    script:
    """
    samtools index ${idSample}.bam
    """
}

if (!save_bam_mapped) tsv_bam_indexed.close()

(tsv_bam_indexed, tsv_bam_indexed_sample) = tsv_bam_indexed.into(2)

// Creating a TSV file to restart from this step
tsv_bam_indexed.map { idPatient, idSample ->
    gender = genderMap[idPatient]
    status = statusMap[idPatient, idSample]
    bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
    bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
}.collectFile(
    name: 'mapped.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
)

tsv_bam_indexed_sample
    .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
        ["mapped_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
}
// STEP 2: MARKING DUPLICATES

process MarkDuplicates {

    container pipeline_container
    memory '70G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'
    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}.bam.metrics") "Reports/${idSample}/MarkDuplicates/${it}"
            else "Preprocessing/${idSample}/DuplicatesMarked/${it}"
        }

    input:
        set idPatient, idSample, file("${idSample}.bam") from bam_mapped_merged

    output:
        set idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bam.bai") into bam_duplicates_marked
        set idPatient, idSample into tsv_bam_duplicates_marked
        file ("${idSample}.bam.metrics") optional true into duplicates_marked_report

    when: !(params.skip_markduplicates)

    script:
    markdup_java_options = task.memory.toGiga() > 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    metrics = 'markduplicates' in skipQC ? '' : "-M ${idSample}.bam.metrics"
    if (params.use_gatk_spark)
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicatesSpark \
        -I ${idSample}.bam \
        -O ${idSample}.md.bam \
        ${metrics} \
        --tmp-dir . \
        --create-output-bam-index true \
        --spark-master local[${task.cpus}]
    """
    else
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicates \
        --INPUT ${idSample}.bam \
        --METRICS_FILE ${idSample}.bam.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${idSample}.md.bam
    
    mv ${idSample}.md.bai ${idSample}.md.bam.bai
    """
}

(tsv_bam_duplicates_marked, tsv_bam_duplicates_marked_sample) = tsv_bam_duplicates_marked.into(2)

// Creating a TSV file to restart from this step
tsv_bam_duplicates_marked.map { idPatient, idSample ->
    gender = genderMap[idPatient]
    status = statusMap[idPatient, idSample]
    bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
    bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
}.collectFile(
    name: 'duplicates_marked_no_table.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
)

tsv_bam_duplicates_marked_sample
    .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
        ["duplicates_marked_no_table_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
}

if ('markduplicates' in skipQC) duplicates_marked_report.close()

if (step == 'preparerecalibration') bam_duplicates_marked = inputSample

bam_duplicates_marked = bam_duplicates_marked.dump(tag:'MD BAM')
duplicates_marked_report = duplicates_marked_report.dump(tag:'MD Report')

if (params.skip_markduplicates) bam_duplicates_marked = bam_mapped_merged_indexed

(bamMD, bamMDToJoin, bam_duplicates_marked) = bam_duplicates_marked.into(3)

bamBaseRecalibrator = bamMD.combine(intBaseRecalibrator)

bamBaseRecalibrator = bamBaseRecalibrator.dump(tag:'BAM FOR BASERECALIBRATOR')

// STEP 2': SENTIEON DEDUP

process Sentieon_Dedup {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}_*.txt" && 'sentieon' in skipQC) null
            else if (it == "${idSample}_*.txt") "Reports/${idSample}/Sentieon/${it}"
            else "Preprocessing/${idSample}/DedupedSentieon/${it}"
        }

    input:
        set idPatient, idSample, file(bam), file(bai) from bam_sentieon_mapped_merged_indexed
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, file("${idSample}.deduped.bam"), file("${idSample}.deduped.bam.bai") into bam_sentieon_dedup

    when: params.sentieon

    script:
    """
    sentieon driver \
        -t ${task.cpus} \
        -i ${bam} \
        -r ${fasta} \
        --algo GCBias --summary ${idSample}_gc_summary.txt ${idSample}_gc_metric.txt \
        --algo MeanQualityByCycle ${idSample}_mq_metric.txt \
        --algo QualDistribution ${idSample}_qd_metric.txt \
        --algo InsertSizeMetricAlgo ${idSample}_is_metric.txt  \
        --algo AlignmentStat ${idSample}_aln_metric.txt

    sentieon driver \
        -t ${task.cpus} \
        -i ${bam} \
        --algo LocusCollector \
        --fun score_info ${idSample}_score.gz

    sentieon driver \
        -t ${task.cpus} \
        -i ${bam} \
        --algo Dedup \
        --rmdup \
        --score_info ${idSample}_score.gz  \
        --metrics ${idSample}_dedup_metric.txt ${idSample}.deduped.bam
    """
}

// STEP 3: CREATING RECALIBRATION TABLES

process BaseRecalibrator {

    container pipeline_container
    memory '70G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'
    tag "${idPatient}-${idSample}-${intervalBed.baseName}"

    input:
        set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamBaseRecalibrator
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(fasta) from ch_fasta
        file(dict) from ch_dict
        file(fastaFai) from ch_fai
        file(knownIndels) from ch_known_indels
        file(knownIndelsIndex) from ch_known_indels_tbi

    output:
        set idPatient, idSample, file("${prefix}${idSample}.recal.table") into tableGatherBQSRReports
        set idPatient, idSample into recalTableTSVnoInt

    when: params.known_indels

    script:
    dbsnpOptions = params.dbsnp ? "--known-sites ${dbsnp}" : ""
    knownOptions = params.known_indels ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    // TODO: --use-original-qualities ???
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        BaseRecalibrator \
        -I ${bam} \
        -O ${prefix}${idSample}.recal.table \
        --tmp-dir . \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        ${knownOptions} \
        --verbosity INFO
    """
}

if (!params.no_intervals) tableGatherBQSRReports = tableGatherBQSRReports.groupTuple(by:[0, 1])

tableGatherBQSRReports = tableGatherBQSRReports.dump(tag:'BQSR REPORTS')

if (params.no_intervals) {
    (tableGatherBQSRReports, tableGatherBQSRReportsNoInt) = tableGatherBQSRReports.into(2)
    recalTable = tableGatherBQSRReportsNoInt
} else recalTableTSVnoInt.close()

// STEP 3.5: MERGING RECALIBRATION TABLES

process GatherBQSRReports {

    container pipeline_container
    memory '70G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}.recal.table" && !params.skip_markduplicates) "Preprocessing/${idSample}/DuplicatesMarked/${it}"
            else "Preprocessing/${idSample}/Mapped/${it}"
        }

    input:
        set idPatient, idSample, file(recal) from tableGatherBQSRReports

    output:
        set idPatient, idSample, file("${idSample}.recal.table") into recalTable
        file("${idSample}.recal.table") into baseRecalibratorReport
        set idPatient, idSample into recalTableTSV

    when: !(params.no_intervals)

    script:
    input = recal.collect{"-I ${it}"}.join(' ')
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GatherBQSRReports \
        ${input} \
        -O ${idSample}.recal.table \
    """
}

if ('baserecalibrator' in skipQC) baseRecalibratorReport.close()

recalTable = recalTable.dump(tag:'RECAL TABLE')

(recalTableTSV, recalTableSampleTSV) = recalTableTSV.mix(recalTableTSVnoInt).into(2)

// Create TSV files to restart from this step
if (params.skip_markduplicates) {
    recalTableTSV.map { idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
        recalTable = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.recal.table"
        "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"
    }.collectFile(
        name: 'mapped_no_duplicates_marked.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
    )

    recalTableSampleTSV
        .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV/") {
            idPatient, idSample ->
            status = statusMap[idPatient, idSample]
            gender = genderMap[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
            recalTable = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.recal.table"
            ["mapped_no_duplicates_marked_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"]
    }
} else {
    recalTableTSV.map { idPatient, idSample ->
    status = statusMap[idPatient, idSample]
    gender = genderMap[idPatient]
    bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
    bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
    recalTable = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.recal.table"

        "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"
    }.collectFile(
        name: 'duplicates_marked.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
    )

    recalTableSampleTSV
        .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV/") {
            idPatient, idSample ->
            status = statusMap[idPatient, idSample]
            gender = genderMap[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
            recalTable = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.recal.table"
            ["duplicates_marked_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"]
    }
}

bamApplyBQSR = bamMDToJoin.join(recalTable, by:[0,1])

if (step == 'recalibrate') bamApplyBQSR = inputSample

bamApplyBQSR = bamApplyBQSR.dump(tag:'BAM + BAI + RECAL TABLE')

bamApplyBQSR = bamApplyBQSR.combine(intApplyBQSR)

bamApplyBQSR = bamApplyBQSR.dump(tag:'BAM + BAI + RECAL TABLE + INT')

// STEP 4: RECALIBRATING

process ApplyBQSR {

    container pipeline_container
    memory '70G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'

    tag "${idPatient}-${idSample}-${intervalBed.baseName}"

    input:
        set idPatient, idSample, file(bam), file(bai), file(recalibrationReport), file(intervalBed) from bamApplyBQSR
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, file("${prefix}${idSample}.recal.bam") into bam_recalibrated_to_merge

    script:
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        ApplyBQSR \
        -R ${fasta} \
        --input ${bam} \
        --output ${prefix}${idSample}.recal.bam \
        ${intervalsOptions} \
        --bqsr-recal-file ${recalibrationReport}
    """
}

(bam_recalibrated_to_merge, bam_recalibrated_to_index) = bam_recalibrated_to_merge.groupTuple(by:[0, 1]).into(2)

// STEP 4': SENTIEON BQSR

bam_sentieon_dedup = bam_sentieon_dedup.dump(tag:'deduped.bam')

process Sentieon_BQSR {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}_recal_result.csv" && 'sentieon' in skipQC) "Reports/${idSample}/Sentieon/${it}"
            else "Preprocessing/${idSample}/RecalSentieon/${it}"
        }

    input:
        set idPatient, idSample, file(bam), file(bai) from bam_sentieon_dedup
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(fasta) from ch_fasta
        file(dict) from ch_dict
        file(fastaFai) from ch_fai
        file(knownIndels) from ch_known_indels
        file(knownIndelsIndex) from ch_known_indels_tbi

    output:
        set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bam_sentieon_recal
        set idPatient, idSample, file(bam), file(bai), file("${idSample}.recal.table") into bam_sentieon_deduped_table
        set idPatient, idSample into tsv_sentieon

    when: params.sentieon

    script:
    known = knownIndels.collect{"--known-sites ${it}"}.join(' ')
    """
    sentieon driver  \
        -t ${task.cpus} \
        -r ${fasta} \
        -i ${idSample}.deduped.bam \
        --algo QualCal \
        -k ${dbsnp} \
        ${idSample}.recal.table

    sentieon driver \
        -t ${task.cpus} \
        -r ${fasta} \
        -i ${idSample}.deduped.bam \
        -q ${idSample}.recal.table \
        --algo QualCal \
        -k ${dbsnp} \
        ${idSample}.table.post \
        --algo ReadWriter ${idSample}.recal.bam

    sentieon driver \
        -t ${task.cpus} \
        --algo QualCal \
        --plot \
        --before ${idSample}.recal.table \
        --after ${idSample}.table.post \
        ${idSample}_recal_result.csv
    """
}

(tsv_sentieon_deduped, tsv_sentieon_deduped_sample, tsv_sentieon_recal, tsv_sentieon_recal_sample) = tsv_sentieon.into(4)

// Creating a TSV file to restart from this step
tsv_sentieon_deduped.map { idPatient, idSample ->
    gender = genderMap[idPatient]
    status = statusMap[idPatient, idSample]
    bam = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam"
    bai = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam.bai"
    table = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.table"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${table}\n"
}.collectFile(
    name: 'sentieon_deduped.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
)

tsv_sentieon_deduped_sample
    .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam.bai"
        table = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.table"
        ["sentieon_deduped_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${table}\n"]
}

// Creating a TSV file to restart from this step
tsv_sentieon_recal.map { idPatient, idSample ->
    gender = genderMap[idPatient]
    status = statusMap[idPatient, idSample]
    bam = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam"
    bai = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam.bai"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
}.collectFile(
    name: 'sentieon_recalibrated.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
)

tsv_sentieon_recal_sample
    .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam.bai"
        ["sentieon_recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
}

// STEP 4.5: MERGING THE RECALIBRATED BAM FILES

process MergeBamRecal {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idPatient}-${idSample}"

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(bam) from bam_recalibrated_to_merge

    output:
        set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bam_recalibrated
        set idPatient, idSample, file("${idSample}.recal.bam") into bam_recalibrated_qc
        set idPatient, idSample into tsv_bam_recalibrated

    when: !(params.no_intervals)

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.recal.bam ${bam}
    samtools index ${idSample}.recal.bam
    """
}

// STEP 4.5': INDEXING THE RECALIBRATED BAM FILES

process IndexBamRecal {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idPatient}-${idSample}"

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file("${idSample}.recal.bam") from bam_recalibrated_to_index

    output:
        set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bam_recalibrated_indexed
        set idPatient, idSample, file("${idSample}.recal.bam") into bam_recalibrated_no_int_qc
        set idPatient, idSample into tsv_bam_recalibrated_no_int

    when: params.no_intervals

    script:
    """
    samtools index ${idSample}.recal.bam
    """
}

bam_recalibrated = bam_recalibrated.mix(bam_recalibrated_indexed)
bam_recalibrated_qc = bam_recalibrated_qc.mix(bam_recalibrated_no_int_qc)
tsv_bam_recalibrated = tsv_bam_recalibrated.mix(tsv_bam_recalibrated_no_int)

(bam_recalibrated_bamqc, bam_recalibrated_samtools_stats) = bam_recalibrated_qc.into(2)
(tsv_bam_recalibrated, tsv_bam_recalibrated_sample) = tsv_bam_recalibrated.into(2)

// Creating a TSV file to restart from this step
tsv_bam_recalibrated.map { idPatient, idSample ->
    gender = genderMap[idPatient]
    status = statusMap[idPatient, idSample]
    bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
    bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
}.collectFile(
    name: 'recalibrated.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
)

tsv_bam_recalibrated_sample
    .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") {
        idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
        ["recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
}

// STEP 5: QC

process SamtoolsStats {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idPatient}-${idSample}"

    publishDir "${params.outdir}/Reports/${idSample}/SamToolsStats", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(bam) from bam_recalibrated_samtools_stats

    output:
        file ("${bam}.samtools.stats.out") into samtoolsStatsReport

    when: !('samtools' in skipQC)

    script:
    """
    samtools stats ${bam} > ${bam}.samtools.stats.out
    """
}

samtoolsStatsReport = samtoolsStatsReport.dump(tag:'SAMTools')

bamBamQC = bamMappedBamQC.mix(bam_recalibrated_bamqc)

process BamQC {

    container pipeline_container
    memory '70G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'

    tag "${idPatient}-${idSample}"

    publishDir "${params.outdir}/Reports/${idSample}/bamQC", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(bam) from bamBamQC
        file(targetBED) from ch_target_bed

    output:
        file("${bam.baseName}") into bamQCReport

    when: !('bamqc' in skipQC)

    script:
    use_bed = params.target_bed ? "-gff ${targetBED}" : ''
    """
    qualimap --java-mem-size=${task.memory.toGiga()}G \
        bamqc \
        -bam ${bam} \
        --paint-chromosome-limits \
        --genome-gc-distr HUMAN \
        $use_bed \
        -nt ${task.cpus} \
        -skip-duplicated \
        --skip-dup-mode 0 \
        -outdir ${bam.baseName} \
        -outformat HTML
    """
}

bamQCReport = bamQCReport.dump(tag:'BamQC')

/*
================================================================================
                            GERMLINE VARIANT CALLING
================================================================================
*/

// When using sentieon for mapping, Channel bam_recalibrated is bam_sentieon_recal
if (params.sentieon && step == 'mapping') bam_recalibrated = bam_sentieon_recal

// When no knownIndels for mapping, Channel bam_recalibrated is bam_duplicates_marked
if (!params.known_indels && step in ['mapping', 'preparerecalibration']) bam_recalibrated = bam_duplicates_marked

// When starting with variant calling, Channel bam_recalibrated is inputSample
if (step == 'variantcalling') bam_recalibrated = inputSample

bam_recalibrated = bam_recalibrated.dump(tag:'BAM for Variant Calling')

// Here we have a recalibrated bam set
// The TSV file is formatted like: "idPatient status idSample bamFile baiFile"
// Manta will be run in Germline mode, or in Tumor mode depending on status
// HaplotypeCaller, TIDDIT and Strelka will be run for Normal and Tumor samples

(bamMantaSingle, bamStrelkaSingle, bamTIDDIT, bamFreebayesSingleNoIntervals, bamHaplotypeCallerNoIntervals, bamRecalAll) = bam_recalibrated.into(6)

(bam_sentieon_DNAseq, bam_sentieon_DNAscope, bam_sentieon_all) = bam_sentieon_deduped_table.into(3)

// To speed Variant Callers up we are chopping the reference into smaller pieces
// Do variant calling by this intervals, and re-merge the VCFs

bamHaplotypeCaller = bamHaplotypeCallerNoIntervals.combine(intHaplotypeCaller)
bamFreebayesSingle = bamFreebayesSingleNoIntervals.combine(intFreebayesSingle)

// STEP GATK HAPLOTYPECALLER.1

process HaplotypeCaller {

    container pipeline_container
    memory '70G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-mediium'

    tag "${idSample}-${intervalBed.baseName}"

    input:
        set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamHaplotypeCaller
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf") into gvcfHaplotypeCaller
        set idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf") into gvcfGenotypeGVCFs

    when: 'haplotypecaller' in tools

    script:
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    dbsnpOptions = params.dbsnp ? "--D ${dbsnp}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        -O ${intervalBed.baseName}_${idSample}.g.vcf \
        -ERC GVCF
    """
}

gvcfHaplotypeCaller = gvcfHaplotypeCaller.groupTuple(by:[0, 1, 2])

if (!params.generate_gvcf) gvcfHaplotypeCaller.close()
else gvcfHaplotypeCaller = gvcfHaplotypeCaller.dump(tag:'GVCF HaplotypeCaller')

// STEP GATK HAPLOTYPECALLER.2

process GenotypeGVCFs {

    tag "${idSample}-${intervalBed.baseName}"
    memory '70G'
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-mediium'
    input:
        set idPatient, idSample, file(intervalBed), file(gvcf) from gvcfGenotypeGVCFs
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
    set val("HaplotypeCaller"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into vcfGenotypeGVCFs

    when: 'haplotypecaller' in tools

    script:
    // Using -L is important for speed and we have to index the interval files also
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    dbsnpOptions = params.dbsnp ? "--D ${dbsnp}" : ""
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        IndexFeatureFile \
        -I ${gvcf}

    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        -V ${gvcf} \
        -O ${intervalBed.baseName}_${idSample}.vcf
    """
}

vcfGenotypeGVCFs = vcfGenotypeGVCFs.groupTuple(by:[0, 1, 2])

// STEP SENTIEON DNAseq

process Sentieon_DNAseq {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idSample}"

    input:
        set idPatient, idSample, file(bam), file(bai), file(recal) from bam_sentieon_DNAseq
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
    set val("SentieonDNAseq"), idPatient, idSample, file("DNAseq_${idSample}.vcf") into vcf_sentieon_DNAseq

    when: 'dnaseq' in tools && params.sentieon

    script:
    """
    sentieon driver \
        -t ${task.cpus} \
        -r ${fasta} \
        -i ${bam} \
        -q ${recal} \
        --algo Haplotyper \
        -d ${dbsnp} \
        DNAseq_${idSample}.vcf
    """
}

vcf_sentieon_DNAseq = vcf_sentieon_DNAseq.dump(tag:'sentieon DNAseq')

// STEP SENTIEON DNAscope

process Sentieon_DNAscope {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idSample}"

    input:
        set idPatient, idSample, file(bam), file(bai), file(recal) from bam_sentieon_DNAscope
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
    set val("SentieonDNAscope"), idPatient, idSample, file("DNAscope_${idSample}.vcf") into vcf_sentieon_DNAscope
    set val("SentieonDNAscope"), idPatient, idSample, file("DNAscope_SV_${idSample}.vcf") into vcf_sentieon_DNAscope_SV

    when: 'dnascope' in tools && params.sentieon

    script:
    """
    sentieon driver \
        -t ${task.cpus} \
        -r ${fasta} \
        -i ${bam} \
        -q ${recal} \
        --algo DNAscope \
        -d ${dbsnp} \
        DNAscope_${idSample}.vcf

    sentieon driver \
        -t ${task.cpus} \
        -r ${fasta}\
        -i ${bam} \
        -q ${recal} \
        --algo DNAscope \
        --var_type bnd \
        -d ${dbsnp} \
        DNAscope_${idSample}.temp.vcf

    sentieon driver \
        -t ${task.cpus} \
        -r ${fasta}\
        -q ${recal} \
        --algo SVSolver \
        -v DNAscope_${idSample}.temp.vcf \
        DNAscope_SV_${idSample}.vcf
    """
}

vcf_sentieon_DNAscope = vcf_sentieon_DNAscope.dump(tag:'sentieon DNAscope')
vcf_sentieon_DNAscope_SV = vcf_sentieon_DNAscope_SV.dump(tag:'sentieon DNAscope SV')

// STEP STRELKA.1 - SINGLE MODE

process StrelkaSingle {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Strelka", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(bam), file(bai) from bamStrelkaSingle
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
        set val("Strelka"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelkaSingle

    when: 'strelka' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configureStrelkaGermlineWorkflow.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        ${options} \
        --runDir Strelka

    python Strelka/runWorkflow.py -m local -j ${task.cpus}

    mv Strelka/results/variants/genome.*.vcf.gz \
        Strelka_${idSample}_genome.vcf.gz
    mv Strelka/results/variants/genome.*.vcf.gz.tbi \
        Strelka_${idSample}_genome.vcf.gz.tbi
    mv Strelka/results/variants/variants.vcf.gz \
        Strelka_${idSample}_variants.vcf.gz
    mv Strelka/results/variants/variants.vcf.gz.tbi \
        Strelka_${idSample}_variants.vcf.gz.tbi
    """
}

vcfStrelkaSingle = vcfStrelkaSingle.dump(tag:'Strelka - Single Mode')

// STEP MANTA.1 - SINGLE MODE

process MantaSingle {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Manta", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(bam), file(bai) from bamMantaSingle
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
        set val("Manta"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfMantaSingle

    when: 'manta' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    status = statusMap[idPatient, idSample]
    inputbam = status == 0 ? "--bam" : "--tumorBam"
    vcftype = status == 0 ? "diploid" : "tumor"
    """
    ${beforeScript}
    configManta.py \
        ${inputbam} ${bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSample}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSample}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSample}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${idSample}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
    """
}

vcfMantaSingle = vcfMantaSingle.dump(tag:'Single Manta')

// STEP TIDDIT

process TIDDIT {

    tag "${idSample}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "TIDDIT_${idSample}.vcf") "VariantCalling/${idSample}/TIDDIT/${it}"
            else "Reports/${idSample}/TIDDIT/${it}"
        }

    input:
        set idPatient, idSample, file(bam), file(bai) from bamTIDDIT
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set val("TIDDIT"), idPatient, idSample, file("*.vcf.gz"), file("*.tbi") into vcfTIDDIT
        set file("TIDDIT_${idSample}.old.vcf"), file("TIDDIT_${idSample}.ploidy.tab"), file("TIDDIT_${idSample}.signals.tab"), file("TIDDIT_${idSample}.wig"), file("TIDDIT_${idSample}.gc.wig") into tidditOut

    when: 'tiddit' in tools

    script:
    """
    tiddit --sv -o TIDDIT_${idSample} --bam ${bam} --ref ${fasta}

    mv TIDDIT_${idSample}.vcf TIDDIT_${idSample}.old.vcf

    grep -E "#|PASS" TIDDIT_${idSample}.old.vcf > TIDDIT_${idSample}.vcf

    bgzip --threads ${task.cpus} -c TIDDIT_${idSample}.vcf > TIDDIT_${idSample}.vcf.gz

    tabix TIDDIT_${idSample}.vcf.gz
    """
}

vcfTIDDIT = vcfTIDDIT.dump(tag:'TIDDIT')

// STEP FREEBAYES SINGLE MODE

process FreebayesSingle {

    tag "${idSample}-${intervalBed.baseName}"

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    
    input:
        set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamFreebayesSingle
        file(fasta) from ch_fasta
        file(fastaFai) from ch_software_versions_yaml
    
    output:
        set val("FreeBayes"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into vcfFreebayesSingle
    
    when: 'freebayes' in tools

    script:
    intervalsOptions = params.no_intervals ? "" : "-t ${intervalBed}"
    """
    freebayes \
        -f ${fasta} \
        --min-alternate-fraction 0.1 \
        --min-mapping-quality 1 \
        ${intervalsOptions} \
        ${bam} > ${intervalBed.baseName}_${idSample}.vcf
    """
}

vcfFreebayesSingle = vcfFreebayesSingle.groupTuple(by: [0,1,2])


/*
================================================================================
                             SOMATIC VARIANT CALLING
================================================================================
*/

// Ascat, pileup, pileups with no intervals, recalibrated BAMs
(bamAscat, bamMpileup, bamMpileupNoInt, bamRecalAll) = bamRecalAll.into(4)

// separate BAM by status
bamNormal = Channel.create()
bamTumor = Channel.create()

bamRecalAll.choice(bamTumor, bamNormal) {statusMap[it[0], it[1]] == 0 ? 1 : 0}

// Mutect2Single, Mutect2, MSIsensorSingle
(singleBamTumor, singleBamMsisensor, bamTumor) = bamTumor.into(3)

// Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
// Remapping channel to remove common key idPatient
pairBam = bamNormal.cross(bamTumor).map {
    normal, tumor ->
    [normal[0], normal[1], normal[2], normal[3], tumor[1], tumor[2], tumor[3]]
}
pairBam = pairBam.dump(tag:'BAM Somatic Pair')

// Manta, Strelka, MSIsensor, Mutect2
(pairBamManta, pairBamStrelka, pairBamStrelkaBP, pairBamMsisensor, pairBamCNVkit, pairBam) = pairBam.into(6)

// Add the intervals
intervalPairBam = pairBam.combine(bedIntervalsPair)
intervalBam = singleBamTumor.combine(bedIntervalsSingle)

// intervals for Mutect2 calls, FreeBayes and pileups for Mutect2 filtering
(pairBamMutect2, pairBamFreeBayes, bamPileupSummaries) = intervalPairBam.into(3)
(singleBamMutect2, bamPileupSummariesSingle) = intervalBam.into(2)

// intervals for MPileup
bamMpileup = bamMpileup.combine(intMpileup)

// Making Pair Bam for Sention

// separate BAM by status
bam_sention_normal = Channel.create()
bam_sentieon_tumor = Channel.create()

bam_sentieon_all.choice(bam_sentieon_tumor, bam_sention_normal) {statusMap[it[0], it[1]] == 0 ? 1 : 0}

// Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
// Remapping channel to remove common key idPatient

bam_pair_sentieon_TNscope = bam_sention_normal.cross(bam_sentieon_tumor).map {
    normal, tumor ->
    [normal[0], normal[1], normal[2], normal[3], normal[4], tumor[1], tumor[2], tumor[3], tumor[4]]
}

// STEP FREEBAYES

process FreeBayes {

    tag "${idSampleTumor}_vs_${idSampleNormal}-${intervalBed.baseName}"

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamFreeBayes
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set val("FreeBayes"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into vcfFreeBayes

    when: 'freebayes' in tools

    script:
    intervalsOptions = params.no_intervals ? "" : "-t ${intervalBed}"
    """
    freebayes \
        -f ${fasta} \
        --pooled-continuous \
        --pooled-discrete \
        --genotype-qualities \
        --report-genotype-likelihood-max \
        --allele-balance-priors-off \
        --min-alternate-fraction 0.03 \
        --min-repeat-entropy 1 \
        --min-alternate-count 2 \
        ${intervalsOptions} \
        ${bamTumor} \
        ${bamNormal} > ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
    """
}

vcfFreeBayes = vcfFreeBayes.groupTuple(by:[0,1,2])

// STEP GATK MUTECT2.1 - RAW CALLS

process Mutect2 {

    tag "${idSampleTumor}_vs_${idSampleNormal}-${intervalBed.baseName}"
    memory '70G'
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-mediium'

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamMutect2
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(germlineResource) from ch_germline_resource
        file(germlineResourceIndex) from ch_germline_resource_tbi
        file(intervals) from ch_intervals
        file(pon) from ch_pon
        file(ponIndex) from ch_pon_tbi

    output:
        set val("Mutect2"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2PairOutput
        set idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf.stats") optional true into intervalStatsFilesPair
        set idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf.stats"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") optional true into mutect2StatsPair

    when: 'mutect2' in tools

    script:
    // please make a panel-of-normals, using at least 40 samples
    // https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
    PON = params.pon ? "--panel-of-normals ${pon}" : ""
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    softClippedOption = params.ignore_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
    """
    # Get raw calls
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      Mutect2 \
      -R ${fasta}\
      -I ${bamTumor} -tumor ${idSampleTumor} \
      -I ${bamNormal} -normal ${idSampleNormal} \
      ${intervalsOptions} \
      ${softClippedOption} \
      --germline-resource ${germlineResource} \
      ${PON} \
      -O ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
    """
}

mutect2PairOutput = mutect2PairOutput.groupTuple(by:[0,1,2])
mutect2StatsPair = mutect2StatsPair.groupTuple(by:[0,1])

process Mutect2Single {

    tag "${idSampleTumor}-${intervalBed.baseName}"
    memory '70G'
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    input:
        set idPatient, idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from singleBamMutect2
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(germlineResource) from ch_germline_resource
        file(germlineResourceIndex) from ch_germline_resource_tbi
        file(intervals) from ch_intervals
        file(pon) from ch_pon
        file(ponIndex) from ch_pon_tbi

    output:
        set val("Mutect2"), idPatient, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}.vcf") into mutect2SingleOutput
        set idPatient, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}.vcf.stats") optional true into intervalStatsFilesSingle
        set idPatient, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}.vcf.stats"), file("${intervalBed.baseName}_${idSampleTumor}.vcf") optional true into mutect2StatsSingle

    when: 'mutect2' in tools

    script:
    // please make a panel-of-normals, using at least 40 samples
    // https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
    PON = params.pon ? "--panel-of-normals ${pon}" : ""
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    softClippedOption = params.ignore_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
    """
    # Get raw calls
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      Mutect2 \
      -R ${fasta}\
      -I ${bamTumor}  -tumor ${idSampleTumor} \
      ${intervalsOptions} \
      ${softClippedOption} \
      --germline-resource ${germlineResource} \
      ${PON} \
      -O ${intervalBed.baseName}_${idSampleTumor}.vcf
    """
}

mutect2SingleOutput = mutect2SingleOutput.groupTuple(by:[0,1,2])
mutect2StatsSingle = mutect2StatsSingle.groupTuple(by:[0,1])

// STEP GATK MUTECT2.2 - MERGING STATS

mutect2Stats = mutect2StatsSingle.mix(mutect2StatsPair)
mutect2Stats = mutect2Stats.dump(tag:'Mutect2 Stats')

process MergeMutect2Stats {

    tag "${idSample}"
    memory '70G'
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-mediium'
    publishDir "${params.outdir}/VariantCalling/${idSample}/Mutect2", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(statsFiles), file(vcf) from mutect2Stats // Actual stats files and corresponding VCF chunks
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(germlineResource) from ch_germline_resource
        file(germlineResourceIndex) from ch_germline_resource_tbi
        file(intervals) from ch_intervals

    output:
        set idPatient, idSample, file("${idSample}.vcf.gz.stats") into mergedStatsFile

    when: 'mutect2' in tools

    script:
               stats = statsFiles.collect{ "-stats ${it} " }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        MergeMutectStats \
        ${stats} \
        -O ${idSample}.vcf.gz.stats
    """
}

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller

// STEP MERGING VCF - FREEBAYES & GATK HAPLOTYPECALLER

vcfConcatenateVCFs = vcfFreeBayes.mix(vcfFreebayesSingle, vcfGenotypeGVCFs, gvcfHaplotypeCaller)
vcfConcatenateVCFs = vcfConcatenateVCFs.dump(tag:'VCF to merge')

process ConcatVCF {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'

    tag "${variantCaller}-${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/${"$variantCaller"}", mode: params.publish_dir_mode

    input:
        set variantCaller, idPatient, idSample, file(vcf) from vcfConcatenateVCFs
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        set variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenated

    when: ('haplotypecaller' in tools || 'freebayes' in tools)

    script:
    if (variantCaller == 'HaplotypeCallerGVCF')
        outputFile = "HaplotypeCaller_${idSample}.g.vcf"
    else
        outputFile = "${variantCaller}_${idSample}.vcf"
    options = params.target_bed ? "-t ${targetBED}" : ""
    intervalsOptions = params.no_intervals ? "-n" : ""
    """
    concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
    """
}

vcfConcatenated = vcfConcatenated.dump(tag:'VCF')

// STEP MERGING VCF - MUTECT2

mutect2Output = mutect2SingleOutput.mix(mutect2PairOutput)
mutect2Output = mutect2Output.dump(tag:'Mutect2 output VCF to merge')

process ConcatVCF_Mutect2 {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'

    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Mutect2", mode: params.publish_dir_mode

    input:
        set variantCaller, idPatient, idSample, file(vcf) from mutect2Output
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        set variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenatedForFilter

    when: 'mutect2' in tools

    script:
    outputFile = "Mutect2_unfiltered_${idSample}.vcf"
    options = params.target_bed ? "-t ${targetBED}" : ""
    intervalsOptions = params.no_intervals ? "-n" : ""
    """
    concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
    """
}

vcfConcatenatedForFilter = vcfConcatenatedForFilter.dump(tag:'Mutect2 unfiltered VCF')

// STEP GATK MUTECT2.3 - GENERATING PILEUP SUMMARIES

bamPileupSummariesSingle = bamPileupSummariesSingle.join(intervalStatsFilesSingle, by:[0,1])

bamPileupSummaries = bamPileupSummaries.map{
  idPatient, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor, intervalBed ->
  [idPatient, idSampleTumor + "_vs_" + idSampleNormal, bamTumor, baiTumor, intervalBed]
}.join(intervalStatsFilesPair, by:[0,1])

bamPileupSummaries = bamPileupSummaries.mix(bamPileupSummariesSingle)
bamPileupSummaries = bamPileupSummaries.dump(tag:'Mutect2 Pileup Summaries')

process PileupSummariesForMutect2 {

    container pipeline_container
    memory '70G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-mediium'

    input:
        set idPatient, idSample, file(bamTumor), file(baiTumor), file(intervalBed), file(statsFile) from bamPileupSummaries
        file(germlineResource) from ch_germline_resource
        file(germlineResourceIndex) from ch_germline_resource_tbi

    output:
        set idPatient, idSample, file("${intervalBed.baseName}_${idSample}_pileupsummaries.table") into pileupSummaries

    when: 'mutect2' in tools

    script:
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        GetPileupSummaries \
        -I ${bamTumor} \
        -V ${germlineResource} \
        ${intervalsOptions} \
        -O ${intervalBed.baseName}_${idSample}_pileupsummaries.table
    """
}

pileupSummaries = pileupSummaries.groupTuple(by:[0,1])

// STEP GATK MUTECT2.4 - MERGING PILEUP SUMMARIES

process MergePileupSummaries {

    container pipeline_container
   memory '60G' 
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-mediium'

    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Mutect2", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(pileupSums) from pileupSummaries
        file(dict) from ch_dict

    output:
        set idPatient, idSample, file("${idSample}_pileupsummaries.table") into mergedPileupFile

    when: 'mutect2' in tools

    script:
    allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        GatherPileupSummaries \
        --sequence-dictionary ${dict} \
        ${allPileups} \
        -O ${idSample}_pileupsummaries.table
    """
}

// STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

process CalculateContamination {

    container pipeline_container
    memory '16G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Mutect2", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(mergedPileup) from mergedPileupFile

     output:
        set idPatient, val("${idSample}"), file("${idSample}_contamination.table") into contaminationTable

    when: 'mutect2' in tools

    script:   
    """
    # calculate contamination
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CalculateContamination \
        -I ${idSample}_pileupsummaries.table \
        -O ${idSample}_contamination.table
    """
}


// STEP GATK MUTECT2.6 - FILTERING CALLS

mutect2CallsToFilter = vcfConcatenatedForFilter.map{
    variantCaller, idPatient, idSample, vcf, tbi ->
    [idPatient, idSample, vcf, tbi]
}.join(mergedStatsFile, by:[0,1]).join(contaminationTable, by:[0,1])

process FilterMutect2Calls {

    container pipeline_container
    memory '16G'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Mutect2", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(unfiltered), file(unfilteredIndex), file(stats), file(contaminationTable) from mutect2CallsToFilter
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(germlineResource) from ch_germline_resource
        file(germlineResourceIndex) from ch_germline_resource_tbi
        file(intervals) from ch_intervals

    output:
        set val("Mutect2"), idPatient, idSample, file("Mutect2_filtered_${idSample}.vcf.gz"), file("Mutect2_filtered_${idSample}.vcf.gz.tbi"), file("Mutect2_filtered_${idSample}.vcf.gz.filteringStats.tsv") into filteredMutect2Output

    when: 'mutect2' in tools

    script:
    """
    # do the actual filtering
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        FilterMutectCalls \
        -V ${unfiltered} \
        --contamination-table ${contaminationTable} \
        --stats ${stats} \
        -R ${fasta} \
        -O Mutect2_filtered_${idSample}.vcf.gz
    """
}

// STEP SENTIEON TNSCOPE

process Sentieon_TNscope {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idSampleTumor}_vs_${idSampleNormal}"

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), file(recalNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(recalTumor) from bam_pair_sentieon_TNscope
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(pon) from ch_pon
        file(ponIndex) from ch_pon_tbi

    output:
        set val("SentieonTNscope"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf") into vcf_sentieon_TNscope

    when: 'tnscope' in tools && params.sentieon

    script:
    PON = params.pon ? "--pon ${pon}" : ""
    """
    sentieon driver \
        -t ${task.cpus} \
        -r ${fasta} \
        -i ${bamTumor} \
        -q ${recalTumor} \
        -i ${bamNormal} \
        -q ${recalNormal} \
        --algo TNscope \
        --tumor_sample ${idSampleTumor} \
        --normal_sample ${idSampleNormal} \
        --dbsnp ${dbsnp} \
        ${PON} \
        TNscope_${idSampleTumor}_vs_${idSampleNormal}.vcf
    """
}

vcf_sentieon_TNscope = vcf_sentieon_TNscope.dump(tag:'Sentieon TNscope')

vcf_sentieon = vcf_sentieon_DNAseq.mix(vcf_sentieon_DNAscope, vcf_sentieon_DNAscope_SV, vcf_sentieon_TNscope)

process CompressSentieonVCF {

    tag "${idSample} - ${vcf}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    publishDir "${params.outdir}/VariantCalling/${idSample}/${variantCaller}", mode: params.publish_dir_mode

    input:
        set variantCaller, idPatient, idSample, file(vcf) from vcf_sentieon

    output:
        set variantCaller, idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcf_sentieon_compressed

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}

vcf_sentieon_compressed = vcf_sentieon_compressed.dump(tag:'Sentieon VCF indexed')

// STEP STRELKA.2 - SOMATIC PAIR

process Strelka {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSampleTumor}_vs_${idSampleNormal}"

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Strelka", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamStrelka
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
        set val("Strelka"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelka

    when: 'strelka' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configureStrelkaSomaticWorkflow.py \
        --tumor ${bamTumor} \
        --normal ${bamNormal} \
        --referenceFasta ${fasta} \
        ${options} \
        --runDir Strelka

    python Strelka/runWorkflow.py -m local -j ${task.cpus}

    mv Strelka/results/variants/somatic.indels.vcf.gz \
        Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
    mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
        Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
    mv Strelka/results/variants/somatic.snvs.vcf.gz \
        Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
    mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
        Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
    """
}

vcfStrelka = vcfStrelka.dump(tag:'Strelka')

// STEP MANTA.2 - SOMATIC PAIR

process Manta {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSampleTumor}_vs_${idSampleNormal}"

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Manta", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamManta
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
        set val("Manta"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfManta
        set idPatient, idSampleNormal, idSampleTumor, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

    when: 'manta' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configManta.py \
        --normalBam ${bamNormal} \
        --tumorBam ${bamTumor} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz.tbi
    mv Manta/results/variants/somaticSV.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz
    mv Manta/results/variants/somaticSV.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz.tbi
    """
}

vcfManta = vcfManta.dump(tag:'Manta')

// Remmaping channels to match input for StrelkaBP
pairBamStrelkaBP = pairBamStrelkaBP.map {
    idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor ->
    [idPatientNormal, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor]
}.join(mantaToStrelka, by:[0,1,2]).map {
    idPatientNormal, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor, mantaCSI, mantaCSIi ->
    [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor, mantaCSI, mantaCSIi]
}

// STEP STRELKA.3 - SOMATIC PAIR - BEST PRACTICES

process StrelkaBP {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idSampleTumor}_vs_${idSampleNormal}"

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Strelka", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(mantaCSI), file(mantaCSIi) from pairBamStrelkaBP
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
        set val("Strelka"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelkaBP

    when: 'strelka' in tools && 'manta' in tools && !params.no_strelka_bp

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configureStrelkaSomaticWorkflow.py \
        --tumor ${bamTumor} \
        --normal ${bamNormal} \
        --referenceFasta ${fasta} \
        --indelCandidates ${mantaCSI} \
        ${options} \
        --runDir Strelka

    python Strelka/runWorkflow.py -m local -j ${task.cpus}

    mv Strelka/results/variants/somatic.indels.vcf.gz \
        StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
    mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
        StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
    mv Strelka/results/variants/somatic.snvs.vcf.gz \
        StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
    mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
        StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
    """
}

vcfStrelkaBP = vcfStrelkaBP.dump(tag:'Strelka BP')

// STEP CNVkit

process CNVkit {

    tag "${idSampleTumor}_vs_${idSampleNormal}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/CNVkit", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamCNVkit
        file(targetBED) from ch_target_bed
        file(fasta) from ch_fasta

    output:
        set idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}*"), file("${idSampleNormal}*") into cnvkitOut

    when: 'cnvkit' in tools && params.target_bed

    script:
    """
    cnvkit.py \
      batch \
      ${bamTumor} \
      --normal ${bamNormal} \
      --targets ${targetBED} \
      --fasta ${fasta} \
      --output-reference output_reference.cnn \
      --output-dir ./ \
      --diagram \
      --scatter 
    """
}

// STEP MSISENSOR.1 - SCAN

// Scan reference genome for microsatellites
process MSIsensor_scan {

    container msisensor_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    
    tag "${fasta}"
    
    input:
    file(fasta) from ch_fasta
    file(fastaFai) from ch_fai

    output:
    file "microsatellites.list" into msi_scan_ch

    when: 'msisensor' in tools

    script:
    """
    msisensor-pro scan -d ${fasta} -o microsatellites.list
    """
}

// STEP MSISENSOR.2 - SCORE

// Score the normal vs somatic pair of bams

process MSIsensor_msi {

    container msisensor_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
        
    tag "${idSampleTumor}_vs_${idSampleNormal}"
    
    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/MSIsensor", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamMsisensor
        file msiSites from msi_scan_ch

    output:
        set val("Msisensor"), idPatient, file("${idSampleTumor}_vs_${idSampleNormal}_msisensor"), file("${idSampleTumor}_vs_${idSampleNormal}_msisensor_dis"), file("${idSampleTumor}_vs_${idSampleNormal}_msisensor_germline"), file("${idSampleTumor}_vs_${idSampleNormal}_msisensor_somatic") into msisensor_out_ch

    when: 'msisensor' in tools

    script:
    """
    msisensor-pro msi -d ${msiSites} \
                      -b 4 \
                      -n ${bamNormal} \
                      -t ${bamTumor} \
                      -o ${idSampleTumor}_vs_${idSampleNormal}_msisensor
    """
}

// TODO: Add the msisensor-pro baseline step
process MSIsensor_msiSingle {

    container msisensor_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
        
    tag "${idSampleTumor}"
    
    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/MSIsensor", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleTumor, file(bamTumor), file(baiTumor) from singleBamMsisensor
        file msiSites from msi_scan_ch

    output:
        set val("MsisensorSingle"), idPatient, file("${idSampleTumor}_msisensor"), file("${idSampleTumor}_msisensor_dis"), file("${idSampleTumor}_msisensor_all") into msisensor_out_ch_single

    when: 'msisensor' in tools

    script:
    """
    msisensor-pro pro -d ${msiSites} \
                      -b 4 \
                      -t ${bamTumor} \
                      -o ${idSampleTumor}_msisensor
    """
}

// STEP ASCAT.1 - ALLELECOUNTER

// Run commands and code from Malin Larsson
// Based on Jesper Eisfeldt's code
process AlleleCounter {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSample}"

    input:
        set idPatient, idSample, file(bam), file(bai) from bamAscat
        file(acLoci) from ch_ac_loci
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, file("${idSample}.alleleCount") into alleleCounterOut

    when: 'ascat' in tools

    script:
    """
    alleleCounter \
        -l ${acLoci} \
        -r ${fasta} \
        -b ${bam} \
        -o ${idSample}.alleleCount;
    """
}

alleleCountOutNormal = Channel.create()
alleleCountOutTumor = Channel.create()

alleleCounterOut
    .choice(alleleCountOutTumor, alleleCountOutNormal) {statusMap[it[0], it[1]] == 0 ? 1 : 0}

alleleCounterOut = alleleCountOutNormal.combine(alleleCountOutTumor, by:0)

alleleCounterOut = alleleCounterOut.map {
    idPatientNormal, idSampleNormal, alleleCountOutNormal,
    idSampleTumor, alleleCountOutTumor ->
    [idPatientNormal, idSampleNormal, idSampleTumor, alleleCountOutNormal, alleleCountOutTumor]
}

// STEP ASCAT.2 - CONVERTALLELECOUNTS

// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process ConvertAlleleCounts {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSampleTumor}_vs_${idSampleNormal}"

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/ASCAT", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleNormal, idSampleTumor, file(alleleCountNormal), file(alleleCountTumor) from alleleCounterOut

    output:
        set idPatient, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convertAlleleCountsOut

    when: 'ascat' in tools

    script:
    gender = genderMap[idPatient]
    """
    convertAlleleCounts.r ${idSampleTumor} ${alleleCountTumor} ${idSampleNormal} ${alleleCountNormal} ${gender}
    """
}

// STEP ASCAT.3 - ASCAT

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process Ascat {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSampleTumor}_vs_${idSampleNormal}"

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/ASCAT", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleNormal, idSampleTumor, file(bafNormal), file(logrNormal), file(bafTumor), file(logrTumor) from convertAlleleCountsOut
        file(acLociGC) from ch_ac_loci_gc

    output:
        set val("ASCAT"), idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.*.{png,txt}") into ascatOut

    when: 'ascat' in tools

    script:
    gender = genderMap[idPatient]
    purity_ploidy = (params.ascat_purity && params.ascat_ploidy) ? "--purity ${params.ascat_purity} --ploidy ${params.ascat_ploidy}" : ""
    """
    for f in *BAF *LogR; do sed 's/chr//g' \$f > tmpFile; mv tmpFile \$f;done
    run_ascat.r \
        --tumorbaf ${bafTumor} \
        --tumorlogr ${logrTumor} \
        --normalbaf ${bafNormal} \
        --normallogr ${logrNormal} \
        --tumorname ${idSampleTumor} \
        --basedir ${projectDir} \
        --gcfile ${acLociGC} \
        --gender ${gender} \
        ${purity_ploidy}
    """
}

ascatOut.dump(tag:'ASCAT')

// STEP MPILEUP.1

process Mpileup {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'

    tag "${idSample}-${intervalBed.baseName}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup" ? "VariantCalling/${idSample}/Control-FREEC/${it}" : null }

    input:
        set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamMpileup
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, file("${prefix}${idSample}.pileup") into mpileupMerge
        set idPatient, idSample into tsv_mpileup

    when: 'controlfreec' in tools || 'mpileup' in tools

    script:
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-l ${intervalBed}"

    """
    # Control-FREEC reads uncompresses the zipped file TWICE in single-threaded mode.
    # we are therefore not using compressed pileups here
    samtools mpileup \
        -f ${fasta} ${bam} \
        ${intervalsOptions} > ${prefix}${idSample}.pileup
    """
}

(tsv_mpileup, tsv_mpileup_sample) = tsv_mpileup.groupTuple(by:[0, 1]).into(2)

// Creating a TSV file to restart from this step
tsv_mpileup.map { idPatient, idSample ->
    gender = genderMap[idPatient]
    status = statusMap[idPatient, idSample]
    mpileup = "${params.outdir}/VariantCalling/${idSample}/Control-FREEC/${idSample}.pileup"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${mpileup}\n"
}.collectFile(
    name: 'control-freec_mpileup.tsv', sort: true, storeDir: "${params.outdir}/VariantCalling/TSV"
)

tsv_mpileup_sample
    .collectFile(storeDir: "${params.outdir}/VariantCalling/TSV") {
        idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        mpileup = "${params.outdir}/VariantCalling/${idSample}/Control-FREEC/${idSample}.pileup"
        ["control-freec_mpileup_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${mpileup}\n"]
}

if (!params.no_intervals) {
    mpileupMerge = mpileupMerge.groupTuple(by:[0, 1])
    mpileupNoInt = Channel.empty()
} else {
    (mpileupMerge, mpileupNoInt) = mpileupMerge.into(2)
    mpileupMerge.close()
}

// STEP MPILEUP.2 - MERGE
process MergeMpileup {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    tag "${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup" ? "VariantCalling/${idSample}/Control-FREEC/${it}" : null }

    input:
        set idPatient, idSample, file(mpileup) from mpileupMerge

    output:
        set idPatient, idSample, file("${idSample}.pileup") into mpileupOut

    when: !(params.no_intervals) && 'controlfreec' in tools || 'mpileup' in tools

    script:
    """
    for i in `ls -1v *.pileup`;
        do cat \$i >> ${idSample}.pileup
    done
    """
}

mpileupOut = mpileupOut.mix(mpileupNoInt)
mpileupOut = mpileupOut.dump(tag:'mpileup')

mpileupOutNormal = Channel.create()
mpileupOutTumor = Channel.create()

if (step == 'controlfreec') mpileupOut = inputSample

mpileupOut
    .choice(mpileupOutTumor, mpileupOutNormal) {statusMap[it[0], it[1]] == 0 ? 1 : 0}

(mpileupOutSingle,mpileupOutTumor) = mpileupOutTumor.into(2)

mpileupOut = mpileupOutNormal.combine(mpileupOutTumor, by:0)

mpileupOut = mpileupOut.map {
    idPatientNormal, idSampleNormal, mpileupOutNormal,
    idSampleTumor, mpileupOutTumor ->
    [idPatientNormal, idSampleNormal, idSampleTumor, mpileupOutNormal, mpileupOutTumor]
}

// STEP CONTROLFREEC.1 - CONTROLFREEC

process ControlFREEC {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSampleTumor}_vs_${idSampleNormal}"

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Control-FREEC", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleNormal, idSampleTumor, file(mpileupNormal), file(mpileupTumor) from mpileupOut
        file(chrDir) from ch_chr_dir
        file(mappability) from ch_mappability
        file(chrLength) from ch_chr_length
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
        set idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.pileup_CNVs"), file("${idSampleTumor}.pileup_ratio.txt"), file("${idSampleTumor}.pileup_BAF.txt") into controlFreecViz
        set file("*.pileup*"), file("${idSampleTumor}_vs_${idSampleNormal}.config.txt") into controlFreecOut

    when: 'controlfreec' in tools

    script:
    config = "${idSampleTumor}_vs_${idSampleNormal}.config.txt"
    gender = genderMap[idPatient]
    // Window has higher priority than coefficientOfVariation if both given
    window = params.cf_window ? "window = ${params.cf_window}" : ""
    coeffvar = params.cf_coeff ? "coefficientOfVariation = ${params.cf_coeff}" : ""
    use_bed = params.target_bed ? "captureRegions = ${targetBED}" : ""
    // This parameter makes Control-FREEC unstable (still in Beta according to the developers)
    // so we disable it by setting it to its default value (it is disabled by default)
    //min_subclone = params.target_bed ? "30" : "20"
    min_subclone = 100
    readCountThreshold = params.target_bed ? "50" : "10"
    breakPointThreshold = params.target_bed ? "1.2" : "0.8"
    breakPointType = params.target_bed ? "4" : "2"
    mappabilitystr = params.mappability ? "gemMappabilityFile = \${PWD}/${mappability}" : ""
    contamination_adjustment = params.cf_contamination_adjustment ? "contaminationAdjustment = TRUE" : ""
    contamination_value = params.cf_contamination ? "contamination = ${params.cf_contamination}" : ""
    """
    touch ${config}
    echo "[general]" >> ${config}
    echo "BedGraphOutput = TRUE" >> ${config}
    echo "chrFiles = \${PWD}/${chrDir.fileName}" >> ${config}
    echo "chrLenFile = \${PWD}/${chrLength.fileName}" >> ${config}
    echo "forceGCcontentNormalization = 1" >> ${config}
    echo "maxThreads = ${task.cpus}" >> ${config}
    echo "minimalSubclonePresence = ${min_subclone}" >> ${config}
    echo "ploidy = ${params.cf_ploidy}" >> ${config}
    echo "sex = ${gender}" >> ${config}
    echo "readCountThreshold = ${readCountThreshold}" >> ${config}
    echo "breakPointThreshold = ${breakPointThreshold}" >> ${config}
    echo "breakPointType = ${breakPointType}" >> ${config}
    echo "${window}" >> ${config}
    echo "${coeffvar}" >> ${config}
    echo "${mappabilitystr}" >> ${config}
    echo "${contamination_adjustment}" >> ${config}
    echo "${contamination_value}" >> ${config}
    echo "" >> ${config}
    
    echo "[control]" >> ${config}
    echo "inputFormat = pileup" >> ${config}
    echo "mateFile = \${PWD}/${mpileupNormal}" >> ${config}
    echo "mateOrientation = FR" >> ${config}
    echo "" >> ${config}

    echo "[sample]" >> ${config}
    echo "inputFormat = pileup" >> ${config}
    echo "mateFile = \${PWD}/${mpileupTumor}" >> ${config}
    echo "mateOrientation = FR" >> ${config}
    echo "" >> ${config}

    echo "[BAF]" >> ${config}
    echo "SNPfile = ${dbsnp.fileName}" >> ${config}
    echo "" >> ${config}

    echo "[target]" >> ${config}
    echo "${use_bed}" >> ${config}

    freec -conf ${config}
    """
}

controlFreecOut.dump(tag:'ControlFREEC')

process ControlFREECSingle {
    container pipeline_container

    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSampleTumor}"

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Control-FREEC", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleTumor, file(mpileupTumor) from mpileupOutSingle
        file(chrDir) from ch_chr_dir
        file(mappability) from ch_mappability
        file(chrLength) from ch_chr_length
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
        set idPatient, idSampleTumor, file("${idSampleTumor}.pileup_CNVs"), file("${idSampleTumor}.pileup_ratio.txt"), file("${idSampleTumor}.pileup_BAF.txt") into controlFreecVizSingle
        set file("*.pileup*"), file("${idSampleTumor}.config.txt") into controlFreecOutSingle

    when: 'controlfreec' in tools

    script:
    config = "${idSampleTumor}.config.txt"
    gender = genderMap[idPatient]
    // Window has higher priority than coefficientOfVariation if both given
    window = params.cf_window ? "window = ${params.cf_window}" : ""
    coeffvar = params.cf_coeff ? "coefficientOfVariation = ${params.cf_coeff}" : ""
    use_bed = params.target_bed ? "captureRegions = ${targetBED}" : ""
    // This parameter makes Control-FREEC unstable (still in Beta according to the developers)
    // so we disable it by setting it to its default value (it is disabled by default)
    //min_subclone = params.target_bed ? "30" : "20"
    min_subclone = 100
    readCountThreshold = params.target_bed ? "50" : "10"
    breakPointThreshold = params.target_bed ? "1.2" : "0.8"
    breakPointType = params.target_bed ? "4" : "2"
    mappabilitystr = params.mappability ? "gemMappabilityFile = \${PWD}/${mappability}" : ""
    contamination_adjustment = params.cf_contamination_adjustment ? "contaminationAdjustment = TRUE" : ""
    contamination_value = params.cf_contamination ? "contamination = ${params.cf_contamination}" : ""
    """
    touch ${config}
    echo "[general]" >> ${config}
    echo "BedGraphOutput = TRUE" >> ${config}
    echo "chrFiles = \${PWD}/${chrDir.fileName}" >> ${config}
    echo "chrLenFile = \${PWD}/${chrLength.fileName}" >> ${config}
    echo "forceGCcontentNormalization = 1" >> ${config}
    echo "maxThreads = ${task.cpus}" >> ${config}
    echo "minimalSubclonePresence = ${min_subclone}" >> ${config}
    echo "ploidy = ${params.cf_ploidy}" >> ${config}
    echo "sex = ${gender}" >> ${config}
    echo "readCountThreshold = ${readCountThreshold}" >> ${config}
    echo "breakPointThreshold = ${breakPointThreshold}" >> ${config}
    echo "breakPointType = ${breakPointType}" >> ${config}
    echo "${window}" >> ${config}
    echo "${coeffvar}" >> ${config}
    echo "${mappabilitystr}" >> ${config}
    echo "${contamination_adjustment}" >> ${config}
    echo "${contamination_value}" >> ${config}
    echo "" >> ${config}

    echo "[sample]" >> ${config}
    echo "inputFormat = pileup" >> ${config}
    echo "mateFile = \${PWD}/${mpileupTumor}" >> ${config}
    echo "mateOrientation = FR" >> ${config}
    echo "" >> ${config}

    echo "[BAF]" >> ${config}
    echo "SNPfile = ${dbsnp.fileName}" >> ${config}
    echo "" >> ${config}

    echo "[target]" >> ${config}
    echo "${use_bed}" >> ${config}

    freec -conf ${config}
    """
}

controlFreecOutSingle.dump(tag:'ControlFREECSingle')

// STEP CONTROLFREEC.3 - VISUALIZATION

process ControlFreecViz {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSampleTumor}_vs_${idSampleNormal}"

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Control-FREEC", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleNormal, idSampleTumor, file(cnvTumor), file(ratioTumor), file(bafTumor) from controlFreecViz

    output:
        set file("*.txt"), file("*.png"), file("*.bed") into controlFreecVizOut

    when: 'controlfreec' in tools

    script:
    """
    echo "############### Calculating significance values for TUMOR CNVs #############"
    cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/assess_significance.R | R --slave --args ${cnvTumor} ${ratioTumor}

    echo "############### Creating graph for TUMOR ratios ###############"
    cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/makeGraph.R | R --slave --args 2 ${ratioTumor} ${bafTumor}

    echo "############### Creating BED files for TUMOR ##############"
    perl /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/freec2bed.pl -f ${ratioTumor} > ${idSampleTumor}.bed
    """
}

controlFreecVizOut.dump(tag:'ControlFreecViz')

process ControlFreecVizSingle {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "${idSampleTumor}"

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Control-FREEC", mode: params.publish_dir_mode

    input:
        set idPatient, idSampleTumor, file(cnvTumor), file(ratioTumor), file(bafTumor) from controlFreecVizSingle

    output:
        set file("*.txt"), file("*.png"), file("*.bed") into controlFreecVizOutSingle

    when: 'controlfreec' in tools

    script:
    """
    echo "############### Calculating significance values for TUMOR CNVs #############"
    cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/assess_significance.R | R --slave --args ${cnvTumor} ${ratioTumor}

    echo "############### Creating graph for TUMOR ratios ###############"
    cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/makeGraph.R | R --slave --args 2 ${ratioTumor} ${bafTumor}

    echo "############### Creating BED files for TUMOR ##############"
    perl /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/freec2bed.pl -f ${ratioTumor} > ${idSampleTumor}.bed
    """
}

controlFreecVizOutSingle.dump(tag:'ControlFreecVizSingle')

// Remapping channels for QC and annotation

(vcfStrelkaIndels, vcfStrelkaSNVS) = vcfStrelka.into(2)
(vcfStrelkaBPIndels, vcfStrelkaBPSNVS) = vcfStrelkaBP.into(2)
(vcfMantaSomaticSV, vcfMantaDiploidSV) = vcfManta.into(2)

vcfKeep = Channel.empty().mix(
    filteredMutect2Output.map{
        variantCaller, idPatient, idSample, vcf, tbi, tsv ->
        [variantCaller, idSample, vcf]
    },
    vcfConcatenated.map{
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf]
    },
    vcf_sentieon_compressed.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf]
    },
    vcfStrelkaSingle.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf[1]]
    },
    vcfMantaSingle.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf[2]]
    },
    vcfMantaDiploidSV.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf[2]]
    },
    vcfMantaSomaticSV.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf[3]]
    },
    vcfStrelkaIndels.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf[0]]
    },
    vcfStrelkaSNVS.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf[1]]
    },
    vcfStrelkaBPIndels.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf[0]]
    },
    vcfStrelkaBPSNVS.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf[1]]
    },
    vcfTIDDIT.map {
        variantCaller, idPatient, idSample, vcf, tbi ->
        [variantCaller, idSample, vcf]
    })

(vcfBCFtools, vcfVCFtools, vcfAnnotation) = vcfKeep.into(3)

// STEP VCF.QC

process BcftoolsStats {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    tag "${variantCaller} - ${vcf}"

    publishDir "${params.outdir}/Reports/${idSample}/BCFToolsStats", mode: params.publish_dir_mode

    input:
        set variantCaller, idSample, file(vcf) from vcfBCFtools

    output:
        file ("*.bcf.tools.stats.out") into bcftoolsReport

    when: !('bcftools' in skipQC)

    script:
    """
    bcftools stats ${vcf} > ${reduceVCF(vcf.fileName)}.bcf.tools.stats.out
    """
}

bcftoolsReport = bcftoolsReport.dump(tag:'BCFTools')

process Vcftools {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    tag "${variantCaller} - ${vcf}"

    publishDir "${params.outdir}/Reports/${idSample}/VCFTools", mode: params.publish_dir_mode

    input:
        set variantCaller, idSample, file(vcf) from vcfVCFtools

    output:
        file ("${reduceVCF(vcf.fileName)}.*") into vcftoolsReport

    when: !('vcftools' in skipQC)

    script:
    """
    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-count \
    --out ${reduceVCF(vcf.fileName)}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-qual \
    --out ${reduceVCF(vcf.fileName)}

    vcftools \
    --gzvcf ${vcf} \
    --FILTER-summary \
    --out ${reduceVCF(vcf.fileName)}
    """
}

vcftoolsReport = vcftoolsReport.dump(tag:'VCFTools')

/*
================================================================================
                                   ANNOTATION
================================================================================
*/

if (step == 'annotate') {
    vcfToAnnotate = Channel.create()
    vcfNoAnnotate = Channel.create()

    if (tsvPath == []) {
    // Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
    // Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
    // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,SentieonDNAseq,SentieonDNAscope,SentieonTNscope,Strelka,TIDDIT}/*.vcf.gz
    // Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
    // The small snippet `vcf.minus(vcf.fileName)[-2]` catches idSample
    // This field is used to output final annotated VCFs in the correct directory
      Channel.empty().mix(
        Channel.fromPath("${params.outdir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
          .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
          .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Mutect2/*.vcf.gz")
          .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAseq/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonDNAseq', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAscope/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonDNAscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonTNscope/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonTNscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
          .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/TIDDIT/*.vcf.gz")
          .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
      ).choice(vcfToAnnotate, vcfNoAnnotate) {
        annotate_tools == [] || (annotate_tools != [] && it[0] in annotate_tools) ? 0 : 1
      }
    } else if (annotate_tools == []) {
    // Annotate user-submitted VCFs
    // If user-submitted, Sarek assume that the idSample should be assumed automatically
      vcfToAnnotate = Channel.fromPath(tsvPath)
        .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
    } else exit 1, "specify only tools or files to annotate, not both"

    vcfNoAnnotate.close()
    vcfAnnotation = vcfAnnotation.mix(vcfToAnnotate)
}

// as now have the list of VCFs to annotate, the first step is to annotate with allele frequencies, if there are any

(vcfSnpeff, vcfVep) = vcfAnnotation.into(2)

vcfVep = vcfVep.map {
  variantCaller, idSample, vcf ->
  [variantCaller, idSample, vcf, null]
}

// STEP SNPEFF

process Snpeff {

    tag "${idSample} - ${variantCaller} - ${vcf}"
    memory '70G'
    container snpeff_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-mediium'
    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_snpEff.ann.vcf") null
        else "Reports/${idSample}/snpEff/${it}"
    }

    input:
        set variantCaller, idSample, file(vcf) from vcfSnpeff
        file(dataDir) from ch_snpeff_cache
        val snpeffDb from ch_snpeff_db

    output:
        set file("${reducedVCF}_snpEff.genes.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv") into snpeffReport
        set variantCaller, idSample, file("${reducedVCF}_snpEff.ann.vcf") into snpeffVCF

    when: 'snpeff' in tools || 'merge' in tools

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    cache = (params.snpeff_cache && params.annotation_cache) ? "-dataDir \${PWD}/${dataDir}" : ""
    """
    snpEff -Xmx${task.memory.toGiga()}g \
        ${snpeffDb} \
        -csvStats ${reducedVCF}_snpEff.csv \
        -nodownload \
        ${cache} \
        -canon \
        -v \
        ${vcf} \
        > ${reducedVCF}_snpEff.ann.vcf

    mv snpEff_summary.html ${reducedVCF}_snpEff.html
    """
}

snpeffReport = snpeffReport.dump(tag:'snpEff report')

// STEP COMPRESS AND INDEX VCF.1 - SNPEFF

process CompressVCFsnpEff {

    tag "${idSample} - ${vcf}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir "${params.outdir}/Annotation/${idSample}/snpEff", mode: params.publish_dir_mode

    input:
        set variantCaller, idSample, file(vcf) from snpeffVCF

    output:
        set variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into (compressVCFsnpEffOut)

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}

compressVCFsnpEffOut = compressVCFsnpEffOut.dump(tag:'VCF')

// STEP VEP.1

process VEP {

    container vep_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'

    tag "${idSample} - ${variantCaller} - ${vcf}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
        else null
    }

    input:
        set variantCaller, idSample, file(vcf), file(idx) from vcfVep
        file(dataDir) from ch_vep_cache
        val cache_version from ch_vep_cache_version
        file(cadd_InDels) from ch_cadd_indels
        file(cadd_InDels_tbi) from ch_cadd_indels_tbi
        file(cadd_WG_SNVs) from ch_cadd_wg_snvs
        file(cadd_WG_SNVs_tbi) from ch_cadd_wg_snvs_tbi

    output:
        set variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf") into vepVCF
        file("${reducedVCF}_VEP.summary.html") into vepReport

    when: 'vep' in tools

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome

    dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
    cadd = (params.cadd_cache && params.cadd_wg_snvs && params.cadd_indels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${params.species} \
        ${cadd} \
        ${genesplicer} \
        --cache \
        --cache_version ${cache_version} \
        --dir_cache ${dir_cache} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${reducedVCF}_VEP.summary.html \
        --total_length \
        --vcf

    rm -rf ${reducedVCF}
    """
}

vepReport = vepReport.dump(tag:'VEP')

// STEP VEP.2 - VEP AFTER SNPEFF

process VEPmerge {

    container vep_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'
    tag "${idSample} - ${variantCaller} - ${vcf}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
        else null
    }

    input:
        set variantCaller, idSample, file(vcf), file(idx) from compressVCFsnpEffOut
        file(dataDir) from ch_vep_cache
        val cache_version from ch_vep_cache_version
        file(cadd_InDels) from ch_cadd_indels
        file(cadd_InDels_tbi) from ch_cadd_indels_tbi
        file(cadd_WG_SNVs) from ch_cadd_wg_snvs
        file(cadd_WG_SNVs_tbi) from ch_cadd_wg_snvs_tbi

    output:
        set variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf") into vepVCFmerge
        file("${reducedVCF}_VEP.summary.html") into vepReportMerge

    when: 'merge' in tools

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
    dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
    cadd = (params.cadd_cache && params.cadd_wg_snvs && params.cadd_indels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${params.species} \
        ${cadd} \
        ${genesplicer} \
        --cache \
        --cache_version ${cache_version} \
        --dir_cache ${dir_cache} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${reducedVCF}_VEP.summary.html \
        --total_length \
        --vcf

    rm -rf ${reducedVCF}
    """
}

vepReportMerge = vepReportMerge.dump(tag:'VEP')

vcfCompressVCFvep = vepVCF.mix(vepVCFmerge)

// STEP COMPRESS AND INDEX VCF.2 - VEP

process CompressVCFvep {
    
    tag "${idSample} - ${vcf}"
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir "${params.outdir}/Annotation/${idSample}/VEP", mode: params.publish_dir_mode

    input:
        set variantCaller, idSample, file(vcf) from vcfCompressVCFvep

    output:
        set variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into compressVCFOutVEP

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}

compressVCFOutVEP = compressVCFOutVEP.dump(tag:'VCF')

/*
================================================================================
                                     MultiQC
================================================================================
*/

// STEP MULTIQC

process MultiQC {

    publishDir "${params.outdir}/Reports/MultiQC", mode: params.publish_dir_mode
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    input:
        file (multiqcConfig) from ch_multiqc_config
        file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
        file (versions) from ch_software_versions_yaml.collect()
        file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
        file ('bamQC/*') from bamQCReport.collect().ifEmpty([])
        file ('BCFToolsStats/*') from bcftoolsReport.collect().ifEmpty([])
        file ('FastQC/*') from fastQCReport.collect().ifEmpty([])
        file ('TrimmedFastQC/*') from trimGaloreReport.collect().ifEmpty([])
        file ('MarkDuplicates/*') from duplicates_marked_report.collect().ifEmpty([])
        file ('DuplicatesMarked/*.recal.table') from baseRecalibratorReport.collect().ifEmpty([])
        file ('SamToolsStats/*') from samtoolsStatsReport.collect().ifEmpty([])
        file ('snpEff/*') from snpeffReport.collect().ifEmpty([])
        file ('VCFTools/*') from vcftoolsReport.collect().ifEmpty([])

    output:
        file "*multiqc_report.html" into ch_multiqc_report
        file "*_data"
        file "multiqc_plots"

    when: !('multiqc' in skipQC)

    script:
    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f ${rtitle} ${rfilename} ${custom_config_file} .
    """
}

ch_multiqc_report.dump(tag:'MultiQC')

// Output Description HTML
process Output_documentation {

    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode
    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    input:
        file output_docs from ch_output_docs
        file images from ch_output_docs_images

    output:
        file "results_description.html"

    when: !('documentation' in skipQC)

    script:
    """
    python ${workflow.launchDir/params.config_files}/bin/markdown_to_html.py $output_docs -o results_description.html
    """
}

process copy_workfiles {

    container pipeline_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-large'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
        file(multiqc_reportz) from ch_multiqc_report


    script:
    """
    cp -r ${workflow.workDir} ${workflow.launchDir}/out
    cp -r ${workflow.launchDir}/.nextflow.log ${workflow.launchDir}/out
    """
}

// Completion e-mail notification
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/sarek] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/sarek] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/sarek] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/sarek] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/sarek] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/sarek] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset  = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/sarek]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/sarek]${c_red} Pipeline completed with errors${c_reset}-"
    }
}

//workflow.onError {//
    // Print unexpected parameters - easiest is to just rerun validation
    //NfcoreSchema.validateParameters(params, json_schema, log)
//}

def checkHostname() {
    def c_reset       = params.monochrome_logs ? '' : "\033[0m"
    def c_white       = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red         = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}

/*
================================================================================
                                 sarek functions
================================================================================
*/

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

// Define list of available tools to annotate
def defineAnnoList() {
    return [
        'haplotypecaller',
        'manta',
        'mutect2',
        'strelka',
        'tiddit'
    ]
}

// Define list of skipable QC tools
def defineSkipQClist() {
    return [
        'bamqc',
        'baserecalibrator',
        'bcftools',
        'documentation',
        'fastqc',
        'markduplicates',
        'multiqc',
        'samtools',
        'sentieon',
        'vcftools',
        'versions'
    ]
}

// Define list of available step
def defineStepList() {
    return [
        'annotate',
        'controlfreec',
        'mapping',
        'preparerecalibration',
        'recalibrate',
        'variantcalling'
    ]
}

// Define list of available tools
def defineToolList() {
    return [
        'ascat',
        'cnvkit',
        'controlfreec',
        'dnascope',
        'dnaseq',
        'freebayes',
        'haplotypecaller',
        'manta',
        'merge',
        'mpileup',
        'mutect2',
        'snpeff',
        'strelka',
        'tiddit',
        'tnscope',
        'vep',
        'msisensor'
    ]
}

// Channeling the TSV file containing BAM.
// Format is: "subject gender status sample bam bai"
def extractBam(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 6)
            def idPatient = row[0]
            def gender    = row[1]
            def status    = returnStatus(row[2].toInteger())
            def idSample  = row[3]
            def bamFile   = returnFile(row[4])
            def baiFile   = returnFile(row[5])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

            return [idPatient, gender, status, idSample, bamFile, baiFile]
        }
}

// Create a channel of germline FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and _R2_ in their names.
def extractFastqFromDir(pattern) {
    def fastq = Channel.create()
    // a temporary channel does all the work
    Channel
        .fromPath(pattern, type: 'dir')
        .ifEmpty { error "No directories found matching pattern '${pattern}'" }
        .subscribe onNext: { sampleDir ->
            // the last name of the sampleDir is assumed to be a unique sample id
            sampleId = sampleDir.getFileName().toString()

            for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
                assert path1.getName().contains('_R1_')
                path2 = file(path1.toString().replace('_R1_', '_R2_'))
                if (!path2.exists()) error "Path '${path2}' not found"
                (flowcell, lane) = flowcellLaneFromFastq(path1)
                String random = org.apache.commons.lang.RandomStringUtils.random(8, true, true) // random string to avoid duplicate names
                patient = sampleId
                gender = 'ZZ'  // unused
                status = 0  // normal (not tumor)
                rgId = "${flowcell}.${sampleId}.${lane}.${random}"
                result = [patient, gender, status, sampleId, rgId, path1, path2]
                fastq.bind(result)
            }
    }, onComplete: { fastq.close() }
    fastq
}

// Extract gender and status from Channel
def extractInfos(channel) {
    def genderMap = [:]
    def statusMap = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        genderMap[idPatient] = gender
        statusMap[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    [genderMap, statusMap, channel]
}

// Channeling the TSV file containing FASTQ or BAM
// Format is: "subject gender status sample lane fastq1 fastq2"
// or: "subject gender status sample lane bam"
def extractFastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def idRun      = row[4]
            def file1      = returnFile(row[5])
            def file2      = "null"
            if (hasExtension(file1, "fastq.gz") || hasExtension(file1, "fq.gz") || hasExtension(file1, "fastq") || hasExtension(file1, "fq")) {
                checkNumberOfItem(row, 7)
                file2 = returnFile(row[6])
                if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")  && !hasExtension(file2, "fastq") && !hasExtension(file2, "fq")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
                if (hasExtension(file1, "fastq") || hasExtension(file1, "fq") || hasExtension(file2, "fastq") || hasExtension(file2, "fq")) {
                    exit 1, "We do recommend to use gziped fastq file to help you reduce your data footprint."
                }
            }
            else if (hasExtension(file1, "bam")) checkNumberOfItem(row, 6)
            else "No recognisable extention for input file: ${file1}"

            [idPatient, gender, status, idSample, idRun, file1, file2]
        }
}

// Channeling the TSV file containing mpileup
// Format is: "subject gender status sample pileup"
def extractPileup(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 5)
            def idPatient = row[0]
            def gender    = row[1]
            def status    = returnStatus(row[2].toInteger())
            def idSample  = row[3]
            def mpileup   = returnFile(row[4])

            if (!hasExtension(mpileup, "pileup")) exit 1, "File: ${mpileup} has the wrong extension. See --help for more information"

            return [idPatient, gender, status, idSample, mpileup]
        }
}

// Channeling the TSV file containing Recalibration Tables.
// Format is: "subject gender status sample bam bai recalTable"
def extractRecal(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 7)
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def bamFile    = returnFile(row[4])
            def baiFile    = returnFile(row[5])
            def recalTable = returnFile(row[6])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
            if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"

            [idPatient, gender, status, idSample, bamFile, baiFile, recalTable]
    }
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String fcid
    int lane
    if (fields.size() == 7) {
        // CASAVA 1.8+ format
        fcid = fields[2]
        lane = fields[3].toInteger()
    } else if (fields.size() == 5) {
        fcid = fields[0]
        lane = fields[1].toInteger()
    }
    [fcid, lane]
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Remove .ann .gz and .vcf extension from a VCF file
def reduceVCF(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}
