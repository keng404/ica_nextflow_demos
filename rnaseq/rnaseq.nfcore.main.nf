#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rnaseq
========================================================================================
 nf-core/rnaseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rnaseq
----------------------------------------------------------------------------------------
*/
// Global default params, used in configs

// Pipeline Options
// Workflow flags
params.genome = 'GRCh37'
params.reads = "data/*{1,2}.fastq.gz"
params.singleEnd = false
params.publish_dir_mode = 'copy'
params.config_files = ''

// References
params.salmon_index = false
params.transcript_fasta = false
params.splicesites = false
params.saveReference = false
params.gencode = false
params.compressedReference = false

params.gff =  false
params.hisat2_index = false

// Strandedness
params.forwardStranded = false
params.reverseStranded = false
params.unStranded = false

// Trimming
params.skipTrimming = false
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0
params.trim_nextseq = 0
params.pico = false
params.saveTrimmed = false

// Ribosomal RNA removal
params.removeRiboRNA = false
params.save_nonrRNA_reads = false

// Alignment
params.aligner = 'star'
params.pseudo_aligner = false
params.stringTieIgnoreGTF = false
params.seq_center = false
params.saveAlignedIntermediates = false
params.skipAlignment = false
params.saveUnaligned = false

// Read Counting
params.fc_extra_attributes = 'gene_name'
params.fc_group_features = 'gene_id'
params.fc_count_type = 'exon'
params.fc_group_features_type = 'gene_biotype'
params.sampleLevel = false
params.skipBiotypeQC = false

// QC
params.skipQC = false
params.skipFastQC = false
params.skipPreseq = false
params.skipDupRadar = false
params.skipQualimap = false
params.skipRseQC = false
params.skipEdgeR = false
params.skipMultiQC = false

// Defaults
params.project = false
params.markdup_java_options = '"-Xms4000m -Xmx7g"' //Established values for markDuplicate memory consumption, see issue PR #689 (in Sarek) for details
params.hisat_build_memory = 200 // Required amount of memory in GB to build HISAT2 index with splice sites
params.readPaths = null
params.star_memory = false // Cluster specific param required for hebbe
params.rRNA_database_manifest = "${params.config_files}/assets/rrna-db-defaults.txt"

// Boilerplate options
params.clusterOptions = false
params.outdir = "${workflow.launchDir}/out"
params.name = false
params.multiqc_config = "${params.config_files}/assets/multiqc_config.yaml"
params.email = false
params.email_on_fail = false
params.max_multiqc_email_size = 25.MB
params.plaintext_email = false
params.monochrome_logs = false
params.help = false
params.tracedir = "${params.outdir}/pipeline_info"
params.awsqueue = false
params.awsregion = 'eu-west-1'
params.igenomesIgnore = false
params.custom_config_version = 'master'
params.custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
params.hostnames = false
params.config_profile_description = false
params.config_profile_contact = false
params.config_profile_url = false

///// Function to ensure that resource requirements don't go beyond
// a maximum limit
def check_max(obj, type) {
    if (type == 'memory') {
        try {
            if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                return params.max_memory as nextflow.util.MemoryUnit
            else
                return obj
        } catch (all) {
            println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
            return obj
        }
    } else if (type == 'time') {
        try {
            if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                return params.max_time as nextflow.util.Duration
            else
                return obj
        } catch (all) {
            println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
            return obj
        }
    } else if (type == 'cpus') {
        try {
            return Math.min( obj, params.max_cpus as int )
        } catch (all) {
            println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
            return obj
        }
    }
}
// Container slug. Stable releases should specify release tag!
// Developmental code should specify :dev
// process.container = 'nfcore/rnaseq:1.4.2'

// Defaults only, expecting to be overwritten
params.max_memory = 1024.GB
params.max_cpus = 16
params.max_time = 240.h
params.igenomes_base = 's3://ngi-igenomes/igenomes/'


// Avoid this error:
// WARNING: Your kernel does not support swap limit capabilities or the cgroup is not mounted. Memory limited without swap.
// Testing this in nf-core after discussion here https://github.com/nf-core/tools/pull/351, once this is established and works well, nextflow might implement this behavior as new default.
// docker.runOptions = '-u \$(id -u):\$(id -g)'

// Load igenomes.config if required

/*
 * -------------------------------------------------
 *  Nextflow config file for iGenomes paths
 * -------------------------------------------------
 * Defines reference genomes, using iGenome paths
 * Can be used by any config that customises the base
 * path using $params.igenomes_base / --igenomes_base
 */
// illumina iGenomes reference file paths
params.genomes = [:]
// GRCh37
params.genomes['GRCh37'] = [:]
params.genomes['GRCh37'].fasta       = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['GRCh37'].bwa         = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"
params.genomes['GRCh37'].bowtie2     = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/"
params.genomes['GRCh37'].star        = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/"
params.genomes['GRCh37'].bismark     = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/BismarkIndex/"
params.genomes['GRCh37'].gtf         = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
params.genomes['GRCh37'].bed12       = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed"
params.genomes['GRCh37'].mito_name   = "MT"
params.genomes['GRCh37'].macs_gsize  = "2.7e9"
params.genomes['GRCh37'].blacklist   = "${params.config_files}/assets/blacklists/GRCh37-blacklist.bed"

//GRCh38
params.genomes['GRCh38'] = [:]
params.genomes['GRCh38'].fasta       = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['GRCh38'].bwa         = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
params.genomes['GRCh38'].bowtie2     = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/"
params.genomes['GRCh38'].star        = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/"
params.genomes['GRCh38'].bismark     = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/BismarkIndex/"
params.genomes['GRCh38'].gtf         = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
params.genomes['GRCh38'].bed12       = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.bed"
params.genomes['GRCh38'].mito_name   = "chrM"
params.genomes['GRCh38'].macs_gsize  = "2.7e9"
params.genomes['GRCh38'].blacklist   = "${params.config_files}/assets/blacklists/hg38-blacklist.bed"

// GRCm38
params.genomes['GRCm38'] = [:]
params.genomes['GRCm38'].fasta       = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['GRCm38'].bwa         = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa"
params.genomes['GRCm38'].bowtie2     = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/"
params.genomes['GRCm38'].star        = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/STARIndex/"
params.genomes['GRCm38'].bismark     = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/BismarkIndex/"
params.genomes['GRCm38'].gtf         = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"
params.genomes['GRCm38'].bed12       = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.bed"
params.genomes['GRCm38'].mito_name   = "MT"
params.genomes['GRCm38'].macs_gsize  = "1.87e9"
params.genomes['GRCm38'].blacklist   = "${params.config_files}/assets/blacklists/GRCm38-blacklist.bed"

//'TAIR10' 
params.genomes['TAIR10'] = [:]
params.genomes['TAIR10'].fasta       = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['TAIR10'].bwa         = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BWAIndex/genome.fa"
params.genomes['TAIR10'].bowtie2     = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/"
params.genomes['TAIR10'].star        = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/STARIndex/"
params.genomes['TAIR10'].bismark     = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BismarkIndex/"
params.genomes['TAIR10'].gtf         = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.gtf"
params.genomes['TAIR10'].bed12       = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.bed"
params.genomes['TAIR10'].mito_name   = "Mt"

// 'EB2' 
params.genomes['EB2'] = [:]
params.genomes['EB2'].fasta       = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['EB2'].bwa         = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/BWAIndex/genome.fa"
params.genomes['EB2'].bowtie2     = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/Bowtie2Index/"
params.genomes['EB2'].star        = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/STARIndex/"
params.genomes['EB2'].bismark     = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/BismarkIndex/"
params.genomes['EB2'].gtf         = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Annotation/Genes/genes.gtf"
params.genomes['EB2'].bed12       = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Annotation/Genes/genes.bed"

//'UMD3.1' 
params.genomes['UMD3.1'] = [:]
params.genomes['UMD3.1'].fasta       = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['UMD3.1'].bwa         = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/BWAIndex/genome.fa"
params.genomes['UMD3.1'].bowtie2     = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/Bowtie2Index/"
params.genomes['UMD3.1'].star        = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/STARIndex/"
params.genomes['UMD3.1'].bismark     = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/BismarkIndex/"
params.genomes['UMD3.1'].gtf         = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Annotation/Genes/genes.gtf"
params.genomes['UMD3.1'].bed12       = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Annotation/Genes/genes.bed"
params.genomes['UMD3.1'].mito_name   = "MT"

// 'WBcel235' 
params.genomes['WBcel235'] = [:]
params.genomes['WBcel235'].fasta       = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['WBcel235'].bwa         = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BWAIndex/genome.fa"
params.genomes['WBcel235'].bowtie2     = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/Bowtie2Index/"
params.genomes['WBcel235'].star        = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/STARIndex/"
params.genomes['WBcel235'].bismark     = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BismarkIndex/"
params.genomes['WBcel235'].gtf         = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
params.genomes['WBcel235'].bed12       = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.bed"
params.genomes['WBcel235'].mito_name   = "MtDNA"
params.genomes['WBcel235'].macs_gsize  = "9e7"

//'CanFam3.1'
params.genomes['CanFam3.1'] = [:]
params.genomes['CanFam3.1'].fasta       = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['CanFam3.1'].bwa         = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/BWAIndex/genome.fa"
params.genomes['CanFam3.1'].bowtie2     = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/Bowtie2Index/"
params.genomes['CanFam3.1'].star        = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/STARIndex/"
params.genomes['CanFam3.1'].bismark     = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/BismarkIndex/"
params.genomes['CanFam3.1'].gtf         = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Annotation/Genes/genes.gtf"
params.genomes['CanFam3.1'].bed12       = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Annotation/Genes/genes.bed"
 params.genomes['CanFam3.1'].mito_name   = "MT"

//'GRCz10' 
params.genomes['GRCz10'] = [:]
params.genomes['GRCz10'].fasta       = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['GRCz10'].bwa         = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/BWAIndex/genome.fa"
params.genomes['GRCz10'].bowtie2     = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/Bowtie2Index/"
params.genomes['GRCz10'].star        = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/STARIndex/"
params.genomes['GRCz10'].bismark     = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/BismarkIndex/"
params.genomes['GRCz10'].gtf         = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Annotation/Genes/genes.gtf"
params.genomes['GRCz10'].bed12       = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Annotation/Genes/genes.bed"
params.genomes['GRCz10'].mito_name   = "MT"

// 'BDGP6' 
params.genomes['BDGP6'] = [:]
params.genomes['BDGP6'].fasta       = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['BDGP6'].bwa         = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/BWAIndex/genome.fa"
params.genomes['BDGP6'].bowtie2     = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/Bowtie2Index/"
params.genomes['BDGP6'].star        = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/STARIndex/"
params.genomes['BDGP6'].bismark     = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/BismarkIndex/"
params.genomes['BDGP6'].gtf         = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Annotation/Genes/genes.gtf"
params.genomes['BDGP6'].bed12       = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Annotation/Genes/genes.bed"
params.genomes['BDGP6'].mito_name   = "M"
params.genomes['BDGP6'].macs_gsize  = "1.2e8"

// 'EquCab2' 
params.genomes['EquCab2'] = [:]
params.genomes['EquCab2'].fasta       = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['EquCab2'].bwa         = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/BWAIndex/genome.fa"
params.genomes['EquCab2'].bowtie2     = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/Bowtie2Index/"
params.genomes['EquCab2'].star        = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/STARIndex/"
params.genomes['EquCab2'].bismark     = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/BismarkIndex/"
params.genomes['EquCab2'].gtf         = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Annotation/Genes/genes.gtf"
params.genomes['EquCab2'].bed12       = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Annotation/Genes/genes.bed"
params.genomes['EquCab2'].mito_name   = "MT"

//'EB1' 
params.genomes['EB1'] = [:]
params.genomes['EB1'].fasta       = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['EB1'].bwa         = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/BWAIndex/genome.fa"
params.genomes['EB1'].bowtie2     = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/Bowtie2Index/"
params.genomes['EB1'].star        = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/STARIndex/"
params.genomes['EB1'].bismark     = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/BismarkIndex/"
params.genomes['EB1'].gtf         = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/Genes/genes.gtf"
params.genomes['EB1'].bed12       = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/Genes/genes.bed"

// 'Galgal4'
params.genomes['Galgal4'] = [:]
params.genomes['Galgal4'].fasta       = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Galgal4'].bwa         = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/BWAIndex/genome.fa"
params.genomes['Galgal4'].bowtie2     = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/Bowtie2Index/"
params.genomes['Galgal4'].star        = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/STARIndex/"
params.genomes['Galgal4'].bismark     = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/BismarkIndex/"
params.genomes['Galgal4'].gtf         = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Annotation/Genes/genes.gtf"
params.genomes['Galgal4'].bed12       = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Annotation/Genes/genes.bed"
params.genomes['Galgal4'].mito_name   = "MT"

// 'Gm01' 
params.genomes['Gm01'] = [:]
params.genomes['Gm01'].fasta       = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Gm01'].bwa         = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/BWAIndex/genome.fa"
params.genomes['Gm01'].bowtie2     = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/Bowtie2Index/"
params.genomes['Gm01'].star        = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/STARIndex/"
params.genomes['Gm01'].bismark     = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/BismarkIndex/"
params.genomes['Gm01'].gtf         = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Annotation/Genes/genes.gtf"
params.genomes['Gm01'].bed12       = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Annotation/Genes/genes.bed"

// 'Mmul_1' 
params.genomes['Mmul_1'] = [:]
params.genomes['Mmul_1'].fasta       = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Mmul_1'].bwa         = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/BWAIndex/genome.fa"
params.genomes['Mmul_1'].bowtie2     = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/Bowtie2Index/"
params.genomes['Mmul_1'].star        = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/STARIndex/"
params.genomes['Mmul_1'].bismark     = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/BismarkIndex/"
params.genomes['Mmul_1'].gtf         = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Annotation/Genes/genes.gtf"
params.genomes['Mmul_1'].bed12       = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Annotation/Genes/genes.bed"
params.genomes['Mmul_1'].mito_name   = "MT"

// 'IRGSP-1.0' 
params.genomes['IRGSP-1.0'] = {}
params.genomes['IRGSP-1.0'].fasta       = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['IRGSP-1.0'].bwa         = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/BWAIndex/genome.fa"
params.genomes['IRGSP-1.0'].bowtie2     = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/Bowtie2Index/"
params.genomes['IRGSP-1.0'].star        = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/STARIndex/"
params.genomes['IRGSP-1.0'].bismark     = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/BismarkIndex/"
params.genomes['IRGSP-1.0'].gtf         = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Annotation/Genes/genes.gtf"
params.genomes['IRGSP-1.0'].bed12       = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Annotation/Genes/genes.bed"
params.genomes['IRGSP-1.0'].mito_name   = "Mt"

// 'CHIMP2.1.4' 
params.genomes['CHIMP2.1.4'] = [:]
params.genomes['CHIMP2.1.4'].fasta       = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['CHIMP2.1.4'].bwa         = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/BWAIndex/genome.fa"
params.genomes['CHIMP2.1.4'].bowtie2     = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/Bowtie2Index/"
params.genomes['CHIMP2.1.4'].star        = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/STARIndex/"
params.genomes['CHIMP2.1.4'].ismark     = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/BismarkIndex/"
params.genomes['CHIMP2.1.4'].gtf         = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/Genes/genes.gtf"
params.genomes['CHIMP2.1.4'].bed12       = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/Genes/genes.bed"
params.genomes['CHIMP2.1.4'].mito_name   = "MT"

//' Rnor_6.0' 
params.genomes['Rnor_6.0'] = [:]
params.genomes['Rnor_6.0'].fasta       = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Rnor_6.0'].bwa         = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/BWAIndex/genome.fa"
params.genomes['Rnor_6.0'].bowtie2     = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/Bowtie2Index/"
params.genomes['Rnor_6.0'].star        = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/STARIndex/"
params.genomes['Rnor_6.0'].bismark     = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/BismarkIndex/"
params.genomes['Rnor_6.0'].gtf         = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.gtf"
params.genomes['Rnor_6.0'].bed12       = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.bed"
params.genomes['Rnor_6.0'].mito_name   = "MT"

// 'R64-1-1' 
params.genomes['R64-1-1'] = [:]
params.genomes['R64-1-1'].fasta       = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['R64-1-1'].bwa         = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/BWAIndex/genome.fa"
params.genomes['R64-1-1'].bowtie2     = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/Bowtie2Index/"
params.genomes['R64-1-1'].star        = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/STARIndex/"
params.genomes['R64-1-1'].bismark     = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/BismarkIndex/"
params.genomes['R64-1-1'].gtf         = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Annotation/Genes/genes.gtf"
params.genomes['R64-1-1'].bed12       = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Annotation/Genes/genes.bed"
params.genomes['R64-1-1'].mito_name   = "MT"
params.genomes['R64-1-1'].macs_gsize  = "1.2e7"

// 'EF2' 
params.genomes['EF2'] = [:]
params.genomes['EF2'].fasta       = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['EF2'].bwa         = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/BWAIndex/genome.fa"
params.genomes['EF2'].bowtie2     = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/Bowtie2Index/"
params.genomes['EF2'].star        = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/STARIndex/"
params.genomes['EF2'].bismark     = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/BismarkIndex/"
params.genomes['EF2'].gtf         = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Genes/genes.gtf"
params.genomes['EF2'].bed12       = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Genes/genes.bed"
params.genomes['EF2'].mito_name   = "MT"
params.genomes['EF2'].macs_gsize  = "1.21e7"

// 'Sbi1' 
params.genomes['Sbi1'] = [:]
params.genomes['Sbi1'].fasta       = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Sbi1'].bwa         = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/BWAIndex/genome.fa"
params.genomes['Sbi1'].bowtie2     = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/Bowtie2Index/"
params.genomes['Sbi1'].star        = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/STARIndex/"
params.genomes['Sbi1'].bismark     = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/BismarkIndex/"
params.genomes['Sbi1'].gtf         = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Annotation/Genes/genes.gtf"
params.genomes['Sbi1'].bed12       = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Annotation/Genes/genes.bed"

// 'Sscrofa10.2' 
params.genomes['Sscrofa10.2'] = [:]
params.genomes['Sscrofa10.2'].fasta       = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Sscrofa10.2'].bwa         = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/BWAIndex/genome.fa"
params.genomes['Sscrofa10.2'].bowtie2     = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/Bowtie2Index/"
params.genomes['Sscrofa10.2'].star        = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/STARIndex/"
params.genomes['Sscrofa10.2'].bismark     = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/BismarkIndex/"
params.genomes['Sscrofa10.2'].gtf         = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/Genes/genes.gtf"
params.genomes['Sscrofa10.2'].bed12       = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/Genes/genes.bed"
params.genomes['Sscrofa10.2'].mito_name   = "MT"

// 'AGPv3' 
params.genomes['AGPv3'] = [:]
params.genomes['AGPv3'].fasta       = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['AGPv3'].bwa         = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/BWAIndex/genome.fa"
params.genomes['AGPv3'].bowtie2     = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/Bowtie2Index/"
params.genomes['AGPv3'].star        = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/STARIndex/"
params.genomes['AGPv3'].bismark     = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/BismarkIndex/"
params.genomes['AGPv3'].gtf         = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Annotation/Genes/genes.gtf"
params.genomes['AGPv3'].bed12       = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Annotation/Genes/genes.bed"
params.genomes['AGPv3'].mito_name   = "Mt"

// 'hg38' 
params.genomes['hg38'] = [:]
params.genomes['hg38'].fasta       = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['hg38'].bwa         = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
params.genomes['hg38'].bowtie2     = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/"
params.genomes['hg38'].star        = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/STARIndex/"
params.genomes['hg38'].bismark     = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/BismarkIndex/"
params.genomes['hg38'].gtf         = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
params.genomes['hg38'].bed12       = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.bed"
params.genomes['hg38'].mito_name   = "chrM"
params.genomes['hg38'].macs_gsize  = "2.7e9"
params.genomes['hg38'].blacklist   = "${params.config_files}/assets/blacklists/hg38-blacklist.bed"

// 'hg19' 
params.genomes['hg19'] = [:]
params.genomes['hg19'].fasta       = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['hg19'].bwa         = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
params.genomes['hg19'].bowtie2     = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/"
params.genomes['hg19'].star        = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/"
params.genomes['hg19'].bismark     = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/BismarkIndex/"
params.genomes['hg19'].gtf         = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
params.genomes['hg19'].bed12       = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed"
params.genomes['hg19'].mito_name   = "chrM"
params.genomes['hg19'].macs_gsize  = "2.7e9"
params.genomes['hg19'].blacklist   = "${params.config_files}/assets/blacklists/hg19-blacklist.bed"

// 'mm10' 
params.genomes['mm10'] = [:]
params.genomes['mm10'].fasta       = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['mm10'].bwa         = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa"
params.genomes['mm10'].bowtie2     = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/"
params.genomes['mm10'].star        = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/STARIndex/"
params.genomes['mm10'].bismark     = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/BismarkIndex/"
params.genomes['mm10'].gtf         = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
params.genomes['mm10'].bed12       = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.bed"
params.genomes['mm10'].mito_name   = "chrM"
params.genomes['mm10'].macs_gsize  = "1.87e9"
params.genomes['mm10'].blacklist   = "${params.config_files}/assets/blacklists/mm10-blacklist.bed"

// 'bosTau8' 
params.genomes['bosTau8'] = [:]
params.genomes['bosTau8'].fasta       = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['bosTau8'].bwa         = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/BWAIndex/genome.fa"
params.genomes['bosTau8'].bowtie2     = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/Bowtie2Index/"
params.genomes['bosTau8'].star        = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/STARIndex/"
params.genomes['bosTau8'].bismark     = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/BismarkIndex/"
params.genomes['bosTau8'].gtf         = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.gtf"
params.genomes['bosTau8'].bed12       = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.bed"
params.genomes['bosTau8'].mito_name   = "chrM"

// 'ce10' 
params.genomes['ce10'] = [:]
params.genomes['ce10'].fasta       = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['ce10'].bwa         = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/BWAIndex/genome.fa"
params.genomes['ce10'].bowtie2     = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/Bowtie2Index/"
params.genomes['ce10'].star        = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/STARIndex/"
params.genomes['ce10'].bismark     = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/BismarkIndex/"
params.genomes['ce10'].gtf         = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Annotation/Genes/genes.gtf"
params.genomes['ce10'].bed12       = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Annotation/Genes/genes.bed"
params.genomes['ce10'].mito_name   = "chrM"
params.genomes['ce10'].macs_gsize  = "9e7"

//  'canFam3' 
params.genomes['canFam3'] = [:]
params.genomes['canFam3'].fasta       = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['canFam3'].bwa         = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/BWAIndex/genome.fa"
params.genomes['canFam3'].bowtie2     = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/Bowtie2Index/"
params.genomes['canFam3'].star        = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/STARIndex/"
params.genomes['canFam3'].bismark     = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/BismarkIndex/"
params.genomes['canFam3'].gtf         = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Annotation/Genes/genes.gtf"
params.genomes['canFam3'].bed12       = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Annotation/Genes/genes.bed"
params.genomes['canFam3'].mito_name   = "chrM"

// 'danRer10' 
params.genomes['danRer10'] = [:]
params.genomes['danRer10'].fasta       = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['danRer10'].bwa         = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/BWAIndex/genome.fa"
params.genomes['danRer10'].bowtie2     = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/Bowtie2Index/"
params.genomes['danRer10'].star        = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/STARIndex/"
params.genomes['danRer10'].bismark     = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/BismarkIndex/"
params.genomes['danRer10'].gtf         = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Annotation/Genes/genes.gtf"
params.genomes['danRer10'].bed12       = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Annotation/Genes/genes.bed"
params.genomes['danRer10'].mito_name   = "chrM"

// 'dm6' 
params.genomes['dm6'] = [:]
params.genomes['dm6'].fasta       = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['dm6'].bwa         = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex/genome.fa"
params.genomes['dm6'].bowtie2     = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/"
params.genomes['dm6'].star        = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/STARIndex/"
params.genomes['dm6'].bismark     = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/BismarkIndex/"
params.genomes['dm6'].gtf         = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf"
params.genomes['dm6'].bed12       = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.bed"
params.genomes['dm6'].mito_name   = "chrM"
params.genomes['dm6'].macs_gsize  = "1.2e8"

// 'equCab2' 
params.genomes['equCab2'] = [:]
params.genomes['equCab2'].fasta       = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['equCab2'].bwa         = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/BWAIndex/genome.fa"
params.genomes['equCab2'].bowtie2     = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/Bowtie2Index/"
params.genomes['equCab2'].star        = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/STARIndex/"
params.genomes['equCab2'].bismark     = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/BismarkIndex/"
params.genomes['equCab2'].gtf         = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Annotation/Genes/genes.gtf"
params.genomes['equCab2'].bed12       = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Annotation/Genes/genes.bed"
params.genomes['equCab2'].mito_name   = "chrM"

// 'galGal4' 
params.genomes['galGal4'] = [:]
params.genomes['galGal4'].fasta       = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['galGal4'].bwa         = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/BWAIndex/genome.fa"
params.genomes['galGal4'].bowtie2     = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/Bowtie2Index/"
params.genomes['galGal4'].star        = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/STARIndex/"
params.genomes['galGal4'].bismark     = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/BismarkIndex/"
params.genomes['galGal4'].gtf         = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Annotation/Genes/genes.gtf"
params.genomes['galGal4'].bed12       = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Annotation/Genes/genes.bed"
params.genomes['galGal4'].mito_name   = "chrM"

// 'panTro4' 
params.genomes['panTro4'] = [:]
params.genomes['panTro4'].fasta       = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['panTro4'].bwa         = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/BWAIndex/genome.fa"
params.genomes['panTro4'].bowtie2     = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/Bowtie2Index/"
params.genomes['panTro4'].star        = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/STARIndex/"
params.genomes['panTro4'].bismark     = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/BismarkIndex/"
params.genomes['panTro4'].gtf         = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Annotation/Genes/genes.gtf"
params.genomes['panTro4'].bed12       = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Annotation/Genes/genes.bed"
params.genomes['panTro4'].mito_name   = "chrM"

// 'rn6' 
params.genomes['rn6'] = [:]
params.genomes['rn6'].fasta       = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['rn6'].bwa         = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/BWAIndex/genome.fa"
params.genomes['rn6'].bowtie2     = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/Bowtie2Index/"
params.genomes['rn6'].star        = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/STARIndex/"
params.genomes['rn6'].bismark     = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/BismarkIndex/"
params.genomes['rn6'].gtf         = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf"
params.genomes['rn6'].bed12       = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.bed"
params.genomes['rn6'].mito_name   = "chrM"

// 'sacCer3' 
params.genomes['sacCer3'] = [:]
params.genomes['sacCer3'].fasta       = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['sacCer3'].bwa         = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BWAIndex/genome.fa"
params.genomes['sacCer3'].bowtie2     = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/"
params.genomes['sacCer3'].star        = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/STARIndex/"
params.genomes['sacCer3'].bismark     = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BismarkIndex/"
params.genomes['sacCer3'].mito_name   = "chrM"
params.genomes['sacCer3'].macs_gsize  = "1.2e7"

// 'susScr3' 
params.genomes['susScr3'] = [:]
params.genomes['susScr3'].fasta       = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['susScr3'].bwa         = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/BWAIndex/genome.fa"
params.genomes['susScr3'].bowtie2     = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/Bowtie2Index/"
params.genomes['susScr3'].star        = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/STARIndex/"
params.genomes['susScr3'].bismark     = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/BismarkIndex/"
params.genomes['susScr3'].gtf         = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Annotation/Genes/genes.gtf"
params.genomes['susScr3'].bed12       = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Annotation/Genes/genes.bed"
params.genomes['susScr3'].mito_name   = "chrM"

// Capture exit codes from upstream processes when piping
// process.shell = ['/bin/bash', '-euo', 'pipefail']

timeline = [:]
timeline.enabled = true
timeline.file = "${params.tracedir}/execution_timeline.html"

report = [:]
report.enabled = true
report.file = "${params.tracedir}/execution_report.html"

trace = [:]
trace.enabled = true
trace.file = "${params.tracedir}/execution_trace.txt"

dag = [:]
dag.enabled = true
dag.file = "${params.tracedir}/pipeline_dag.svg"



docker_container = 'nfcore/rnaseq:1.4.2'
instance_size = 'himem-medium'
input_filter_string = params.input_filter_string

input = ''
input = Channel.fromPath(params.reads)
input1 = file(params.reads)
input_basename = input1.getBaseName()
input = "${workflow.launchDir}/${input_basename}/${input_filter_string}"

def helpMessage() {
    log.info nfcoreHeader()
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nf-core/rnaseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile docker
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
    Generic:
      --singleEnd                   Specifies that the input is single-end reads
    References:                     If not specified in the configuration file or you wish to overwrite any of the references.
      --genome                      Name of iGenomes reference
      --star_index                  Path to STAR index
      --hisat2_index                Path to HiSAT2 index
      --salmon_index                Path to Salmon index
      --fasta                       Path to genome fasta file
      --transcript_fasta            Path to transcript fasta file
      --splicesites                 Path to splice sites file for building HiSat2 index
      --gtf                         Path to GTF file
      --gff                         Path to GFF3 file
      --bed12                       Path to bed12 file
      --saveReference               Save the generated reference files to the results directory
      --gencode                     Use fc_group_features_type = 'gene_type' and pass '--gencode' flag to Salmon
    Strandedness:
      --forwardStranded             The library is forward stranded
      --reverseStranded             The library is reverse stranded
      --unStranded                  The default behaviour
    Trimming:
      --skipTrimming                Skip Trim Galore step
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]   Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed
      --trim_nextseq [int]          Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails
      --pico                        Sets trimming and standedness settings for the SMARTer Stranded Total RNA-Seq Kit - Pico Input kit. Equivalent to: --forwardStranded --clip_r1 3 --three_prime_clip_r2 3
      --saveTrimmed                 Save trimmed FastQ file intermediates
    Ribosomal RNA removal:
      --removeRiboRNA               Removes ribosomal RNA using SortMeRNA
      --save_nonrRNA_reads          Save FastQ file intermediates after removing rRNA
      --rRNA_database_manifest      Path to file that contains file paths for rRNA databases, optional
    Alignment:
      --aligner                     Specifies the aligner to use (available are: 'hisat2', 'star')
      --pseudo_aligner              Specifies the pseudo aligner to use (available are: 'salmon'). Runs in addition to `--aligner`
      --stringTieIgnoreGTF          Perform reference-guided de novo assembly of transcripts using StringTie i.e. dont restrict to those in GTF file
      --seq_center                  Add sequencing center in @RG line of output BAM header
      --saveAlignedIntermediates    Save the BAM files from the aligment step - not done by default
      --saveUnaligned               Save unaligned reads from either STAR, HISAT2 or Salmon to extra output files.
      --skipAlignment               Skip alignment altogether (usually in favor of pseudoalignment)
    Read Counting:
      --fc_extra_attributes         Define which extra parameters should also be included in featureCounts (default: 'gene_name')
      --fc_group_features           Define the attribute type used to group features. (default: 'gene_id')
      --fc_count_type               Define the type used to assign reads. (default: 'exon')
      --fc_group_features_type      Define the type attribute used to group features based on the group attribute (default: 'gene_biotype')
    QC:
      --skipQC                      Skip all QC steps apart from MultiQC
      --skipFastQC                  Skip FastQC
      --skipPreseq                  Skip Preseq
      --skipDupRadar                Skip dupRadar (and Picard MarkDuplicates)
      --skipQualimap                Skip Qualimap
      --skipBiotypeQC               Skip Biotype QC
      --skipRseQC                   Skip RSeQC
      --skipEdgeR                   Skip edgeR MDS plot and heatmap
      --skipMultiQC                 Skip MultiQC
    Other options
      --sampleLevel                 Used to turn off the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail               Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size      Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Check if genome exists in the config
try {
    params.genomes[params.genome]
} catch(MissingPropertyException){
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. "
}


// Reference index path configuration
// Define these here - after the profiles are loaded with the iGenomes paths
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
// Check if genome exists in the config 

if (!params.fasta) {
     exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. "
}
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
if(!params.gtf){
    params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false
}
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
if(params.aligner != 'star'){
    params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
}
ch_mdsplot_header = Channel.fromPath("${workflow.launchDir/params.config_files}/assets/mdsplot_header.txt", checkIfExists: true)
ch_heatmap_header = Channel.fromPath("${workflow.launchDir/params.config_files}/assets/heatmap_header.txt", checkIfExists: true)
ch_biotypes_header = Channel.fromPath("${workflow.launchDir/params.config_files}/assets/biotypes_header.txt", checkIfExists: true)
Channel.fromPath("${workflow.launchDir/params.config_files}/assets/where_are_my_files.txt", checkIfExists: true)
       .into{ch_where_trim_galore; ch_where_star; ch_where_hisat2; ch_where_hisat2_sort}

// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
forwardStranded = params.forwardStranded
reverseStranded = params.reverseStranded
unStranded = params.unStranded

// Preset trimming options
if (params.pico) {
    clip_r1 = 3
    clip_r2 = 0
    three_prime_clip_r1 = 0
    three_prime_clip_r2 = 3
    forwardStranded = true
    reverseStranded = false
    unStranded = false
}

// Get rRNA databases
// Default is set to bundled DB list in `assets/rrna-db-defaults.txt`

rRNA_database = file(params.rRNA_database_manifest)
if (rRNA_database.isEmpty()) {exit 1, "File ${rRNA_database.getName()} is empty!"}
Channel
    .from( rRNA_database.readLines() )
    .map { row -> file(row) }
    .set { sortmerna_fasta }

// Validate inputs
if (params.aligner != 'star' && params.aligner != 'hisat2') {
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2'"
}
if (params.pseudo_aligner && params.pseudo_aligner != 'salmon') {
    exit 1, "Invalid pseudo aligner option: ${params.pseudo_aligner}. Valid options: 'salmon'"
}

if (params.star_index && params.aligner == 'star' && !params.skipAlignment) {
  if (hasExtension(params.star_index, 'gz')) {
    star_index_gz = Channel
        .fromPath(params.star_index, checkIfExists: true)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
  } else{
    star_index = Channel
        .fromPath(params.star_index, checkIfExists: true)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
  }
}
else if (params.hisat2_index && params.aligner == 'hisat2' && !params.skipAlignment) {
  if (hasExtension(params.hisat2_index, 'gz')) {
    hs2_indices_gz = Channel
        .fromPath("${params.hisat2_index}", checkIfExists: true)
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
  } else {
    hs2_indices = Channel
        .fromPath("${params.hisat2_index}*", checkIfExists: true)
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
  }
}
else if (params.fasta && !params.skipAlignment) {
  if (hasExtension(params.fasta, 'gz')) {
    Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
        .set { genome_fasta_gz }
  } else {
    Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
        .into { ch_fasta_for_star_index; ch_fasta_for_hisat_index }
  }

} else if (params.skipAlignment) {
  println "Skipping alignment ..."
}
else {
    exit 1, "No reference genome files specified!"
}

if (params.aligner == 'hisat2' && params.splicesites) {
    Channel
        .fromPath(params.bed12, checkIfExists: true)
        .ifEmpty { exit 1, "HISAT2 splice sites file not found: $alignment_splicesites" }
        .into { indexing_splicesites; alignment_splicesites }
}

// Separately check for whether salmon needs a genome fasta to extract
// transcripts from, or can use a transcript fasta directly
if (params.pseudo_aligner == 'salmon') {
    if (params.salmon_index) {
      if (hasExtension(params.salmon_index, 'gz')) {
        salmon_index_gz = Channel
            .fromPath(params.salmon_index, checkIfExists: true)
            .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
      } else {
        salmon_index = Channel
            .fromPath(params.salmon_index, checkIfExists: true)
            .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
      }
    } else if (params.transcript_fasta) {
      if (hasExtension(params.transcript_fasta, 'gz')) {
        transcript_fasta_gz = Channel
            .fromPath(params.transcript_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Transcript fasta file not found: ${params.transcript_fasta}" }
      } else {
        ch_fasta_for_salmon_index = Channel
            .fromPath(params.transcript_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Transcript fasta file not found: ${params.transcript_fasta}" }
      }
    } else if (params.fasta && (params.gff || params.gtf)) {
      log.info "Extracting transcript fastas from genome fasta + gtf/gff"
      if (hasExtension(params.fasta, 'gz')) {
        Channel.fromPath(params.fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
            .set { genome_fasta_gz }
      } else {
        Channel.fromPath(params.fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
            .set { ch_fasta_for_salmon_transcripts }
      }
    } else {
      exit 1, "To use with `--pseudo_aligner 'salmon'`, must provide either --transcript_fasta or both --fasta and --gtf"
    }
}

if (params.gtf) {
  if (params.gff) {
      // Prefer gtf over gff
      log.info "Both GTF and GFF have been provided: Using GTF as priority."
  }
  if (hasExtension(params.gtf, 'gz')) {
  gtf_gz = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
  } else {
    Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_makeHisatSplicesites; gtf_makeHISATindex; gtf_makeSalmonIndex; gtf_makeBED12;
                gtf_star; gtf_dupradar; gtf_qualimap;  gtf_featureCounts; gtf_stringtieFPKM; gtf_salmon; gtf_salmon_merge }

  }
  } else if (params.gff) {
  if (hasExtension(params.gff, 'gz')) {
    gff_gz = Channel.fromPath(params.gff, checkIfExists: true)
                  .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }

  } else {
    gffFile = Channel.fromPath(params.gff, checkIfExists: true)
                  .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }

  }
} else {
    exit 1, "No GTF or GFF3 annotation specified!"
}

if (params.bed12) {
    bed12 = Channel
        .fromPath(params.bed12, checkIfExists: true)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .set { bed_rseqc }
}

if (params.gencode) {
  biotype = "gene_type"
} else {
  biotype = params.fc_group_features_type
}

if (params.skipAlignment && !params.pseudo_aligner) {
  exit 1, "--skipAlignment specified without --pseudo_aligner .. did you mean to specify --pseudo_aligner salmon"
}

if (workflow.profile == 'uppmax' || workflow.profile == 'uppmax-devel') {
    if (!params.project) exit 1, "No UPPMAX project ID found! Use --project"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}

if (workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)
ch_output_docs = file("${params.config_files}/docs/output.md", checkIfExists: true)

/*
 * Create a channel for input read files
 */
if (params.readPaths) {
    if (params.singleEnd) {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_trimgalore }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_trimgalore }
    }
} else {
    Channel
        .fromFilePairs( input, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${input}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_trimgalore }
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Reads'] = input
summary['Data Type'] = params.singleEnd ? 'Single-End' : 'Paired-End'
if (params.genome) summary['Genome'] = params.genome
if (params.pico) summary['Library Prep'] = "SMARTer Stranded Total RNA-Seq Kit - Pico Input"
summary['Strandedness'] = (unStranded ? 'None' : forwardStranded ? 'Forward' : reverseStranded ? 'Reverse' : 'None')
summary['Trimming'] = "5'R1: $clip_r1 / 5'R2: $clip_r2 / 3'R1: $three_prime_clip_r1 / 3'R2: $three_prime_clip_r2 / NextSeq Trim: $params.trim_nextseq"
if (params.aligner == 'star') {
    summary['Aligner'] = "STAR"
    if (params.star_index)summary['STAR Index'] = params.star_index
    else if (params.fasta)summary['Fasta Ref'] = params.fasta
} else if (params.aligner == 'hisat2') {
    summary['Aligner'] = "HISAT2"
    if (params.hisat2_index)summary['HISAT2 Index'] = params.hisat2_index
    else if (params.fasta)summary['Fasta Ref'] = params.fasta
    if (params.splicesites)summary['Splice Sites'] = params.splicesites
}
if (params.pseudo_aligner == 'salmon') {
    summary['Pseudo Aligner'] = "Salmon"
    if (params.transcript_fasta)summary['Transcript Fasta'] = params.transcript_fasta
}
if (params.gtf) summary['GTF Annotation'] = params.gtf
if (params.gff) summary['GFF3 Annotation'] = params.gff
if (params.bed12) summary['BED Annotation'] = params.bed12
if (params.gencode) summary['GENCODE'] = params.gencode
if (params.stringTieIgnoreGTF) summary['StringTie Ignore GTF'] = params.stringTieIgnoreGTF
summary['Remove Ribosomal RNA'] = params.removeRiboRNA
if (params.fc_group_features_type) summary['Biotype GTF field'] = biotype
summary['Save prefs'] = "Ref Genome: "+(params.saveReference ? 'Yes' : 'No')+" / Trimmed FastQ: "+(params.saveTrimmed ? 'Yes' : 'No')+" / Alignment intermediates: "+(params.saveAlignedIntermediates ? 'Yes' : 'No')
summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir'] = params.outdir
summary['Launch dir'] = workflow.launchDir
summary['Working dir'] = workflow.workDir
summary['Script dir'] = workflow.projectDir
summary['User'] = workflow.userName
if (workflow.profile == 'awsbatch') {
  summary['AWS Region']     = params.awsregion
  summary['AWS Queue']      = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
  summary['E-mail Address']    = params.email
  summary['E-mail on failure'] = params.email_on_fail
  summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-rnaseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/rnaseq Workflow Summary'
    section_href: 'https://github.com/nf-core/rnaseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * Parse software version numbers
 */
process get_software_versions {

    container docker_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf(".csv") > 0) filename
            else null
        }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version &> v_ngi_rnaseq.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    sortmerna --version &> v_sortmerna.txt
    STAR --version &> v_star.txt
    hisat2 --version &> v_hisat2.txt
    stringtie --version &> v_stringtie.txt
    preseq &> v_preseq.txt
    read_duplication.py --version &> v_rseqc.txt
    echo \$(bamCoverage --version 2>&1) > v_deeptools.txt
    featureCounts -v &> v_featurecounts.txt
    salmon --version &> v_salmon.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    samtools --version &> v_samtools.txt
    multiqc --version &> v_multiqc.txt
    Rscript -e "library(edgeR); write(x=as.character(packageVersion('edgeR')), file='v_edgeR.txt')"
    Rscript -e "library(dupRadar); write(x=as.character(packageVersion('dupRadar')), file='v_dupRadar.txt')"
    unset DISPLAY && qualimap rnaseq  > v_qualimap.txt 2>&1 || true
    python ${workflow.launchDir/params.config_files}/bin/scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

compressedReference = hasExtension(params.fasta, 'gz') || hasExtension(params.transcript_fasta, 'gz') || hasExtension(params.star_index, 'gz') || hasExtension(params.hisat2_index, 'gz')

if (compressedReference) {
  // This complex logic is to prevent accessing the genome_fasta_gz variable if
  // necessary indices for STAR, HiSAT2, Salmon already exist, or if
  // params.transcript_fasta is provided as then the transcript sequences dont need to be extracted.
  need_star_index = params.aligner == 'star' && !params.star_index
  need_hisat2_index = params.aligner == 'hisat2' && !params.hisat2_index
  need_aligner_index = need_hisat2_index || need_star_index
  alignment_no_indices = !params.skipAlignment && need_aligner_index
  pseudoalignment_no_indices = params.pseudo_aligner == "salmon" && !(params.transcript_fasta || params.salmon_index)
  if (params.fasta && (alignment_no_indices || pseudoalignment_no_indices)) {
    process gunzip_genome_fasta {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from genome_fasta_gz

        output:
        file "${gz.baseName}" into ch_fasta_for_star_index, ch_fasta_for_hisat_index, ch_fasta_for_salmon_transcripts

        script:
        """
        gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
  }
  if (params.gtf) {
    process gunzip_gtf {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from gtf_gz

        output:
        file "${gz.baseName}" into gtf_makeSTARindex, gtf_makeHisatSplicesites, gtf_makeHISATindex, gtf_makeSalmonIndex, gtf_makeBED12,
                                        gtf_star, gtf_dupradar, gtf_featureCounts, gtf_stringtieFPKM, gtf_salmon, gtf_salmon_merge, gtf_qualimap

        script:
        """
        gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
  }
  if (params.gff && !params.gtf) {
    process gunzip_gff {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from gff_gz

        output:
        file "${gz.baseName}" into gffFile

        script:
        """
        gunzip --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
  }
  if (params.transcript_fasta && params.pseudo_aligner == 'salmon' && !params.salmon_index) {
    process gunzip_transcript_fasta {
        container docker_container
       
        
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_transcriptome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from transcript_fasta_gz

        output:
        file "${gz.baseName}" into ch_fasta_for_salmon_index

        script:
        """
        gunzip --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
  }
  if (params.bed12) {
    process gunzip_bed12 {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from bed12_gz

        output:
        file "${gz.baseName}" into bed_rseqc

        script:
        """
        gunzip --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
  }
  if (!params.skipAlignment && params.star_index && params.aligner == "star") {
    process gunzip_star_index {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome/star" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from star_index_gz

        output:
        file "${gz.simpleName}" into star_index

        script:
        // Use tar as the star indices are a folder, not a file
        """
        tar -xzvf ${gz}
        """
    }
  }
  if (!params.skipAlignment && params.hisat2_index && params.aligner == 'hisat2') {
    process gunzip_hisat_index {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome/hisat2" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from hs2_indices_gz

        output:
        file "*.ht2*" into hs2_indices

        script:
        // Use tar as the hisat2 indices are a folder, not a file
        """
        tar -xzvf ${gz}
        """
    }
  }
  if (params.salmon_index && params.pseudo_aligner == 'salmon') {
    process gunzip_salmon_index {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gz"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_transcriptome/hisat2" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gz from salmon_index_gz

        output:
        file "${gz.simpleName}" into salmon_index

        script:
        // Use tar as the hisat2 indices are a folder, not a file
        """
        tar -xzvf ${gz}
        """
    }
  }
}

/*
 * PREPROCESSING - Convert GFF3 to GTF
 */
if (params.gff && !params.gtf) {
    process convertGFFtoGTF {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gff"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gff from gffFile

        output:
        file "${gff.baseName}.gtf" into gtf_makeSTARindex, gtf_makeHisatSplicesites, gtf_makeHISATindex, gtf_makeSalmonIndex, gtf_makeBED12,
                                        gtf_star, gtf_dupradar, gtf_featureCounts, gtf_stringtieFPKM, gtf_salmon, gtf_salmon_merge, gtf_qualimap

        script:
        """
        gffread $gff --keep-exon-attrs -F -T -o ${gff.baseName}.gtf
        """
    }
}

/*
 * PREPROCESSING - Build BED12 file
 */
if (!params.bed12) {
    process makeBED12 {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeBED12

        output:
        file "${gtf.baseName}.bed" into bed_rseqc

        script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
        """
        perl ${workflow.launchDir/params.config_files}/bin/gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}

/*
 * PREPROCESSING - Build STAR index
 */
if (!params.skipAlignment) {
  if (params.aligner == 'star' && !params.star_index && params.fasta) {
      process makeSTARindex {
       
        
        container docker_container
        memory '16G'
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
          
          tag "$fasta"
          publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                     saveAs: { params.saveReference ? it : null }, mode: 'copy'

          input:
          file fasta from ch_fasta_for_star_index
          file gtf from gtf_makeSTARindex

          output:
          file "star" into star_index

          script:
          def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
          """
          mkdir star
          STAR \\
              --runMode genomeGenerate \\
              --runThreadN ${task.cpus} \\
              --sjdbGTFfile $gtf \\
              --genomeDir star/ \\
              --genomeFastaFiles $fasta \\
              $avail_mem
          """
      }
  }

  /*
   * PREPROCESSING - Build HISAT2 splice sites file
   */
  if (params.aligner == 'hisat2' && !params.splicesites) {
      process makeHisatSplicesites {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
          tag "$gtf"
          publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                     saveAs: { params.saveReference ? it : null }, mode: 'copy'

          input:
          file gtf from gtf_makeHisatSplicesites

          output:
          file "${gtf.baseName}.hisat2_splice_sites.txt" into indexing_splicesites, alignment_splicesites

          script:
          """
          python ${workflow.launchDir/params.config_files}/bin/hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
          """
      }
  }

  /*
   * PREPROCESSING - Build HISAT2 index
   */
  if (params.aligner == 'hisat2' && !params.hisat2_index && params.fasta) {
      process makeHISATindex {
       
        
        container docker_container
        memory '200G'
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
          tag "$fasta"
          publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                     saveAs: { params.saveReference ? it : null }, mode: 'copy'

          input:
          file fasta from ch_fasta_for_hisat_index
          file indexing_splicesites from indexing_splicesites
          file gtf from gtf_makeHISATindex

          output:
          file "${fasta.baseName}.*.ht2*" into hs2_indices

          script:
          if (!task.memory) {
              log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
              avail_mem = 0
          } else {
              log.info "[HISAT2 index build] Available memory: ${task.memory}"
              avail_mem = task.memory.toGiga()
          }
          if (avail_mem > params.hisat_build_memory) {
              log.info "[HISAT2 index build] Over ${params.hisat_build_memory} GB available, so using splice sites and exons in HISAT2 index"
              extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt"
              ss = "--ss $indexing_splicesites"
              exon = "--exon ${gtf.baseName}.hisat2_exons.txt"
          } else {
              log.info "[HISAT2 index build] Less than ${params.hisat_build_memory} GB available, so NOT using splice sites and exons in HISAT2 index."
              log.info "[HISAT2 index build] Use --hisat_build_memory [small number] to skip this check."
              extract_exons = ''
              ss = ''
              exon = ''
          }
          """
          $extract_exons
          hisat2-build -p ${task.cpus} $ss $exon $fasta ${fasta.baseName}.hisat2_index
          """
      }
  }
}


/*
 * PREPROCESSING - Create Salmon transcriptome index
 */
if (params.pseudo_aligner == 'salmon' && !params.salmon_index) {
    if (!params.transcript_fasta) {
        process transcriptsToFasta {
            container docker_container
           
            
            pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
            tag "$fasta"
            publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                               saveAs: { params.saveReference ? it : null }, mode: 'copy'


            input:
            file fasta from ch_fasta_for_salmon_transcripts
            file gtf from gtf_makeSalmonIndex

            output:
            file "*.fa" into ch_fasta_for_salmon_index

            script:
	          // filter_gtf_for_genes_in_genome.py is bundled in this package, in rnaseq/bin
            """
            python ${workflow.launchDir/params.config_files}/bin/filter_gtf_for_genes_in_genome.py --gtf $gtf --fasta $fasta -o ${gtf.baseName}__in__${fasta.baseName}.gtf
            gffread -F -w transcripts.fa -g $fasta ${gtf.baseName}__in__${fasta.baseName}.gtf
            """
        }
    }
    process makeSalmonIndex {
        container docker_container
       
        
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
        
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                           saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_salmon_index

        output:
        file 'salmon_index' into salmon_index

        script:
        def gencode = params.gencode  ? "--gencode" : ""
        """
        salmon index --threads $task.cpus -t $fasta $gencode -i salmon_index
        """
    }
}

/*
 * STEP 1 - FastQC
 */
process fastqc {

    container docker_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    tag "$name"
    
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }

    when:
    !params.skipQC && !params.skipFastQC

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
if (!params.skipTrimming) {
    process trim_galore {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
        
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
                else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
                else null
            }

        input:
        set val(name), file(reads) from raw_reads_trimgalore
        file wherearemyfiles from ch_where_trim_galore.collect()

        output:
        set val(name), file("*fq.gz") into trimgalore_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
        file "where_are_my_files.txt"


        script:
        c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
        c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
        tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
        tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
        nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $c_r1 $tpc_r1 $nextseq $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $nextseq $reads
            """
        }
    }
}else{
   raw_reads_trimgalore
       .set {trimgalore_reads}
   trimgalore_results = Channel.empty()
}


/*
 * STEP 2+ - SortMeRNA - remove rRNA sequences on request
 */
if (!params.removeRiboRNA) {
    trimgalore_reads
        .into { trimmed_reads_alignment; trimmed_reads_salmon }
    sortmerna_logs = Channel.empty()
} else {
    process sortmerna_index {
        container docker_container
       
        
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
        
        tag "${fasta.baseName}"

        input:
        file(fasta) from sortmerna_fasta

        output:
        val("${fasta.baseName}") into sortmerna_db_name
        file("$fasta") into sortmerna_db_fasta
        file("${fasta.baseName}*") into sortmerna_db

        script:
        """
        indexdb_rna --ref $fasta,${fasta.baseName} -m 3072 -v
        """
    }

    process sortmerna {
        container docker_container
       
        
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    
        tag "$name"
        publishDir "${params.outdir}/SortMeRNA", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_rRNA_report.txt") > 0) "logs/$filename"
                else if (params.saveNonRiboRNAReads) "reads/$filename"
                else null
            }

        input:
        set val(name), file(reads) from trimgalore_reads
        val(db_name) from sortmerna_db_name.collect()
        file(db_fasta) from sortmerna_db_fasta.collect()
        file(db) from sortmerna_db.collect()

        output:
        set val(name), file("*.fq.gz") into trimmed_reads_alignment, trimmed_reads_salmon
        file "*_rRNA_report.txt" into sortmerna_logs


        script:
        //concatenate reference files: ${db_fasta},${db_name}:${db_fasta},${db_name}:...
        def Refs = ''
        for (i=0; i<db_fasta.size(); i++) { Refs+= ":${db_fasta[i]},${db_name[i]}" }
        Refs = Refs.substring(1)

        if (params.singleEnd) {
            """
            gzip -d --force < ${reads} > all-reads.fastq
            sortmerna --ref ${Refs} \
                --reads all-reads.fastq \
                --num_alignments 1 \
                -a ${task.cpus} \
                --fastx \
                --aligned rRNA-reads \
                --other non-rRNA-reads \
                --log -v
            gzip --force < non-rRNA-reads.fastq > ${name}.fq.gz
            mv rRNA-reads.log ${name}_rRNA_report.txt
            """
        } else {
            """
            gzip -d --force < ${reads[0]} > reads-fw.fq
            gzip -d --force < ${reads[1]} > reads-rv.fq
            merge-paired-reads.sh reads-fw.fq reads-rv.fq all-reads.fastq
            sortmerna --ref ${Refs} \
                --reads all-reads.fastq \
                --num_alignments 1 \
                -a ${task.cpus} \
                --fastx --paired_in \
                --aligned rRNA-reads \
                --other non-rRNA-reads \
                --log -v
            unmerge-paired-reads.sh non-rRNA-reads.fastq non-rRNA-reads-fw.fq non-rRNA-reads-rv.fq
            gzip < non-rRNA-reads-fw.fq > ${name}-fw.fq.gz
            gzip < non-rRNA-reads-rv.fq > ${name}-rv.fq.gz
            mv rRNA-reads.log ${name}_rRNA_report.txt
            """
        }
    }
}    

/*
 * STEP 3 - align with STAR
 */
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skipped_poor_alignment = []
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if (percent_aligned.toFloat() <= '5'.toFloat()) {
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        skipped_poor_alignment << logname
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}
if (!params.skipAlignment) {
  if (params.aligner == 'star') {
      hisat_stdout = Channel.from(false)
      process star {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: instance_size
         
          tag "$name"
          publishDir "${params.outdir}/STAR", mode: 'copy',
              saveAs: {filename ->
                  if (filename.indexOf(".bam") == -1) "logs/$filename"
                  else if (params.saveUnaligned && filename != "where_are_my_files.txt" && 'Unmapped' in filename) unmapped/filename
                  else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                  else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                  else null
              }

          input:
          set val(name), file(reads) from trimmed_reads_alignment
          file index from star_index.collect()
          file gtf from gtf_star.collect()
          file wherearemyfiles from ch_where_star.collect()

          output:
          set file("*Log.final.out"), file ('*.bam') into star_aligned
          file "*.out" into alignment_logs
          file "*SJ.out.tab"
          file "*Log.out" into star_log
          file "where_are_my_files.txt"
          file "*Unmapped*" optional true
          file "${prefix}Aligned.sortedByCoord.out.bam.bai" into bam_index_rseqc, bam_index_genebody

          script:
          prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
          def star_mem = task.memory ?: params.star_memory ?: false
          def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''
          seq_center = params.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$params.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
          unaligned = params.saveUnaligned ? "--outReadsUnmapped Fastx" : ''
          """
          STAR --genomeDir $index \\
              --sjdbGTFfile $gtf \\
              --readFilesIn $reads  \\
              --runThreadN ${task.cpus} \\
              --twopassMode Basic \\
              --outWigType bedGraph \\
              --outSAMtype BAM SortedByCoordinate $avail_mem \\
              --readFilesCommand zcat \\
              --runDirPerm All_RWX $unaligned \\
              --outFileNamePrefix $prefix $seq_center
          samtools index ${prefix}Aligned.sortedByCoord.out.bam
          """
      }
      // Filter removes all 'aligned' channels that fail the check
      star_aligned
          .filter { logs, bams -> check_log(logs) }
          .flatMap {  logs, bams -> bams }
      .into { bam_count; bam_rseqc; bam_qualimap; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM; bam_forSubsamp; bam_skipSubsamp  }
  }


  /*
   * STEP 3 - align with HISAT2
   */
  if (params.aligner == 'hisat2') {
      star_log = Channel.from(false)
      process hisat2Align {
        container docker_container
       
        
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
        
          tag "$name"
          publishDir "${params.outdir}/HISAT2", mode: 'copy',
              saveAs: {filename ->
                  if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
                  else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                  else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                  else null
              }

          input:
          set val(name), file(reads) from trimmed_reads_alignment
          file hs2_indices from hs2_indices.collect()
          file alignment_splicesites from alignment_splicesites.collect()
          file wherearemyfiles from ch_where_hisat2.collect()

          output:
          file "${prefix}.bam" into hisat2_bam
          file "${prefix}.hisat2_summary.txt" into alignment_logs
          file "where_are_my_files.txt"
          file "unmapped.hisat2*" optional true

          script:
          index_base = hs2_indices[0].toString() - ~/.\d.ht2l?/
          prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
          seq_center = params.seq_center ? "--rg-id ${prefix} --rg CN:${params.seq_center.replaceAll('\\s','_')} SM:$prefix" : "--rg-id ${prefix} --rg SM:$prefix"
          def rnastrandness = ''
          if (forwardStranded && !unStranded) {
              rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
          } else if (reverseStranded && !unStranded) {
              rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
          }

          if (params.singleEnd) {
              unaligned = params.saveUnaligned ? "--un-gz unmapped.hisat2.gz" : ''
              """
              hisat2 -x $index_base \\
                     -U $reads \\
                     $rnastrandness \\
                     --known-splicesite-infile $alignment_splicesites \\
                     -p ${task.cpus} $unaligned\\
                     --met-stderr \\
                     --new-summary \\
                     --dta \\
                     --summary-file ${prefix}.hisat2_summary.txt $seq_center \\
                     | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
              """
          } else {
              unaligned = params.saveUnaligned ? "--un-conc-gz unmapped.hisat2.gz" : ''
              """
              hisat2 -x $index_base \\
                     -1 ${reads[0]} \\
                     -2 ${reads[1]} \\
                     $rnastrandness \\
                     --known-splicesite-infile $alignment_splicesites \\
                     --no-mixed \\
                     --no-discordant \\
                     -p ${task.cpus} $unaligned\\
                     --met-stderr \\
                     --new-summary \\
                     --summary-file ${prefix}.hisat2_summary.txt $seq_center \\
                     | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
              """
          }
      }

      process hisat2_sortOutput {
       
        
        container docker_container
        memory '70G'
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
          
          tag "${hisat2_bam.baseName}"
          publishDir "${params.outdir}/HISAT2", mode: 'copy',
              saveAs: { filename ->
                  if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                  else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") "aligned_sorted/$filename"
                  else null
              }

          input:
          file hisat2_bam
          file wherearemyfiles from ch_where_hisat2_sort.collect()

          output:
          file "${hisat2_bam.baseName}.sorted.bam" into bam_count, bam_rseqc, bam_qualimap, bam_preseq, bam_markduplicates, bam_featurecounts, bam_stringtieFPKM,bam_forSubsamp, bam_skipSubsamp
          file "${hisat2_bam.baseName}.sorted.bam.bai" into bam_index_rseqc, bam_index_genebody
          file "where_are_my_files.txt"

          script:
          def suff_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
          def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / task.cpus}" : ''
          """
          samtools sort \\
              $hisat2_bam \\
              -@ ${task.cpus} ${avail_mem} \\
              -o ${hisat2_bam.baseName}.sorted.bam
          samtools index ${hisat2_bam.baseName}.sorted.bam
          """
      }
  }


  /*
   * STEP 4 - RSeQC analysis
   */
  process rseqc {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      
      tag "${bam_rseqc.baseName - '.sorted'}"
      publishDir "${params.outdir}/rseqc" , mode: 'copy',
          saveAs: {filename ->
                   if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
              else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
              else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
              else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
              else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
              else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
              else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
              else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
              else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
              else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
              else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
              else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
              else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
              else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
              else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
              else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
              else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
              else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
              else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
              else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
              else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
              else filename
          }

      when:
      !params.skipQC && !params.skipRseQC

      input:
      file bam_rseqc
      file index from bam_index_rseqc
      file bed12 from bed_rseqc.collect()

      output:
      file "*.{txt,pdf,r,xls}" into rseqc_results

      script:
      """
      infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
      junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
      bam_stat.py -i $bam_rseqc 2> ${bam_rseqc.baseName}.bam_stat.txt
      junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 2> ${bam_rseqc.baseName}.junction_annotation_log.txt
      inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
      read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.read_distribution.txt
      read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication
      """
  }

  /*
   * STEP 5 - preseq analysis
   */
  process preseq {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      tag "${bam_preseq.baseName - '.sorted'}"
      publishDir "${params.outdir}/preseq", mode: 'copy'

      when:
      !params.skipQC && !params.skipPreseq

      input:
      file bam_preseq

      output:
      file "${bam_preseq.baseName}.ccurve.txt" into preseq_results

      script:
      """
      preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq.baseName}.ccurve.txt
      """
  }

  /*
   * STEP 6 - Mark duplicates
   */
  process markDuplicates {
       
        
        container docker_container
        memory '70G'
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      tag "${bam.baseName - '.sorted'}"
      publishDir "${params.outdir}/markDuplicates", mode: 'copy',
          saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

      when:
      !params.skipQC && !params.skipDupRadar

      input:
      file bam from bam_markduplicates

      output:
      file "${bam.baseName}.markDups.bam" into bam_md
      file "${bam.baseName}.markDups_metrics.txt" into picard_results
      file "${bam.baseName}.markDups.bam.bai"

      script:
      markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""
      """
      picard ${markdup_java_options} MarkDuplicates \\
          INPUT=$bam \\
          OUTPUT=${bam.baseName}.markDups.bam \\
          METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
          REMOVE_DUPLICATES=false \\
          ASSUME_SORTED=true \\
          PROGRAM_RECORD_ID='null' \\
          VALIDATION_STRINGENCY=LENIENT
      samtools index ${bam.baseName}.markDups.bam
      """
  }

  /*
   * STEP 7 - Qualimap
   */
  process qualimap {
       
        
        container docker_container
        memory '70G'
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    
      tag "${bam.baseName}"
      publishDir "${params.outdir}/qualimap", mode: 'copy'

      when:
      !params.skipQC && !params.skipQualimap

      input:
      file bam from bam_qualimap
      file gtf from gtf_qualimap.collect()

      output:
      file "${bam.baseName}" into qualimap_results

      script:
      def qualimap_direction = 'non-strand-specific'
      if (forwardStranded) {
          qualimap_direction = 'strand-specific-forward'
      }else if (reverseStranded) {
          qualimap_direction = 'strand-specific-reverse'
      }
      def paired = params.singleEnd ? '' : '-pe'
      memory = task.memory.toGiga() + "G"
      """
      unset DISPLAY
      qualimap --java-mem-size=${memory} rnaseq $qualimap_direction $paired -s -bam $bam -gtf $gtf -outdir ${bam.baseName}
      """
  }

  /*
   * STEP 8 - dupRadar
   */
  process dupradar {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      
      tag "${bam_md.baseName - '.sorted.markDups'}"
      publishDir "${params.outdir}/dupradar", mode: 'copy',
          saveAs: {filename ->
              if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
              else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
              else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
              else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
              else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
              else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
              else "$filename"
          }

      when:
      !params.skipQC && !params.skipDupRadar

      input:
      file bam_md
      file gtf from gtf_dupradar.collect()

      output:
      file "*.{pdf,txt}" into dupradar_results

      script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
      def dupradar_direction = 0
      if (forwardStranded && !unStranded) {
          dupradar_direction = 1
      } else if (reverseStranded && !unStranded) {
          dupradar_direction = 2
      }
      def paired = params.singleEnd ? 'single' :  'paired'
      """
      Rscript ${workflow.launchDir/params.config_files}/bin/dupRadar.r $bam_md $gtf $dupradar_direction $paired ${task.cpus}
      """
  }

  /*
   * STEP 9 - Feature counts
   */
  process featureCounts {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      
      tag "${bam_featurecounts.baseName - '.sorted'}"
      publishDir "${params.outdir}/featureCounts", mode: 'copy',
          saveAs: {filename ->
              if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
              else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
              else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
              else "$filename"
          }

      input:
      file bam_featurecounts
      file gtf from gtf_featureCounts.collect()
      file biotypes_header from ch_biotypes_header.collect()

      output:
      file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
      file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
      file "${bam_featurecounts.baseName}_biotype_counts*mqc.{txt,tsv}" optional true into featureCounts_biotype

      script:
      def featureCounts_direction = 0
      def extraAttributes = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
      if (forwardStranded && !unStranded) {
          featureCounts_direction = 1
      } else if (reverseStranded && !unStranded) {
          featureCounts_direction = 2
      }
      // Try to get real sample name
      sample_name = bam_featurecounts.baseName - 'Aligned.sortedByCoord.out' - '_subsamp.sorted'
      biotype_qc = params.skipBiotypeQC ? '' : "featureCounts -a $gtf -g $biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts"
      mod_biotype = params.skipBiotypeQC ? '' : "cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${bam_featurecounts.baseName}_biotype_counts_mqc.txt && python ${workflow.launchDir/params.config_files}/bin/mqc_features_stat.py ${bam_featurecounts.baseName}_biotype_counts_mqc.txt -s $sample_name -f rRNA -o ${bam_featurecounts.baseName}_biotype_counts_gs_mqc.tsv"
      """
      featureCounts -a $gtf -g ${params.fc_group_features} -t ${params.fc_count_type} -o ${bam_featurecounts.baseName}_gene.featureCounts.txt $extraAttributes -p -s $featureCounts_direction $bam_featurecounts
      $biotype_qc
      $mod_biotype
      """
  }



  /*
   * STEP 10 - Merge featurecounts
   */
  process merge_featureCounts {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      
      tag "${input_files[0].baseName - '.sorted'}"
      publishDir "${params.outdir}/featureCounts", mode: 'copy'

      input:
      file input_files from featureCounts_to_merge.collect()

      output:
      file 'merged_gene_counts.txt' into featurecounts_merged

      script:
      // Redirection (the `<()`) for the win!
      // Geneid in 1st column and gene_name in 7th
      gene_ids = "<(tail -n +2 ${input_files[0]} | cut -f1,7 )"
      counts = input_files.collect{filename ->
        // Remove first line and take third column
        "<(tail -n +2 ${filename} | sed 's:.sorted.bam::' | cut -f8)"}.join(" ")
      """
      paste $gene_ids $counts > merged_gene_counts.txt
      """
  }

  /*
   * STEP 12 - stringtie FPKM
   */
  process stringtieFPKM {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      tag "${bam_stringtieFPKM.baseName - '.sorted'}"
      publishDir "${params.outdir}/stringtieFPKM", mode: 'copy',
          saveAs: {filename ->
              if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
              else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
              else if (filename.indexOf("ballgown") > 0) "ballgown/$filename"
              else "$filename"
          }

      input:
      file bam_stringtieFPKM
      file gtf from gtf_stringtieFPKM.collect()

      output:
      file "${bam_stringtieFPKM.baseName}_transcripts.gtf"
      file "${bam_stringtieFPKM.baseName}.gene_abund.txt"
      file "${bam_stringtieFPKM}.cov_refs.gtf"
      file "${bam_stringtieFPKM.baseName}_ballgown"

      script:
      def st_direction = ''
      if (forwardStranded && !unStranded) {
          st_direction = "--fr"
      } else if (reverseStranded && !unStranded) {
          st_direction = "--rf"
      }
      def ignore_gtf = params.stringTieIgnoreGTF ? "" : "-e"
      """
      stringtie $bam_stringtieFPKM \\
          $st_direction \\
          -o ${bam_stringtieFPKM.baseName}_transcripts.gtf \\
          -v \\
          -G $gtf \\
          -A ${bam_stringtieFPKM.baseName}.gene_abund.txt \\
          -C ${bam_stringtieFPKM}.cov_refs.gtf \\
          -b ${bam_stringtieFPKM.baseName}_ballgown \\
          $ignore_gtf
      """
  }

  /*
   * STEP 13 - edgeR MDS and heatmap
   */
  process sample_correlation {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      
      tag "${input_files[0].toString() - '.sorted_gene.featureCounts.txt' - 'Aligned'}"
      publishDir "${params.outdir}/sample_correlation", mode: 'copy'

      when:
      !params.skipQC && !params.skipEdgeR

      input:
      file input_files from geneCounts.collect()
      val num_bams from bam_count.count()
      file mdsplot_header from ch_mdsplot_header
      file heatmap_header from ch_heatmap_header

      output:
      file "*.{txt,pdf,csv}" into sample_correlation_results

      when:
      num_bams > 2 && (!params.sampleLevel)

      script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
      """
      Rscript ${workflow.launchDir/params.config_files}/bin/edgeR_heatmap_MDS.r $input_files
      cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
      mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
      cat $heatmap_header log2CPM_sample_correlation_mqc.csv >> tmp_file
      mv tmp_file log2CPM_sample_correlation_mqc.csv
      """
  }


} else {
  star_log = Channel.from(false)
  hisat_stdout = Channel.from(false)
  alignment_logs = Channel.from(false)
  rseqc_results = Channel.from(false)
  qualimap_results = Channel.from(false)
  sample_correlation_results = Channel.from(false)
  featureCounts_logs = Channel.from(false)
  dupradar_results = Channel.from(false)
  preseq_results = Channel.from(false)
  featureCounts_biotype = Channel.from(false)
}


/*
 * STEP 11 - Transcriptome quantification with Salmon
 */
if (params.pseudo_aligner == 'salmon') {
    process salmon {
        container docker_container
       
        
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
        
        tag "$sample"
        publishDir "${params.outdir}/salmon", mode: 'copy'

        input:
        set sample, file(reads) from trimmed_reads_salmon
        file index from salmon_index.collect()
        file gtf from gtf_salmon.collect()

        output:
        file "${sample}/" into salmon_logs
        set val(sample), file("${sample}/") into salmon_tximport, salmon_parsegtf

        script:
        def rnastrandness = params.singleEnd ? 'U' : 'IU'
        if (forwardStranded && !unStranded) {
            rnastrandness = params.singleEnd ? 'SF' : 'ISF'
        } else if (reverseStranded && !unStranded) {
            rnastrandness = params.singleEnd ? 'SR' : 'ISR'
        }
        def endedness = params.singleEnd ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
        unmapped = params.saveUnaligned ? "--writeUnmappedNames" : ''
        """
        salmon quant --validateMappings \\
                        --seqBias --useVBOpt --gcBias \\
                        --geneMap ${gtf} \\
                        --threads ${task.cpus} \\
                        --libType=${rnastrandness} \\
                        --index ${index} \\
                        $endedness $unmapped\\
                        -o ${sample}
        """
        }


    process salmon_tx2gene {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      
      publishDir "${params.outdir}/salmon", mode: 'copy'

      input:
      file ("salmon/*") from salmon_parsegtf.collect()
      file gtf from gtf_salmon_merge

      output:
      file "tx2gene.csv" into salmon_tx2gene, salmon_merge_tx2gene

      script:
      """
      python ${workflow.launchDir/params.config_files}/bin/parse_gtf.py --gtf $gtf --salmon salmon --id ${params.fc_group_features} --extra ${params.fc_extra_attributes} -o tx2gene.csv
      """
    }

    process salmon_tximport {
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: instance_size
      
      publishDir "${params.outdir}/salmon", mode: 'copy'

      input:
      set val(name), file ("salmon/*") from salmon_tximport
      file tx2gene from salmon_tx2gene.collect()

      output:
      file "${name}_salmon_gene_tpm.csv" into salmon_gene_tpm
      file "${name}_salmon_gene_counts.csv" into salmon_gene_counts
      file "${name}_salmon_transcript_tpm.csv" into salmon_transcript_tpm
      file "${name}_salmon_transcript_counts.csv" into salmon_transcript_counts

      script:
      """
      Rscript ${workflow.launchDir/params.config_files}/bin/tximport.r NULL salmon ${name}
      """
    }

    process salmon_merge {
       
        
        container docker_container
        pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
      publishDir "${params.outdir}/salmon", mode: 'copy'

      input:
      file gene_tpm_files from salmon_gene_tpm.collect()
      file gene_count_files from salmon_gene_counts.collect()
      file transcript_tpm_files from salmon_transcript_tpm.collect()
      file transcript_count_files from salmon_transcript_counts.collect()
      file tx2gene from salmon_merge_tx2gene

      output:
      file "salmon_merged*.csv" into salmon_merged_ch
      file "*.rds"

      script:
      // First field is the gene/transcript ID
      gene_ids = "<(cut -f1 -d, ${gene_tpm_files[0]} | tail -n +2 | cat <(echo '${params.fc_group_features}') - )"
      transcript_ids = "<(cut -f1 -d, ${transcript_tpm_files[0]} | tail -n +2 | cat <(echo 'transcript_id') - )"

      // Second field is counts/TPM
      gene_tpm = gene_tpm_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
      gene_counts = gene_count_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
      transcript_tpm = transcript_tpm_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
      transcript_counts = transcript_count_files.collect{f -> "<(cut -d, -f2 ${f})"}.join(" ")
      """
      paste -d, $gene_ids $gene_tpm > salmon_merged_gene_tpm.csv
      paste -d, $gene_ids $gene_counts > salmon_merged_gene_counts.csv
      paste -d, $transcript_ids $transcript_tpm > salmon_merged_transcript_tpm.csv
      paste -d, $transcript_ids $transcript_counts > salmon_merged_transcript_counts.csv
      Rscript ${workflow.launchDir/params.config_files}/bin/se.r NULL salmon_merged_gene_counts.csv salmon_merged_gene_tpm.csv
      Rscript ${workflow.launchDir/params.config_files}/bin/se.r NULL salmon_merged_transcript_counts.csv salmon_merged_transcript_tpm.csv
      """
    }
} else {
    salmon_logs = Channel.empty()
}


/*
 * STEP 14 - MultiQC
 */
process multiqc {

    container docker_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-medium'
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skipMultiQC

    input:
    file multiqc_config from ch_multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimgalore/*') from trimgalore_results.collect().ifEmpty([])
    file ('alignment/*') from alignment_logs.collect().ifEmpty([])
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    file ('qualimap/*') from qualimap_results.collect().ifEmpty([])
    file ('preseq/*') from preseq_results.collect().ifEmpty([])
    file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    file ('featureCounts/*') from featureCounts_logs.collect().ifEmpty([])
    file ('featureCounts_biotype/*') from featureCounts_biotype.collect()
    file ('salmon/*') from salmon_logs.collect().ifEmpty([])
    file ('sample_correlation_results/*') from sample_correlation_results.collect().ifEmpty([]) // If the Edge-R is not run create an Empty array
    file ('sortmerna/*') from sortmerna_logs.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    file "multiqc*"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m picard -m preseq -m rseqc -m featureCounts -m hisat2 -m star -m cutadapt -m sortmerna -m fastqc -m qualimap -m salmon
    """
}

/*
 * STEP 15 - Output Description HTML
 */
process output_documentation {

    container docker_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    Rscript ${workflow.launchDir/params.config_files}/bin/markdown_to_html.r $output_docs results_description.html
    """
}

/*
 * STEP 16 copy intermediate work files and logs
 */

process copy_workfiles {

    container docker_container
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
        file(multiqc_reportz) from multiqc_report


    script:
    """
    cp -r ${workflow.workDir} ${workflow.launchDir}/out
    cp -r ${workflow.launchDir}/.nextflow.log ${workflow.launchDir}/out
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/rnaseq] Successful: $workflow.runName"
    if (skipped_poor_alignment.size() > 0) {
        subject = "[nf-core/rnaseq] Partially Successful (${skipped_poor_alignment.size()} skipped): $workflow.runName"
    }
    if (!workflow.success) {
      subject = "[nf-core/rnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
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
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['skipped_poor_alignment'] = skipped_poor_alignment
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success && !params.skipMultiQC) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/rnaseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/rnaseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
          if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/rnaseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, email_address ].execute() << email_txt
          log.info "[nf-core/rnaseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = file("${output_d}/pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = file("${output_d}/pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (skipped_poor_alignment.size() > 0) {
        log.info "${c_purple}[nf-core/rnaseq]${c_red} WARNING - ${skipped_poor_alignment.size()} samples skipped due to poor mapping percentages!${c_reset}"
    }
    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[nf-core/rnaseq]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/rnaseq]${c_red} Pipeline completed with errors${c_reset}"
    }

}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/rnaseq v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}