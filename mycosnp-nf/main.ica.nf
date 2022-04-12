#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/mycosnp
========================================================================================
    Github : https://github.com/CDCgov/mycosnp-nf
    Wiki   : https://github.com/CDCgov/mycosnp-nf/wiki
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl = 2
/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/
params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.input = null
params.add_sra_file = null
params.add_vcf_file = null
params.genome = null
params.igenomes_base = 's3://ngi-igenomes/igenomes'
params.igenomes_ignore = false
params.multiqc_config = null
params.multiqc_title = null
params.max_multiqc_email_size = '25.MB'
params.outdir = './results'
params.tracedir = "${params.outdir}/pipeline_info"
params.email = null
params.email_on_fail = null
params.plaintext_email = false
params.monochrome_logs = false
params.help = false
params.validate_params = true
params.show_hidden_params = false
params.schema_ignore_params = 'genomes'
params.enable_conda = false
params.custom_config_version = 'master'
params.custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
params.config_profile_description = null
params.config_profile_contact = null
params.config_profile_url = null
params.config_profile_name = null
params.max_memory = '6.GB'
params.max_cpus = "4"
params.max_time = '240.h'
params.save_reference = true
params.save_alignment = true
params.sample_ploidy = "1"
params.coverage = "70"
params.rate = ''
params.gvcfs_filter = QD"
params.vcftools_filter = --min_GQ"
params.max_amb_samples = "10000000"
params.max_perc_amb_samples = "10"
params.publish_dir_mode = 'copy'
params.rapidnj = true
params.fasttree = true
params.iqtree = false
params.raxmlng = false
params.save_debug = false
params.mask = true
params.tmpdir = "$projectDir/tmp"
params.skip_samples = ""
params.skip_samples_file = null
params.skip_combined_analysis = false
params.skip_phylogeny = false
params.ref_dir = null
params.ref_masked_fasta = null
params.ref_fai = null
params.ref_bwa = null
params.ref_dict = null
params.TMPDIR = "$params.tmpdir"
params.params.enable_conda = true
params.params.genomes = "[:]"
params.genomes = [:]
params.genomes['GRCh37'] = [:]
params.genomes['GRCh37'].fasta = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['GRCh37'].bwa = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"
params.genomes['GRCh37'].bowtie2 = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/"
params.genomes['GRCh37'].star = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/"
params.genomes['GRCh37'].bismark = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Sequence/BismarkIndex/"
params.genomes['GRCh37'].gtf = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
params.genomes['GRCh37'].bed12 = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed"
params.genomes['GRCh37'].readme = "${params.igenomes_base}/Homo_sapiens/Ensembl/GRCh37/Annotation/README.txt"
params.genomes['GRCh37'].mito_name = "MT"
params.genomes['GRCh37'].macs_gsize = "2.7e9"
params.genomes['GRCh37'].blacklist = "${projectDir}/assets/blacklists/GRCh37-blacklist.bed"
params.genomes['GRCh38'] = [:]
params.genomes['GRCh38'].fasta = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['GRCh38'].bwa = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
params.genomes['GRCh38'].bowtie2 = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/"
params.genomes['GRCh38'].star = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/"
params.genomes['GRCh38'].bismark = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/BismarkIndex/"
params.genomes['GRCh38'].gtf = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
params.genomes['GRCh38'].bed12 = "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.bed"
params.genomes['GRCh38'].mito_name = "chrM"
params.genomes['GRCh38'].macs_gsize = "2.7e9"
params.genomes['GRCh38'].blacklist = "${projectDir}/assets/blacklists/hg38-blacklist.bed"
params.genomes['GRCm38'] = [:]
params.genomes['GRCm38'].fasta = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['GRCm38'].bwa = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa"
params.genomes['GRCm38'].bowtie2 = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/"
params.genomes['GRCm38'].star = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/STARIndex/"
params.genomes['GRCm38'].bismark = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Sequence/BismarkIndex/"
params.genomes['GRCm38'].gtf = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"
params.genomes['GRCm38'].bed12 = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.bed"
params.genomes['GRCm38'].readme = "${params.igenomes_base}/Mus_musculus/Ensembl/GRCm38/Annotation/README.txt"
params.genomes['GRCm38'].mito_name = "MT"
params.genomes['GRCm38'].macs_gsize = "1.87e9"
params.genomes['GRCm38'].blacklist = "${projectDir}/assets/blacklists/GRCm38-blacklist.bed"
params.genomes['TAIR10'] = [:]
params.genomes['TAIR10'].fasta = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['TAIR10'].bwa = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BWAIndex/genome.fa"
params.genomes['TAIR10'].bowtie2 = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/"
params.genomes['TAIR10'].star = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/STARIndex/"
params.genomes['TAIR10'].bismark = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BismarkIndex/"
params.genomes['TAIR10'].gtf = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.gtf"
params.genomes['TAIR10'].bed12 = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.bed"
params.genomes['TAIR10'].readme = "${params.igenomes_base}/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/README.txt"
params.genomes['TAIR10'].mito_name = "Mt"
params.genomes['EB2'] = [:]
params.genomes['EB2'].fasta = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['EB2'].bwa = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/BWAIndex/genome.fa"
params.genomes['EB2'].bowtie2 = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/Bowtie2Index/"
params.genomes['EB2'].star = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/STARIndex/"
params.genomes['EB2'].bismark = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Sequence/BismarkIndex/"
params.genomes['EB2'].gtf = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Annotation/Genes/genes.gtf"
params.genomes['EB2'].bed12 = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Annotation/Genes/genes.bed"
params.genomes['EB2'].readme = "${params.igenomes_base}/Bacillus_subtilis_168/Ensembl/EB2/Annotation/README.txt"
params.genomes['UMD3.1'] = [:]
params.genomes['UMD3.1'].fasta = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['UMD3.1'].bwa = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/BWAIndex/genome.fa"
params.genomes['UMD3.1'].bowtie2 = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/Bowtie2Index/"
params.genomes['UMD3.1'].star = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/STARIndex/"
params.genomes['UMD3.1'].bismark = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Sequence/BismarkIndex/"
params.genomes['UMD3.1'].gtf = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Annotation/Genes/genes.gtf"
params.genomes['UMD3.1'].bed12 = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Annotation/Genes/genes.bed"
params.genomes['UMD3.1'].readme = "${params.igenomes_base}/Bos_taurus/Ensembl/UMD3.1/Annotation/README.txt"
params.genomes['UMD3.1'].mito_name = "MT"
params.genomes['WBcel235'] = [:]
params.genomes['WBcel235'].fasta = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['WBcel235'].bwa = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BWAIndex/genome.fa"
params.genomes['WBcel235'].bowtie2 = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/Bowtie2Index/"
params.genomes['WBcel235'].star = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/STARIndex/"
params.genomes['WBcel235'].bismark = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BismarkIndex/"
params.genomes['WBcel235'].gtf = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
params.genomes['WBcel235'].bed12 = "${params.igenomes_base}/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.bed"
params.genomes['WBcel235'].mito_name = "MtDNA"
params.genomes['WBcel235'].macs_gsize = "9e7"
params.genomes['CanFam3.1'] = [:]
params.genomes['CanFam3.1'].fasta = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['CanFam3.1'].bwa = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/BWAIndex/genome.fa"
params.genomes['CanFam3.1'].bowtie2 = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/Bowtie2Index/"
params.genomes['CanFam3.1'].star = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/STARIndex/"
params.genomes['CanFam3.1'].bismark = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Sequence/BismarkIndex/"
params.genomes['CanFam3.1'].gtf = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Annotation/Genes/genes.gtf"
params.genomes['CanFam3.1'].bed12 = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Annotation/Genes/genes.bed"
params.genomes['CanFam3.1'].readme = "${params.igenomes_base}/Canis_familiaris/Ensembl/CanFam3.1/Annotation/README.txt"
params.genomes['CanFam3.1'].mito_name = "MT"
params.genomes['GRCz10'] = [:]
params.genomes['GRCz10'].fasta = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['GRCz10'].bwa = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/BWAIndex/genome.fa"
params.genomes['GRCz10'].bowtie2 = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/Bowtie2Index/"
params.genomes['GRCz10'].star = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/STARIndex/"
params.genomes['GRCz10'].bismark = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Sequence/BismarkIndex/"
params.genomes['GRCz10'].gtf = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Annotation/Genes/genes.gtf"
params.genomes['GRCz10'].bed12 = "${params.igenomes_base}/Danio_rerio/Ensembl/GRCz10/Annotation/Genes/genes.bed"
params.genomes['GRCz10'].mito_name = "MT"
params.genomes['BDGP6'] = [:]
params.genomes['BDGP6'].fasta = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['BDGP6'].bwa = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/BWAIndex/genome.fa"
params.genomes['BDGP6'].bowtie2 = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/Bowtie2Index/"
params.genomes['BDGP6'].star = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/STARIndex/"
params.genomes['BDGP6'].bismark = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Sequence/BismarkIndex/"
params.genomes['BDGP6'].gtf = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Annotation/Genes/genes.gtf"
params.genomes['BDGP6'].bed12 = "${params.igenomes_base}/Drosophila_melanogaster/Ensembl/BDGP6/Annotation/Genes/genes.bed"
params.genomes['BDGP6'].mito_name = "M"
params.genomes['BDGP6'].macs_gsize = "1.2e8"
params.genomes['EquCab2'] = [:]
params.genomes['EquCab2'].fasta = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['EquCab2'].bwa = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/BWAIndex/genome.fa"
params.genomes['EquCab2'].bowtie2 = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/Bowtie2Index/"
params.genomes['EquCab2'].star = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/STARIndex/"
params.genomes['EquCab2'].bismark = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Sequence/BismarkIndex/"
params.genomes['EquCab2'].gtf = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Annotation/Genes/genes.gtf"
params.genomes['EquCab2'].bed12 = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Annotation/Genes/genes.bed"
params.genomes['EquCab2'].readme = "${params.igenomes_base}/Equus_caballus/Ensembl/EquCab2/Annotation/README.txt"
params.genomes['EquCab2'].mito_name = "MT"
params.genomes['EB1'] = [:]
params.genomes['EB1'].fasta = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['EB1'].bwa = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/BWAIndex/genome.fa"
params.genomes['EB1'].bowtie2 = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/Bowtie2Index/"
params.genomes['EB1'].star = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/STARIndex/"
params.genomes['EB1'].bismark = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/BismarkIndex/"
params.genomes['EB1'].gtf = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/Genes/genes.gtf"
params.genomes['EB1'].bed12 = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/Genes/genes.bed"
params.genomes['EB1'].readme = "${params.igenomes_base}/Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/README.txt"
params.genomes['Galgal4'] = [:]
params.genomes['Galgal4'].fasta = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Galgal4'].bwa = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/BWAIndex/genome.fa"
params.genomes['Galgal4'].bowtie2 = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/Bowtie2Index/"
params.genomes['Galgal4'].star = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/STARIndex/"
params.genomes['Galgal4'].bismark = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Sequence/BismarkIndex/"
params.genomes['Galgal4'].gtf = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Annotation/Genes/genes.gtf"
params.genomes['Galgal4'].bed12 = "${params.igenomes_base}/Gallus_gallus/Ensembl/Galgal4/Annotation/Genes/genes.bed"
params.genomes['Galgal4'].mito_name = "MT"
params.genomes['Gm01'] = [:]
params.genomes['Gm01'].fasta = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Gm01'].bwa = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/BWAIndex/genome.fa"
params.genomes['Gm01'].bowtie2 = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/Bowtie2Index/"
params.genomes['Gm01'].star = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/STARIndex/"
params.genomes['Gm01'].bismark = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Sequence/BismarkIndex/"
params.genomes['Gm01'].gtf = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Annotation/Genes/genes.gtf"
params.genomes['Gm01'].bed12 = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Annotation/Genes/genes.bed"
params.genomes['Gm01'].readme = "${params.igenomes_base}/Glycine_max/Ensembl/Gm01/Annotation/README.txt"
params.genomes['Mmul_1'] = [:]
params.genomes['Mmul_1'].fasta = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Mmul_1'].bwa = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/BWAIndex/genome.fa"
params.genomes['Mmul_1'].bowtie2 = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/Bowtie2Index/"
params.genomes['Mmul_1'].star = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/STARIndex/"
params.genomes['Mmul_1'].bismark = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Sequence/BismarkIndex/"
params.genomes['Mmul_1'].gtf = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Annotation/Genes/genes.gtf"
params.genomes['Mmul_1'].bed12 = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Annotation/Genes/genes.bed"
params.genomes['Mmul_1'].readme = "${params.igenomes_base}/Macaca_mulatta/Ensembl/Mmul_1/Annotation/README.txt"
params.genomes['Mmul_1'].mito_name = "MT"
params.genomes['IRGSP-1.0'] = [:]
params.genomes['IRGSP-1.0'].fasta = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['IRGSP-1.0'].bwa = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/BWAIndex/genome.fa"
params.genomes['IRGSP-1.0'].bowtie2 = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/Bowtie2Index/"
params.genomes['IRGSP-1.0'].star = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/STARIndex/"
params.genomes['IRGSP-1.0'].bismark = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/BismarkIndex/"
params.genomes['IRGSP-1.0'].gtf = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Annotation/Genes/genes.gtf"
params.genomes['IRGSP-1.0'].bed12 = "${params.igenomes_base}/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Annotation/Genes/genes.bed"
params.genomes['IRGSP-1.0'].mito_name = "Mt"
params.genomes['CHIMP2.1.4'] = [:]
params.genomes['CHIMP2.1.4'].fasta = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['CHIMP2.1.4'].bwa = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/BWAIndex/genome.fa"
params.genomes['CHIMP2.1.4'].bowtie2 = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/Bowtie2Index/"
params.genomes['CHIMP2.1.4'].star = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/STARIndex/"
params.genomes['CHIMP2.1.4'].bismark = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/BismarkIndex/"
params.genomes['CHIMP2.1.4'].gtf = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/Genes/genes.gtf"
params.genomes['CHIMP2.1.4'].bed12 = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/Genes/genes.bed"
params.genomes['CHIMP2.1.4'].readme = "${params.igenomes_base}/Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/README.txt"
params.genomes['CHIMP2.1.4'].mito_name = "MT"
params.genomes['Rnor_5.0'] = [:]
params.genomes['Rnor_5.0'].fasta = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Rnor_5.0'].bwa = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/BWAIndex/genome.fa"
params.genomes['Rnor_5.0'].bowtie2 = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/Bowtie2Index/"
params.genomes['Rnor_5.0'].star = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/STARIndex/"
params.genomes['Rnor_5.0'].bismark = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/BismarkIndex/"
params.genomes['Rnor_5.0'].gtf = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_5.0/Annotation/Genes/genes.gtf"
params.genomes['Rnor_5.0'].bed12 = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_5.0/Annotation/Genes/genes.bed"
params.genomes['Rnor_5.0'].mito_name = "MT"
params.genomes['Rnor_6.0'] = [:]
params.genomes['Rnor_6.0'].fasta = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Rnor_6.0'].bwa = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/BWAIndex/genome.fa"
params.genomes['Rnor_6.0'].bowtie2 = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/Bowtie2Index/"
params.genomes['Rnor_6.0'].star = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/STARIndex/"
params.genomes['Rnor_6.0'].bismark = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/BismarkIndex/"
params.genomes['Rnor_6.0'].gtf = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.gtf"
params.genomes['Rnor_6.0'].bed12 = "${params.igenomes_base}/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.bed"
params.genomes['Rnor_6.0'].mito_name = "MT"
params.genomes['R64-1-1'] = [:]
params.genomes['R64-1-1'].fasta = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['R64-1-1'].bwa = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/BWAIndex/genome.fa"
params.genomes['R64-1-1'].bowtie2 = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/Bowtie2Index/"
params.genomes['R64-1-1'].star = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/STARIndex/"
params.genomes['R64-1-1'].bismark = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/BismarkIndex/"
params.genomes['R64-1-1'].gtf = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Annotation/Genes/genes.gtf"
params.genomes['R64-1-1'].bed12 = "${params.igenomes_base}/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Annotation/Genes/genes.bed"
params.genomes['R64-1-1'].mito_name = "MT"
params.genomes['R64-1-1'].macs_gsize = "1.2e7"
params.genomes['EF2'] = [:]
params.genomes['EF2'].fasta = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['EF2'].bwa = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/BWAIndex/genome.fa"
params.genomes['EF2'].bowtie2 = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/Bowtie2Index/"
params.genomes['EF2'].star = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/STARIndex/"
params.genomes['EF2'].bismark = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/BismarkIndex/"
params.genomes['EF2'].gtf = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Genes/genes.gtf"
params.genomes['EF2'].bed12 = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Genes/genes.bed"
params.genomes['EF2'].readme = "${params.igenomes_base}/Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/README.txt"
params.genomes['EF2'].mito_name = "MT"
params.genomes['EF2'].macs_gsize = "1.21e7"
params.genomes['Sbi1'] = [:]
params.genomes['Sbi1'].fasta = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Sbi1'].bwa = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/BWAIndex/genome.fa"
params.genomes['Sbi1'].bowtie2 = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/Bowtie2Index/"
params.genomes['Sbi1'].star = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/STARIndex/"
params.genomes['Sbi1'].bismark = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Sequence/BismarkIndex/"
params.genomes['Sbi1'].gtf = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Annotation/Genes/genes.gtf"
params.genomes['Sbi1'].bed12 = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Annotation/Genes/genes.bed"
params.genomes['Sbi1'].readme = "${params.igenomes_base}/Sorghum_bicolor/Ensembl/Sbi1/Annotation/README.txt"
params.genomes['Sscrofa10.2'] = [:]
params.genomes['Sscrofa10.2'].fasta = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['Sscrofa10.2'].bwa = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/BWAIndex/genome.fa"
params.genomes['Sscrofa10.2'].bowtie2 = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/Bowtie2Index/"
params.genomes['Sscrofa10.2'].star = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/STARIndex/"
params.genomes['Sscrofa10.2'].bismark = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/BismarkIndex/"
params.genomes['Sscrofa10.2'].gtf = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/Genes/genes.gtf"
params.genomes['Sscrofa10.2'].bed12 = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/Genes/genes.bed"
params.genomes['Sscrofa10.2'].readme = "${params.igenomes_base}/Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/README.txt"
params.genomes['Sscrofa10.2'].mito_name = "MT"
params.genomes['AGPv3'] = [:]
params.genomes['AGPv3'].fasta = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['AGPv3'].bwa = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/BWAIndex/genome.fa"
params.genomes['AGPv3'].bowtie2 = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/Bowtie2Index/"
params.genomes['AGPv3'].star = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/STARIndex/"
params.genomes['AGPv3'].bismark = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Sequence/BismarkIndex/"
params.genomes['AGPv3'].gtf = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Annotation/Genes/genes.gtf"
params.genomes['AGPv3'].bed12 = "${params.igenomes_base}/Zea_mays/Ensembl/AGPv3/Annotation/Genes/genes.bed"
params.genomes['AGPv3'].mito_name = "Mt"
params.genomes['hg38'] = [:]
params.genomes['hg38'].fasta = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['hg38'].bwa = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
params.genomes['hg38'].bowtie2 = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/"
params.genomes['hg38'].star = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/STARIndex/"
params.genomes['hg38'].bismark = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Sequence/BismarkIndex/"
params.genomes['hg38'].gtf = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
params.genomes['hg38'].bed12 = "${params.igenomes_base}/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.bed"
params.genomes['hg38'].mito_name = "chrM"
params.genomes['hg38'].macs_gsize = "2.7e9"
params.genomes['hg38'].blacklist = "${projectDir}/assets/blacklists/hg38-blacklist.bed"
params.genomes['hg19'] = [:]
params.genomes['hg19'].fasta = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['hg19'].bwa = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
params.genomes['hg19'].bowtie2 = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/"
params.genomes['hg19'].star = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/"
params.genomes['hg19'].bismark = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Sequence/BismarkIndex/"
params.genomes['hg19'].gtf = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
params.genomes['hg19'].bed12 = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed"
params.genomes['hg19'].readme = "${params.igenomes_base}/Homo_sapiens/UCSC/hg19/Annotation/README.txt"
params.genomes['hg19'].mito_name = "chrM"
params.genomes['hg19'].macs_gsize = "2.7e9"
params.genomes['hg19'].blacklist = "${projectDir}/assets/blacklists/hg19-blacklist.bed"
params.genomes['mm10'] = [:]
params.genomes['mm10'].fasta = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['mm10'].bwa = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa"
params.genomes['mm10'].bowtie2 = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/"
params.genomes['mm10'].star = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/STARIndex/"
params.genomes['mm10'].bismark = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Sequence/BismarkIndex/"
params.genomes['mm10'].gtf = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
params.genomes['mm10'].bed12 = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.bed"
params.genomes['mm10'].readme = "${params.igenomes_base}/Mus_musculus/UCSC/mm10/Annotation/README.txt"
params.genomes['mm10'].mito_name = "chrM"
params.genomes['mm10'].macs_gsize = "1.87e9"
params.genomes['mm10'].blacklist = "${projectDir}/assets/blacklists/mm10-blacklist.bed"
params.genomes['bosTau8'] = [:]
params.genomes['bosTau8'].fasta = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['bosTau8'].bwa = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/BWAIndex/genome.fa"
params.genomes['bosTau8'].bowtie2 = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/Bowtie2Index/"
params.genomes['bosTau8'].star = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/STARIndex/"
params.genomes['bosTau8'].bismark = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Sequence/BismarkIndex/"
params.genomes['bosTau8'].gtf = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.gtf"
params.genomes['bosTau8'].bed12 = "${params.igenomes_base}/Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.bed"
params.genomes['bosTau8'].mito_name = "chrM"
params.genomes['ce10'] = [:]
params.genomes['ce10'].fasta = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['ce10'].bwa = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/BWAIndex/genome.fa"
params.genomes['ce10'].bowtie2 = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/Bowtie2Index/"
params.genomes['ce10'].star = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/STARIndex/"
params.genomes['ce10'].bismark = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Sequence/BismarkIndex/"
params.genomes['ce10'].gtf = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Annotation/Genes/genes.gtf"
params.genomes['ce10'].bed12 = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Annotation/Genes/genes.bed"
params.genomes['ce10'].readme = "${params.igenomes_base}/Caenorhabditis_elegans/UCSC/ce10/Annotation/README.txt"
params.genomes['ce10'].mito_name = "chrM"
params.genomes['ce10'].macs_gsize = "9e7"
params.genomes['canFam3'] = [:]
params.genomes['canFam3'].fasta = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['canFam3'].bwa = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/BWAIndex/genome.fa"
params.genomes['canFam3'].bowtie2 = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/Bowtie2Index/"
params.genomes['canFam3'].star = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/STARIndex/"
params.genomes['canFam3'].bismark = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Sequence/BismarkIndex/"
params.genomes['canFam3'].gtf = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Annotation/Genes/genes.gtf"
params.genomes['canFam3'].bed12 = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Annotation/Genes/genes.bed"
params.genomes['canFam3'].readme = "${params.igenomes_base}/Canis_familiaris/UCSC/canFam3/Annotation/README.txt"
params.genomes['canFam3'].mito_name = "chrM"
params.genomes['danRer10'] = [:]
params.genomes['danRer10'].fasta = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['danRer10'].bwa = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/BWAIndex/genome.fa"
params.genomes['danRer10'].bowtie2 = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/Bowtie2Index/"
params.genomes['danRer10'].star = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/STARIndex/"
params.genomes['danRer10'].bismark = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Sequence/BismarkIndex/"
params.genomes['danRer10'].gtf = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Annotation/Genes/genes.gtf"
params.genomes['danRer10'].bed12 = "${params.igenomes_base}/Danio_rerio/UCSC/danRer10/Annotation/Genes/genes.bed"
params.genomes['danRer10'].mito_name = "chrM"
params.genomes['danRer10'].macs_gsize = "1.37e9"
params.genomes['dm6'] = [:]
params.genomes['dm6'].fasta = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['dm6'].bwa = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex/genome.fa"
params.genomes['dm6'].bowtie2 = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/"
params.genomes['dm6'].star = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/STARIndex/"
params.genomes['dm6'].bismark = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Sequence/BismarkIndex/"
params.genomes['dm6'].gtf = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf"
params.genomes['dm6'].bed12 = "${params.igenomes_base}/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.bed"
params.genomes['dm6'].mito_name = "chrM"
params.genomes['dm6'].macs_gsize = "1.2e8"
params.genomes['equCab2'] = [:]
params.genomes['equCab2'].fasta = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['equCab2'].bwa = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/BWAIndex/genome.fa"
params.genomes['equCab2'].bowtie2 = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/Bowtie2Index/"
params.genomes['equCab2'].star = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/STARIndex/"
params.genomes['equCab2'].bismark = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Sequence/BismarkIndex/"
params.genomes['equCab2'].gtf = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Annotation/Genes/genes.gtf"
params.genomes['equCab2'].bed12 = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Annotation/Genes/genes.bed"
params.genomes['equCab2'].readme = "${params.igenomes_base}/Equus_caballus/UCSC/equCab2/Annotation/README.txt"
params.genomes['equCab2'].mito_name = "chrM"
params.genomes['galGal4'] = [:]
params.genomes['galGal4'].fasta = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['galGal4'].bwa = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/BWAIndex/genome.fa"
params.genomes['galGal4'].bowtie2 = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/Bowtie2Index/"
params.genomes['galGal4'].star = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/STARIndex/"
params.genomes['galGal4'].bismark = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Sequence/BismarkIndex/"
params.genomes['galGal4'].gtf = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Annotation/Genes/genes.gtf"
params.genomes['galGal4'].bed12 = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Annotation/Genes/genes.bed"
params.genomes['galGal4'].readme = "${params.igenomes_base}/Gallus_gallus/UCSC/galGal4/Annotation/README.txt"
params.genomes['galGal4'].mito_name = "chrM"
params.genomes['panTro4'] = [:]
params.genomes['panTro4'].fasta = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['panTro4'].bwa = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/BWAIndex/genome.fa"
params.genomes['panTro4'].bowtie2 = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/Bowtie2Index/"
params.genomes['panTro4'].star = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/STARIndex/"
params.genomes['panTro4'].bismark = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Sequence/BismarkIndex/"
params.genomes['panTro4'].gtf = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Annotation/Genes/genes.gtf"
params.genomes['panTro4'].bed12 = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Annotation/Genes/genes.bed"
params.genomes['panTro4'].readme = "${params.igenomes_base}/Pan_troglodytes/UCSC/panTro4/Annotation/README.txt"
params.genomes['panTro4'].mito_name = "chrM"
params.genomes['rn6'] = [:]
params.genomes['rn6'].fasta = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['rn6'].bwa = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/BWAIndex/genome.fa"
params.genomes['rn6'].bowtie2 = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/Bowtie2Index/"
params.genomes['rn6'].star = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/STARIndex/"
params.genomes['rn6'].bismark = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Sequence/BismarkIndex/"
params.genomes['rn6'].gtf = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf"
params.genomes['rn6'].bed12 = "${params.igenomes_base}/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.bed"
params.genomes['rn6'].mito_name = "chrM"
params.genomes['sacCer3'] = [:]
params.genomes['sacCer3'].fasta = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['sacCer3'].bwa = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BWAIndex/genome.fa"
params.genomes['sacCer3'].bowtie2 = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/"
params.genomes['sacCer3'].star = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/STARIndex/"
params.genomes['sacCer3'].bismark = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BismarkIndex/"
params.genomes['sacCer3'].readme = "${params.igenomes_base}/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/README.txt"
params.genomes['sacCer3'].mito_name = "chrM"
params.genomes['sacCer3'].macs_gsize = "1.2e7"
params.genomes['susScr3'] = [:]
params.genomes['susScr3'].fasta = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/WholeGenomeFasta/genome.fa"
params.genomes['susScr3'].bwa = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/BWAIndex/genome.fa"
params.genomes['susScr3'].bowtie2 = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/Bowtie2Index/"
params.genomes['susScr3'].star = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/STARIndex/"
params.genomes['susScr3'].bismark = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Sequence/BismarkIndex/"
params.genomes['susScr3'].gtf = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Annotation/Genes/genes.gtf"
params.genomes['susScr3'].bed12 = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Annotation/Genes/genes.bed"
params.genomes['susScr3'].readme = "${params.igenomes_base}/Sus_scrofa/UCSC/susScr3/Annotation/README.txt"
params.genomes['susScr3'].mito_name = "chrM"
params.ext.skip_samples = "params.skip_samples"
params.ext.skip_samples_file = "params.skip_samples_file"
params.ext.args = "{"
/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/
WorkflowMain.initialise(workflow, params, log)
/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/
include { MYCOSNP } from './workflows/mycosnp'
//
// WORKFLOW: Run main nf-core/mycosnp analysis pipeline
//
workflow NFCORE_MYCOSNP {
    MYCOSNP ()
}
/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/
//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    NFCORE_MYCOSNP ()
}
/*
========================================================================================
    THE END
========================================================================================
*/
workflow.onComplete {
// copy intermediate files + directories
println("Getting intermediate files from ICA")
['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
// return trace files
println("Returning workflow run-metric reports from ICA")
['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null', '|', 'grep','"report"' ,'|','xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()
}
workflow.onError {
// copy intermediate files + directories
println("Getting intermediate files from ICA")
['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
// return trace files
println("Returning workflow run-metric reports from ICA")
['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null', '|', 'grep','"report"' ,'|','xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()
}
