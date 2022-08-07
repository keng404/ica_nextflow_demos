# UPHL-BioNGS_Granduer_V2
Source Git Repository: [https://github.com/UPHL-BioNGS/Grandeur](https://github.com/UPHL-BioNGS/Grandeur)

By: [@erinyoung](https://github.com/erinyoung)

This ICA pipeline was set up to run Granduer on Fastq reads only, and not with some of the other optionality that is avaiable with Granduer.
In order for this pipeline to run successfully you must provide the following files which are labeled with thier associated fields in the XML:

## READS
A direcotry containing paired or single end reads as fastq files.

## PROJECT_DIRS FROM EDITED GIT REPO
The directories found in this Github directory:
```
configs
data
modules
subworkflows
```

You must also include two databases as directories with these exact names:
```
blast_db_refseq
kraken2_db
```

You can get the kraken2_db using these commands:
``` 
mkdir kraken2_db
cd kraken2_db
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz
tar -zxvf minikraken2_v2_8GB_201904.tgz 
```

You can get the blast database using these commands:
```
mkdir blast_db_refseq
cd blast_db_refseq

###NOTE: This is for the complete blast database; I am working on finding some code for just the refseq database####
# get the taxdump file
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
# get the nt files
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz"
# decompress the nt files
for file in *tar.gz ; do tar -zxvf $file ; done
```



