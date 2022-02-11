# How to use the CLI to launch nextflow pipelines

# prequisite commands before you launch your pipeline
- Download CLI (https://illumina.gitbook.io/ica/command-line-interface/cli-releasehistory)) and ensure it's executable ``` chmod u+x icav2 ```
- login ``` icav2 config ```
- enter project context ``` icav2 projects enter ``` or ``` icav2 projects list ```


# Main Concepts:
1) Each pipeline, data (file/folder), and project have pipeline identifiers. These identifiers consist of long alphanumeric strings and you'll frequently use them in the ICAv2 CLI + API
2) The CLI and API is project-centric in that data and pipelines are present in a specific project
3) Storage can be provided as an additional input. This determines the scratch space disk size used for the duration of your pipeline run
4) If creating nextflow pipelines or CWL pipelines in advanced/code mode, an XML file will need to be made when crafting pipelines to make it obvious to users what parameters/options are available to them at the time of run. Examples can be found [here](https://git.illumina.com/keng/ICAv2_demos/blob/main/flow/nextflow_demos/grandeur/grandeur.parameters.v2.xml) and [here](https://git.illumina.com/keng/ICAv2_demos/blob/main/flow/nextflow_demos/peak/peak.parameters.v2.xml)
 
# Example ICAv2 command lines of launching pipeline

```bash
icav2 projectpipelines start nextflow --pipeline-id 3902e2e5-2cdd-4d66-b111-7b5a46a96288 --user-reference my_cli_grandeur_pipeline_45 --parameters local_db_type:'nt' --parameters mash_options:'-v 0 -d 0.5' --parameters prokka_options:'--mincontiglen 500 --compliant --locustag locus_tag' --parameters cg_pipeline_options:'--qual_offset 33 --minLength 1' --parameters kleborate_options:'-all' --parameters seqsero2_options_fasta:'-t 4 -m k' --parameters seqsero2_options_fastq:'-t 2 -m a -b mem' --parameters seqyclean_options:'-minlen 25 -qual' --input fastq_files:fol.7fe68a5a614247be379108d9adf7b017 --input fasta_files:fol.4aa93c4bc4194f8e379708d9adf7b017 --input additional_configs:fol.c6ac9d91c7314786379a08d9adf7b017 --input blast_db:fol.6efad1ef932348a6378d08d9adf7b017 --input kraken2_db:fol.bf3c5b7a361b4611379008d9adf7b017 --storage-size 3fab13dd-46e7-4b54-bb34-b80a01a99379 -c ~/.icav2/isc-balt.config.yaml
```

## Main options in addition to specify data input and pipeline parameters are:
- pipeline-id : alphanumeric identifier of your pipeline, you can obtain this via an ``` icav2 projectpipelines list ```
- user-reference : name your pipeline run to your choosing
- storage-size : determines the scratch space of the disk your pipeline uses [small, medium, large] corresponding to [ 1.2 Tb, 2.4 Tb, 7.2 Tb ] disks
- technical-tags : in addition to user-reference can help you query + monitor specific pipeline runs

Alternatively you can specify ``` --input ``` as ``` --input input_name:data_id,input_name2:data_id2,...,input_nameN:data_idN ``` and ``` --parameters ``` as ``` parameter_name:data_id,parameter_name2:str_2,...,parameter_nameN:str_N ```


### API formatting of data inputs + parameters

```bash
######### Data inputs
[   {   'data_ids': ['fil.7c7f884331ef4c21c1e608d9c3e7d489'],
        'parameter_code': 'api_key_file'},
    {   'data_ids': ['fil.c60d79aca0ab4b179ae308d9e67f1b0e'],
        'parameter_code': 'input_json'}]
############## parameter inputs
[   {   'code': 'bsshruns2icav2tool__run_id',
        'value': 'friday_debug_transfer_testv2'},
    {'code': 'bsshruns2icav2tool__project_name', 'value': 'ken_test'}]
```
