#[ Date 11-17-2022 ]
As of ICAv2.10, it is recommended to bring your own nextflow.config when creating your nf pipeline and refer to other configs (if desired, includeConfig 'conf/my_ica.config') as a way to declare pod annotations [specific to ICA](https://help.ica.illumina.com/project/p-flow/f-pipelines#compute-types). The demos in this repo are still compatible, but the work needed to liftover existing nextflow pipelines is much lower. These demos were made to illustrate work needed to be done when users could not bring their own nextflow.config to ICA.

This repo contains example nextflow scripts and XML files that can be used to create nextflow pipelines in ICAv2.

Please contact keng@illumina.com of there are any questions, comments, suggestions.
