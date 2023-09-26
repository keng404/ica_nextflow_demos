# Background
As of ICAv2.10, it is recommended to bring your own nextflow.config when creating your nf pipeline and refer to other configs (if desired, includeConfig 'conf/my_ica.config') as a way to declare pod annotations [specific to ICA](https://help.ica.illumina.com/project/p-flow/f-pipelines#compute-types). The demos in this repo are still compatible, but the work needed to liftover existing nextflow pipelines is much lower. These demos were made to illustrate work needed to be done when users could not bring their own nextflow.config to ICA.

It is now recommended to look at some lifted over pipelines [here](https://github.com/keng404/ica_nextflow_demos_v2).

For an updated method for lifting over existing Nextflow pipelines to ICA, take a look at [this](https://github.com/keng404/nextflow-to-icav2-config).

## Additonal recommendations are to: 
  1) consider adding cpu and memory declarations to your nextflow process/module or in your configuration file. Because ICA  monitors your nextflow containers, allowing for 1 CPU and 2 Gb of memory to not be used by your process may help ICA monitor your pipeline. You may have to decrease CPU and memory in increments of 1 CPU and 2 Gb of memory if you find that ICA times out (i.e. kubernetes timeout errors or timeout error given by ICA).
  2) If your process/module has input channels that can contain several items, consider using the ```maxForks``` declaration to limit the number of concurrent tasks that get run. As with the point above, you may tune this to allow for ICA to better monitor a pipeline run.

### Notes
One additional recommendation is to ensure whereever possible to have any process/module in your pipeline to have a docker container reference. As of now, ICA can only run nf processes that are docker-enabled/ can run inside of a container

This repo contains example nextflow scripts and XML files that can be used to create nextflow pipelines in ICAv2.

Please contact keng@illumina.com of there are any questions, comments, suggestions.
