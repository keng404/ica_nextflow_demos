# Steps run and scripts used for pipeline port
https://github.com/keng404/nextflow-to-icav2-config.

Documentation/ wiki will be updated in the URL above to reflect 'best practices'/tips in using these scripts.
These scripts will convert your nextflow pipeline, configurations, and create your pipeline in ICA.

It can also give you an ICA CLI template that can be modified to run your pipeline. This documentation will be updated to reflect this

# Concept 0: XML configuration
Fields in this XML file indicate input data/folders and pipeline settings that can be configured in your ICA pipeline.

Parse the files: ```nextflow_schema.json```, ```nextflow.config```, and other relevant configuration files (if your pipeline uses igenomes)

You can find a copy of the XML file for this pipeline [here](https://github.com/keng404/ica_nextflow_demos/blob/master/phoenix/phoenix.pipeline.xml)


# Concept 1: nextflow configuration

1) ```nextflow.config```
contains pipeline settings and references additional configuration files
	- Challenge is ICA parses this file to strip out certain parameters and assumes on the params pipeline settings are in this file
	In doing so it may ignore or not be able to handle other configurations in the original file
	- make sure that docker is enabled and conda is turned-off by default
	- To overcome this challenge, we parse the original nextflow.config file and migrate additional configurations to a ```conf/base.ica.config``` file
	if ```conf/base.ica.config``` does not exist these scripts create this file
2) ```base.config```
usually contains default and process-specific compute configuration (i.e. cpu, memory, docker image, error strategy/retries)
We add pod directives that tell ICA which compute instances to use.

If this file does not exist, these scripts create ths file

# Concept 2: update pipeline code

Your nextflow pipeline code may evolve and the current ICA experience only allows you to update files in our visual editor.
  - You will need to log into the ICA website, navigate to your ICA project, select your pipeline of interest, enter 'edit' mode, find your file(s) of interest and make your edits -- making sure to save as you go
  - Instead these scripts will load your pipeline scripts and underlying assets as files and folders in ICA.
	   - you can use `icav2 projectdata upload` to update files and folders of interest
	- All files and folders will be placed in the `workflow.launchDir`. This means all custom libraries and scripts will be found by ICA.


# Concept 3: minor updates to pipeline code so that it works for ICA

- Finding executable scripts bundled with your nextflow pipeline and ensuring that the path to these files is accurate so that your pipeline can find them
- Convert any underlying process in your pipeline from local executor to docker mode. This is because at run time, your pipeline scripts are run with a minimal ubuntu docker image. You are also non-root. So any underlying assumptions of what binaries/libraries are on your current compute server might not be present in this ubuntu image. So another docker image will be specified and any reference to local executor will be removed for the ICA-confvured version.

# Concept 4: cdc-gov/phoenix pipeline-specific modification(s)

```spades.nf : {PIPELINE_DIR}/modules/local/spades.nf```
  - Added ```spades_status``` and ```spades_fail_status```

```afterSpades.sh : {PIPELINE_DIR}/bin/afterSpades.sh```
  - Modified script to take in 1 argument. Previously this script was relying on an environmental variable to run

```nextflow.config : {PIPELINE_DIR}/nextflow.config```
  - Added ```docker.runOptions = '-m 6G --memory-swap 8G'``` to prevent pipeline from erroring out from ```KRAKEN2_TRIMD (originally {PIPELINE_DIR}/modules/local/kraken2.nf)````
without this option, this process requests more and more memory, creating an exit code 71 on ICA. Error message specificaly refers to out-of-memory error.

A copy of a running pipeline can be found here:
https://github.com/keng404/ica_nextflow_demos/tree/master/phoenix

You may need to repeat some of the steps under this concept in future pipeline versions. 
There also may be additional steps to complete your pipeline port in updated Phoenix pipeline.
