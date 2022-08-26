process DOWNSAMPLE_RATE {
	publishDir  path: { "${params.outdir}/ratedownsample"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    container 'library/ubuntu:20.04'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-small'
    
cpus 1
    
memory '6 GB'
    errorStrategy 'ignore'
    time '1day'
    maxForks 10
    input:
    tuple val(meta), path(reads)
path(reference_fasta)
val(coverage)
val(rate)
    output:
    tuple val(meta), env(SAMPLE_RATE),  env(SAMPLED_NUM_READS),   emit: downsampled_rate
env(SAMPLED_NUM_READS)                                    ,   emit: number_to_sample
script:
"""
REFERENCE_LEN=\$(awk '!/^>/ {len+=length(\$0)} END {print len}' < ${reference_fasta})
READS_LEN=\$(zcat ${reads} | awk '/^@/ {getline; len+=length(\$0)} END {print len}')
if [ ${rate} == 1 ];
then SAMPLE_RATE=1
else
SAMPLE_RATE=\$(echo "${coverage} \${READS_LEN} \${REFERENCE_LEN}" | awk '{x=\$1/(\$2/\$3); x=(1<x?1:x)} END {print x}')
fi
# Calculate number of reads
NUM_READS=\$(zcat ${reads[0]} | awk '/^@/ {lines+=1} END {print lines}')
SAMPLED_NUM_READS=\$(echo "\${NUM_READS} \${SAMPLE_RATE}" | awk '{x=\$1*\$2} END {printf "%.0f", x}')
"""
}
