process DOWNSAMPLE_RATE {
    tag "$meta.id"
    container null
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-small'
    errorStrategy 'ignore'
    time '1day'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    errorStrategy 'ignore'
    time '1day'
    input:
    tuple val(meta), path(reads)
path(reference_fasta)
val(coverage)
    output:
    tuple val(meta), env(SAMPLE_RATE),  env(SAMPLED_NUM_READS),   emit: downsampled_rate
env(SAMPLED_NUM_READS)                                    ,   emit: number_to_sample
script:
"""
REFERENCE_LEN=\$(awk '!/^>/ {len+=length(\$0)} END {print len}' < ${reference_fasta})
READS_LEN=\$(zcat ${reads} | awk '/^@/ {getline; len+=length(\$0)} END {print len}')
SAMPLE_RATE=\$(echo "${coverage} \${READS_LEN} \${REFERENCE_LEN}" | awk '{x=\$1/(\$2/\$3); x=(1<x?1:x)} END {print x}')
# Calculate number of reads
NUM_READS=\$(zcat ${reads[0]} | awk '/^@/ {lines+=1} END {print lines}')
SAMPLED_NUM_READS=\$(echo "\${NUM_READS} \${SAMPLE_RATE}" | awk '{x=\$1*\$2} END {printf "%.0f", x}')
"""
}
