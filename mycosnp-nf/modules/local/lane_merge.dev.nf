process LANE_MERGE {
	publishDir  path: { "${params.outdir}/mergelane"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    tag "$meta.id"
    container 'library/ubuntu:20.04'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    
cpus 6
    
memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    maxForks 10
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("combined/*.fastq.gz", includeInputs: true)       , emit: reads
    script:
    numReads = reads.size()
    fileEnding = "fastq"
    if(reads[0].getName().endsWith(".gz"))
    {
        fileEnding = "fastq.gz"
    }
}
    
    """
    echo ${meta.id} $reads
    mkdir combined
    if [[ $fileEnding == "fastq" ]]; then
        if [[ $numReads == 1 ]]; then
            gzip -c ${reads[0]} > combined/${meta.id}.fastq.gz
        elif [[ $numReads == 2 ]]; then
            gzip -c ${reads[0]} > combined/${meta.id}_R1.fastq.gz
            gzip -c ${reads[1]} > combined/${meta.id}_R2.fastq.gz
        elif [[ $numReads == 4 ]]; then
            gzip -c ${reads[0]} ${reads[2]} > combined/${meta.id}_R1.fastq.gz
            gzip -c ${reads[1]} ${reads[3]} > combined/${meta.id}_R2.fastq.gz
        elif [[ $numReads == 6 ]]; then
            gzip -c ${reads[0]} ${reads[2]} ${reads[4]} > combined/${meta.id}_R1.fastq.gz
            gzip -c ${reads[1]} ${reads[3]} ${reads[5]} > combined/${meta.id}_R2.fastq.gz
        fi
    else
        if [[ $numReads == 1 ]]; then
            cd combined; ln -s ../${reads[0]} ${meta.id}.$fileEnding; cd ..
        elif [[ $numReads == 2 ]]; then
            cd combined; ln -s ../${reads[0]} ${meta.id}_R1.$fileEnding; cd ..
            cd combined; ln -s ../${reads[1]} ${meta.id}_R2.$fileEnding; cd ..
        elif [[ $numReads == 4 ]]; then
            cat ${reads[0]} ${reads[2]} > combined/${meta.id}_R1.$fileEnding
            cat ${reads[1]} ${reads[3]} > combined/${meta.id}_R2.$fileEnding
        elif [[ $numReads == 6 ]]; then
            cat ${reads[0]} ${reads[2]} ${reads[4]} > combined/${meta.id}_R1.$fileEnding
            cat ${reads[1]} ${reads[3]} ${reads[5]} > combined/${meta.id}_R2.$fileEnding
        fi
    fi
    """
}
