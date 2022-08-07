process bwa {
  tag "${sample}"
  label "medcpus"
  maxForks 10
  pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-medium'
  cpus 4
  errorStrategy 'ignore'

  publishDir = [ path: params.outdir, mode: 'copy' ]

  container  'staphb/bwa:0.7.17'

  input:
  tuple val(sample), file(reads), file(contig)

  output:
  tuple val(sample), file("bwa/${sample}.sam")             , emit: sam
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}", emit: log

  shell:
  '''
    mkdir -p bwa logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    bwa index !{contig} 2>> $err_file >> $log_file

    bwa mem !{params.bwa_options} \
      -t !{task.cpus} \
      !{contig} \
      !{reads} \
      2>> $err_file \
      > bwa/!{sample}.sam
  '''
}
