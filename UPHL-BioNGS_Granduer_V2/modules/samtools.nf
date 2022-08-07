process sort {
  tag "${sample}"
  label "medcpus"
  maxForks 10
  pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-medium'
  cpus 6
  errorStrategy 'ignore'

  publishDir = [ path: params.outdir, mode: 'copy' ]

  container  'staphb/samtools:1.15'

  input:
  tuple val(sample), file(sam)

  output:
  tuple val(sample), file("aligned/${sample}.sorted.bam*")                   , emit: bam
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p aligned logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    samtools --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    samtools sort !{params.samtools_sort_options} \
      --threads !{task.cpus} \
      !{sam} \
      -o aligned/!{sample}.sorted.bam \
      --write-index \
      2>> $err_file >> $log_file
  '''
}
