process prokka {
  tag "${sample}"
  label "maxcpus"
  maxForks 10
  pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'

  errorStrategy 'ignore'

  publishDir = [ path: params.outdir, mode: 'copy' ]

  container  'staphb/prokka:1.14.5'

  when:
  params.phylogenetic_processes =~ /prokka/

  input:
  tuple val(sample), file(contigs), val(genus), val(species)

  output:
  path "prokka/${sample}/*"                                            , emit: prokka_files
  path "prokka/${sample}/${sample}.txt"                                , emit: for_multiqc
  path "gff/${sample}.gff"                                             , emit: gffs
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}", emit: log

  shell:
  '''
    mkdir -p prokka gff logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    prokka -v >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    prokka !{params.prokka_options} \
      --cpu !{task.cpus} \
      --outdir prokka/!{sample} \
      --prefix !{sample} \
      --genus !{genus} \
      --species !{species} \
      --force !{contigs} 2>> $err_file | tee -a $log_file

    cp prokka/!{sample}/!{sample}.gff gff/!{sample}.gff
  '''
}
