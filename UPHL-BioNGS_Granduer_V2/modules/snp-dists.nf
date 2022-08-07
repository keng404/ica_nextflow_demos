process snp_dists {
  tag "SNP matrix"
  maxForks 10
  pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-medium'

  errorStrategy 'ignore'

  publishDir = [ path: params.outdir, mode: 'copy' ]

  container  'staphb/snp-dists:0.8.2'

  when:
  params.phylogenetic_processes =~ /snpdists/

  input:
  file(contigs)

  output:
  path "snp-dists/snp_matrix.txt"                                            , emit: snp_matrix
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}", emit: log_files

  shell:
  '''
    mkdir -p snp-dists logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    snp-dists -v >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    snp-dists !{params.snp_dists_options} \
      !{contigs} \
      2>> $err_file \
      > snp-dists/snp_matrix.txt
  '''
}
