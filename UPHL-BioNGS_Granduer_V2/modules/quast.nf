process quast {
  tag "${sample}"
  maxForks 10
  pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-medium'
  cpus 6
  errorStrategy 'ignore'

  publishDir = [ path: params.outdir, mode: 'copy' ]

  container 'staphb/quast:5.0.2'

  when:
  params.contig_processes =~ /quast/

  input:
  tuple val(sample), file(contigs)

  output:
  path "quast/${sample}"                                                , emit: files
  path "quast/${sample}_quast_report.tsv"     , optional: true          , emit: for_multiqc
  path "quast/${sample}/transposed_report.tsv", optional: true          , emit: collect
  tuple val(sample), env(gc)                                            , emit: gc
  tuple val(sample), env(num_contigs)                                   , emit: contigs
  tuple val(sample), env(n50)                                           , emit: nfifty
  tuple val(sample), env(length)                                        , emit: length
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    quast.py --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    quast.py !{params.quast_options} \
      !{contigs} \
      --output-dir quast/!{sample} \
      --threads !{task.cpus} \
      2>> $err_file | tee -a $log_file

    gc=$(grep "GC (" quast/!{sample}/report.txt | awk '{print $3}' )
    num_contigs=$(grep "contigs" quast/!{sample}/report.txt | grep -v "(" | awk '{print $3}' )
    n50=$(grep "N50" quast/!{sample}/report.txt | awk '{print $2}' )
    length=$(grep "Total length" quast/!{sample}/report.txt | grep -v "(" | awk '{print $3}' )
    if [ -z "$gc" ] ; then gc='NA' ; fi
    if [ -z "$num_contigs" ] ; then num_contigs='NA' ; fi
    if [ -z "$n50" ] ; then n50='NA' ; fi
    if [ -z "$length" ] ; then length='NA' ; fi

    if [ -f "quast/!{sample}/report.tsv" ] ; then cp quast/!{sample}/report.tsv quast/!{sample}_quast_report.tsv ; fi
  '''
}
