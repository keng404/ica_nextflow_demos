process gunzip_blast_db {
  maxForks 10
  pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-medium'
  cpus 3
  memory '5 GB'

  input:
  val(blastdb)
  val(output_prefix)

  script:
  """
    gunzip -c ${workflow.launchDir/blastdb} > ${workflow.launchDir/output_prefix}
  """
}