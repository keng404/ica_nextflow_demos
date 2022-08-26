process MULTIQC_TSV_FROM_LIST {
	publishDir  path: { "${params.outdir}/listfromtsvmultiqc"}, mode: "copy", saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    executor 'local'
    memory 100.MB
    container library/ubuntu:20.04
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-small'
    
cpus 6
    
memory '48 GB'
    errorStrategy 'ignore'
    time '1day'
    input:
    val tsv_data   // [ ['foo', 1], ['bar', 1] ]
    val header     // [ 'name', 'number' ]
    val out_prefix
    output:
    path "*.tsv"
    when:
    task.ext.when == null || task.ext.when
    exec:
    // Generate file contents
    def contents = ""
    if (tsv_data.size() > 0) {
        contents += "${header.join('\t')}\n"
        contents += tsv_data.join('\n')
    }
    // Write to file
    def mqc_file = task.workDir.resolve("${out_prefix}_mqc.tsv")
    mqc_file.text = contents
}
