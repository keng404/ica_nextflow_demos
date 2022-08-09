params.file = 'test.txt'
params.outdir = './out/fastq'
def sra_ids = []
new File(params.file).withReader { reader ->
   while (line = reader.readLine()) {
        sra_ids << line
   }
}

process downloadSRA{
    container 'keng404/sra-tools:0.0.3'
    publishDir params.outdir, mode: 'copy'
    pod annotation: 'scheduler.illumina.com/presetSize' , value: 'standard-large'
    cpus 7
    memory '30 GB'
    maxForks 10
    time '2h'
    input:
        val(ids) from sra_ids

    output:
        file("*fastq.gz") into reads
    
    script:
    """
        fasterq-dump ${ids}
        pigz *fastq
    """
}


workflow.onError {
// copy intermediate files + directories
println("Getting intermediate files from ICA")
['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
// return trace files
println("Returning workflow run-metric reports from ICA")
['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null', '|', 'grep','"report"' ,'|','xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()
}
