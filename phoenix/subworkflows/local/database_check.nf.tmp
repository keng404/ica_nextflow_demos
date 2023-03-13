/*
========================================================================================
    Processes
========================================================================================
*/
process database_check {
    db_ch = Channel.fromPath(${params.databases}, checkIfExists: true )
    input:
    path(db_path)
    output:
    script:
    """
        database_checker.sh ${db_path}
    """
}
workflow database_check {
  database_check(params.databases)
}
workflow.onError{ 
// copy intermediate files + directories
println("Getting intermediate files from ICA")
['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
// return trace files
println("Returning workflow run-metric reports from ICA")
['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null', '|', 'grep','"report"' ,'|','xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()
}
