// Check mandatory parameters
if (params.reads) { raw_reads = Channel.fromPath(params.reads) } else { exit 1, 'Please move your FASTQ files to the "FASTQs" folder!' }
// Add modules needed
include { BBMAP_BBDUK } from '../modules/nf-core/modules/bbmap/bbduk/main'
include { FASTP } from '../modules/nf-core/modules/fastp/main'
workflow RAW_TO_TRIM {
    
    SEQKIT_PAIR ( raw_reads )
    BBMAP_BBDUK ( SEQKIT_PAIR.out.reads, phiX)
    FASTP ( BBMAP_BBDUK.out.reads, true, true )
}
workflow.onError{ 
// copy intermediate files + directories
println("Getting intermediate files from ICA")
['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
// return trace files
println("Returning workflow run-metric reports from ICA")
['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null', '|', 'grep','"report"' ,'|','xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()
}
