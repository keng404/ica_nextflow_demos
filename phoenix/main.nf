#!/usr/bin/env nextflow
/*
========================================================================================
    CDCgov/phoenix
========================================================================================
    Github : https://github.com/CDCgov/phoenix
    Slack  : https://staph-b-dev.slack.com/channels/phoenix-dev
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl = 2
/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/
WorkflowMain.initialise(workflow, params, log)
/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/
include { PHOENIX_EXTERNAL       } from './workflows/phoenix'
include { PHOENIX_EXQC           } from './workflows/cdc_phoenix'
include { SRA_PHOENIX            } from './workflows/sra_phoenix'
include { SCAFFOLD_EXTERNAL      } from './workflows/scaffolds'
//
// WORKFLOW: Run main cdcgov/phoenix analysis pipeline
//
workflow PHOENIX {
    main:
        PHOENIX_EXTERNAL ()
    emit:
        scaffolds        = PHOENIX_EXTERNAL.out.scaffolds
        trimmed_reads    = PHOENIX_EXTERNAL.out.trimmed_reads
        mlst             = PHOENIX_EXTERNAL.out.mlst
        amrfinder_report = PHOENIX_EXTERNAL.out.amrfinder_report
        gamma_ar         = PHOENIX_EXTERNAL.out.gamma_ar
}
//
// WORKFLOW: Run internal version of cdcgov/phoenix analysis pipeline that includes BUSCO, SRST2 and KRAKEN_ASMBLED
//
workflow CDC_PHOENIX {
    main:
        PHOENIX_EXQC ()
    emit:
        scaffolds        = PHOENIX_EXQC.out.scaffolds
        trimmed_reads    = PHOENIX_EXQC.out.trimmed_reads
        paired_trmd_json = PHOENIX_EXQC.out.paired_trmd_json
        amrfinder_report = PHOENIX_EXQC.out.amrfinder_report
        gamma_ar         = PHOENIX_EXQC.out.gamma_ar
        summary_report   = PHOENIX_EXQC.out.summary_report
}
/*
========================================================================================
    RUN SRA WORKFLOWS
========================================================================================
*/
//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow SRA_LIST {
    main:
        SRA_PHOENIX ()
    emit:
        scaffolds        = SRA_PHOENIX.out.scaffolds
        trimmed_reads    = SRA_PHOENIX.out.trimmed_reads
        amrfinder_report = SRA_PHOENIX.out.amrfinder_report
}
/*
========================================================================================
    RUN SCAFFOLD WORKFLOW
========================================================================================
*/
//
// WORKFLOW: Entry point to analyze scaffold file(s) and run everything after Spades
//
workflow SCAFFOLD_LIST {
    main: SCAFFOLD_EXTERNAL ()
}
/*
========================================================================================
    THE END
========================================================================================
*/
workflow{
    PHOENIX()
}
workflow.onError{ 
// copy intermediate files + directories
println("Getting intermediate files from ICA")
['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
// return trace files
println("Returning workflow run-metric reports from ICA")
['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null', '|', 'grep','"report"' ,'|','xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()
}
