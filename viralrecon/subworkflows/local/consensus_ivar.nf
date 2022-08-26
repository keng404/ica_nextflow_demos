params = [:]
process = [:]
process.ext = [:]
params.consensus_caller = 'ivar'
params.save_mpileup = false
def consensus_caller = 'ivar'
//
// Consensus calling with iVar and downstream processing QC
//
include { IVAR_CONSENSUS } from '../../modules/nf-core/modules/ivar/consensus/main'
include { CONSENSUS_QC   } from './consensus_qc'
workflow CONSENSUS_IVAR {
    take:
    bam           // channel: [ val(meta), [ bam ] ]
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    nextclade_db  // channel: /path/to/nextclade_db/
    main:
    ch_versions = Channel.empty()
    //
    // Call consensus sequence with iVar
    //
if (!params.skip_consensus && params.consensus_caller == 'ivar') {
    process.ext.args = '-t 0.75 -q 20 -m 10 -n N'
    process.ext.args2 = '--count-orphans --no-BAQ --max-depth 0 --min-BQ 0 -aa'
}
    IVAR_CONSENSUS (
        bam,
        fasta,
        params.save_mpileup
    )
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())
    //
    // Consensus sequence QC
    //
    CONSENSUS_QC (
        IVAR_CONSENSUS.out.fasta,
        fasta,
        gff,
        nextclade_db
    )
    ch_versions = ch_versions.mix(CONSENSUS_QC.out.versions.first())
    emit:
    consensus        = IVAR_CONSENSUS.out.fasta          // channel: [ val(meta), [ fasta ] ]
    consensus_qual   = IVAR_CONSENSUS.out.qual           // channel: [ val(meta), [ qual.txt ] ]
    quast_results    = CONSENSUS_QC.out.quast_results    // channel: [ val(meta), [ results ] ]
    quast_tsv        = CONSENSUS_QC.out.quast_tsv        // channel: [ val(meta), [ tsv ] ]
    pangolin_report  = CONSENSUS_QC.out.pangolin_report  // channel: [ val(meta), [ csv ] ]
    nextclade_report = CONSENSUS_QC.out.nextclade_report // channel: [ val(meta), [ csv ] ]
    bases_tsv        = CONSENSUS_QC.out.bases_tsv        // channel: [ val(meta), [ tsv ] ]
    bases_pdf        = CONSENSUS_QC.out.bases_pdf        // channel: [ val(meta), [ pdf ] ]
    versions         = ch_versions                       // channel: [ versions.yml ]
}
