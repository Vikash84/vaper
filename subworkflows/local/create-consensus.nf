//
// Check input samplesheet and get read channels
//

include { BWA_MEM as BWA_MEM               } from '../../modules/local/bwa_mem'
include { IVAR_CONSENSUS as IVAR_CONSENSUS } from '../../modules/local/ivar_consensus'


workflow CREATE_CONSENSUS {
    take:
    ref_list // channel: [meta, ref]
    refs_tar  // channel: [refs_tar]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run BWA MEM
    //
    BWA_MEM (
        ref_list,
        refs_tar
    )

    //
    // MODULE: Run Ivar
    //
    IVAR_CONSENSUS (
        BWA_MEM.out.bam
    )
}
