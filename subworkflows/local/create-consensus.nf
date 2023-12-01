//
// Check input samplesheet and get read channels
//

include { BWA_MEM          } from '../../modules/local/bwa_mem'
include { SAMTOOLSTATS2TBL } from '../../modules/local/samtoolstats2tbl'
include { MAPPED_FASTQ     } from '../../modules/local/get_mapped_fastq'
include { IVAR_CONSENSUS   } from '../../modules/local/ivar_consensus'


workflow CREATE_CONSENSUS {
    take:
    ref_list // channel: [meta, ref, reads]
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
    // MODULE: Summarize Samtool stats
    //
    SAMTOOLSTATS2TBL (
        BWA_MEM.out.stats
    )

    //
    // MODULE: Convert mapped reads to FASTQ
    //
    MAPPED_FASTQ (
        BWA_MEM.out.read_list.map{meta, read_list -> [meta, read_list] }.join(ref_list, by: 0)
    )

    //
    // MODULE: Run Ivar
    //
    IVAR_CONSENSUS (
        BWA_MEM.out.bam
    )

    emit:
    samtoolstats2tbl = SAMTOOLSTATS2TBL.out.tbl   // channel: [ val(meta), path(stats_table)) ]
    assembly_stats = IVAR_CONSENSUS.out.stats
}
