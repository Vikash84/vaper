//
// Check input samplesheet and get read channels
//

include { BWA_MEM          } from '../../modules/local/bwa_mem'
include { SAMTOOLSTATS2TBL } from '../../modules/local/samtoolstats2tbl'
include { MAPPED_FASTQ     } from '../../modules/local/get_mapped_fastq'
include { IVAR_CONSENSUS   } from '../../modules/local/ivar_consensus'
include { MAFFT            } from '../../modules/local/mafft'
include { SNPSITES         } from '../../modules/local/snp-sites'


workflow CREATE_CONSENSUS {
    take:
    ref_list // channel: [meta, ref, reads]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run BWA MEM
    //
    BWA_MEM (
        ref_list,
        params.refs
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
        BWA_MEM.out.read_list.map{meta, ref, read_list -> [meta, ref, read_list] }.join(ref_list, by: [0,1])
    )

    //
    // MODULE: Run Ivar
    //
    IVAR_CONSENSUS (
        BWA_MEM.out.bam
    )

    //
    // MODULE: Align reference to consensus using MAFFT
    //
    MAFFT (
        IVAR_CONSENSUS.out.consensus,
        params.refs
    )

    //
    // MODULE: Get variant sites compared to reference
    //
    SNPSITES (
        MAFFT.out.aln
    )

    emit:
    samtoolstats2tbl = SAMTOOLSTATS2TBL.out.tbl // channel: [ val(meta), val(ref), path(stats)) ]
    assembly_stats = IVAR_CONSENSUS.out.stats   // channel: [ val(meta), val(ref), path(stats)) ]
}
