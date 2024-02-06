//
// Check input samplesheet and get read channels
//

include { BWA_MEM          } from '../../modules/local/bwa_mem'
include { SAMTOOLSTATS2TBL } from '../../modules/local/samtoolstats2tbl'
include { MAPPED_FASTQ     } from '../../modules/local/get_mapped_fastq'
include { IVAR_CONSENSUS   } from '../../modules/local/ivar_consensus'
include { MAFFT            } from '../../modules/local/mafft'
include { SNPSITES         } from '../../modules/local/snp-sites'


workflow ASSEMBLE {
    take:
    ref_list // channel: [sample_meta, ref_id, ref_path, reads]

    main:

    ch_versions = Channel.empty()
    
    /* 
    =============================================================================================================================
        ALIGN FULL READ SET TO REFERENCES
    =============================================================================================================================
    */ 

    // MODULE: Run BWA MEM
    BWA_MEM (
        ref_list    
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // MODULE: Summarize Samtool stats
    SAMTOOLSTATS2TBL (
        BWA_MEM.out.stats
    )

    // MODULE: Convert mapped reads to FASTQ
    MAPPED_FASTQ (
        BWA_MEM.out.read_list.map{meta, ref, read_list -> [ meta, ref, read_list ] }.join(ref_list, by: [0,1])
    )

    /* 
    =============================================================================================================================
        CREATE CONSENSUS ASSEMBLY
    =============================================================================================================================
    */ 

    // MODULE: Run Ivar
    IVAR_CONSENSUS (
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())

    // MODULE: Align reference to consensus using MAFFT
    MAFFT (
        IVAR_CONSENSUS.out.consensus.join(ref_list.map{ meta, ref_id, ref_path, reads -> [ meta, ref_id, ref_path ] }, by: [0,1])
    )
    ch_versions = ch_versions.mix(MAFFT.out.versions.first())

    // MODULE: Get variant sites compared to reference
    SNPSITES (
        MAFFT.out.aln
    )
    ch_versions = ch_versions.mix(SNPSITES.out.versions.first())

    emit:
    samtoolstats2tbl = SAMTOOLSTATS2TBL.out.tbl // channel: [ val(meta), val(ref), path(stats)) ]
    assembly_stats   = IVAR_CONSENSUS.out.stats // channel: [ val(meta), val(ref), path(stats)) ]
    assembly_vcf     = SNPSITES.out.vcf         // channel: [ val(meta), val(ref), path(vcf)) ]
}
