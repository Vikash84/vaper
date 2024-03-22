//
// Check input samplesheet and get read channels
//

include { BWA_MEM          } from '../../modules/local/bwa_mem'
include { SAMTOOLSTATS2TBL } from '../../modules/local/samtoolstats2tbl'
include { MAPPED_FASTQ     } from '../../modules/local/get_mapped_fastq'
include { IVAR_CONSENSUS   } from '../../modules/local/ivar_consensus'
include { NEXTCLADE_CONFIG } from '../../modules/local/nexclade/config/main'
include { NEXTCLADE_RUN    } from '../../modules/local/nexclade/run/main'


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

    /* 
    =============================================================================================================================
        ASSEMBLY QC METRICS
    =============================================================================================================================
    */
    // MODULE: Configure Nextclade QC file
    NEXTCLADE_CONFIG (
        ref_list.map{ meta, ref_id, ref_path, reads -> [ ref_id, ref_path ] }.unique(),
        file("${baseDir}/assets/nextclade-template.json", checkIfExists: true)
    )
    // MODULE: Run Nextclade
    NEXTCLADE_RUN (
        IVAR_CONSENSUS
            .out
            .consensus
            .combine(NEXTCLADE_CONFIG.out.config.map{ ref_id, ref, config -> [ ref, ref_id, config ] }, by: 1)
            .map{ ref_id, meta, consensus, ref, config -> [ meta, ref_id, consensus, ref, config ] }
    )

    emit:
    samtoolstats2tbl = SAMTOOLSTATS2TBL.out.tbl // channel: [ val(meta), val(ref), path(stats)) ]
    nextflow         = NEXTCLADE_RUN.out.tsv    // channel: [ val(meta), val(ref), path(stats)) ]
}
