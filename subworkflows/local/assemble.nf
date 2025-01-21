//
// Check input samplesheet and get read channels
//

include { BWA_MEM        } from '../../modules/local/bwa_mem'
include { IVAR_CONSENSUS } from '../../modules/local/ivar_consensus'
include { IRMA           } from '../../modules/local/irma'
include { CONDENSE       } from '../../modules/local/condense'
include { BAM_STATS      } from '../../modules/local/bam_stats'
include { MAPPED_FASTQ   } from '../../modules/local/get_mapped_fastq'
include { NEXTCLADE      } from '../../modules/local/nextclade'

workflow ASSEMBLE {
    take:
    ch_ref_list // channel: [meta, ref_id, ref_path, reads]

    main:

    ch_versions = Channel.empty()
    ch_bam      = Channel.empty()

    /* 
    =============================================================================================================================
        ASSEMBLY OPTION 1: IRMA
    =============================================================================================================================
    */
    
    if (params.cons_assembler == "irma"){
        // MODULE: Run IRMA
        IRMA (
            ch_ref_list.groupTuple(by:0).map{ meta, ref_ids, ref_paths, reads -> [ meta, ref_paths, reads.get(0) ]},
            file("${baseDir}/assets/IRMA_MODULE/", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(IRMA.out.versions)

        /*
        Troubleshooting
        */
        IRMA
            .out
            .bam
            .transpose()
            .map{ meta, bam -> [ meta, bam.getSimpleName(), bam ] }
            .set{ ch_bam }
        IRMA
            .out
            .consensus
            .transpose()
            .map{ meta, consensus -> [meta, consensus.getSimpleName().replace(meta.id+'_', ''), consensus] }
            .set{ ch_consensus }
    }
    
    /* 
    =============================================================================================================================
        ASSEMBLY OPTION 2: IVAR
    =============================================================================================================================
    */

    if (params.cons_assembler == "ivar"){
        // MODULE: Run BWA MEM
        BWA_MEM (
            ch_ref_list
        )
        ch_versions = ch_versions.mix(BWA_MEM.out.versions)
        BWA_MEM.out.bam.set{ ch_bam }

        // MODULE: Run Ivar
        IVAR_CONSENSUS (
            ch_bam
        )
        ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions)
        IVAR_CONSENSUS.out.consensus.set{ch_consensus}
    }

    /* 
    =============================================================================================================================
        GET MAPPING STATS & EXPORT MAPPED READS TO FASTQ
    =============================================================================================================================
    */
    // Module: Get BAM stats
    BAM_STATS (
        ch_bam
    )
    ch_versions = ch_versions.mix(BAM_STATS.out.versions)

    // MODULE: Convert mapped reads to FASTQ
    MAPPED_FASTQ (
        BAM_STATS.out.read_list.join(ch_ref_list, by: [0,1])
    )
    ch_versions = ch_versions.mix(MAPPED_FASTQ.out.versions)

    /* 
    =============================================================================================================================
        CONDENSE DUPLICATE ASSEMBLIES
    =============================================================================================================================
    */
    CONDENSE(
        ch_consensus
            .groupTuple(by: 0)
            .join(BAM_STATS.out.coverage.groupTuple(by: 0), by: 0)
            .map{ meta, ref_id1, assemblies, ref_id2, stats -> [ meta, assemblies, stats ] }

    )
    ch_versions = ch_versions.mix(CONDENSE.out.versions)
    CONDENSE
        .out
        .assembly
        .transpose()
        .map{ meta, assembly -> [ meta, assembly.getSimpleName().replaceAll("${meta.id}_", ''), assembly ] }
        .set{ ch_consensus }

    /* 
    =============================================================================================================================
        ASSEMBLY QC METRICS
    =============================================================================================================================
    */
    // MODULE: Run Nextclade
    NEXTCLADE (
        ch_consensus
            .combine(ch_ref_list.map{ meta, ref_id, ref_path, reads -> [ ref_path, ref_id ] }.unique(), by: 1)
    )
    ch_versions = ch_versions.mix(NEXTCLADE.out.versions)

    emit:
    bam_stats = BAM_STATS.out.stats // channel: [ val(meta), val(ref), path(stats)) ]
    nextclade = NEXTCLADE.out.tsv   // channel: [ val(meta), val(ref), path(stats)) ]
    consensus = ch_consensus        // channel: [ val(meta), val(ref_id), path(assembly) ]
    versions  = ch_versions
}
