//
// Check input samplesheet and get read channels
//
include { FORMAT_REFS                           } from '../../modules/local/format_refs'
include { SHOVILL                               } from '../../modules/nf-core/shovill/main'
include { MINIMAP2_ALIGN                        } from '../../modules/nf-core/minimap2/align/main'
include { SM_SKETCH_REF                         } from '../../modules/local/sourmash_sketch_ref'
include { SOURMASH_SKETCH as SM_SKETCH_SAMPLE   } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_GATHER as SM_GATHER_SELECT   } from '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_GATHER as SM_GATHER_SAMPLE   } from '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_METAGENOME as SM_META_SAMPLE } from '../../modules/nf-core/sourmash/metagenome/main'
include { SUMMARIZE_TAXA                        } from '../../modules/local/summarize_taxa'
include { SM2REFS                               } from '../../modules/local/sm2refs'
include { NCBI_DATASETS                         } from '../../modules/local/ncbi-datasets'


workflow CLASSIFY {
    take:
    ch_reads // channel: [ val(meta), path(reads), path/val(reference ) ]
    ch_refs  // channel: [ val(meta), path(reference) ]

    main:

    ch_versions = Channel.empty()
    
    ch_man_refs = ch_reads.map{ meta, reads, ref -> [ meta, ref] }
    ch_reads    = ch_reads.map{ meta, reads, ref -> [ meta, reads ] }

    /* 
    =============================================================================================================================
        CREATE SOURMASH SKETCHES
    =============================================================================================================================
    */

    // MODULE: Sample sketches
    SM_SKETCH_SAMPLE (
        ch_reads.map{ meta, reads -> [ meta, reads.get(0) ] } // only uses forward read
    )
    ch_versions = ch_versions.mix(SM_SKETCH_SAMPLE.out.versions)

    /* 
    =============================================================================================================================
        DETERMINE DOMINANT VIRAL TAXA
    =============================================================================================================================
    */

    // MODULE: Classify viral species in each sample using Sourmash
    SM_GATHER_SAMPLE (
        SM_SKETCH_SAMPLE.out.signatures,
        params.sm_db,
        false,
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(SM_GATHER_SAMPLE.out.versions)

    // MODULE: Summarize taxa from sourmash gather
    SM_META_SAMPLE (
        SM_GATHER_SAMPLE.out.result,
        file(params.sm_taxa, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(SM_META_SAMPLE.out.versions)

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: ACCURATE
    =============================================================================================================================
    */
    if (params.ref_mode == "accurate"){
        ch_refs
            .map{ meta, assembly -> assembly }
            .collect()
            .map{ assembly -> [ assembly ] }
            .set{ ch_ref_assemblies }
        // MODULE: Prepare references
        FORMAT_REFS (
            ch_ref_assemblies
                .flatten()
                .map{ assembly -> file(assembly).getName() }
                .collectFile(name: 'ref-list.txt',newLine: true)
                .combine(ch_ref_assemblies)
        )

        // MODULE: Run Shovill
        SHOVILL (
            ch_reads
        )
        ch_versions = ch_versions.mix(SHOVILL.out.versions)

        // MODULE: Map contigs to the references
        MINIMAP2_ALIGN (
            SHOVILL.out.contigs,
            FORMAT_REFS.out.refs.map{ refs -> [ "reference", refs ] },
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        // Set reference list & reference composition summary
        MINIMAP2_ALIGN.out.paf.set{ ch_ref_list }
        FORMAT_REFS.out.refs_comp.set{ ch_refs_comp }

    }

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: FAST
    =============================================================================================================================
    */
    if (params.ref_mode == "fast"){

        // MODULE: Reference sketches
        SM_SKETCH_REF (
            ch_refs.map{meta, assembly -> assembly }.collect()
        )
        ch_versions = ch_versions.mix(SM_SKETCH_REF.out.versions)

        // MODULE: Run Sourmash gather against the reference pool using the forward reads
        SM_GATHER_SELECT (
            SM_SKETCH_SAMPLE.out.signatures,
            SM_SKETCH_REF.out.signatures,
            false,
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(SM_GATHER_SELECT.out.versions)

        // Set reference list & empty reference composition summary
        SM_GATHER_SELECT.out.result.set{ ch_ref_list }
        ch_refs_comp = []

    }

    /* 
    =============================================================================================================================
        SUMMARIZE TAXONOMY
    =============================================================================================================================
    */
    // combine reference mapping, sample taxonomy, and reference taxonomy
    ch_ref_list
        .join(SM_GATHER_SAMPLE.out.result)
        .join(SM_META_SAMPLE.out.result)
        .set{ ch_taxa_sample }

    SUMMARIZE_TAXA(
        ch_taxa_sample,
        ch_refs_comp
    )
    ch_versions = ch_versions.mix(SUMMARIZE_TAXA.out.versions)

    // Update reference list
    SUMMARIZE_TAXA
        .out
        .ref_list
        .splitCsv(header: false, elem: 1)
        .transpose()
        .map{ meta, ref -> [ file(ref).getSimpleName(), meta ] }
        .combine(ch_refs.map{ ref_meta, assembly -> [ ref_meta.id, assembly ] }, by: 0)
        .map{ ref_id, meta, ref -> [ meta, ref_id, ref ] }
        .set{ ch_ref_list }

    // Create Sourmash summary channel
    SUMMARIZE_TAXA
        .out
        .sm_summary
        .set{ ch_sm_summary }

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: Kitchen Sink
        - Downloads NCBI assembly for all Sourmash hits and adds them to the reference list
    =============================================================================================================================
    */
    if (params.ref_kitchensink){

        SM2REFS (
        ch_sm_summary
        )

        SM2REFS
            .out
            .refs
            .splitCsv(header: true, elem: 1)
            .map{ meta, result -> [ meta, result.taxa, result.accession ] }
            .set{ ch_sm_refs }

        NCBI_DATASETS (
            ch_sm_refs.map{ meta, ref_id, accession -> [ ref_id, accession ] }.unique()
        )
        ch_versions = ch_versions.mix(NCBI_DATASETS.out.versions)

        ch_sm_refs
            .map{ meta, ref_id, accession -> [ ref_id, meta ] }
            .combine(NCBI_DATASETS.out.assembly, by: 0)
            .map{ ref_id, meta, assembly -> [ meta, assembly ] }
            .concat(ch_ref_list)
            .set{ ch_ref_list }

    }

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: Manual
        - Adds list of references provided in the samplesheet
    =============================================================================================================================
    */
    ch_man_refs
        .filter{ meta, ref -> ref }
        .map{ meta, ref -> [ meta, ref.tokenize(';') ] }
        .transpose()
        .map{ meta, ref -> [ meta, file(ref).exists() ? file(ref).getSimpleName() : ref, file(ref).exists() ? file(ref) : null ] }
        .set{ ch_man_refs }
    ch_man_refs
        .filter{ meta, ref_id, ref -> ! ref }
        .map{ meta, ref_id, ref -> [ ref_id, meta ] }
        .combine(ch_refs.map{ meta, assembly -> [ meta.id, assembly ] }, by: 0)
        .map{ ref_id, meta, ref -> [ meta, ref_id, ref ] }
        .concat(ch_man_refs.filter{ meta, ref_id, ref -> ref })
        .concat( ch_ref_list )
        .unique()
        .set{ ch_ref_list }

    emit:
    ref_list   = ch_ref_list   // channel: [ val(sample_meta), val(ref_id), path(ref_path) ]
    sm_summary = ch_sm_summary // channel: [ val(meta), val(result) ]
    versions   = ch_versions   // channel: [ versions.yml ]
}