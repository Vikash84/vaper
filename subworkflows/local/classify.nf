//
// Check input samplesheet and get read channels
//
include { FORMAT_REFS                         } from '../../modules/local/format_refs'
include { SHOVILL                             } from '../../modules/nf-core/shovill/main'
include { MINIMAP2_ALIGN                      } from '../../modules/nf-core/minimap2/align/main'
include { SOURMASH_SKETCH as SM_SKETCH_REF    } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_SKETCH as SM_SKETCH_SAMPLE } from '../../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_GATHER as SM_GATHER_SELECT } from '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_GATHER as SM_GATHER_SAMPLE } from '../../modules/nf-core/sourmash/gather/main'
include { SUMMARIZE_TAXA                      } from '../../modules/local/summarize_taxa'
include { SM2REFS                             } from '../../modules/local/sm2refs'
include { NCBI_DATASETS                       } from '../../modules/local/ncbi-datasets'


workflow CLASSIFY {
    take:
    ch_reads // channel: [ val(meta), path(reads) ]
    ch_refs  // channel: [ path(reference) ]

    main:

    ch_versions = Channel.empty()

    /* 
    =============================================================================================================================
        CREATE SOURMASH SKETCHES
    =============================================================================================================================
    */

    // MODULE: Sample sketches
    SM_SKETCH_SAMPLE (
        ch_reads.map{ meta, reads -> [ meta, reads.get(0) ] } // only uses forward read
    )

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

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: ACCURATE
    =============================================================================================================================
    */
    if (params.mode == "accurate"){
        ch_refs
            .map{ meta, refs -> refs }
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

        // MODULE: Map contigs to the references
        MINIMAP2_ALIGN (
            SHOVILL.out.contigs,
            FORMAT_REFS.out.refs.map{ refs -> [ "reference", refs ] }.first(),
            false,
            false,
            false
        )

        // Set ch_seq from reads (fast mode) to contigs (accurate mode)
        SHOVILL.out.contigs.set{ ch_seq }

        // Set reference list & reference composition summary
        MINIMAP2_ALIGN.out.paf.set{ ch_ref_list }
        FORMAT_REFS.out.refs_comp.set{ ch_refs_comp }

    }

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: FAST
    =============================================================================================================================
    */
    if (params.mode == "fast"){

        // MODULE: Reference sketches
        SM_SKETCH_REF (
            ch_refs
        )

        // MODULE: Run Sourmash gather against the reference pool using the forward reads
        SM_GATHER_SELECT (
            SM_SKETCH_SAMPLE.out.signatures,
            SM_SKETCH_REF.out.signatures.map{ meta, ref -> ref }.collect(),
            false,
            false,
            false,
            false
        )

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
        .set{ ch_taxa_sample }

    SUMMARIZE_TAXA(
        ch_taxa_sample,
        ch_refs_comp.first()
    )

    // Update reference list
    SUMMARIZE_TAXA
        .out
        .ref_list
        .splitCsv(header: false, elem: 1)
        .map{ meta, ref -> [ meta, file(ref.get(0)).baseName ] }
        .combine(ch_refs.map{ meta, ref -> [ meta.id, params.mode == "accurate" ? file(ref).baseName : meta.id, ref ] }, by: 1)
        .map{ index, meta, ref_id, ref -> [ meta, ref_id, ref ] }
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
    if (params.kitchen_sink){

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

        ch_sm_refs
            .map{ meta, ref_id, accession -> [ ref_id, meta ] }
            .combine(NCBI_DATASETS.out.assembly, by: 0)
            .map{ ref_id, meta, assembly -> [ meta, ref_id, assembly ] }
            .concat(ch_ref_list)
            .set{ ch_ref_list }

    }
        
    emit:
    ref_list   = ch_ref_list   // channel: [ val(sample_meta), val(ref_id), path(ref_path) ]
    sm_summary = ch_sm_summary // channel: [ val(meta), val(result) ]
    //versions = SHOVILL.out.versions // channel: [ versions.yml ]
}