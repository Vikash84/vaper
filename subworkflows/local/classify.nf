//
// Check input samplesheet and get read channels
//
include { FORMAT_REFS         } from '../../modules/local/format_refs'
include { SHOVILL             } from '../../modules/nf-core/shovill/main'
include { MINIMAP2_ALIGN      } from '../../modules/local/minimap2_align'
include { SM_SKETCH_REF       } from '../../modules/local/sourmash_sketch_ref'
include { SOURMASH_SKETCH     } from '../../modules/nf-core/sourmash/sketch/main'
include { SM_GATHER_SELECT    } from '../../modules/local/sourmash_select'
include { SOURMASH_GATHER     } from '../../modules/nf-core/sourmash/gather/main'
include { SOURMASH_METAGENOME } from '../../modules/local/sourmash_metagenome'
include { KRONA_KTIMPORTTEXT  } from '../../modules/nf-core/krona/ktimporttext/main'
include { REF_SELECT_ACCURATE } from '../../modules/local/ref-select_accurate'
include { SUMMARIZE_TAXA      } from '../../modules/local/summarize_taxa'
include { TAR2REFS            } from '../../modules/local/tar2refs'
include { SM2REFS             } from '../../modules/local/sm2refs'
include { NCBI_DATASETS       } from '../../modules/local/ncbi-datasets'
include { FINALIZE_REFS       } from '../../modules/local/finalize-refs'


workflow CLASSIFY {
    take:
    ch_reads     // channel: [ val(meta), path(reads) ]
    ch_refs      // channel: [ path(refs.fa.gz), path(refs-comp.txt.gz), path(refs.tar.gz)  ]
    ch_refs_man  // channel: [ val(meta), path/val(ref)  ]

    main:

    ch_versions  = Channel.empty()
    ch_ref_list  = Channel.empty()

    // Parse reference versions
    ch_refs.map{ refs, comp, tar, sheet -> refs  }.set{ ch_refs_fmt  }
    ch_refs.map{ refs, comp, tar, sheet -> comp  }.set{ ch_refs_comp }
    ch_refs.map{ refs, comp, tar, sheet -> tar   }.set{ ch_refs_tar  }
    ch_refs.map{ refs, comp, tar, sheet -> sheet }.set{ ch_refsheet }

    /* 
    =============================================================================================================================
        CREATE SOURMASH SKETCHES
    =============================================================================================================================
    */
    // MODULE: Sample sketches
    SOURMASH_SKETCH (
        ch_reads.map{ meta, reads -> [ meta, reads.get(0) ] } // only uses forward read
    )
    ch_versions = ch_versions.mix(SOURMASH_SKETCH.out.versions)

    /* 
    =============================================================================================================================
        DETERMINE DOMINANT VIRAL TAXA
    =============================================================================================================================
    */
    // MODULE: Classify viral species in each sample using Sourmash
    SOURMASH_GATHER (
        SOURMASH_SKETCH.out.signatures,
        params.sm_db,
        false,
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(SOURMASH_GATHER.out.versions)

    // MODULE: Summarize taxa from sourmash gather
    SOURMASH_METAGENOME (
        SOURMASH_GATHER.out.result,
        file(params.sm_taxa, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(SOURMASH_METAGENOME.out.versions)

    // MODULE: Create Krona plot from Sourmash metagenome output
    KRONA_KTIMPORTTEXT (
        SOURMASH_METAGENOME.out.krona
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)

    // Combine reference mapping, sample taxonomy, and reference taxonomy
    SUMMARIZE_TAXA(
        SOURMASH_METAGENOME.out.result
    )
    ch_versions = ch_versions.mix(SUMMARIZE_TAXA.out.versions)
    // Acount for samples that had no metagenome
    SUMMARIZE_TAXA
        .out
        .sm_summary
        .join(ch_reads.map{ meta, reads -> meta }, by: 0, remainder: true)
        .map{ meta, summary -> if(summary){
                [ meta, summary ]
            }else{
                fwork = file(workflow.workDir).resolve("${meta.id}-sm_summary.txt")
                fwork.text = '100% unclassified'
                [ meta, fwork ]
            }
        }.set{ ch_sm_summary }

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: Manually Provided
        - Adds list of references provided in the samplesheet
    =============================================================================================================================
    */
    ch_refs_man
        .filter{ meta, ref -> ref }
        .transpose()
        .map{ meta, ref -> [ meta, file(ref).exists() ? getRefName(file(ref)) : "${ref}.fa.gz", file(ref).exists() ? file(ref) : null ] }
        .branch{ meta, ref_id, ref ->
            internal: ! ref
            external: ref
        }
        .set{ ch_man_refs }
    /* 
    =============================================================================================================================
        REFERENCE SELECTION: Accurate Mode
    =============================================================================================================================
    */
    if (params.ref_mode == "accurate"){
        // MODULE: Run Shovill
        SHOVILL (
            ch_reads
        )
        ch_versions = ch_versions.mix(SHOVILL.out.versions)

        // MODULE: Map contigs to the references
        MINIMAP2_ALIGN (
            SHOVILL
                .out
                .contigs
                .combine( ch_refs_fmt )
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        // MODULE: Select references using accurate approach
        REF_SELECT_ACCURATE (
            MINIMAP2_ALIGN.out.paf.combine( ch_refs_comp )
        )
        REF_SELECT_ACCURATE
            .out
            .ref_list
            .splitCsv(header: false, elem: 1)
            .transpose()
            .concat(ch_man_refs.internal.map{ meta, ref_id, ref -> [ meta, ref_id ] })
            .set{ ch_ref_list }
    }

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: Fast Mode
    =============================================================================================================================
    */
    if (params.ref_mode == "fast"){

        // MODULE: Reference sketches
        SM_SKETCH_REF (
            ch_refs_tar
        )
        ch_versions = ch_versions.mix(SM_SKETCH_REF.out.versions)

        // MODULE: Run Sourmash gather against the reference pool using the forward reads
        SM_GATHER_SELECT (
            SOURMASH_SKETCH.out.signatures.combine( SM_SKETCH_REF.out.signatures.map{ [ it ] } )
        )
        ch_versions = ch_versions.mix(SM_GATHER_SELECT.out.versions)

        // Set reference list
        SM_GATHER_SELECT
            .out
            .result
            .splitCsv(header: true, decompress: true)
            .filter{ meta, data -> data.f_match >= params.ref_genfrac && data.average_abund >= params.qc_depth }
            .map{ meta, data -> [ meta, data.filename.replaceAll('.sig$', '') ] }
            .concat(ch_man_refs.internal.map{ meta, ref_id, ref -> [ meta, ref_id.replaceAll('.sig$', '') ] })
            .set{ ch_ref_list }
    }

    /* 
    =============================================================================================================================
        REFERENCE SELECTION: Kitchen Sink Mode
        - Downloads NCBI assembly for all Sourmash hits and adds them to the reference list
    =============================================================================================================================
    */
    if (params.ref_mode == "kitchen-sink"){
        SOURMASH_GATHER
            .out
            .result
            .splitCsv(header: true, decompress: true)
            .map{ meta, data -> [ meta, ( data.name =~ /(GCA_|GFA)[0-9]+\.[0-9]+/ ).collect{ it[0] }[0] ] }
            .set{ ch_sm_acc }
        NCBI_DATASETS (
            ch_sm_acc.map{ meta, accession -> accession }.unique()
        )
        ch_versions = ch_versions.mix(NCBI_DATASETS.out.versions)
        ch_sm_acc
            .map{ meta, accession -> [ accession, meta ] }
            .combine(NCBI_DATASETS.out.assembly.transpose(), by: 0)
            .map{ accession, meta, assembly -> [ meta, getRefName(assembly), assembly ] }
            .set{ ch_ref_list }
        
        Channel.of('assembly,ref_info')
            .concat (NCBI_DATASETS.out.refsheet.splitCsv(header: true).map{ it.values().join(',') } )
            .collectFile(name: 'refsheet.csv', newLine: true, sort: 'index')
            .set{ ch_refsheet }
    }

    /*
    =============================================================================================================================
        REFERENCE SELECTION: Final Clean Up
        - Add reference paths for 'Accurate' and 'Fast' modes
        - Add manually supplied references
    =============================================================================================================================
    */
    if( params.ref_mode == 'fast' || params.ref_mode == 'accurate' ){
        // MODULE: Gather references from the tar file
        TAR2REFS (
            ch_ref_list
                .map{ meta, ref_id -> ref_id }
                .unique()
                .filter{ it != 'none_selected' }
                .collect()
                .map{ [ it ] }
                .combine(ch_refs_tar)
        )
        // Combine reference files with samples
        ch_ref_list
            .map{ meta, ref -> [ getRefName(file(ref)), meta ] }
            .combine(TAR2REFS.out.refs.flatten().map{ [ getRefName(it), it ] }, by: 0)
            .map{ ref_id, meta, ref -> [ meta, ref_id, ref ] }
            .set{ ch_ref_list }        
    }
    // Combine any manually supplied reference paths
    ch_ref_list.concat( ch_man_refs.external ).set{ ch_ref_list }
    // Perform final check on all references
    FINALIZE_REFS (
        ch_ref_list
            .map{ meta, ref_id, ref -> [ ref_id, ref ] }
            .unique()
    ).refs.set{ ch_refs_final }
    // Update finalized refs
    ch_ref_list
        .map{ meta, ref_id, ref -> [ ref_id, meta ] }
        .combine( ch_refs_final, by: 0 )
        .map{ ref_id, meta, ref -> [ meta, ref_id, ref ] }
        .set{ ch_ref_list }

    emit:
    ref_list   = ch_ref_list    // channel: [ val(sample_meta), val(ref_id), path(ref_path) ]
    refsheet   = ch_refsheet    // channel: [ path(refinfo.csv) ]
    sm_summary = ch_sm_summary  // channel: [ val(meta), val(result) ]
    versions   = ch_versions    // channel: [ versions.yml ]
}

// Remove file extension for references - allows reference names to have periods in them (GCA_7839185.1.fa.gz)
def getRefName(filename){
    def refname = filename
    def tokens  = filename.getName().tokenize('\\.')
    if(tokens.size() > 1){
        refname = filename.toString().endsWith('.gz') ? tokens[0..-3].join('.') : tokens[0..-2].join('.')
    }
    
    return refname
}