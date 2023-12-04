//
// Check input samplesheet and get read channels
//
include { PREPARE_REFS   } from '../../modules/local/prepare_refs'
include { SHOVILL        } from '../../modules/nf-core/shovill/main'
include { KRAKEN2        } from '../../modules/local/kraken2'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { SUMMARIZE_TAXA } from '../../modules/local/summarize_taxa'


workflow CLASSIFY_VIRUSES {
    take:
    reads // channel: [meta, reads]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Prepapre references
    //
    PREPARE_REFS (
        params.refs
    )

    //
    // MODULE: Run Shovill
    //
    SHOVILL (
        reads
    )
    ch_versions = ch_versions.mix(SHOVILL.out.versions.first())

    //
    // MODULE: Run Kraken2
    //
    KRAKEN2 (
        SHOVILL.out.contigs,
        params.k2db
    )

    //
    // MODULE: Map contigs to the references
    //
    MINIMAP2_ALIGN (
        SHOVILL.out.contigs,
        PREPARE_REFS.out.refs.map{ refs -> ["reference",refs] },
        false,
        false,
        false
    )

    //
    // MODULE: Summarize taxonomy
    //
    KRAKEN2
        .out
        .output
        .map{ meta, output, cov -> [meta, output, cov] }
        .set{ k2_output }
    MINIMAP2_ALIGN
        .out
        .paf
        .map{ meta, paf -> [meta, paf] }
        .join(k2_output, by: 0)
        .set{ taxa_files }

    SUMMARIZE_TAXA (
        taxa_files,
        params.ncbi_assembly_stats
    )

    SUMMARIZE_TAXA
        .out
        .ref_list
        .splitCsv(header: false)
        .map{ meta, ref -> [meta, ref.get(0)] }
        .combine(reads.map{ meta, reads -> [meta, reads] }, by: 0)
        .set{ ref_list }

    emit:
    ref_list   // channel: [ val(meta), val(reference), path(reads) ]
    k2_summary = SUMMARIZE_TAXA.out.k2_summary
    versions = SHOVILL.out.versions // channel: [ versions.yml ]
}