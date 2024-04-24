//
// Check input samplesheet and get read channels
//

include { PAIR    } from '../../modules/local/val_pair'
include { METRICS } from '../../modules/local/val_metrics'
include { JOIN    } from '../../modules/local/val_join'
include { REPORT  } from '../../modules/local/val_report'


workflow VALIDATE {
    take:
    ch_samplesheet // channel: [ val(meta), path(truth), val(inter_group), val(intra_group) ]
    ch_consensus   // channel: [ val(meta), path(consensus) ]

    main:

    ch_versions = Channel.empty()

    /* 
    =============================================================================================================================
        CREATE DATASETS
    =============================================================================================================================
    */

    // Group consensus assemblies by meta.id - accounts for segmented genomes and/or co-infections
    ch_consensus.groupTuple(by: 0).set{ ch_consensus }

    // Accuracy dataset
    ch_samplesheet
        .map{ meta, truth, inter_group, intra_group -> [ meta, truth ] }
        .filter{ meta, truth -> truth }
        .combine(ch_consensus, by: 0)
        .map{ meta, truth, consensus -> [ meta.id, truth, consensus, "accuracy", null ] }
        .set{ ch_accuracy }

    // Inter-Assay Reproducibility
    ch_samplesheet
        .map{ meta, truth, inter_group, intra_group -> [ meta, inter_group ] }
        .filter{ meta, inter_group -> inter_group }
        .combine(ch_consensus, by: 0)
        .set{ ch_inter_group }
    ch_inter_group
        .combine(ch_inter_group, by: 1)
        .filter{ group, meta1, consensus1, meta2, consensus2 -> meta1 != meta2 }
        .map{ group, meta1, consensus1, meta2, consensus2 -> [ group, consensus1, consensus2 ] }
        .groupTuple(by: 0)
        .map{ group, consensus1, consensus2 -> [ group, consensus1.get(0), consensus2.get(0), "precision", "inter" ] }
        .set{ ch_inter_group }

    // Intra-Assay Reproducibility
    ch_samplesheet
        .map{ meta, truth, inter_group, intra_group -> [ meta, intra_group ] }
        .filter{ meta, intra_group -> intra_group }
        .combine(ch_consensus, by: 0)
        .set{ ch_intra_group }
    ch_intra_group
        .combine(ch_intra_group, by: 1)
        .filter{ group, meta1, consensus1, meta2, consensus2 -> meta1 != meta2 }
        .map{ group, meta1, consensus1, meta2, consensus2 -> [ group, consensus1, consensus2 ] }
        .groupTuple(by: 0)
        .map{ group, consensus1, consensus2 -> [ group, consensus1.get(0), consensus2.get(0), "precision", "intra" ] }
        .set{ ch_intra_group }

    ch_accuracy.concat(ch_inter_group).concat(ch_intra_group).set{ ch_datasets }

    PAIR (
        ch_datasets
    )

    METRICS (
        PAIR.out.fasta.transpose()
    )

    METRICS
        .out
        .result
        .groupTuple(by: [0,1])
        .join(PAIR.out.pairs.groupTuple(by: [0,1]), by: [0,1])
        .set{ ch_results }
    JOIN (
        ch_results
    )

    REPORT (
        JOIN.out.results.flatten().collect()
    )

    //emit:
    //samtoolstats2tbl = BAM_STATS.out.stats    // channel: [ val(meta), val(ref), path(stats)) ]
    //nextclade        = NEXTCLADE_RUN.out.tsv // channel: [ val(meta), val(ref), path(stats)) ]
}
