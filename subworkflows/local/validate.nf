//
// Check input samplesheet and get read channels
//

include { PAIR    } from '../../modules/local/val_pair'
include { GATHER  } from '../../modules/local/val_gather'
include { REPORT  } from '../../modules/local/val_report'


workflow VALIDATE {
    take:
    ch_samplesheet // channel: [ val(meta), path(truth), val(inter_group), val(intra_group) ]
    ch_consensus   // channel: [ val(meta), path(consensus) ]
    ch_summary     // channel: [ path(summary) ]

    main:

    ch_versions = Channel.empty()

    /* 
    =============================================================================================================================
        FILTER SAMPLES PASSING QC
    =============================================================================================================================
    */

    ch_summary
        .splitCsv(header: true)
        .filter{ it.ASSEMBLY_QC == "PASS" }
        .map{ [[id: it.ID, single_end: false], it.ID+"_"+it.REFERENCE+".fa"] }
        .join(ch_consensus.map{meta, con -> [ meta, con.name, con ]}, by: [0,1])
        .map{ meta, name, consensus -> [ meta, consensus ] }
        .set{ ch_consensus }  

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

    /* 
    =============================================================================================================================
        DETERMINE PAIRWISE COMPARISONS
    =============================================================================================================================
    */
    PAIR (
        ch_datasets
    )

    /* 
    =============================================================================================================================
        GATHER METRICS
    =============================================================================================================================
    */
    GATHER (
        PAIR.out.fasta.transpose()
    )

    /* 
    =============================================================================================================================
        JOIN ALL METRICS/PAIRS INTO SINGLE CHANNEL
    =============================================================================================================================
    */
    // Accuracy
    GATHER
        .out
        .result
        .filter{ metric, type, results -> metric == "accuracy" }
        .map{ metric, type, results -> results }
        .splitText()
        .filter(line -> line != "Sample,Truth,Total,Compared,Correct,Incorrect,Accuracy\n")
        .collectFile(name: "accuracy_results.csv")
        .set{ ch_acc_res }
    
    PAIR
        .out
        .pairs
        .filter{ metric, type, pairs -> metric == "accuracy" }
        .map{ metric, type, pairs -> pairs }
        .splitText()
        .collectFile(name: "accuracy_pairs.tsv")
        .set{ ch_acc_pair }
    
    // Intra-assay Precision
    GATHER
        .out
        .result
        .filter{ metric, type, results -> metric == "precision" && type == "inter" }
        .map{ metric, type, results -> results }
        .splitText()
        .filter(line -> line != "Sample1,Sample2,Total,Agreements,Disagreements,Precision\n")
        .collectFile(name: "precision_inter_results.csv")
        .set{ ch_prec_inter_res }
        
    
    PAIR
        .out
        .pairs
        .filter{ metric, type, pairs -> metric == "precision" && type == "inter" }
        .map{ metric, type, pairs -> pairs }
        .splitText()
        .collectFile(name: "precision_inter_pairs.tsv")
        .set{ ch_prec_inter_pair }
    
    // Intra-assay Precision
    GATHER
        .out
        .result
        .filter{ metric, type, results -> metric == "precision" && type == "intra" }
        .map{ metric, type, results -> results }
        .splitText()
        .filter(line -> line != "Sample1,Sample2,Total,Agreements,Disagreements,Precision\n")
        .collectFile(name: "precision_intra_results.csv")
        .set{ ch_prec_intra_res }
    
    PAIR
        .out
        .pairs
        .filter{ metric, type, pairs -> metric == "precision" && type == "intra" }
        .map{ metric, type, pairs -> pairs }
        .splitText()
        .collectFile(name: "precision_intra_pairs.tsv")
        .set{ ch_prec_intra_pair }

    // Combine all
    ch_acc_res
        .concat(ch_acc_pair)
        .concat(ch_prec_inter_res)
        .concat(ch_prec_inter_pair)
        .concat(ch_prec_intra_res)
        .concat(ch_prec_intra_pair)
        .flatten()
        .collect()
        .set{ ch_results }

    /* 
    =============================================================================================================================
        GENERATE REPORT
    =============================================================================================================================
    */
    REPORT (
       ch_results
    )

    //emit:
    //samtoolstats2tbl = BAM_STATS.out.stats    // channel: [ val(meta), val(ref), path(stats)) ]
    //nextclade        = NEXTCLADE_RUN.out.tsv // channel: [ val(meta), val(ref), path(stats)) ]
}
