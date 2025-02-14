#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
=============================================================================================================================
    PIPELINE PARAMETERS
=============================================================================================================================
*/
// Define results file
def result_file = file(params.outdir).resolve('metrics.csv')

// Check parameters
if( ! params.input ){ exit 1, "Please supply a samplesheet using the '--input' paramater." }
if( ! file(params.input).exists() ){ exit 1, "${params.input} does not exist!" }

/*
=============================================================================================================================
    MAIN WORKFLOW
=============================================================================================================================
*/
workflow {

    /*
    =============================================================================================================================
        LOAD SAMPLESHEET
    =============================================================================================================================
    */ 
    // Create channel from samplesheet
    Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map{ checkSampleSheet(it) } // Check samplesheet
        .set{ ch_input }
    // Get sequence lengths
    ch_input
        .map{ [ it.sample, it, it.assembly ] }
        .splitFasta(record: [id: true, seqString: true])
        .map{ sample, it, fa -> it + [ length: fa.seqString.length() ] }
        .set{ ch_input }
    // Create pairwise comparisons between replicates
    ch_input
        .map{ [ it.sample, [ it ] ] }
        .groupTuple(by:0)
        .map{ sample, data -> createPairs(data) }
        .transpose()
        .map{ sample, data -> data }
        .set{ ch_pairs }

    /*
    =============================================================================================================================
        COMPARE REPLICATE ASSEMBLIES
        - Gather initial metrics using Nextclade
        - Capture any Nextclade errors
        - Account for unaligned sequences on the terminal ends (Nextlade ignores these)
        - Calculate number of reproducibility
    =============================================================================================================================
    */
    // Gather initial metrics with Nextclade
    NEXTCLADE(
        ch_pairs.map{ [ it.sample, it.reps, it.assemblies ] }
    )
    // Load Nextclade results
    NEXTCLADE
        .out
        .json
        .splitJson()
        .set{ ch_nc }
    // Capture Nextclade errors
    ch_nc
        .filter{ sample, reps, item -> item.key == 'errors' }
        .map{ sample, reps, item -> [ sample, reps, item.value[0] ? item.value[0].errors[0] : null ]  }
        .set{ ch_ncErrors }
    // Extract Nextclade results and gather assembly metrics
    ch_nc
        .filter{ sample, reps, item -> item.key == 'results' }
        .join( ch_pairs.map{ [ it.sample, it.reps, it.lengths ] }, by: [0,1] )
        .join( ch_ncErrors, by: [0,1] )
        .map{ sample, reps, results, lengths, ncErrors -> [ sample: sample, reps: reps, ncErrors: ncErrors ] + gatherMetrics(results + [ lengths: lengths ]) }
        .set{ ch_assembly_metrics }
    
    /*
    =============================================================================================================================
        MERGE METRIC CHANNELS
    =============================================================================================================================
    */
    // Get a list of linked hashmap keys across all samples and fill any missing keys in each sample with 'null'
    ch_assembly_metrics
        .map{ it.keySet() }
        .collect()
        .map{ [ it.flatten().unique() ] }
        .combine( ch_assembly_metrics )
        .map{ cols, metrics -> fillMissingKeys(cols, metrics) }
        .set{ ch_assembly_metrics }
    // Export the metrics to a file
    ch_assembly_metrics
        .take(1)
        .map{ it.keySet().join(',') }
        .concat(ch_assembly_metrics.map{ it.values().collect{ v -> "\"${v}\"" }.join(',') })
        .collectFile(name: 'result.csv', sort: 'index', newLine: true)
        .subscribe{ file(it).copyTo(result_file) }
}

// Actions upon workflow completion
workflow.onComplete {
    println "\nWorkflow complete!\nResults can be found at ${result_file}"
}

/*
=============================================================================================================================
    FUNCTIONS
=============================================================================================================================
*/

def createPairs(data){
    def sample = data.sample[0][0]
    // Assign replicate numbers
    def reps = [:]
    data.eachWithIndex{ it, idx -> reps[ idx + 1 ] = it}
    def repKeys = reps.keySet().toList()
    if(repKeys.size() == 1){ exit 1, "${sample} does not have any replicates!" }
    // Generate non-repeating pairwise comparisons
    def pairs = []
    repKeys.eachWithIndex { rep1, index1 ->
        repKeys[(index1 + 1)..<repKeys.size()].each { rep2 ->
            pairs << [ sample: sample, 
                       reps: [ rep1, rep2 ], 
                       assemblies: [ reps[rep1].assembly[0], reps[rep2].assembly[0] ],
                       lengths: [ reps[rep1].length[0], reps[rep2].length[0] ] ]
        }
    }
    return [ sample, pairs ]

}
// Function for checking and formatting the input samplesheet
def checkSampleSheet(LinkedHashMap row){
    // Check for missing columns
    def missing_cols = [ "sample", "assembly", ].findAll { !row.containsKey(it) }
    if( missing_cols ){ exit 1, "The supplied samplesheet does not contain one or more required column: ${missing_cols.join(', ')}" }
    
    // Check if input files exist
    def assembly = file(row.assembly, checkIfExists: true)

    return [ sample: row.sample, assembly: assembly ]
}

// Function for gathering assembly metrics
def gatherMetrics(LinkedHashMap row){
    // Define metrics
    def metrics            = row.value[0] ? row.value[0] : null

    /*
    ==============================================================
        Default Nextclade Metrics
    ==============================================================
    */        
    def totalDeletions     = metrics ? metrics.totalDeletions : null
    def totalNonACGTNs     = metrics ? metrics.totalNonACGTNs : null

    /*
    ==============================================================
        Modified Nextclade Metrics
    ==============================================================
    */
    // Substitutions
    def ignoredSubstitutions = 0
    if( metrics ){ metrics.substitutions.each{ ignoredSubstitutions = ignoredSubstitutions + it.refNuc.count( 'N' ) + it.refNuc.count( '-' ) } }
    def totalSubstitutions = metrics ? metrics.totalSubstitutions - ignoredSubstitutions : null
    // Missing data - Nextclade does not count sites missing on the termini
    def termMissing        = metrics ? ( metrics.nucleotideComposition.'-' ? metrics.nucleotideComposition.'-' : 0 ) : null
    def totalMissing       = metrics ? metrics.totalMissing + termMissing : null
    // False insertions due to differences between the reference used
    def fauxInsertions     = metrics ? 0 : null
    if( metrics ){ metrics.insertions.each{ fauxInsertions = fauxInsertions + it.ins.count( 'N' ) + it.ins.count( '-' ) } } 
    def totalInsertions    = metrics ? metrics.totalInsertions - fauxInsertions : null

    /*
    ==============================================================
        New metrics
    ==============================================================
    */
    def alignLen           = metrics ? row.lengths[0] - totalMissing - ignoredSubstitutions : null
    def disagreement       = metrics ? totalSubstitutions + totalDeletions + totalInsertions + totalNonACGTNs : null
    def agreement          = metrics ? alignLen - disagreement : null
    def precision          = metrics ? 100 * agreement / alignLen : null

    /*
    ==============================================================
        Final Metrics
    ==============================================================
    */    
    return [ lengths: row.lengths,
             substitutions: totalSubstitutions,
             ignoredSubstitutions: ignoredSubstitutions,
             deletions: totalDeletions, 
             insertions: totalInsertions,
             fauxInsertions: fauxInsertions,
             missing: totalMissing,
             terminiMissing: termMissing,
             nonACGTNs: totalNonACGTNs,
             disagreements: disagreement,
             aligned: alignLen,
             precision: precision ]
}

// Function for filling in missing linked hashmap keys with null
def fillMissingKeys(cols, row){
    def missing_keys = cols.findAll { !row.containsKey(it) }
    if(missing_keys.size() > 0){
        def makeup_keys = [:]
        missing_keys.each{ makeup_keys[it] = null }
        row = row + makeup_keys
    }

    return row
}

/*
=============================================================================================================================
    PROCESSES
=============================================================================================================================
*/

process NEXTCLADE {
    container 'docker.io/nextstrain/nextclade:3.9.1'
    publishDir "${params.outdir}/nextclade/", mode: 'symlink'

    input:
    tuple val(sample), val(reps), path(assemblies)

    output:
    tuple val(sample), val(reps), path("${prefix}.aligned.fasta"), emit: fasta
    tuple val(sample), val(reps), path("${prefix}.json"),  emit: json
    
    script:
    prefix = "${sample}_${reps.join("-")}"
    """
    nextclade run -r ${assemblies[0]} -O ./ -n ${prefix} ${assemblies[1]}
    """
}