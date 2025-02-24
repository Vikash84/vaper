//
// Check input samplesheet and get read channels
//

include { FASTERQDUMP } from '../../modules/local/fasterqdump'
include { REFS2TAR    } from '../../modules/local/refs2tar'
include { FORMAT_REFS } from '../../modules/local/format_refs'

workflow PREPARE {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    refs        // file: /path/to/ref-list.csv or ref-list.tar.gz

    main:
    ch_versions = Channel.empty()

    /* 
    =============================================================================================================================
        SAMPLESHEET
    =============================================================================================================================
    */
    // Create samplesheet channel
    Channel
        .fromPath(samplesheet)
        .splitCsv ( header:true, sep:',', quote: '"' )
        .map { create_sample_channel(it) }
        .set { ch_reads }

    // Fix duplicate sample names
    ch_reads
        .map{ it -> [ it.meta, it ] }
        .groupTuple(by: 0)
        .map{ meta, its -> its.eachWithIndex{ it, idx -> it.meta.id = "${it.meta.id}_T${idx + 1}"
                                                           it }
                             its }
        .flatten()
        .set{ ch_reads }

    // MODULE: Download reads from SRA
    FASTERQDUMP(
        ch_reads.filter{ it -> it.sra }.map{ it -> [ it.meta, it.sra ] }
    )
    ch_versions = ch_versions.mix(FASTERQDUMP.out.versions)
    // Update reads channel with SRA reads
    ch_reads
        .map{ it -> [ it.meta, it ] }
        .join(FASTERQDUMP.out.reads.map{ meta, fastq_1, fastq_2 -> [ meta, [ fastq_1, fastq_2 ] ] }, by: 0, remainder: true)
        .map{ meta, it, reads -> it + [ sra_reads: reads ] }
        .map{ it -> it.fastq_12 = it.sra_reads ? it.sra_reads : it.fastq_12
                    it }
        .set{ ch_reads }

    /* 
    =============================================================================================================================
        REFERENCES
    =============================================================================================================================
    */
    // Get all references into compressed format
    if (! refs.toString().endsWith('.tar.gz') ){
        // MODULE: Compress references into tar.gz file
        REFS2TAR (
            Channel
                .fromPath(refs)
                .map{ [ it, it ] }
                .splitCsv( header: true, elem: 1 )
                .map{ refsheet, data -> [ refsheet, it.assembly ] }
        )
        REFS2TAR.out.tar.set{ ch_refs_tar }
    }else{
        Channel
            .fromPath(refs)
            .set{ ch_refs_tar }
    }
    // Combine all references into a single file
    FORMAT_REFS (
        ch_refs_tar
    )

    emit:
    reads    = ch_reads                       // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ], path(fastq_l), val(sra), val/path(references) ]
    refs     = FORMAT_REFS.out.refs           // channel: [ path(refs.fa.gz), path(refs-comp.txt.gz), path(refs.tar.gz), path(refsheet.csv.gz) ]
    versions = ch_versions                    // channel: [ versions.yml ]
}

/*
=============================================================================================================================
    SUBWORKFLOW FUNCTIONS
=============================================================================================================================
*/
//// Function for checking files
def check_file(f, patterns) {
    // Check that the file exists
    if( ! f.exists() ){ exit 1, "ERROR: ${f} does not exist!" }
    // Check that the file is compressed
    if( f.getExtension() == ".gz" ){ exit 1, "ERROR: ${f} is not gz compressed!" }
    // Check the file extension
    f_ext = f.getBaseName().tokenize('.')[-1]
    if( patterns.any{ ext -> ext == f_ext } ){ exit 1, "ERROR: ${f} does not match one of the expected file extensions (${patterns})" }
}

// Function to prepare the samplesheet
def create_sample_channel(LinkedHashMap row) {
    //// Sample Name
    // Check that the sample name was provided
    if( ! row.sample ){exit 1, "ERROR: Sample name not provided!\n${row}"}
    sample = row.sample
    // Check name length
    if(sample.length() > 20){ println "TIP: Short sample names are recommended (${sample})"}

    //// Reads
    fastq_1 = row.fastq_1 ? file(row.fastq_1) : null
    fastq_2 = row.fastq_2 ? file(row.fastq_2) : null
    fastq_l = row.fastq_l ? file(row.fastq_l) : null
    sra     = row.sra ? row.sra : null
    // Check that at least one read type was provided
    if( ! (fastq_1 && fastq_2) && ! fastq_l && ! sra ){ exit 1, "ERROR: Reads must be provided via the 'fastq_1' and 'fastq_2', 'fastq_l', or 'sra' columns for ${sample}." }
    // Check that only one read type was provided
    if( ((fastq_1 && fastq_2) && fastq_l) || ((fastq_1 && fastq_2) && sra) || (fastq_l && sra) ){ exit 1, "ERROR: Multiple read inputs provided - choose one" }
    // Check read formats
    if( fastq_1 ){ check_file( fastq_1, [".fq",".fastq"] ) }
    if( fastq_2 ){ check_file( fastq_2, [".fq",".fastq"] ) }
    if( fastq_l ){ check_file( fastq_l, [".fq",".fastq"] ) }
    if( sra ){ if( ! sra ==~ /^SRR\d{8}$/ ){ exit 1, "ERROR: ${sra} does not look like a SRA accession." } }

    //// Optional fields
    reference = row.reference ? row.reference.tokenize(';') : []

    //// Build output
    result = [ meta: [ id: sample, single_end: false ], fastq_12: [ fastq_1, fastq_2 ], fastq_l: fastq_l, sra: sra, reference: reference ]
    
    return result
}
