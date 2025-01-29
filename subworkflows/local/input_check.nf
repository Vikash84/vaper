//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { REFS2TAR          } from '../../modules/local/refs2tar'
include { FORMAT_REFS       } from '../../modules/local/format_refs'



workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    refs        // file: /path/to/ref-list.csv or ref-list.tar.gz

    main:
    /* 
    =============================================================================================================================
        SAMPLESHEET
    =============================================================================================================================
    */
    // Create samplesheet channel
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { ch_reads }

    /* 
    =============================================================================================================================
        REFERENCES
    =============================================================================================================================
    */
    // Get all references into compressed format
    if (! refs.toString().endsWith('.tar.gz')){
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
    // Further format references
    FORMAT_REFS (
        ch_refs_tar
    )

    emit:
    reads    = ch_reads                       // channel: [ val(meta), [ reads ] ]
    refs     = FORMAT_REFS.out.refs           // channel: [ path(refs.fa.gz), path(refs-comp.txt.gz), path(refs.tar.gz), path(refsheet.csv.gz) ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // define manual references
    ref = row.reference ? row.reference : []

    // define validation fields
    truth = row.truth ? file(row.truth, checkIfExists: true) : null
    inter_group = row.inter_group ? row.inter_group : null
    intra_group = row.intra_group ? row.intra_group : null

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ], ref, truth, inter_group, intra_group ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ], ref, truth, inter_group, intra_group ]
    }    
    return fastq_meta
}
