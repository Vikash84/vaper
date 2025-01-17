//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { TAR2REFS          } from '../../modules/local/tar2refs'


workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    refs        // file: /path/to/ref-list.csv or ref-list.tar.gz

    main:
    // Create samplesheet channel
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { ch_reads }

    // Create reference assembly channel
    if (refs.toString().endsWith('.tar.gz')){
        // MODULE: Extract reference set from a tar.gz compressed file.
        TAR2REFS (
            refs
        )
        TAR2REFS
            .out
            .refs
            .flatten()
            .map{ create_ref_channel([ it.getName(), it ]) }
            .set{ ch_refs }
        TAR2REFS
            .out
            .refsheet
            .set{ ch_refsheet }
    }else {
        Channel
            .fromPath(refs)
            .splitCsv( header: true )
            .map{ it -> create_ref_channel( [ file(it.assembly).getName(), it.assembly ] ) }
            .set{ ch_refs }
        Channel
            .fromPath(refs)
            .set{ ch_refsheet }
    }

    emit:
    reads    = ch_reads                       // channel: [ val(meta), [ reads ] ]
    refs     = ch_refs                        // channel: [ val(meta), path(refs), path(refsheet) ]
    refsheet = ch_refsheet                    // channel: [ path(refsheet) ] 
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // define manual references
    ref = row.reference ? row.reference : null

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

// Function to get list of [ meta, assembly ]
def create_ref_channel(row) {
    // extract elements
    def name     = row[0]
    def ref      = file(row[1], checkIfExists: true)
    // create meta map - mimics the samplesheet
    def meta        = [:]
    meta.id         = file(name).getSimpleName()
    meta.single_end = true

    // check reference assemblies
    if ( ref.getExtension() != 'gz' ) {
        exit 1, "ERROR: Reference assemblies must be gzip compressed. Please fix ${ref}"
    }
    if ( ref.baseName.toString().chars().filter(it -> it == '.').count() > 2 ) {
        exit 1, "ERROR: Reference assemblies cannot contain periods in their name. Please fix ${ref}"
    }

    // add path of reference assembly to the meta map
    def refs = [ meta, ref ]

    return refs
}
