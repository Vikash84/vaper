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
        .set { reads }

    // Create reference assembly channel
    if (params.refs.endsWith('.tar.gz')){
        // MODULE: Extract reference set from a tar.gz compressed file.
        TAR2REFS (
            file(params.refs)
        )

        TAR2REFS
            .out
            .refsheet
            .splitCsv(header: true)
            .map{ tuple(it.assembly, it.taxa, it.segment)  }
            .join(TAR2REFS.out.refs.flatten().map{ assembly -> [ file(assembly).getName(), assembly ] }, by: 0)
            .map{ name, taxa, segment, assembly -> [ taxa: taxa, segment: segment, assembly: assembly ] }
            .map{ it -> create_ref_channel(it) }
            .set{ refs }
    }else {
        Channel
        .fromPath(refs)
        .splitCsv( header: true, sep:',' )
        .map{ it -> create_ref_channel(it) }
        .set{ refs }
    }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    refs                                      // channel: [ val(meta), path(refs) ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

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
        fastq_meta = [ meta, [ file(row.fastq_1) ], truth, inter_group, intra_group ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ], truth, inter_group, intra_group ]
    }    
    return fastq_meta
}

// Function to get list of [ meta, assembly ]
def create_ref_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.taxa
    meta.single_end = true

    // check reference assemblies
    def ref = file(row.assembly, checkIfExists: true)
    if ( ref.getExtension() != 'gz' ) {
        exit 1, "ERROR: Reference assemblies must be gzip compressed. Please fix ${ref}"
    }

    // add path of reference assembly to the meta map
    def refs = [ meta, row.segment, ref ]

    return refs
}
