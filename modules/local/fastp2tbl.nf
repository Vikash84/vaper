process FASTP2TBL {
    label 'process_low'

    input:
    tuple val(meta), path(fastp_json)

    output:
    path "*.csv", emit: tbl

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # convert Fastp read summary to table
    fastp2tbl.sh ${fastp_json} > ${prefix}.fastp2tbl.csv
    """
}
