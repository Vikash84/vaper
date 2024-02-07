process SAMTOOLSTATS2TBL {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/jdj0303/waphl-viral-base:1.0.0"

    input:
    tuple val(meta), val(ref_id), path(stats)

    output:
    tuple val(meta), val(ref_id), path("*.csv"), emit: tbl

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # convert Fastp read summary to table
    samtoolstats2tbl.sh ${stats} > ${prefix}.${ref_id}.samtoolstats2tbl.csv
    """
}
