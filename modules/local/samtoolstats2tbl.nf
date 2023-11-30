process SAMTOOLSTATS2TBL {
    label 'process_low'

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*.csv"), emit: tbl

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # convert Fastp read summary to table
    samtoolstats2tbl.sh ${stats} > ${prefix}.samtoolstats2tbl.csv
    """
}
