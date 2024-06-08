process SUMMARYLINE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), val(ref_id), path(samtoolstats2tbl), path(nextclade), path(fastp2tbl), val(sm_summary), path(refsheet)

    output:
    tuple val(meta), path("*.summaryline.csv"), emit: summaryline

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # create summaryline
    summaryline.R "${fastp2tbl}" "${sm_summary}" "${samtoolstats2tbl}" "${nextclade}" "${prefix}" "${ref_id}" "${refsheet}"
    # rename using prefix and reference
    mv summaryline.csv "${prefix}-${ref_id}.summaryline.csv"
    """
}
