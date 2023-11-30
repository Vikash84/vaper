process SUMMARYLINE {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/jdj0303/waphl-viral-base:1.0.0"

    input:
    tuple val(meta), val(ref), path(k2_summary), path(samtoolstats2tbl), path(assembly_stats), path(fastp2tbl)

    output:
    tuple val(meta), path("*.summaryline.csv"), emit: summaryline

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # create summaryline
    summaryline.R "${fastp2tbl}" "${k2_summary}" "${samtoolstats2tbl}" "${assembly_stats}" "${prefix}" "${ref}"
    # rename using prefix and reference
    mv summaryline.csv ${prefix}-${ref}.summaryline.csv
    """
}
