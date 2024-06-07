process COMBINE_SUMMARYLINES {
    label 'process_low'

    input:
    path summaries

    output:
    path "combined-summary.csv", emit: summary

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # combine summaries
    combine-summary.R "${params.qc_depth}" "${params.qc_genfrac}"
    """
}
