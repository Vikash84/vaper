process COMBINE_SUMMARYLINES {
    label 'process_low'

    input:
    path summaries

    output:
    path "combined-summary.csv", emit: summary
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # combine summaries
    combine-summary.R "${params.qc_depth}" "${params.qc_genfrac}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combine-summary.R: \$(combine-summary.R version)
    END_VERSIONS
    """
}
