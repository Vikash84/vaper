process REPORT {
    label 'process_low'

    input:
    path results

    output:
    path "*",  emit: report, includeInputs: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    val_report.R
    """

}
