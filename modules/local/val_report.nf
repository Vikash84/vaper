process REPORT {
    label 'process_low'

    input:
    path results

    output:
    path "validation-report.csv",  emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    val_report.R
    """

}
