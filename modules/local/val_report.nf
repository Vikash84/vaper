process REPORT {
    label 'process_low'

    container "docker.io/jdj0303/waphl-viral-base:1.0.0"

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
