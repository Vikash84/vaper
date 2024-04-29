process JOIN {
    label 'process_low'

    container "docker.io/jdj0303/waphl-viral-base:1.0.0"

    input:
    tuple val(metric), val(precision_type), path(results), path(pairs)
 
    output:
    path "${metric}_${precision_type}*",  emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    [[ "${metric}" == "accuracy" ]] && header="Sample,Truth,TP,FP,FN,Result" && echo \${header} > "${metric}_${precision_type}_results.csv"
    [[ "${metric}" == "precision" ]] && header="Sample1,Sample2,TP,PP,Result" && echo \${header} > "${metric}_${precision_type}_results.csv"
    cat ${results} | grep -v "\${header}" | tr '\t' ',' >> "${metric}_${precision_type}_results.csv"
    cat ${pairs} | tr '\t' ',' > "${metric}_${precision_type}_pairs.csv"
    """

}
