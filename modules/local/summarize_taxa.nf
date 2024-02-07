process SUMMARIZE_TAXA {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/jdj0303/waphl-viral-base:1.0.0"

    input:
    tuple val(meta), path(ref_info), path(sm_gather, stageAs: "sm_gather.csv.gz")
    path refs_comp

    output:
    tuple val(meta), path("*.ref-summary.csv"), emit: ref_summary, optional: true
    tuple val(meta), path("*.ref-list.csv"),    emit: ref_list
    tuple val(meta), path("*.taxa-summary.csv"),  emit: sm_summary

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    #---- REFERENCE SELECTION: FAST ----#
    if [ "${params.mode}" == "fast" ]
    then
        zcat ${ref_info} | awk -F ',' -v d=${params.avg_depth} -v g=${params.gen_frac} 'NR > 1 && \$3 >= g && \$6*2 >= d {print \$9}' > ${prefix}.ref-list.csv
        cp ${ref_info} ${prefix}.ref-summary.csv
    fi

    #---- REFERENCE SELECTION: ACCURATE ----#
    if [ "${params.mode}" == "accurate" ] && [ -s ${ref_info} ]
    then
        ref-select_accurate.R ${ref_info} ${refs_comp} ${prefix} ${params.gen_frac}
    else
        touch ${prefix}.ref-list.csv        
    fi

    #---- TAXA SUMMARY ----#
    # subset taxa with >= 10% abundance and 1X average coverage
    zcat ${sm_gather} | awk -F ',' 'NR > 1 && \$3 >= 0.1 && \$6*2 >= 1 {print \$10" ("int(100*\$3)"%/"int(\$6*2)"X)"}' | tr -d '"' | tr '\n' ';' > ${prefix}.taxa-summary.csv
    """
}