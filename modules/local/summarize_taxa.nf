process SUMMARIZE_TAXA {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/jdj0303/vaper-base:beta"

    input:
    tuple val(meta), path(ref_info), path(sm_gather, stageAs: "sm_gather.csv.gz")
    path refs_comp

    output:
    tuple val(meta), path("*.ref-summary.csv"),  emit: ref_summary, optional: true
    tuple val(meta), path("*.ref-list.csv"),     emit: ref_list
    tuple val(meta), env(sm_summary),          emit: sm_summary
    tuple val(meta), path("*.jpg"),              emit: plots, optional: true


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
        ref-select_accurate.R ${ref_info} ${refs_comp} "${prefix}" "${params.gen_frac}" "${params.cov_plot ? 'TRUE' : 'FALSE' }"
    else
        touch ${prefix}.ref-list.csv        
    fi

    #---- TAXA SUMMARY ----#
    zcat sm_gather.csv.gz > sm_gather.csv 
    # summarize taxa at >= 1X coverage and >= 1% relative abundance
    sm_summary.R sm_gather.csv ${prefix}
    sm_summary=\$(cat ${prefix}.taxa-summary.csv)
    """
}
