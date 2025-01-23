process SUMMARIZE_TAXA {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(ref_info), path(sm_gather, stageAs: "sm_gather.csv.gz"), path(sm_meta)
    path refs_comp

    output:
    tuple val(meta), path("*.ref-summary.csv"),          emit: ref_summary, optional: true
    tuple val(meta), path("*.ref-list.csv"),             emit: ref_list
    tuple val(meta), path("*.taxa-summary.csv"), emit: sm_summary
    tuple val(meta), path("*.jpg"),                      emit: plots, optional: true
    path "versions.yml",                                 emit: versions



    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    #---- REFERENCE SELECTION: FAST ----#
    if [ "${params.ref_mode}" == "fast" ]
    then
        zcat ${ref_info} | awk -F ',' -v d=${params.qc_depth} -v g=${params.ref_genfrac} 'NR > 1 && \$3 >= g && \$6*2 >= d {print \$9}' > ${prefix}.ref-list.csv
        zcat ${ref_info} > ${prefix}.ref-summary.csv

    #---- REFERENCE SELECTION: ACCURATE ----#
    elif [ "${params.ref_mode}" == "accurate" ] && [ -s ${ref_info} ]
    then
        zcat ${refs_comp} > refs-comp.txt
        ref-select_accurate.R ${ref_info} refs-comp.txt "${prefix}" "${params.ref_genfrac}" "${params.ref_covplot ? 'TRUE' : 'FALSE' }"
    else
        echo 'none_selected' > ${prefix}.ref-list.csv
    fi

    #---- TAXA SUMMARY ----#
    # summarize taxa at >= 1X coverage and >= 1% relative abundance
    sm_summary.R ${sm_meta} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sm_summary.R: \$(combine-summary.R version)
    END_VERSIONS
    """
}
