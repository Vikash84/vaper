process REF_SELECT_ACCURATE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(ref_info), path(refs_comp)
    
    output:
    tuple val(meta), path("*.ref-summary.csv"), emit: ref_summary, optional: true
    tuple val(meta), path("*.ref-list.csv"),    emit: ref_list
    tuple val(meta), path("*.jpg"),             emit: plots, optional: true

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    if [ -s ${ref_info} ]
    then
        zcat ${refs_comp} > refs-comp.txt
        ref-select_accurate.R ${ref_info} refs-comp.txt "${prefix}" "${params.ref_genfrac}" "${params.ref_covplot ? 'TRUE' : 'FALSE' }"
    else
        echo 'none_selected' > ${prefix}.ref-list.csv
    fi
    """
}
