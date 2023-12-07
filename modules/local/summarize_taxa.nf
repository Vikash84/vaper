process SUMMARIZE_TAXA {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/jdj0303/waphl-viral-base:1.0.0"

    input:
    tuple val(meta), path(paf), path(k2_output), path(contig_cov)
    path refs_comp
    path ncbi_assembly_stats

    output:
    tuple val(meta), path("*.k2-summary.csv"),  emit: k2_summary
    tuple val(meta), path("*.ref-summary.csv"), emit: ref_summary, optional: true
    tuple val(meta), path("*.ref-list.csv"),    emit: ref_list

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # summarize kraken2 output
    k2_summary.R ${k2_output} ${contig_cov} ${prefix} ${ncbi_assembly_stats}
    # select references for consensus assembly generation
    ## check if the alignment file is empty
    if [ -s ${paf} ]
    then
        refs_summary.R ${paf} ${refs_comp} ${prefix} ${params.gen_frac}
    else
        echo "none_selected" > ${prefix}.ref-list.csv
    fi

    """
}