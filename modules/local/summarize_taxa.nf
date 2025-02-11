process SUMMARIZE_TAXA {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(sm_meta)
    
    output:
    tuple val(meta), path("*.taxa-summary.csv"), emit: sm_summary
    tuple val(meta), path("*.jpg"), emit: plot,  optional: true
    path "versions.yml",                         emit: versions



    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    #---- TAXA SUMMARY ----#
    # summarize taxa at >= 1X coverage and >= 1% relative abundance
    sm_summary.R "${sm_meta}" ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sm_summary.R: \$(sm_summary.R version)
    END_VERSIONS
    """
}
