process SUMMARYLINE {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/jdj0303/waphl-viral-base:1.0.0"

    input:
    tuple val(meta), val(ref_id), path(samtoolstats2tbl), path(assembly_stats), path(vcf), path(fastp2tbl), path(sm_summary)

    output:
    tuple val(meta), path("*.summaryline.csv"), emit: summaryline

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # determine number of SNVs
    if [ -f "${vcf}" ]
    then
        snvs=\$(cat ${vcf} | grep -v "#" | wc -l)
    else
        snvs="NA"
    fi

    # create summaryline
    summaryline.R "${fastp2tbl}" "${sm_summary}" "${samtoolstats2tbl}" "${assembly_stats}" "${prefix}" "${ref_id}"  "\$(echo \$snvs)"
    # rename using prefix and reference
    mv summaryline.csv "${prefix}-${ref_id}.summaryline.csv"
    """
}
