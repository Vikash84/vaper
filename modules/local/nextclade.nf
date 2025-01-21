process NEXTCLADE {
    input:
    tuple val(ref_id), val(meta), path(consensus), path(ref)

    output:
    tuple val(meta), val(ref_id), path("${prefix}.metrics.tsv"), emit: tsv
    path "versions.yml"                                        , emit: versions
    
    script:
    prefix = "${meta.id}-${ref_id}"
    """
    # Run nextclade
    nextclade \\
        run \\
        -r ${ref} \\
        -O ./ \\
        ${consensus}

    # add extra metrics
    echo -e "ASSEMBLY_LENGTH\tREF_LENGTH\tASSEMBLY_TERMINAL_GAPS" > extra_metrics.tsv
    cat ${consensus} | grep -v ">" | tr -d '\t\n\r ' | wc -c > ASSEMBLY_LENGTH
    cat ${ref} | grep -v ">" | tr -d '\t\n\r ' | wc -c > ASSEMBLY_LENGTH
    cat ${prefix}.aligned.fasta | grep -oE '^[-]+|[-]+\$' | tr -d '\n\r\t ' | wc -c > TERMINAL_GAPS || true
    paste ASSEMBLY_LENGTH ASSEMBLY_LENGTH TERMINAL_GAPS >> extra_metrics.tsv
    paste nextclade.tsv extra_metrics.tsv > ${prefix}.metrics.tsv

    # version info
    echo -e "\\"${task.process}\\":\\n    nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')" > versions.yml
    """
}
