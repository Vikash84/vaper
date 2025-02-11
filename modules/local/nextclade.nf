process NEXTCLADE {
    tag "${prefix}"
    
    input:
    tuple val(ref_id), val(meta), path(consensus), path(ref)

    output:
    tuple val(meta), val(ref_id), env(ASSEMBLY_LENGTH), env(REF_LENGTH), path("nextclade.json"), emit: json
    path "versions.yml",                                                                         emit: versions
    
    script:
    prefix = "${meta.id}_${ref_id}"
    """
    # Run nextclade
    nextclade \\
        run \\
        -r ${ref} \\
        -O ./ \\
        ${consensus}

    # Get reference and assembly length
    ASSEMBLY_LENGTH=\$(zcat ${consensus} | grep -v ">" | tr -d '\t\n\r ' | wc -c)
    REF_LENGTH=\$(zcat ${ref} | grep -v ">" | tr -d '\t\n\r ' | wc -c)

    # version info
    echo -e "\\"${task.process}\\":\\n    nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')" > versions.yml
    """
}
