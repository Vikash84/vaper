process NEXTCLADE_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade%3A3.1.0--h9ee0642_0' :
        'biocontainers/nextclade:3.1.0--h9ee0642_0' }"

    input:
    tuple val(meta), val(ref_id), path(consensus), path(ref)
    path db_template

    output:
    tuple val(meta), val(ref_id), path("${prefix}-${ref_id}.len.tsv"), emit: tsv
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # concat reference contigs
    echo ">${ref_id}" > ref.fa
    cat ${ref} | grep -v ">" | tr -d '\n\t\r ' >> ref.fa

    nextclade \\
        run \\
        $args \\
        --jobs $task.cpus \\
        -r ref.fa \\
        -p ${db_template} \\
        --output-all ./ \\
        --output-basename "${prefix}-${ref_id}" \\
        ${consensus}

    # add assembly length to output
    echo "ASSEMBLY_LENGTH" > LENGTH && cat ${consensus} | grep -v ">" | tr -d '\t\n\r ' | wc -c >> LENGTH
    paste ${prefix}-${ref_id}.tsv LENGTH > ${prefix}-${ref_id}.len.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """
}
