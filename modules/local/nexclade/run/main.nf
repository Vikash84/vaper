process NEXTCLADE_RUN {
    tag "${prefix}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade%3A3.1.0--h9ee0642_0' :
        'biocontainers/nextclade:3.1.0--h9ee0642_0' }"

    input:
    tuple val(meta), val(ref_id), path(consensus), path(ref), path(config)

    output:
    tuple val(meta), val(ref_id), path("${prefix}.metrics.tsv"), emit: tsv
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}-${ref_id}"
    """
    nextclade \\
        run \\
        $args \\
        --jobs $task.cpus \\
        -r ${ref} \\
        -p ${config} \\
        --output-all ./ \\
        --output-basename "${prefix}" \\
        ${consensus}

    # add extra metrics
    echo "ASSEMBLY_LENGTH" > CON_LENGTH && cat ${consensus} | grep -v ">" | tr -d '\t\n\r ' | wc -c >> CON_LENGTH
    echo "REF_LENGTH" > REF_LENGTH && cat ${ref} | grep -v ">" | tr -d '\t\n\r ' | wc -c >> REF_LENGTH
    echo "ASSEMBLY_TERMINAL_GAPS" > TERMINAL_GAPS && cat ${prefix}.aligned.fasta | grep -oE '^[-]+|[-]+\$' | tr -d '\n\r\t ' | wc -c >> TERMINAL_GAPS || true
    paste ${prefix}.tsv CON_LENGTH REF_LENGTH TERMINAL_GAPS > ${prefix}.metrics.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """
}
