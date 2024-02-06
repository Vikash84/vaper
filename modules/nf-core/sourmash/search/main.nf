process SOURMASH_SEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.4--hdfd78af_0':
        'biocontainers/sourmash:4.8.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(signature)
    path(database)

    output:
    tuple val(meta), path('*.csv.gz')             , optional:true, emit: result
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sourmash search \\
        $args \\
        --output ${prefix}.csv.gz \\
        ${signature} \\
        ${database}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}
