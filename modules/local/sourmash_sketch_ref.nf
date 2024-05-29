process SM_SKETCH_REF {
    label 'process_single'
    stageInMode 'copy'
    stageOutMode 'copy'

    conda "bioconda::sourmash=4.8.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.4--hdfd78af_0':
        'biocontainers/sourmash:4.8.4--hdfd78af_0' }"

    input:
    path sequence

    output:
    path "*.sig",        emit: signatures
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // required defaults for the tool to run, but can be overridden
    def args = task.ext.args ?: "dna --param-string 'scaled=1000,k=21,k=31,k=51,abund'"
    """
    sourmash sketch \\
        $args \\
        $sequence

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}