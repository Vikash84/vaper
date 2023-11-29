process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kraken2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    tuple val(meta), path(contigs)
    path  db

    output:
    tuple val(meta), path('*report.txt'),                            emit: report
    tuple val(meta), path('*output.txt'), path('*cov.csv'), emit: output
    path "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # run Kraken2
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.k2.report.txt \\
        --output ${prefix}.k2.output.txt \\
        --use-names \\
        $args \\
        $contigs

    # create summary of contig stats from Shovill output
    cat $contigs | grep ">" | tr -d '>' | cut -f 1,3 -d ' ' | sed 's/ ...=/,/g' > ${prefix}.cov.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
