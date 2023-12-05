process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kraken2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    tuple val(meta), path(contigs)
    path  db_tar

    output:
    tuple val(meta), path('*report.txt'), emit: report
    tuple val(meta), path('*output.txt'), emit: output
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # extract Kraken2 database
    mkdir k2_db
    tar -xzvhf ${db_tar} -C k2_db
    # run Kraken2
    kraken2 \\
        --db k2_db \\
        --threads $task.cpus \\
        --report ${prefix}.k2.report.txt \\
        --output ${prefix}.k2.output.txt \\
        --use-names \\
        $args \\
        $contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
