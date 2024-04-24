process METRICS {
    label 'process_low'

    conda "bioconda::mafft=7.520"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.520--h031d066_3':
        'biocontainers/mafft:7.520--h031d066_3' }"

    input:
    tuple val(metric), val(precision_type), path(fasta)

    output:
    tuple val(metric), val(precision_type), path("*.csv"),  emit: result
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    mafft --auto ${fasta} > ${fasta.baseName}.aln
    validate.sh ${fasta.baseName}.aln "${metric}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """

}
