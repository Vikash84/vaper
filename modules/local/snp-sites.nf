process SNPSITES {
    label 'process_medium'

    conda "bioconda::snp-sites=2.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snp-sites:2.5.1--hed695b0_0' :
        'biocontainers/snp-sites:2.5.1--hed695b0_0' }"

    input:
    tuple val(meta), val(ref), path(aln)

    output:
    tuple val(meta), val(ref), path("*.vcf"), emit: vcf
    path "versions.yml" ,                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    snp-sites \\
        $args \\
        -v \\
        -o ${prefix}-${ref}.vcf \\
        $aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsites: \$(snp-sites -V 2>&1 | sed 's/snp-sites //')
    END_VERSIONS
    """
}
