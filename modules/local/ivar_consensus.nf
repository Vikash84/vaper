process IVAR_CONSENSUS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ivar"
    container "staphb/ivar:1.3.1-titan"

    input:
    tuple val(meta), path(bam), val(ref)

    output:
    tuple val(meta), path('*.fa'), emit: consensus
    tuple val(meta), path('*.csv'), emit: stats
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # setup for pipe
    set -euxo pipefail

    # create mpilup and call consensus
    samtools mpileup -aa -A -Q 0 -d 0 ${bam} | ivar consensus -p ${prefix}-${ref} -m 10 -n N -t 0.5

    # gather stats
    assembly-stats.sh ${prefix}-${ref}.fa > ${prefix}-${ref}.assembly-stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools 2>&1 | grep "Version" | cut -f 2 -d ' ')
        ivar: \$(ivar version | head -n 1)
    END_VERSIONS
    """
}
