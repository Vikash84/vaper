process IVAR_CONSENSUS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ivar"
    container "staphb/ivar:1.4.2"

    input:
    tuple val(meta), val(ref), path(bam)

    output:
    tuple val(meta), val(ref), path('*.fa'), emit: consensus
    tuple val(meta), val(ref), path('*.csv'), emit: stats
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
    samtools mpileup -aa -A -Q 0 -d 0 ${bam} | \\
       ivar consensus \\
       -p ${prefix}-${ref} \\
       -m ${params.ivar_m} \\
       -n ${params.ivar_n} \\
       -t ${params.ivar_t} \\
       -q ${params.ivar_q} \\
       ${args}

    # gather stats
    assembly-stats.sh ${prefix}-${ref}.fa > ${prefix}-${ref}.assembly-stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools 2>&1 | grep "Version" | cut -f 2 -d ' ')
        ivar: \$(ivar version | head -n 1)
    END_VERSIONS
    """
}
