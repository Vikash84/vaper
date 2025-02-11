process EXPORT_ALN {
    tag "${prefix}"
    label 'process_high'
    stageInMode 'copy'

    input:
    tuple val(meta), val(ref_id), path(ref), path(bam)

    output:
    tuple val(meta), path("*", includeInputs: true)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}-${ref_id}"
    """
    # make sure reference is not compressed
    gzip -d * || true
    # index bam and reference
    samtools index ${bam}
    samtools faidx ${ ref.name.replaceAll('.gz$', '') }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
