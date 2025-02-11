process SM2REFS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(sm_csv)

    output:
    tuple val(meta), path("${prefix}.sm-refs.csv"), emit: refs, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${meta.id}"
    """
    zcat ${}
    """
}
