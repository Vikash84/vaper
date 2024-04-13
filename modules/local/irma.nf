process IRMA {
    tag "${prefix}"
    label 'process_high'

    container "docker.io/jdj0303/waphl-viral-base:1.0.0"

    input:
    tuple val(meta), path(refs), path(reads)
    path irma

    output:
    tuple val(meta), path('results/*.fasta'),  emit: consensus
    tuple val(meta), path('results/*.bam'),    emit: bam
    //path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}"

    """
    # combine references into single file
    cat ${refs} > flu-amd/IRMA_RES/modules/ORG/reference/consensus.fasta
    # run IRMA
    flu-amd/IRMA ORG ${reads[0]} ${reads[1]} results
    """
}
