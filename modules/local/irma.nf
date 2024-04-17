process IRMA {
    tag "${prefix}"
    label 'process_high'
    stageInMode 'copy'
    
    container "docker.io/jdj0303/vaper-base:beta"

    input:
    tuple val(meta), path(refs), path(reads)
    path module_template

    output:
    tuple val(meta), path('results/*.fasta'),  emit: consensus
    tuple val(meta), path('results/*.bam'),    emit: bam
    tuple val(meta), path('results/logs/'),    emit: logs
    tuple val(meta), path('results/figures/'), emit: figures
    //path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}"

    """
    # set permissions
    echo "Checking Permissions:"
    chown -R \$(id -u):\$(id -g) flu-amd/ || true
    chmod -R +wrx flu-amd/ || true
    ls -la flu-amd/
    echo "\n"
    # combine references into single file
    cat ${refs} > flu-amd/IRMA_RES/modules/ORG/reference/consensus.fasta
    # run IRMA
    flu-amd/IRMA ORG ${reads[0]} ${reads[1]} results
    """
}
