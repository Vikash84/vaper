process PAIR {
    label 'process_low'
    tag "${id}"

    conda "bioconda::mash=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(id), path(seqs1, stageAs: "seqs1/*"), path(seqs2, stageAs: "seqs2/*"), val(metric), val(precision_type)

    output:
    tuple val(metric), val(precision_type), path("*.pair.fa"),   emit: fasta
    tuple val(metric), val(precision_type), path("*.pairs.txt"), emit: pairs
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    # combine seq1 and seq2 into single files
    cat ${seqs1} > seq1.fa
    cat ${seqs2} > seq2.fa

    # create contig pairs
    val_pair.sh seq1.fa seq2.fa "${id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """

}
