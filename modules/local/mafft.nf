process MAFFT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mafft=7.520"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.520--h031d066_3':
        'biocontainers/mafft:7.520--h031d066_3' }"

    input:
    tuple val(meta), val(ref_id), path(consensus), path(ref)

    output:
    tuple val(meta), val(ref_id), path("comb.aln"), emit: aln
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    """
    # concat reference contigs and create multi-fasta with the consensus
    echo ">Reference" > ref.fa && cat ${ref} | grep -v ">" | tr -d '\t\n\r ' >> ref.fa && echo "\n" >> ref.fa
    cat ref.fa ${consensus} > comb.fa

    mafft \\
        --thread ${task.cpus} \\
        ${args} \\
        comb.fa \\
        > comb.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """

}
