process MAPPED_FASTQ {
    tag "${meta.id}-${ref}"
    label 'process_single'

    conda "bioconda::seqtk=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    tuple val(meta), path(filter_list), val(ref), path(reads)

    output:
    path "*.gz"         , emit: reads
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = "fa"
    """
    # forward reads
    seqtk \\
        subseq \\
        $args \\
        ${reads[0]} \\
        $filter_list | \\
        gzip --no-name > ${prefix}-${ref}_R1.fastq.gz

    # reverse reads
    seqtk \\
        subseq \\
        $args \\
        ${reads[0]} \\
        $filter_list | \\
        gzip --no-name > ${prefix}-${ref}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
