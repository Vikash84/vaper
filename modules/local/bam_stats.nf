process BAM_STATS {
    tag "${prefix}"
    label 'process_high'

    conda "bioconda::bwa"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0' :
        'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0' }"

    input:
    tuple val(meta), val(ref_id), path(bam)

    output:
    tuple val(meta), val(ref_id), path('*.coverage.txt'),         emit: coverage
    tuple val(meta), val(ref_id), path('*.samtoolstats2tbl.csv'), emit: stats
    tuple val(meta), val(ref_id), path("*.read-list.txt"),        emit: read_list
    //path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}-${ref_id}"
    """
    # setup for pipe
    set -euxo pipefail

    # gather read stats
    samtools coverage ${bam} > ${prefix}.coverage.txt
    samtools stats --threads ${task.cpus} ${bam} > ${prefix}.stats.txt
    samtoolstats2tbl.sh ${prefix}.stats.txt > ${prefix}.samtoolstats2tbl.csv

    # get list of read headers for fastq extraction
    samtools view ${bam} | cut -f 1 | sort | uniq > ${prefix}.read-list.txt
    """
}
