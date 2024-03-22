process BWA_MEM {
    tag "${prefix}"
    label 'process_high'

    conda "bioconda::bwa"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0' :
        'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0' }"

    input:
    tuple val(meta), val(ref_id), path(ref), path(reads)

    output:
    tuple val(meta), val(ref_id), path('*.bam'),           emit: bam
    tuple val(meta), val(ref_id), path('*.coverage.txt'),  emit: coverage
    tuple val(meta), val(ref_id), path('*.stats.txt'),     emit: stats
    tuple val(meta), val(ref_id), path("*.read-list.txt"), emit: read_list
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}-${ref_id}"

    """
    # setup for pipe
    set -euxo pipefail

    # index the reference
    bwa index ${ref}

    # run bwa mem, select only mapped reads, convert to .bam, and sort
    bwa mem -t ${task.cpus} ${ref} ${reads[0]} ${reads[1]} | samtools view -b -F 4 - | samtools sort - > ${prefix}.bam

    # gather read stats
    samtools coverage ${prefix}.bam > ${prefix}.coverage.txt
    samtools stats --threads ${task.cpus} ${prefix}.bam > ${prefix}.stats.txt

    # get list of read headers for fastq extraction
    samtools view ${prefix}.bam | cut -f 1 | sort | uniq > ${prefix}.read-list.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem: \$(bwa 2>&1 | grep "Version" | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
