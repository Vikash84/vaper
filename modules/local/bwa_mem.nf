process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bwa"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0' :
        'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0' }"

    input:
    tuple val(meta), val(ref), path(reads)
    path  refs_tar

    output:
    tuple val(meta), path('*.bam'), val(ref), emit: bam
    tuple val(meta), path('*.txt'),           emit: coverage_detail
    tuple val(meta), path('*.txt'),           emit: coverage_summary


    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # setup for pipe
    set -euxo pipefail

    # extract references
    mkdir refs
    gzip -d ${refs_tar}
    tar -xvhf *.tar -C refs

    # index the reference
    bwa index refs/*/${ref}

    # run bwa mem, select only mapped reads, convert to .bam, and sort
    bwa mem -t 3 refs/*/${ref} ${reads[0]} ${reads[1]} | samtools view -b -F 4 - | samtools sort - > ${prefix}-${ref}.bam

    # gather read stats
    samtools coverage ${prefix}-${ref}.bam > ${prefix}-${ref}.coverage-detail.txt
    tail -n+2 ${prefix}-${ref}.coverage-detail.txt | awk 'BEGIN {OFS = ","} {tot_reads+=$4; cov_sum+=$6; depth_sum+=$7; baseq_sum+=$8; mapq_sum+=$9; n++ } END{ if (n > 0) print tot_reads,cov_sum/n,depth_sum/n,baseq_sum/n,mapq_sum/n;}' > ${prefix}-${ref}.coverage-summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem: \$(bwa 2>&1 | grep "Version" | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
