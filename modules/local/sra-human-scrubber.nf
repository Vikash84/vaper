process SRA_HUMAN_SCRUBBER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::sra-human-scrubber"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber%3A2.2.1--hdfd78af_0' :
        'ncbi/sra-human-scrubber' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.scrubbed.fastq.gz'), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    prefix = task.ext.prefix ?: "${meta.id}"
    version = "2.2.1"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # run sra-human-scrubber
    ## forward reads
    scrub.sh \
      -i ${reads[0]} \
      -o ${prefix}_R1.scrubbed.fastq \
      -x -r

    ## reverse reads
    scrub.sh \
      -i ${reads[1]} \
      -o ${prefix}_R2.scrubbed.fastq \
      -x -r

    ## compress
    gzip *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-human-scrubber: ${version}
    END_VERSIONS
    """
}
