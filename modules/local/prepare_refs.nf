process PREPARE_REFS {
    tag "$refs_tar"
    label 'process_low'

    conda "bioconda::seqtk"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fusioncatcher-seqtk%3A1.2--hed695b0_1' :
        'staphb/seqtk:1.3' }"

    input:
    path refs_tar

    output:
    path "refs.fa",      emit: refs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    version = "1.0"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # extract references
    mkdir refs
    gzip -d ${refs_tar}
    tar -xvhf *.tar -C refs
    # rename fasta headers and combine into multi-fasta file
    refs=\$(ls -d refs/*/*)
    for ref in \${refs}
    do
        seqtk seq -C \${ref} > tmp.fa
        seqtk rename tmp.fa \${ref##*/} >> refs.fa
    done
    rm tmp.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | grep "Version" | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
