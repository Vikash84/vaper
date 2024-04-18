process FORMAT_REFS {
    label 'process_low'

    conda "bioconda::seqtk"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fusioncatcher-seqtk%3A1.2--hed695b0_1' :
        'staphb/seqtk:1.3' }"

    input:
    tuple path(ref_list), path(assemblies)

    output:
    path "refs.fa",       emit: refs
    path "refs-comp.txt", emit: refs_comp
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    version = "1.0"
    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    # rename fasta headers and combine into multi-fasta file
    for ref in \$(cat ${ref_list})
    do
        seqtk seq -C \${ref} > tmp.fa
        seqtk rename tmp.fa \${ref##*/} >> refs.fa
    done
    rm tmp.fa

    # get contig lengths
    seqtk comp refs.fa > refs-comp.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | grep "Version" | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
